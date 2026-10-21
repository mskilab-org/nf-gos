# CBS model signal and provenance. This file is sourced beside each R driver.

cbs_numeric_field <- function(gr, field, label) {
    if (!methods::is(gr, "GRanges"))
        stop(label, " must be a GRanges object")
    if (!(field %in% names(S4Vectors::mcols(gr))))
        stop(label, " is missing required field '", field, "'")
    value <- S4Vectors::mcols(gr)[[field]]
    if (!is.numeric(value) || length(value) != length(gr))
        stop(label, " field '", field, "' must be numeric")
    value
}

cbs_normal_index <- function(tumor, normal) {
    key <- function(gr) paste(as.character(GenomicRanges::seqnames(gr)),
                              BiocGenerics::start(gr), BiocGenerics::end(gr),
                              as.character(GenomicRanges::strand(gr)), sep = ":")
    tumor.key <- key(tumor)
    normal.key <- key(normal)
    if (anyDuplicated(tumor.key) || anyDuplicated(normal.key))
        stop("Tumor and normal coverage must have unique bin coordinates")
    index <- match(tumor.key, normal.key)
    if (length(tumor) != length(normal) || anyNA(index))
        stop("Tumor and normal coverage must have exactly matching bin coordinates")
    index
}

prepare_cbs_coverage <- function(tumor, normal = NULL, ncn = NULL,
                                 input_field = "foreground") {
    counts <- cbs_numeric_field(tumor, input_field, "Tumor coverage")
    if (!length(counts)) stop("Tumor coverage has no bins")
    weights <- as.numeric(BiocGenerics::width(tumor))
    if (any(weights <= 0)) stop("Coverage bins must have positive widths")
    mode <- if (is.null(normal)) "tumor_only" else "paired"
    if (mode == "paired") {
        normal.counts <- cbs_numeric_field(normal, input_field, "Normal coverage")
        index <- cbs_normal_index(tumor, normal)
        if (!is.numeric(ncn) || length(ncn) != length(normal))
            stop("Paired coverage requires measured normal copy number for every bin")
        normal.counts <- normal.counts[index]
        ncn <- ncn[index]
        denominator <- normal.counts / ncn
        valid <- is.finite(counts) & counts >= 0 &
                 is.finite(normal.counts) & normal.counts > 0 &
                 is.finite(ncn) & ncn > 0 &
                 is.finite(denominator) & denominator > 0
        if (!any(valid)) stop("Paired coverage has no valid model bins")
        tumor.mean <- stats::weighted.mean(counts[valid], weights[valid])
        denominator.mean <- stats::weighted.mean(denominator[valid], weights[valid])
        if (!is.finite(tumor.mean) || tumor.mean <= 0 || !is.finite(denominator.mean))
            stop("Paired coverage requires a finite positive mean tumor signal")
        # Deterministic whole-genome rescale replaces the legacy random subset.
        rescale <- denominator.mean / tumor.mean
        ratio <- rep(NA_real_, length(tumor))
        ratio[valid] <- counts[valid] / denominator[valid] * rescale
        ratio.mean <- stats::weighted.mean(ratio[valid], weights[valid])
        if (any(!is.finite(ratio[valid])) || !is.finite(ratio.mean) || ratio.mean <= 0)
            stop("Paired coverage ratio cannot be represented on a finite positive scale")
        foreground.scale <- tumor.mean / ratio.mean
        foreground <- ratio * foreground.scale
        if (!is.finite(foreground.scale) || foreground.scale <= 0 ||
            any(!is.finite(foreground[valid])))
            stop("Paired coverage foreground cannot be represented on a finite positive scale")
        tumor$norm.counts <- normal.counts
    } else {
        if (!is.null(ncn))
            stop("Tumor-only coverage has no measured normal segmentation; baseline must be 2")
        ncn <- rep(2, length(tumor))
        valid <- is.finite(counts) & counts >= 0 & is.finite(ncn) & ncn > 0
        foreground <- counts
        if (any(!valid)) foreground[!valid] <- NA_real_
        ratio <- foreground
        foreground.scale <- 1
        # Normal counts are an audit field only for paired coverage.
        tumor$norm.counts <- NULL
    }
    tumor$tum.counts <- counts
    tumor$ncn <- ncn
    tumor$ratio <- ratio
    tumor$foreground <- foreground
    list(cov = tumor, mode = mode, input_field = input_field,
         foreground_scale = foreground.scale)
}

cbs_log_signal <- function(cov) {
    signal <- cbs_numeric_field(cov, "ratio", "CBS coverage")
    ncn <- cbs_numeric_field(cov, "ncn", "CBS coverage")
    # Zero is a valid model observation, but not a finite logratio for DNAcopy.
    # Keep historical CBS seg.mean units: foreground differs by one global scalar.
    index <- which(is.finite(signal) & signal > 0 & is.finite(ncn) & ncn > 0)
    if (!length(index)) stop("CBS coverage has no positive finite bins for log segmentation")
    list(index = index, value = log(signal[index]))
}

cbs_segment_ncn <- function(seg, cov) {
    ncn <- cbs_numeric_field(cov, "ncn", "CBS coverage")
    keep <- which(is.finite(ncn) & ncn > 0)
    hits <- GenomicRanges::findOverlaps(seg, cov[keep], ignore.strand = TRUE)
    query <- S4Vectors::queryHits(hits)
    subject <- keep[S4Vectors::subjectHits(hits)]
    weights <- pmin(BiocGenerics::end(seg)[query], BiocGenerics::end(cov)[subject]) -
               pmax(BiocGenerics::start(seg)[query], BiocGenerics::start(cov)[subject]) + 1
    result <- rep(NA_real_, length(seg))
    if (length(query)) {
        totals <- rowsum(cbind(value = ncn[subject] * weights, weight = weights),
                        query, reorder = FALSE)
        result[as.integer(rownames(totals))] <- totals[, "value"] / totals[, "weight"]
    }
    result
}

cbs_validate_info <- function(info, label) {
    if (!is.list(info) || !identical(info$version, 1L) ||
        !is.character(info$id) || length(info$id) != 1L || is.na(info$id) ||
        !grepl("^cbs-v1-[a-f0-9]{32}$", info$id) ||
        !is.character(info$mode) || length(info$mode) != 1L ||
        !(info$mode %in% c("paired", "tumor_only")) ||
        !identical(info$field, "foreground") ||
        !is.character(info$input_field) || length(info$input_field) != 1L ||
        is.na(info$input_field) || !nzchar(info$input_field) ||
        !is.numeric(info$foreground_scale) || length(info$foreground_scale) != 1L ||
        !is.finite(info$foreground_scale) || info$foreground_scale <= 0)
        stop(label, " lacks valid CBS model-coverage provenance; regenerate a complete CBS bundle")
    invisible(info)
}

cbs_coverage_info <- function(cov) {
    if (!methods::is(cov, "GRanges")) stop("CBS coverage must be a GRanges object")
    info <- S4Vectors::metadata(cov)$cbs_coverage
    cbs_validate_info(info, "Coverage")
    foreground <- cbs_numeric_field(cov, "foreground", "CBS coverage")
    ratio <- cbs_numeric_field(cov, "ratio", "CBS coverage")
    counts <- cbs_numeric_field(cov, "tum.counts", "CBS coverage")
    ncn <- cbs_numeric_field(cov, "ncn", "CBS coverage")
    valid <- is.finite(foreground)
    if (!any(valid) || any(is.infinite(foreground)) || any(is.infinite(ratio)) ||
        any(valid != is.finite(ratio)) || any(foreground[valid] < 0) ||
        any(!is.finite(ncn[valid]) | ncn[valid] <= 0) ||
        any(!is.finite(counts[valid]) | counts[valid] < 0) ||
        !isTRUE(all.equal(foreground[valid], ratio[valid] * info$foreground_scale,
                         tolerance = 1e-10, check.attributes = FALSE)))
        stop("CBS coverage has inconsistent foreground, ratio, tumor counts or normal baseline")
    if (info$mode == "paired") {
        normal <- cbs_numeric_field(cov, "norm.counts", "CBS coverage")
        denominator <- normal / ncn
        expected.valid <- is.finite(counts) & counts >= 0 & is.finite(normal) & normal > 0 &
                          is.finite(ncn) & ncn > 0 & is.finite(denominator) & denominator > 0
        weights <- as.numeric(BiocGenerics::width(cov))[valid]
        tumor.mean <- stats::weighted.mean(counts[valid], weights)
        denominator.mean <- stats::weighted.mean(denominator[valid], weights)
        expected.ratio <- counts[valid] / denominator[valid] * denominator.mean / tumor.mean
        if (any(valid != expected.valid) ||
            !isTRUE(all.equal(ratio[valid], expected.ratio, tolerance = 1e-10,
                             check.attributes = FALSE)) ||
            !isTRUE(all.equal(stats::weighted.mean(foreground[valid], weights),
                              tumor.mean, tolerance = 1e-10)))
            stop("Paired CBS coverage has inconsistent ratio or valid-bin weighted tumor depth")
    } else {
        expected.valid <- is.finite(counts) & counts >= 0 & is.finite(ncn) & ncn > 0
        if ("norm.counts" %in% names(S4Vectors::mcols(cov)) || any(!is.finite(ncn) | ncn != 2) ||
            !identical(info$foreground_scale, 1) || any(valid != expected.valid) ||
            !isTRUE(all.equal(foreground[valid], counts[valid], tolerance = 0,
                             check.attributes = FALSE)))
            stop("Tumor-only CBS coverage must preserve the input signal unchanged")
    }
    info
}

cbs_bundle_id <- function(cov, seg, nseg, info) {
    # R serialization embeds the writer's R version, even with version = 2.
    # Hash a version-independent stream instead: UTF-8 strings with NUL delimiters,
    # big-endian vector lengths/doubles, and only model-relevant range fields.
    file <- tempfile("cbs-coverage-", fileext = ".bin")
    connection <- file(file, open = "wb")
    on.exit({
        if (!is.null(connection)) close(connection)
        unlink(file)
    }, add = TRUE)
    write.length <- function(n) writeBin(as.double(n), connection, size = 8L, endian = "big")
    write.text <- function(value) {
        write.length(length(value))
        writeBin(enc2utf8(as.character(value)), connection)
    }
    write.numeric <- function(value) {
        value <- as.double(value)
        value[is.na(value)] <- NA_real_
        value[!is.na(value) & value == 0] <- 0 # Canonicalize signed zero.
        write.length(length(value))
        writeBin(value, connection, size = 8L, endian = "big")
    }
    write.ranges <- function(gr, fields) {
        write.text(as.character(GenomicRanges::seqnames(gr)))
        write.numeric(BiocGenerics::start(gr))
        write.numeric(BiocGenerics::end(gr))
        write.text(as.character(GenomicRanges::strand(gr)))
        for (field in fields) {
            write.text(field)
            if (field %in% names(S4Vectors::mcols(gr)))
                write.numeric(cbs_numeric_field(gr, field, "CBS bundle"))
            else write.length(-1)
        }
    }
    write.text(c("cbs-coverage-canonical-v1", info$mode, "foreground", info$input_field))
    write.numeric(info$foreground_scale)
    write.text("cov")
    write.ranges(cov, c("foreground", "ratio", "tum.counts", "ncn",
                       if (info$mode == "paired") "norm.counts"))
    write.text("seg")
    write.ranges(seg, c("ncn", "seg.mean", "num.mark"))
    write.text("nseg")
    write.ranges(nseg, c("cn", "mean", "sd"))
    close(connection)
    connection <- NULL
    paste0("cbs-v1-", unname(tools::md5sum(file)))
}

finalize_cbs_bundle <- function(cov, seg, nseg, mode, input_field, foreground_scale) {
    info <- list(version = 1L, id = NULL, mode = mode, field = "foreground",
                 input_field = input_field, foreground_scale = foreground_scale)
    info$id <- cbs_bundle_id(cov, seg, nseg, info)
    S4Vectors::metadata(cov)$cbs_coverage <- info
    S4Vectors::metadata(seg)$cbs_coverage <- info
    S4Vectors::metadata(nseg)$cbs_coverage <- info
    validate_cbs_bundle(cov, seg, nseg, mode)
    list(cov = cov, seg = seg, nseg = nseg)
}

validate_cbs_bundle <- function(cov, seg, nseg = NULL, mode = NULL) {
    if (!methods::is(cov, "GRanges")) stop("CBS coverage must be a GRanges object")
    if (!methods::is(seg, "GRanges")) stop("CBS seg must be a GRanges object")
    if (!is.null(nseg) && !methods::is(nseg, "GRanges"))
        stop("CBS nseg must be a GRanges object")
    if (!length(seg)) stop("CBS tumor segmentation must not be empty")
    tagged <- c(!is.null(S4Vectors::metadata(cov)$cbs_coverage),
                !is.null(S4Vectors::metadata(seg)$cbs_coverage),
                !is.null(nseg) && !is.null(S4Vectors::metadata(nseg)$cbs_coverage))
    if (!any(tagged)) {
        # Historical coverage is already the requested foreground signal. Do not
        # normalize it or require fields introduced with the tagged CBS model.
        cbs_numeric_field(cov, "foreground", "Legacy CBS coverage")
        if (identical(mode, "paired") && !length(nseg))
            stop("Paired CBS coverage requires normal segmentation")
        return(invisible(NULL))
    }
    if (!all(tagged)) stop("Mixed CBS tagged and untagged coverage/segmentation provenance")
    info <- cbs_coverage_info(cov)
    if (!is.null(mode) && !identical(mode, info$mode))
        stop("CBS coverage mode does not match the requested sample mode")
    for (name in c("seg", "nseg")) {
        object <- if (name == "seg") seg else nseg
        if (!methods::is(object, "GRanges")) stop("CBS ", name, " must be a GRanges object")
        other <- S4Vectors::metadata(object)$cbs_coverage
        cbs_validate_info(other, paste("CBS", name))
        if (!identical(info, other)) stop("Mixed CBS coverage/segmentation provenance: ", name)
    }
    segment.ncn <- cbs_numeric_field(seg, "ncn", "CBS tumor segmentation")
    if (any(!is.finite(segment.ncn) | segment.ncn <= 0))
        stop("CBS tumor segmentation is missing a finite positive normal baseline")
    if (info$mode == "paired") {
        if (!length(nseg)) stop("Paired CBS coverage requires normal segmentation")
        cbs_numeric_field(nseg, "cn", "CBS normal segmentation")
    } else if (length(nseg)) {
        stop("Tumor-only CBS coverage must have an empty normal segmentation")
    }
    if (!identical(info$id, cbs_bundle_id(cov, seg, nseg, info)))
        stop("CBS bundle contents do not match their provenance ID")
    invisible(info)
}

validate_cbs_graph <- function(graph, cov) {
    if (!methods::is(cov, "GRanges")) stop("CBS coverage must be a GRanges object")
    info <- S4Vectors::metadata(cov)$cbs_coverage
    if (is.null(info))
        cbs_numeric_field(cov, "foreground", "Legacy CBS coverage")
    else
        info <- cbs_coverage_info(cov)
    if (!inherits(graph, "gGraph") && is.list(graph) && !is.null(graph$gg)) {
        outer.info <- graph$meta$cbs_coverage
        if (!is.null(outer.info) && !identical(info, outer.info))
            stop("JaBbA object and coverage have mixed CBS provenance")
        graph <- graph$gg
    }
    graph.info <- tryCatch(graph$meta$cbs_coverage, error = function(e) NULL)
    if (is.null(info) && is.null(graph.info)) return(invisible(NULL))
    if (is.null(info) || is.null(graph.info))
        stop("Graph and coverage have mixed CBS tagged and untagged provenance")
    cbs_validate_info(graph.info, "Graph")
    if (!identical(info, graph.info)) stop("Graph and coverage have mixed CBS provenance")
    invisible(info)
}
