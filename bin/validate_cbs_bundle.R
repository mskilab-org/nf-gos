#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L || !args[[1L]] %in% c("paired", "tumor_only")) {
    stop("Usage: validate_cbs_bundle.R paired|tumor_only cov.rds seg.rds|NULL nseg.rds|NULL [graph_name=graph.rds ...]",
         call. = FALSE)
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(sub("^--file=", "", script_arg[[1L]]))
source(file.path(script_dir, "cbs_coverage.R"))

tryCatch({
    cov <- readRDS(args[[2L]])
    seg <- if (args[[3L]] == "NULL") NULL else readRDS(args[[3L]])
    nseg <- if (args[[4L]] == "NULL") NULL else readRDS(args[[4L]])
    cbs_numeric_field(cov, "foreground", "CBS coverage")
    if (is.null(seg) && is.null(S4Vectors::metadata(cov)$cbs_coverage)) {
        # NI can resume an untagged graph from raw coverage without CBS segments.
        if (!is.null(nseg)) {
            if (!methods::is(nseg, "GRanges")) stop("CBS nseg must be a GRanges object")
            if (!is.null(S4Vectors::metadata(nseg)$cbs_coverage))
                stop("Mixed CBS tagged and untagged coverage/segmentation provenance")
        }
    } else {
        validate_cbs_bundle(cov, seg, nseg, mode = args[[1L]])
    }
    if (length(args) > 4L) {
        suppressPackageStartupMessages(library(gGnome))
        for (argument in args[-seq_len(4L)]) {
            separator <- regexpr("=", argument, fixed = TRUE)[[1L]]
            if (separator < 2L) stop("Graph arguments must be graph_name=graph.rds", call. = FALSE)
            name <- substr(argument, 1L, separator - 1L)
            path <- substring(argument, separator + 1L)
            tryCatch(validate_cbs_graph(readRDS(path), cov), error = function(e) {
                stop(name, ": ", conditionMessage(e), call. = FALSE)
            })
        }
    }
}, error = function(e) {
    stop("CBS input validation: ", conditionMessage(e),
         ". Use wholly untagged legacy inputs, or cov/seg/nseg and graph outputs from the same tagged CBS bundle.",
         call. = FALSE)
})
