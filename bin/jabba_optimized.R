#!/usr/bin/env Rscript

# JaBbA 1.1 delegates its LP/MILP construction to gGnome::balance(). The
# packaged function leaves copied binary rows with inherited infinite bounds
# and uses a fixed M=1000 for every copy-number and indicator constraint. Patch
# the imported function before running JaBbA so CPLEX receives an equivalent,
# tighter formulation.
suppressPackageStartupMessages(library(JaBbA))

script.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script.dir <- dirname(sub("^--file=", "", script.arg))
source(file.path(script.dir, "cbs_coverage.R"), chdir = FALSE)

patch_balance <- function(fun) {
    source <- paste(deparse(fun, width.cutoff = 500L), collapse = "\n")

    copy_anchor <- "gg = gg$copy"
    adaptive_m <- paste(
        copy_anchor,
        "observed.cn = suppressWarnings(max(abs(gg$nodes$dt$cn), na.rm = TRUE))",
        "if (!is.finite(observed.cn)) observed.cn = 0",
        "M = max(50, ceiling(2 * observed.cn))",
        "if (verbose) message(\"Using adaptive big-M = \", M, \" (max observed CN = \", signif(observed.cn, 6), \")\")",
        sep = "\n    "
    )
    if (!grepl(copy_anchor, source, fixed = TRUE)) {
        stop("Unsupported gGnome::balance(): graph-copy anchor not found")
    }
    source <- sub(copy_anchor, adaptive_m, source, fixed = TRUE)

    bounds_anchor <- 'vars[type %in% c("loose.in", "loose.out"), `:=`(lb = 0, ub = Inf)]'
    tight_bounds <- paste(
        bounds_anchor,
        'vars[type %in% c("loose.in", "loose.out"), `:=`(ub = M)]',
        'vars[vtype == "B", `:=`(lb = 0, ub = 1)]',
        sep = "\n    "
    )
    if (!grepl(bounds_anchor, source, fixed = TRUE)) {
        stop("Unsupported gGnome::balance(): variable-bound anchor not found")
    }
    source <- sub(bounds_anchor, tight_bounds, source, fixed = TRUE)

    eval(parse(text = source), envir = environment(fun))
}

patched_balance <- patch_balance(get("balance", asNamespace("gGnome")))
jabba_imports <- parent.env(asNamespace("JaBbA"))
unlockBinding("balance", jabba_imports)
assign("balance", patched_balance, envir = jabba_imports)
lockBinding("balance", jabba_imports)

jba <- system.file("extdata", "jba", package = "JaBbA")
if (!nzchar(jba)) stop("JaBbA command script not found")
# Keep the installed command's option parsing and solver call unchanged. Its
# JaBbA call is intercepted only to enforce the CBS input contract before fit.
coverage.info <- NULL
JaBbA <- function(junctions, coverage, seg = NULL, nseg = NULL, field, ...) {
    if (!identical(field, "foreground")) {
        stop("JaBbA requires --field foreground from the CBS coverage bundle")
    }
    cov <- readRDS(coverage)
    cbs.seg <- if (is.null(seg)) NULL else readRDS(seg)
    cbs.nseg <- if (is.null(nseg)) NULL else readRDS(nseg)
    validate_cbs_bundle(cov, cbs.seg, cbs.nseg)
    coverage.info <<- cbs_coverage_info(cov)
    # An empty tumor-only nseg is a valid bundle member, but JaBbA treats a
    # supplied empty normal segmentation as an all-missing normal baseline.
    if (coverage.info$mode == "tumor_only" && length(cbs.nseg) == 0L) nseg <- NULL
    JaBbA::JaBbA(junctions = junctions, coverage = coverage, seg = seg,
                 nseg = nseg, field = "foreground", ...)
}
# JaBbA itself reads global opt$epgap after fitting, so retain the installed
# command's original global option scope while intercepting its function call.
source(jba, chdir = FALSE)

# The command returns after persisting both legacy JaBbA lists and gGraphs.
# Stamp those files, not just its returned graph, so imported outputs retain
# the exact source bundle identity. Solver intermediates are left untouched.
if (is.null(coverage.info)) stop("JaBbA did not consume a validated CBS bundle")
required.outputs <- file.path(opt$outdir, c("jabba.simple.rds", "jabba.simple.gg.rds"))
if (!all(file.exists(required.outputs))) stop("JaBbA final graph outputs are missing")
output.files <- file.path(opt$outdir, c("jabba.simple.rds", "jabba.simple.gg.rds",
                                          "jabba.raw.rds", "jabba.rds", "jabba.gg.rds",
                                          "karyograph.rds"))
for (output.file in output.files[file.exists(output.files)]) {
    graph <- readRDS(output.file)
    if (inherits(graph, "gGraph")) {
        graph$set(cbs_coverage = coverage.info)
    } else if (is.list(graph)) {
        graph$meta$cbs_coverage <- coverage.info
        if (inherits(graph$gg, "gGraph")) graph$gg$set(cbs_coverage = coverage.info)
    } else {
        stop("Unsupported JaBbA graph output: ", output.file)
    }
    saveRDS(graph, output.file)
}
