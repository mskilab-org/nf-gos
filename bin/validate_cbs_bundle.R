#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L || !args[[1L]] %in% c("paired", "tumor_only")) {
    stop("Usage: validate_cbs_bundle.R paired|tumor_only cov.rds seg.rds nseg.rds [graph_name=graph.rds ...]",
         call. = FALSE)
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(sub("^--file=", "", script_arg[[1L]]))
source(file.path(script_dir, "cbs_coverage.R"))

tryCatch({
    cov <- readRDS(args[[2L]])
    validate_cbs_bundle(cov, readRDS(args[[3L]]), readRDS(args[[4L]]), mode = args[[1L]])
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
    stop("CBS coverage migration: ", conditionMessage(e),
         ". Supply cov/seg/nseg and graph outputs from the same CBS bundle, or remove the legacy imports and regenerate CBS and all graph descendants.",
         call. = FALSE)
})
