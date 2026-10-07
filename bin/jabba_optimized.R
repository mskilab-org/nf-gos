#!/usr/bin/env Rscript

# JaBbA 1.1 delegates its LP/MILP construction to gGnome::balance(). The
# packaged function leaves copied binary rows with inherited infinite bounds
# and uses a fixed M=1000 for every copy-number and indicator constraint. Patch
# the imported function before running JaBbA so CPLEX receives an equivalent,
# tighter formulation.
suppressPackageStartupMessages(library(JaBbA))

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
source(jba, chdir = FALSE)
