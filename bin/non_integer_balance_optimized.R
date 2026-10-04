#!/usr/bin/env Rscript

# Match jabba_optimized.R: patch a private copy, not the installed namespace.
# Adaptive M is a heuristic CN cap, not a proof of formulation equivalence.
suppressPackageStartupMessages(library(gGnome))

patch_non_integer_balance <- function(fun) {
    source <- paste(deparse(fun, width.cutoff = 500L), collapse = "\n")
    copy_anchor <- "gg = gg$copy"
    adaptive_m <- paste(
        copy_anchor,
        "observed.cn = suppressWarnings(max(abs(gg$nodes$dt$cn), na.rm = TRUE))",
        "if (!is.finite(observed.cn)) observed.cn = 0",
        "M = max(50, ceiling(2 * observed.cn))",
        'if (verbose) message("Using adaptive big-M = ", M, " (max observed CN = ", signif(observed.cn, 6), ")")',
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
    # Insert before nonintegral handling: retain its negative fractional lower
    # bounds on REF edges/terminal slack and its chromosome slush variables.
    source <- sub(bounds_anchor, tight_bounds, source, fixed = TRUE)

    solver_anchor <- "sol = Rcplex2("
    if (!grepl(solver_anchor, source, fixed = TRUE)) {
        stop("Unsupported gGnome::balance(): CPLEX call anchor not found")
    }
    source <- sub(solver_anchor, paste0("control$mipemphasis = as.integer(mipemphasis)\n        ", solver_anchor), source, fixed = TRUE)
    # An R-level parameter; no native package rebuild is needed.
    fun <- eval(parse(text = source), envir = environment(fun))
    formals(fun)$mipemphasis <- 0L
    fun
}

configure_non_integer_cplex <- function(threads, mipemphasis, outdir) {
    if (length(threads) != 1L || !is.finite(threads) || threads < 1 || threads != floor(threads)) {
        stop("--threads must be a positive integer")
    }
    if (length(mipemphasis) != 1L || !is.finite(mipemphasis) || !(mipemphasis %in% 0:4)) {
        stop("--mipemphasis must be an integer from 0 through 4")
    }
    # Rcplex2's control list does not support threads. CPLEX reads this file at
    # CPXopenCPLEX. Preserve other settings supplied by an existing parameter file.
    inherited <- Sys.getenv("ILOG_CPLEX_PARAMETER_FILE")
    parameters <- if (nzchar(inherited)) readLines(inherited, warn = FALSE) else {
        # Runtime verified in mskilab/jabba:0.0.8; do not infer from the cplex CLI.
        "CPLEX Parameter File Version 22.1.2.0"
    }
    parameters <- parameters[!grepl("^\\s*(CPXPARAM_Threads|CPX_PARAM_THREADS|CPXPARAM_Emphasis_MIP|CPX_PARAM_MIPEMPHASIS)\\s", parameters)]
    param.file <- file.path(normalizePath(outdir, mustWork = TRUE), "non_integer.cplex.prm")
    writeLines(c(parameters,
                 paste("CPXPARAM_Threads", as.integer(threads), sep = "\t")), param.file)
    Sys.setenv(ILOG_CPLEX_PARAMETER_FILE = param.file)
    message("Using CPLEX threads = ", threads, ", mipemphasis = ", mipemphasis, " via ", param.file)
}

patched_balance <- patch_non_integer_balance(get("balance", asNamespace("gGnome")))
run <- new.env(parent = globalenv())
run$balance <- function(...) {
    opt <- run$opt
    if (!opt$gurobi) {
        configure_non_integer_cplex(opt$threads, opt$mipemphasis, opt$outdir)
    }
    patched_balance(..., mipemphasis = opt$mipemphasis)
}

# Both files are staged as Nextflow path inputs, so changes invalidate task cache.
script.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script.dir <- dirname(sub("^--file=", "", script.arg))
source(file.path(script.dir, "non_integer_balance.R"), local = run, chdir = FALSE)
