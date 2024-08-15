{
    library(optparse)

    message("Loading input arguments")
    if (!exists('opt'))
    {
        option_list = list(
            make_option(c("--complex"), type = "character", help = "output of Events.task - path to .rds gGraph object annotated with complex SV event classes"),
            make_option(c("--homeology"), type = "character", help = "output of homeology_hrd.task - path to homeology results"),
            make_option(c("--homeology_stats"), type = "character", help = "output of homeology_hrd.task - path to homeology results, table of each sequence heatmap component for every junction"),
            make_option(c("--hrdetect_results"), type = "character",
                        help = "output of hrdetect_og.task - preferably run with Strelka2 snv and SvAbA indels"),
            make_option(c("--model"), type = "character", help = "cached path to trained model"),
            make_option(c("--outdir"), type = "character", default = './',
                        help = "output directory"),
            make_option(c("--libdir"), type = "character",
                        help = "Directory containing this R file")
        )
        message("Found input files")

        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)
        message("Parsed input arguments")


        if (is.null(opt$complex) | (is.null(opt$homeology)) | is.null(opt$homeology_stats) |
            is.null(opt$hrdetect_results) | is.null(opt$model))  # was opt$multinomial
            stop(print_help(parseobj))

        print(opt)

        print(.libPaths())
        options(error={expression({traceback(2); quit('no', 1)})})

        message("Establishing record files")
        ## keep record of run
        writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }

    message("Loading Packages")
    suppressWarnings(expr = {
        suppressMessages(expr = {
            suppressPackageStartupMessages(expr = {
                library(gGnome)
                library(dplyr)
                library(khtools)
                library(glmnet)
            })
        })
    })

    message("Done Loading Packages")

    setDTthreads(1)

    system(paste('mkdir -p',  opt$outdir))

##########
##########
########## Oneness Twoness
##########
##########

    ## .NotYetImplemented()
    message("Processing QRP")
    gg = readRDS(opt$complex)
    ## gg = readRDS('/gpfs/commons/groups/imielinski_lab/projects/Starr/BRCA/Flow/testing/events/20210810_qrp/Events/0751/complex.rds')
    gg = gGnome::refresh(gg)
    all.events = gg$meta$events
    ev.types = c("qrppos", "qrpmin", "qrpmix", "qrdup", "qrdel", "tib") ## including the old naming scheme for backwards compatibility
    all.events$type = factor(all.events$type, ev.types)
    ## expl_variables = all.events %>% dcast.data.table("sample" ~ type, fun.aggregate = length, drop = FALSE)

    expl_variables = all.events %>% reshape2::dcast("sample" ~ type, fun.aggregate = length, drop = FALSE)

    if (NROW(all.events) > 0) {
        .new = c(expl_variables$qrppos, expl_variables$qrpmin, expl_variables$qrpmix)
        .old = c(expl_variables$qrdup, expl_variables$qrdel, expl_variables$tib)
        ## .new = expl_variables[, c(qrppos, qrpmin, qrpmix)]
        ## .old = expl_variables[, c(qrdup, qrdel, tib)]

        if (!identical(.new, .old)) {
            if (all(.new == 0) && any(.old) > 0) {
                expl_variables$qrppos = expl_variables$qrdup
                expl_variables$qrpmin = expl_variables$qrdel
                expl_variables$qrpmix = expl_variables$tib
            }  else if (all(.old == 0) && any(.new) > 0) {
                expl_variables$qrdup = expl_variables$qrppos
                expl_variables$qrdel = expl_variables$qrpmin
                expl_variables$tib = expl_variables$qrpmix
            }
            ## expl_variables[, !names(expl_variables) %in% c(""), with = FALSE]
        }
    } else {
        expl_variables$qrppos = 0
        expl_variables$qrpmin = 0
        expl_variables$qrpmix = 0
        expl_variables$qrdup = 0
        expl_variables$qrdel = 0
        expl_variables$tib = 0
    }



    ## .NotYetImplemented()
    message("Processing homeologous dels")
    jhom = fread(opt$homeology)
    ## jhom = fread('/gpfs/commons/groups/imielinski_lab/projects/Starr/BRCA/Flow/homeology/20210408/Homeology/0751/junctions.txt')
    if (NROW(jhom) > 0) {
        bp1 = parse.gr(jhom$bp1)
        bp2 = parse.gr(jhom$bp2)
        bp1 = gr.fix(bp1, bp2)
        bp2 = gr.fix(bp2, bp1)
        jhom$jspan = jJ(grl.pivot(GRangesList(bp1, bp2)))$span
    }
    jhom_stats = fread(opt$homeology_stats)
    ## jhom_stats = fread("/gpfs/commons/groups/imielinski_lab/projects/Starr/BRCA/Flow/homeology/20210408/Homeology/0751/rawstats.txt")
    dels = jhom[!is.na(jhom$del),colnames(jhom),drop=F,with=F]
    if (NROW(jhom_stats)) {
        dels = merge.repl(
            dels,
            jhom_stats[,.(
                hlen = max(max(ifelse(na2false(.SD$r > 0.9), .SD$minpx, 0L)), 0L)
            ), by = edge.id],
            by = "edge.id")
    }
    num_ihdels = NROW(dels[dels$hlen >= 10 & dels$jspan > 1000,])
    expl_variables$ihdel = num_ihdels


    message("Processing HRDetect inputs: del.mh.prop, RS3, RS5, hrd-LOH score, SNV3, SNV8")
    res = readRDS(opt$hrdetect_results)
    ## res = readRDS(mpairs['0751']$hrd_results)
    hrd = res$data_matrix
    if (!identical(type(hrd), "double")) hrd = data.matrix(hrd)
    hrd = as.data.table(hrd)

    hrd = setcols(hrd, c("SV3", "SV5"), c("RS3", "RS5"))

    expl_variables = expl_variables %>% mutate(hrd)
    ## c("tib","ihdel","qrdup","qrdel","RS3","RS5","del.mh.prop","hrd","SNV8","SNV3")
    ## opt$model = "multinomial_glmnet_model.rds"

    if (!grepl("\\/", opt$model))
        mod_path = paste0(opt$libdir, '/', opt$model)
    else
        mod_path = opt$model
    message("Path to model:\n", mod_path)
    mod = readRDS(mod_path)

    tmp = expl_variables[, c("tib","ihdel","qrdup","qrdel","RS3","RS5","del.mh.prop","hrd","SNV8","SNV3"),drop=F]

    ### begin glmnet_utils.R
        ### begin .Rprofile
    options(stringsAsFactors = FALSE)
    options(bitmapType="cairo")
    options(device = grDevices::png)
    options(scipen = 0)

    isTRUE = function(x) {
        identical(x, TRUE)
    }

    isFALSE = function(x) {
        identical(x, FALSE)
    }

    tplot = function(...) {
        this.mar = par()$mar
        this.mai = par()$mai
        this.oma = par()$oma
        on.exit({par(mar = this.mar, mai = this.mai)})
        par(mar = c(0,0,0,0), mai = c(0,0,0,0))
        plot(c(0.25, 0.75), c(0.25, 0.75), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste(...),
              cex = 1.6, col = "black")
        ## text(x = 0.5, y = 0.5, paste("The following is text that'll appear in a plot window.\n",
        ##                          "As you can see, it's in the plot window\n",
        ##                          "One might imagine useful informaiton here"),
        ##       cex = 1.6, col = "black")
    }


    quiet = function(this_expr, do_global = TRUE, set = FALSE) {
        pf = parent.frame()
        fout = file(nullfile(), open = "wt")
        fout2 = file(nullfile(), open = "wt")
        if (!isTRUE(set)) {
            on.exit({
                i <- sink.number(type = "message")
                if (i > 0L)
                    sink(stderr(), type = "message")
                n <- sink.number()
                if (n > 0L)
                    for (i in seq_len(n)) sink()
                gc()
                invisible()
            })
        }
        suppressMessages({
            suppressWarnings({
                suppressPackageStartupMessages({
                    sink(fout, type = "output")
                    sink(fout2, type = "message");
                    tryCatch({
                        if (do_global)
                            eval(this_expr, envir = globalenv())
                        else
                            eval(this_expr, envir = pf)
                    }, error = function(e) invisible())
                    tryCatch(close(fout), error = function(e) invisible())
                    tryCatch(close(fout2), error = function(e) invisible())
                    tryCatch(sink(), error = function(e) invisible())
                    tryCatch(sink(), error = function(e) invisible())
                    invisible()
                })
            })
        })
    }

    unquiet = function() {
        i <- sink.number(type = "message")
        if (i > 0L)
            sink(stderr(), type = "message")
        n <- sink.number()
        if (n > 0L)
            for (i in seq_len(n)) sink()
        gc()
        invisible()
    }

    close_all_connections = function() {
        i <- sink.number(type = "message")
        if (i > 0L)
            sink(stderr(), type = "message")
        n <- sink.number()
        if (n > 0L)
            for (i in seq_len(n)) sink()
        gc()
        ## only connections that are not jupyter related
        cons = as.integer(rownames(subset(
            as.data.frame(showConnections(all = FALSE)),
            !description == "output" & !class == "textConnection" &
            !mode == "wr" & !text == "text" &
            !isopen == "opened" & !`can read` == "no" &
            !`can write` == "yes"
        )))
        # set <- getAllConnections()
        # set <- set[set > 2L]
        for (i in seq_along(cons)) close(getConnection(cons[i]))
        invisible()
    }

    Sys.setenv("PLOTCACHE" = "~/Dropbox/plotcache")
    replotly = function(plotly_obj, plotcache_dir = NULL, verbose = TRUE) {
        require(plotly)
        if (is.null(plotcache_dir)) {
            plotcache_dir = Sys.getenv("PLOTCACHE")
        }
        if (!file.exists(plotcache_dir)) {
            plotcache_dir = "~/Dropbox/plotcache"
        }
        if (!file.exists(plotcache_dir)) {
            plotcache_dir = "~/"
        }
        set.seed(NULL) # prevent overwrites
        if (is.character(plotly_obj) && grepl(".rds$", plotly_obj)) {
            if (file.exists(plotly_obj)) {
                outrds = plotly_obj
                outhtml = sub(".rds$", ".html", outrds)
                print(list(rds = outrds, html = outhtml))
                plotly_obj = readRDS(plotly_obj)
                return(plotly_obj)
            } else {
                NULL
            }
        }
        if (verbose) {
            message("plotting to ", plotcache_dir)
        }
        outfile = tempfile(pattern = "plot", tmpdir = plotcache_dir)
        outrds = paste0(outfile, ".rds")
        outhtml = paste0(outfile, ".html")
        if (!inherits(plotly_obj, c("plotly", "htmlwidget"))) {
            plotly_obj = plotly_build(plotly_obj)
        }
        saveRDS(plotly_obj, outrds, version = 2)
        htmlwidgets::saveWidget(ggpl, outhtml, selfcontained = TRUE)
        return(list(rds = outrds, html = outhtml))
    }



    names2 = function(x) {
        nm = names(x)
        if (is.null(nm))
            return(rep_len("", length(x)))
        else
            return(nm)
    }

    `names2=` = function(x, value, useempty = FALSE) {
        names(x) = if (!is.null(value))
                       rep_len(value, length(x))
                   else {
                       if (useempty)
                           rep_len("", length(x))
                   }
        return(x)
    }

    forceload = function(envir = globalenv(), .force = FALSE) {
        if (!exists("envload_34507213048974")) {
            envload_34507213048974 = new.env(parent = globalenv())
            globasn(envload_34507213048974)
            sesh =  sessionInfo()
            pkgs = c(sesh$basePkgs,
                     names(sesh$otherPkgs),
                     names(sesh$loadedOnly))
            pkvec = rep(FALSE, length(pkgs))
            names(pkvec) = pkgs
            envload_34507213048974$pkvec = pkvec
        }
        force = function(x) x
        ## pkgs = gsub("package:", "", grep('package:', search(), value = TRUE))
        ## pkgs = c(pkgs, names(sessionInfo()$loadedOnly))
        if (!exists("sesh")) sesh =  sessionInfo()
        if (!exists("pkgs")) {
            pkgs = c(sesh$basePkgs,
                     names(sesh$otherPkgs),
                     names(sesh$loadedOnly))
        }
        pkvec = envload_34507213048974$pkvec
        notloaded_firsttime = setdiff(pkgs, names(pkvec))
        pkvec = c(pkvec, setNames(rep_len(FALSE, length(notloaded_firsttime)), notloaded_firsttime))
        if (.force)
            notloaded = names(pkvec)
        else
            notloaded = names(which(pkvec == FALSE))
        if (length(notloaded)) {
            for (pkg in notloaded) {
                tryCatch( {
                    message("force loading ", pkg)
                    ## invisible(eval(as.list((asNamespace(pkg))), envir = envir))
                    ## invisible(eval(eapply(asNamespace(pkg), force, all.names = TRUE), envir = envir))
                    invisible(eval(parse(text = sprintf("as.list((asNamespace(\"%s\")))", pkg)), envir = envir))
                    invisible(eval(parse(text = sprintf("eapply(asNamespace(\"%s\"), force, all.names = TRUE)", pkg)), envir = envir))
                    pkvec[pkg] = TRUE
                }, error = function(e) message("could not force load ", pkg))
            }
            envload_34507213048974$pkvec = pkvec
        } else {
            message("nothing to forceload")
        }
        invisible()
    }

    forcefun = function(envir = globalenv(), evalenvir = globalenv()) {
        funnames = as.character(lsf.str(envir = envir))
        for (fun in funnames) {
            tryCatch( {
                message("force loading ", fun)
                eval(force(get(fun, envir = envir)), envir = evalenvir)
            }, error = function(e) message("could not force load ", fun))
        }
        invisible()
    }



    relib2 = function(lib = 'Flow', force = TRUE, unload = TRUE)
    {
        suppressMessages(forceload(.force = T))
        if (sprintf("package:%s", lib) %in% search())
        {
            expr = sprintf("detach(package:%s, force = force, unload = unload)", lib)
            eval(parse(text = expr))
            ## tryCatch(unload(lib), error = function(e) NULL) ## DO NOT use this line...
            ## it will break re-librarying
        }
        txt = sprintf("library2(%s)", lib)
        eval(parse(text = txt))
        suppressMessages(forceload(.force = T))
        invisible()
    }

    relib3 = function(..., force = TRUE, unload = TRUE)
    {
        if (!exists("envload_34507213048974")) {
            envload_34507213048974 = new.env(parent = globalenv())
            globasn(envload_34507213048974)
            sesh =  sessionInfo()
            pkgs = c(sesh$basePkgs,
                     names(sesh$otherPkgs),
                     names(sesh$loadedOnly))
            pkvec = rep(FALSE, length(pkgs))
            names(pkvec) = pkgs
            envload_34507213048974$pkvec = pkvec
        }
        suppressMessages(forceload(.force = T))
        names2 = function(x) {
            nm = names(x)
            if (is.null(nm))
                return(rep_len("", length(x)))
            else
                return(nm)
        }
        lst.arg = as.list(match.call(expand.dots = F))$`...`
        nm = names2(lst.arg)
        otherarg = lst.arg[nzchar(nm)]
        pkgarg = lst.arg[!nzchar(nm)]
        pkgarg = pkgarg[sapply(pkgarg, function(x) is.call(x) || is.character(x))]
        charvec = as.character(all.vars(match.call()))
        if (length(charvec)) {
            notfound= { set.seed(10); paste0("notfound_", round(runif(1) * 1e9)); }
            vars = mget(charvec, ifnotfound=notfound, mode = "character", inherits = T)
            ## charvec = unlist(strsplit(toString(vars[[1]]), ", "))
            charvec = unique(c(names2(vars[vars == notfound]), unlist(vars[vars != notfound])))
        }
        charvec = c(charvec, unlist(as.vector(sapply(pkgarg,
                                                     function(x) tryCatch(eval(x), error = function(e) NULL)))))
        for (lib in charvec) {
            if (sprintf("package:%s", lib) %in% search())
            {
                expr = sprintf("detach(package:%s, force = force, unload = unload)", lib)
                eval(parse(text = expr))
                ## tryCatch(unload(lib), error = function(e) NULL) ## DO NOT use this line...
                ## it will break re-librarying
            }
            pkvec = envload_34507213048974$pkvec
            if (lib %in% names(pkvec)) {
                pkvec = pkvec[!names(pkvec) %in% lib]
                envload_34507213048974$pkvec = pkvec
            } else {
                pev = packageEvent(lib, "onLoad")
                gh = getHook(pev)
                if (length(gh) == 0 || is.null(gh$forceall12340987)) {
                    setHook(pev,
                            list("forceall12340987" = function(...) forceall(envir = asNamespace(lib))))
                }
            }
            expr = parse(text = sprintf("library(%s)", lib))
            eval(expr, globalenv())
            ## library(lib, character.only = T)
        }
        suppressMessages(forceload(.force = T))
        invisible()
    }

    rereq3 = function(..., force = TRUE, unload = TRUE)
    {
        if (!exists("envload_34507213048974")) {
            envload_34507213048974 = new.env(parent = globalenv())
            globasn(envload_34507213048974)
            sesh =  sessionInfo()
            pkgs = c(sesh$basePkgs,
                     names(sesh$otherPkgs),
                     names(sesh$loadedOnly))
            pkvec = rep(FALSE, length(pkgs))
            names(pkvec) = pkgs
            envload_34507213048974$pkvec = pkvec
        }
        suppressMessages(forceload(.force = T))
        names2 = function(x) {
            nm = names(x)
            if (is.null(nm))
                return(rep_len("", length(x)))
            else
                return(nm)
        }
        lst.arg = as.list(match.call(expand.dots = F))$`...`
        nm = names2(lst.arg)
        otherarg = lst.arg[nzchar(nm)]
        pkgarg = lst.arg[!nzchar(nm)]
        pkgarg = pkgarg[sapply(pkgarg, function(x) is.call(x) || is.character(x))]
        charvec = as.character(all.vars(match.call()))
        if (length(charvec)) {
            notfound= { set.seed(10); paste0("notfound_", round(runif(1) * 1e9)); }
            vars = mget(charvec, ifnotfound=notfound, mode = "character", inherits = T)
            ## charvec = unlist(strsplit(toString(vars[[1]]), ", "))
            charvec = unique(c(names2(vars[vars == notfound]), unlist(vars[vars != notfound])))
        }
        charvec = c(charvec, unlist(as.vector(sapply(pkgarg,
                                                     function(x) tryCatch(eval(x), error = function(e) NULL)))))
        for (lib in charvec) {
            if (sprintf("package:%s", lib) %in% search())
            {
                expr = sprintf("detach(package:%s, force = force, unload = unload)", lib)
                eval(parse(text = expr))
                ## tryCatch(unload(lib), error = function(e) NULL) ## DO NOT use this line...
                ## it will break re-librarying
            }
            pkvec = envload_34507213048974$pkvec
            if (lib %in% names(pkvec)) {
                pkvec = pkvec[!names(pkvec) %in% lib]
                envload_34507213048974$pkvec = pkvec
            } else {
                pev = packageEvent(lib, "onLoad")
                gh = getHook(pev)
                if (length(gh) == 0 || is.null(gh$forceall12340987)) {
                    setHook(pev,
                            list("forceall12340987" = function(...) forceall(envir = asNamespace(lib))))
                }
                ## do.call(require, c(alist(package = lib, character.only = T),
                ##                    otherarg))
                if (NROW(otherarg)) {
                    is.char = sapply(otherarg, is.character)
                    otherarg[is.char] = paste0("\"", otherarg[is.char], "\"")
                    otherargs = paste(paste(names(otherarg), "=", unlist(otherarg)), collapse = ",")
                    eval(parse(text = sprintf("require(%s,%s)", lib, otherargs)), globalenv())
                } else {
                    eval(parse(text = sprintf("require(%s)", lib)), globalenv())
                }
            }
        }
        suppressMessages(forceload(.force = T))
        invisible()
    }

    no.dev = function() {
        evalq({
            for (d in dev.list()) {
                dev.off(d)
            }
        }, envir = globalenv())
        invisible()
    }

    detach2 = function(lib = "Flow", force = TRUE, unload = TRUE) {
        suppressMessages(forceload(.force = T))
        if (sprintf("package:%s", lib) %in% search())
        {
            expr = sprintf("detach(package:%s, force = force, unload = unload)", lib)
            suppressMessages(eval(parse(text = expr)))
            tryCatch(unload(lib), error = function(e) NULL)
        }
        suppressMessages(forceload(.force = T))
        invisible()
    }

    library2 = function(x, ...) {
        suppressMessages(forceload(.force = T))
        arg = as.list(match.call())[["x"]]
        if (is.symbol(arg)) {
            lib = tryCatch(as.character(eval(arg)), error = function(e) arg)
            if (!is.character(lib)) {
                lib = toString(lib)
            }
        } else {
            lib = x
        }
        library(lib, character.only = T, ...)
        suppressMessages(forceload(.force = T))
        invisible()
    }

    library3 = function (...)
    {
        names2 = function(x) {
            nm = names(x)
            if (is.null(nm))
                return(rep_len("", length(x)))
            else
                return(nm)
        }
        suppressMessages(forceload(.force = T))
        lst.arg = as.list(match.call(expand.dots = F))$`...`
        nm = names2(lst.arg)
        otherarg = lst.arg[nzchar(nm)]
        pkgarg = lst.arg[!nzchar(nm)]
        pkgarg = pkgarg[sapply(pkgarg, function(x) is.call(x) || is.character(x))]
        charvec = as.character(all.vars(match.call()))
        if (length(charvec)) {
            notfound= { set.seed(10); paste0("notfound_", round(runif(1) * 1e9)); }
            vars = mget(charvec, ifnotfound=notfound, mode = "character", inherits = T)
            ## charvec = unlist(strsplit(toString(vars[[1]]), ", "))
            charvec = unique(c(names2(vars[vars == notfound]), unlist(vars[vars != notfound])))
        }
        charvec = c(charvec, unlist(as.vector(sapply(pkgarg,
                                                     function(x) tryCatch(eval(x), error = function(e) NULL)))))
        for (lib in charvec) {
            pev = packageEvent(lib, "onLoad")
            gh = getHook(pev)
            if (length(gh) == 0 || is.null(gh$forceall12340987)) {
                setHook(pev,
                        list("forceall12340987" = function(...) forceall(envir = asNamespace(lib))))
            }
            if (NROW(otherarg)) {
                is.char = sapply(otherarg, is.character)
                otherarg[is.char] = paste0("\"", otherarg[is.char], "\"")
                otherargs = paste(paste(names(otherarg), "=", unlist(otherarg)), collapse = ",")
                eval(parse(text = sprintf("library(%s,%s)", lib, otherargs)), globalenv())
            } else {
                eval(parse(text = sprintf("library(%s)", lib)), globalenv())
            }
            ## do.call(library, c(alist(package = lib, character.only = T),
            ##                    otherarg))
        }
        suppressMessages(forceload(.force = T))
        invisible()
    }

    require3 = function (...)
    {
        names2 = function(x) {
            nm = names(x)
            if (is.null(nm))
                return(rep_len("", length(x)))
            else
                return(nm)
        }
        suppressMessages(forceload(.force = T))
        lst.arg = as.list(match.call(expand.dots = F))$`...`
        nm = names2(lst.arg)
        otherarg = lst.arg[nzchar(nm)]
        pkgarg = lst.arg[!nzchar(nm)]
        pkgarg = pkgarg[sapply(pkgarg, function(x) is.call(x) || is.character(x))]
        charvec = as.character(all.vars(match.call()))
        if (length(charvec)) {
            notfound= { set.seed(10); paste0("notfound_", round(runif(1) * 1e9)); }
            vars = mget(charvec, ifnotfound=notfound, mode = "character", inherits = T)
            ## charvec = unlist(strsplit(toString(vars[[1]]), ", "))
            charvec = unique(c(names2(vars[vars == notfound]), unlist(vars[vars != notfound])))
        }
        charvec = c(charvec, unlist(as.vector(sapply(pkgarg,
                                                     function(x) tryCatch(eval(x), error = function(e) NULL)))))
        for (lib in charvec) {
            pev = packageEvent(lib, "onLoad")
            gh = getHook(pev)
            if (length(gh) == 0 || is.null(gh$forceall12340987)) {
                setHook(pev,
                        list("forceall12340987" = function(...) forceall(envir = asNamespace(lib))))
            }
            ## do.call(require, c(alist(package = lib, character.only = T),
            ##                    otherarg))
            if (NROW(otherarg)) {
                is.char = sapply(otherarg, is.character)
                otherarg[is.char] = paste0("\"", otherarg[is.char], "\"")
                otherargs = paste(paste(names(otherarg), "=", unlist(otherarg)), collapse = ",")
                eval(parse(text = sprintf("require(%s,%s)", lib, otherargs)), globalenv())
            } else {
                eval(parse(text = sprintf("require(%s)", lib)), globalenv())
            }
        }
        suppressMessages(forceload(.force = T))
        invisible()
    }


    force2 = function(x)
        tryCatch(x, error = function(e) NULL)


    forceall = function(invisible = TRUE, envir = parent.frame(), evalenvir = parent.frame()) {
        if (!exists("envload_34507213048974")) {
            envload_34507213048974 = new.env(parent = globalenv())
            globasn(envload_34507213048974)
            sesh =  sessionInfo()
            pkgs = c(sesh$basePkgs,
                     names(sesh$otherPkgs),
                     names(sesh$loadedOnly))
            pkvec = rep(FALSE, length(pkgs))
            names(pkvec) = pkgs
            envload_34507213048974$pkvec = pkvec
        }
        pkg = environmentName(envir)
        pkvec = envload_34507213048974$pkvec
        if ( { pkg %in% names(pkvec) && isFALSE(pkvec[pkg]); } ||
             { ! pkg %in% names(pkvec); } ) {
            if (invisible == TRUE)  {
                ## invisible(eval(as.list(envir), envir = evalenvir))
                ## invisible(eval(eapply(envir, force, all.names = TRUE), envir = evalenvir))
                invisible(eval(parse(text = sprintf("as.list(asNamespace(\"%s\"))", pkg)), evalenvir))
                invisible(eval(parse(text = sprintf("eapply(asNamespace(\"%s\"), force, all.names = TRUE)", pkg)), envir = evalenvir))
            } else {
                ## print(eval(as.list(envir), envir = evalenvir))
                ## print(eval(eapply(envir, force, all.names = TRUE), envir = evalenvir))
                print(eval(parse(text = sprintf("as.list(asNamespace(\"%s\"))", pkg)), evalenvir))
                print(eval(parse(text = sprintf("eapply(asNamespace(\"%s\"), force, all.names = TRUE)", pkg)), envir = evalenvir))
            }
            addon = TRUE
            names(addon) = pkg
            envload_34507213048974$pkvec = c(pkvec, addon)
        } else {
            message("nothing to load")
        }
    }

    overwriteR6 = function (newfun, oldfun, r6gen, meth = "public_methods", package = NULL,
                            envir = globalenv())
    {
        meth = ifelse(grepl("^pub", meth), "public_methods", ifelse(grepl("^pri",
                                                                          meth), "private_methods", ifelse(grepl("^act", meth),
                                                                                                           "active", NA_character_)))
        if (is.na(meth))
            stop("method must refer to public, private, or active method")
        if (!is.null(package)) {
            if (is.character(package))
                envpkg = asNamespace(package)
            else if (isNamespace(package))
                envpkg = package
            nmpkg = environmentName(envpkg)
        }
        r6 = get(r6gen)
        tmpfun = r6[[meth]][[oldfun]]
        .newfun = get(newfun)
        environment(.newfun) = environment(tmpfun)
        attributes(.newfun) = attributes(tmpfun)
        r6[[meth]][[oldfun]] = .newfun
        invisible()
    }

    globasn = function (obj, var = NULL, return_obj = TRUE, envir = .GlobalEnv,
                         verbose = TRUE, vareval = F)
    {
        var = as.list(match.call())$var
        if (is.null(var)) {
            globx = as.character(substitute(obj))
        }
        else {
            if (is.name(var)) {
                if (isFALSE(vareval))
                    var = as.character(var)
                else var = eval(var, parent.frame())
            }
            else if (!is.character(var)) {
                stop("var must be coercible to a character")
            }
            if (inherits(var, "character")) {
                globx = var
            }
            else {
                globx = as.character(substitute(var))
            }
        }
        if (verbose)
            message("variable being assigned to ", globx)
        assign(globx, value = obj, envir = envir)
        if (return_obj) {
            invisible(obj)
        }
        else {
            invisible()
        }
    }

    overwritefun = function (newfun, oldfun, package, envir = globalenv())
    {
        if (is.character(newfun) && is.character(oldfun) && missing(package))
            stop("must specify package for oldfun")
        if (!missing(package)) {
            if (is.character(package))
                envpkg = asNamespace(package)
            else if (isNamespace(package))
                envpkg = package
        } else {
            if (missing(package)) {
                envpkg = asNamespace(environment(oldfun))
            }
        }
        if (!is.character(oldfun)) {
            oldfun = deparse(tail(as.list(substitute(oldfun)), 1)[[1]])
        }
        if (!is.character(newfun)) {
            newfunenv = environment(newfun)
            newfun = deparse(tail(as.list(substitute(newfun)), 1)[[1]])
        } else {
            newfunenv = parent.frame()
        }
        nmpkg = environmentName(envpkg)
        tmpfun = get(oldfun, envir = envpkg)
        .newfun = get(newfun, envir = newfunenv)
        envtmpfun = environment(tmpfun)
        if (!is.null(envtmpfun))
            environment(.newfun) = envtmpfun
        attributes(.newfun) = attributes(tmpfun)
        evalq(asn2(oldfun, .newfun, ns = nmpkg), environment(), parent.frame())
        globasn(.newfun, oldfun, vareval = T)
        invisible()
    }

    asn2 = function (x, value, ns, pos = -1, envir = as.environment(pos))
    {
        nf = sys.nframe()
        if (missing(ns)) {
            nm = attr(envir, "name", exact = TRUE)
            if (is.null(nm) || substr(nm, 1L, 8L) != "package:")
                stop("environment specified is not a package")
            ns = asNamespace(substring(nm, 9L))
        }
        else ns = asNamespace(ns)
        ns_name = getNamespaceName(ns)
        if (bindingIsLocked(x, ns)) {
            in_load = Sys.getenv("_R_NS_LOAD_")
            if (nzchar(in_load)) {
                if (in_load != ns_name) {
                    msg = gettextf("changing locked binding for %s in %s whilst loading %s",
                                    sQuote(x), sQuote(ns_name), sQuote(in_load))
                    if (!in_load %in% c("Matrix", "SparseM"))
                        warning(msg, call. = FALSE, domain = NA, immediate. = TRUE)
                }
            }
            else if (nzchar(Sys.getenv("_R_WARN_ON_LOCKED_BINDINGS_"))) {
                warning(gettextf("changing locked binding for %s in %s",
                                 sQuote(x), sQuote(ns_name)), call. = FALSE, domain = NA,
                        immediate. = TRUE)
            }
            unlockBinding(x, ns)
            assign(x, value, envir = ns, inherits = FALSE)
            w = options("warn")
            on.exit(options(w))
            options(warn = -1)
            lockBinding(x, ns)
        }
        else {
            assign(x, value, envir = ns, inherits = FALSE)
        }
        if (!isBaseNamespace(ns)) {
            S3 = .getNamespaceInfo(ns, "S3methods")
            if (!length(S3))
                return(invisible(NULL))
            S3names = S3[, 3L]
            if (x %in% S3names) {
                i = match(x, S3names)
                genfun = get(S3[i, 1L], mode = "function", envir = parent.frame())
                if (.isMethodsDispatchOn() && methods::is(genfun,
                                                          "genericFunction"))
                    genfun = tryCatch(methods::slot(genfun, "default")@methods$ANY,
                                       error = function(e) genfun)
                defenv = if (typeof(genfun) == "closure") {
                              environment(genfun)
                          } else {
                              .BaseNamespaceEnv
                          }
                S3Table = get(".__S3MethodsTable__.", envir = defenv)
                remappedName = paste(S3[i, 1L], S3[i, 2L], sep = ".")
                if (exists(remappedName, envir = S3Table, inherits = FALSE))
                    assign(remappedName, value, S3Table)
            }
        }
        invisible(NULL)
    }

    saveRDS = function (object, file = "", ascii = FALSE,
                         ## version = NULL,
                         version = 2,
                         compress = TRUE,
                         refhook = NULL) {
        if (is.character(file)) {
            if (file == "")
                stop("'file' must be non-empty string")
            if (!dir.exists(dirname(file)))
                system2("mkdir", c("-p", dirname(file)))
            object = object
            mode = if (ascii %in% FALSE)
                        "wb"
                    else "w"
            con = if (is.logical(compress))
                       if (compress)
                           gzfile(file, mode)
                       else file(file, mode)
                   else switch(compress, bzip2 = bzfile(file, mode), xz = xzfile(file,
                                                                                 mode), gzip = gzfile(file, mode), stop("invalid 'compress' argument: ",
                                                                                                                        compress))
            on.exit(close(con))
        }
        else if (inherits(file, "connection")) {
            if (!missing(compress))
                warning("'compress' is ignored unless 'file' is a file name")
            con = file
        }
        else stop("bad 'file' argument")
        .Internal(serializeToConn(object, con, ascii, version, refhook))
    }
    overwritefun('saveRDS','saveRDS', package = "base")

    timestamp = function ()
    {
        return(gsub("[\\:\\-]", "", gsub("\\s", "_", Sys.time())))
    }

    staveRDS = function (object, file, note = NULL, version = 2, ..., verbose = FALSE)
    {

        stamped.file = gsub(".rds$", paste(".", timestamp(), ".rds",
            sep = ""), file, ignore.case = TRUE)
        saveRDS(object, stamped.file, version = version, ...)
        if (file.exists(file)) {
            if (verbose)
                message("Removing existing ", file)
            system(paste("rm", file))
        }
        if (verbose)
            message("Symlinking ", file, " to ", stamped.file)
        system(paste("ln -sfn", normalizePath(stamped.file), file))
        if (!is.null(note)) {
            writeLines(note, paste0(stamped.file, ".readme"))
        }
    }

    ## proceed with caution
    #######################
    #######################
    #######################
    ## if you don't wrap this in a function,
    ## and just run it manually after startup,
    ## everything goes bollocks...
    startup = function() {
        `:::.new` = function (pkg, name) {
            pkg = as.character(substitute(pkg))
            name = as.character(substitute(name))
            pev = packageEvent(pkg, "onLoad")
            gh = getHook(pev)
            if (length(gh) == 0 || is.null(gh$forceall12340987)) {
                setHook(pev,
                        list("forceall12340987" = function(...) forceall(envir = asNamespace(pkg))))
            }
            out = get(name, envir = asNamespace(pkg), inherits = FALSE)
            ## forceall(envir = asNamespace(pkg), evalenvir = globalenv())
            return(out)
        }

        `::.new` = function (pkg, name) {
            pkg = as.character(substitute(pkg))
            name = as.character(substitute(name))
            pev = packageEvent(pkg, "onLoad")
            gh = getHook(pev)
            if (length(gh) == 0 || is.null(gh$forceall12340987)) {
                setHook(pev,
                        list("forceall12340987" = function(...) forceall(envir = asNamespace(pkg))))
            }
            out = getExportedValue(pkg, name)
            ## forceall(envir = asNamespace(pkg), evalenvir = globalenv())
            return(out)
        }

        overwritefun("::.new", "::", package = asNamespace("base"))
        overwritefun(":::.new", ":::", package = asNamespace("base"))
    }

    #######################
    #######################
    #######################

    Sys.setenv(R_DATATABLE_NUM_THREADS = 1)
    Sys.setenv(R_REMOTES_UPGRADE = "never")
    Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
    Sys.setenv("GENCODE_DIR" = "~/DB/GENCODE")

    Sys.setenv("BASH_FUNC_blip()" = "() { echo \"hoohah\"; }")

    Sys.setenv(DEFAULT_GENOME = "~/DB/references/hg19/human_g1k_v37_decoy.chrom.sizes")
    Sys.setenv(DEFAULT_BSGENOME = "~/DB/references/hg19/human_g1k_v37_decoy.chrom.sizes")

    ww = with
    wn = within

    go.R = function() {
        eval(
            quote(
                .libPaths(
                    unique(
                        c(## "/gpfs/commons/groups/imielinski_lab/lib/R-4.0.2_KH",
                            .libPaths()
                        )
                    )
                )
            ), globalenv()
        )
        evalq(
        {
            suppressWarnings({suppressPackageStartupMessages({
            ### begin startup.R
            .libPaths(unique(c("/gpfs/commons/groups/imielinski_lab/lib/R-4.0.2_KH", .libPaths())))

            Sys.setenv(DEFAULT_BSGENOME = "/gpfs/commons/groups/imielinski_lab/DB/GATK/human_g1k_v37_decoy.chrom.sizes")
            Sys.setenv(DEFAULT_GENOME = "/gpfs/commons/groups/imielinski_lab/DB/GATK/human_g1k_v37_decoy.chrom.sizes")

            Sys.setenv(PATH = paste("~/software/bcftools-1.9", Sys.getenv("PATH"), sep = ":"))
            Sys.setenv(BCFTOOLS_PLUGINS = "/gpfs/commons/groups/imielinski_lab/Software/bcftools-1.9/plugins")

            startup();
            library(withr)
            library(skitools)
            library(gGnome)
            library(Flow)
            library(dplyr)
            library(bamUtils)
            library(skidb)
            library(naturalsort)
            library(magrittr)
            library(signature.tools.lib)
            library(tidyr)
            library(MASS)
            library(khtools)
            library(wesanderson)
            forceload(.force = T)

            ## require3(
            ##     withr,
            ##     skitools,
            ##     gGnome,
            ##     Flow,
            ##     dplyr,
            ##     bamUtils,
            ##     skidb,
            ##     naturalsort,
            ##     magrittr,
            ##     signature.tools.lib,
            ##     tidyr,
            ##     MASS,
            ##     khtools,
            ##     wesanderson
            ## )

            ## source("~/lab/home/khadi/git/khscripts/patch_gtrack.R")

            exprs = expression({
                brewer.master <- skitools::brewer.master
                `%$%` <- gUtils::`%$%`
                `%&%` <- gUtils::`%&%`
                ## `%Q%` <- gUtils::`%Q%`
                `%Q%` <- khtools::`%Q%`
                `%+%` <- gUtils::`%+%`
                `%-%` <- gUtils::`%-%`
                `%^%` <- gUtils::`%^%`
                width <- GenomicRanges::width
                reduce <- GenomicRanges::reduce
                between <- data.table::between
                setnames <- data.table::setnames
                select <- dplyr::select
                update <- Flow::update
                cache <- Flow::cache
                set <- data.table::set
                set_names <- rlang::set_names
                matches <- dplyr::matches
                n <- dplyr::n
                last <- dplyr::last
                first <- dplyr::first
                tailf <- khtools::tailf
                coalesce <- khtools::coalesce
                ppdf <- khtools::ppdf
                ppng <- khtools::ppng
                with_libpaths <- withr::with_libpaths
                with_options <- withr::with_options
                melt <- data.table::melt
                overwritefun(khtools::gr.flipstrand, gUtils::gr.flipstrand); gr.flipstrand <- khtools::gr.flipstrand
                .S3method("merge", "Junction", merge.Junction)
                registerS3method("merge", "Junction", merge.Junction, envir = globalenv())
                .S3method("dcast", "data.table", dcast.data.table)
                registerS3method("dcast", "data.table", dcast.data.table, envir = globalenv())
            })

            eix = seq(2, length(exprs[[1]]))
            for (i in eix) {
                print(exprs[[1]][[i]])
                tryCatch({
                    eval(exprs[[1]][[i]], envir = globalenv())
                }, error = function(e) NULL)
            }

            ## tryCatch({
            ##     .S3method("merge", "Junction", merge.Junction)
            ##     registerS3method("merge", "Junction", merge.Junction, envir = parent.frame())
            ## }, error = function(e) NULL)



            ## tryCatch({
            ##     .S3method("dcast", "data.table", dcast.data.table)
            ##     registerS3method("dcast", "data.table", dcast.data.table, envir = parent.frame())
            ## }, error = function(e) NULL)
            ### end startup.R
            })})
        }, globalenv()
        )
    }

    do.dev = function() {
        ## evalq({startup(); library3(devtools)}, globalenv())
        ## eval(quote(.libPaths(unique(c("/gpfs/commons/groups/imielinski_lab/lib/R-4.0.2_KH", .libPaths())))), globalenv())
        eval(quote({
            startup();
            require3(devtools, withr, roxygen2);
            with_libpaths = withr::with_libpaths;
            iinstall = function(...) {
                pf = parent.frame();
                eval(quote({install(dependencies = F, quick = F)}), envir = pf)
            };
            dinstall = function(...) {
                pf = parent.frame();
                eval(quote({document(); install(dependencies = F, quick = F)}), envir = pf)
            };
            bla = ""}), globalenv())
        invisible()
    }

    private_lib = function(suffix = "_KH") {
        libs = .libPaths()
        orig = utils::tail(libs, 1)
        addon = Sys.getenv("R_LIBS")
        addon = normalizePath(unlist(strsplit(addon, "[,|;:]")))
        priv_lib = setdiff(libs, union(addon, orig))
        if (length(priv_lib) == 0) {
            if (dir.exists(suffix))
                priv_lib = suffix
            else
                priv_lib = paste0(.libPaths()[1], suffix)
            expr = parse(text = paste0(".libPaths(c(", paste(paste0("'", unique(c(priv_lib, .libPaths())), "'"), collapse = ","), "))"))
            message("setting library path(s) to: ", paste0(unique(c(priv_lib, .libPaths())), collapse = ", "))
            eval(expr, globalenv())
        } else {
            expr = parse(text = paste0(".libPaths(c(", paste(paste0("'", union(addon, orig), "'"), collapse = ","), "))"))
            message("setting library path(s) to: ", paste0(union(addon, orig)), collapse = ", ")
            eval(expr, globalenv())
        }
        invisible()
    }

    test.start = function() {
        eval(quote({
            startup(); library3(khtools, skitools);
            tailf = khtools::tailf
        }), globalenv())
        ## evalq({
        ##     startup(); library3(khtools, skitools);
        ##     tailf = khtools::tailf
        ## }, globalenv())
    }

    system.time2 = function (expr, gcFirst = FALSE)
    {
        ppt = function(y) {
            if (!is.na(y[4L]))
                y[1L] = y[1L] + y[4L]
            if (!is.na(y[5L]))
                y[2L] = y[2L] + y[5L]
            paste(formatC(y[1L:3L]), collapse = " ")
        }
        if (gcFirst)
            gc(FALSE)
        time = proc.time()
        on.exit(message("Timing stopped at: ", ppt(proc.time() -
                                                   time)))
        expr
        new.time = proc.time()
        on.exit()
        structure(new.time - time, class = "proc_time")
    }
    overwritefun("system.time2", "system.time", "base")

    somejit = function(x, factor = 1e-6) {set.seed(10); jitter(x, factor = factor)}

    ### end .Rprofile
    require3(khtools, glmnet, randomForest, caret, dplyr, signature.tools.lib, pROC)


    ## robust.scale <- function(x, qlow = 0.1, qhigh = 0.9) {
    ##   out = x
    ##   ql = quantile(out, qlow)
    ##   qh = quantile(out, qhigh)
    ##   out = pmax(out, ql)
    ##   out = pmin(out, qh)
    ##   return(out)
    ## }

    robust.scale <- function(x, qlow = 0.1, qhigh = 0.9) {
        MEDIAN =  median(x, na.rm = T)
        (x - MEDIAN) / quantile(x, qlow)
    }



    mystats = function(x, na.rm = TRUE) {
        if (is.null(dim(x))) {
            x = as.data.frame(x)
        }
        fun = function(x) {
            if (is.numeric(x)) {
                c(
                    mean = mean(x, na.rm = na.rm),
                    median = median(x, na.rm = na.rm),
                    sd = sd(x, na.rm = na.rm),
                    mad = mad(x, na.rm = na.rm),
                    var = stats::var(x, na.rm = na.rm),
                    min = min(x, na.rm = na.rm),
                    max = max(x, na.rm = na.rm),
                    q25 = unname(quantile(x, 0.25, na.rm = na.rm)),
                    q75 = unname(quantile(x, 0.75, na.rm = na.rm)),
                    q10 = unname(quantile(x, 0.10, na.rm = na.rm)),
                    q90 = unname(quantile(x, 0.90, na.rm = na.rm))
                )
            } else {
                c(
                    mean = NA_real_,
                    median = NA_real_,
                    sd = NA_real_,
                    mad = NA_real_,
                    var = NA_real_,
                    min = NA_real_,
                    max = NA_real_,
                    q25 = NA_real_,
                    q75 = NA_real_,
                    q10 = NA_real_,
                    q90 = NA_real_
                )
            }
        }
        cn = colnames(x)
        out.lst = lapply(x, fun)
        out.lst
        ## transp(out.lst, c)
    }


    my_preprocess = function(dat) {
        lapply(dat, function(x) {
            if (inherits(x, c("numeric", "integer"))) {
                unlist(mystats(x))
            }
        })
    }

    my_predict = function(mypre, newdat, method = "robust") {
        goods = sapply(mypre, Negate(is.null))
        nm = names(goods[goods])
        qm = qmat(newdat,,nm)
        if (identical(method, "robust")) {
            mp = mapply(function(x,y) {
                (x - y["x.median"]) / (y["x.q90"] - y["x.q10"])
            }, qm, mypre[nm], SIMPLIFY = F)
        }
        do.assign(newdat, setColnames(as.data.frame(mp), nm))
    }


    do.pp = function(dat, vars, method = c("center", "scale"), prefun = NULL, rangeBounds = c(0, 1), verbose = TRUE) {
        require(caret)
        require(dplyr)
        odat = copy3(dat)
        if (!is.null(prefun) && is.function(prefun)) {
            for (v in vars)
                dat[[v]] = prefun(dat[[v]])
        } else {
            prefun = NULL
        }
        if (!is.null(method) && is.character(method)){
            if (identical(method, "robust")) {
                pp = dat %>% select(!!vars) %>% my_preprocess
                dat = my_predict(pp, dat)
            } else {
                pp = preProcess(dat %>% select(!!vars), method = method, rangeBounds = rangeBounds, verbose = verbose)
                dat = predict(pp, dat)
            }
        } else {
            pp = NULL
        }
        return(structure(list(pp = pp, dat = dat, prefun = prefun, vars = vars, method = method, orig_dat = odat), class = "carprep"))
    }

    predict.pp = function(pp.res, newdat, apply.prefun = TRUE) {
        require(caret)
        require(dplyr)
        odat = copy3(newdat)
        if (is.function(apply.prefun) || (!is.null(pp.res$prefun) && isTRUE(apply.prefun) && is.function(pp.res$prefun))) {
          if (is.function(apply.prefun)) pp.res$prefun = apply.prefun
          for (v in pp.res$vars)
            newdat[[v]] = pp.res$prefun(newdat[[v]])
        } else {
          pp.res$prefun = NULL
        }
        if (!is.null(pp.res$pp) && !identical(pp.res$method, "robust")) {
            newdat = predict(pp.res$pp, newdat)
        } else if (identical(pp.res$method, "robust")) {
            newdat = my_predict(pp.res$pp, newdat)
        }
        return(structure(list(pp = pp.res$pp, newdat = newdat, prefun = pp.res$prefun, vars = pp.res$vars, method = pp.res$method, orig_newdat = odat), class = "carpost"))
    }


    ## nysv = c("tib", "ihdel", "qrdup", "qrdel")
    ## uksv = c("RS3", "RS5")
    ## fullmod = c("RS3", "RS5", "del.mh.prop", "hrd", "SNV8", "SNV3")

    ## fitglmnet = function(var, lambda = NULL, family = "multinomial", ycol, dat) {
    ##   if (missing(ycol)) ycol = "fmut"
    ##   if (missing(dat)) {
    ##     message("attempting to grab variable 'fulldat' from env stack")
    ##     fulldat = dg(fulldat)
    ##   } else {
    ##     fulldat = dat
    ##   }
    ##   source("/gpfs/commons/groups/imielinski_lab/home/khadi/git/khscripts/.Rprofile")
    ##   require3(glmnet)
    ##   set.seed(10);
    ##   x = select(fulldat, !!var) %>% asm
    ##   ## x.train = model.matrix(~., (select(traind, !!vars)))
    ##   y = fulldat[[ycol]]
    ##   if (family == "multinomial") {
    ##       y = fulldat[[ycol]]
    ##       modglm = glmnet(x = x,
    ##                       y = y, family = family, type.multinomial = "grouped",
    ##                       lambda = lambda, lower.limits = rep_len(0, NCOL(x)),
    ##                       n_lambda = 1000)
    ##   } else if (family == "binomial") {
    ##       y = refactor(fulldat[[ycol]], "WT") %>% fct_recode("HRD" = "OTHER") %>% relevel("WT")
    ##       modglm = glmnet(x = x,
    ##                       y = y, family = family,
    ##                       lambda = lambda, lower.limits = rep_len(0, NCOL(x)),
    ##                       n_lambda = 1000)
    ##   }
    ##   return(modglm)
    ## }

    fitglmnet = function(var, lambda = NULL, family = "multinomial", ycol, dat) {
      if (missing(ycol)) ycol = "fmut"
      if (missing(dat)) {
        message("attempting to grab variable 'fulldat' from env stack")
        fulldat = dg(fulldat)
      } else {
        fulldat = dat
      }
      ## source("/gpfs/commons/groups/imielinski_lab/home/khadi/git/khscripts/.Rprofile")
      ## require3(glmnet)
      set.seed(10);
      x = select(fulldat, !!var) %>% asm
      ## x.train = model.matrix(~., (select(traind, !!vars)))
      y = fulldat[[ycol]]
      if (family == "multinomial") {
          y = fulldat[[ycol]]
          modglm = cv.glmnet(x = x,
                          y = y, family = family, type.multinomial = "grouped",
                          lower.limits = rep_len(0, NCOL(x)))
      } else if (family == "binomial") {
          y = refactor(fulldat[[ycol]], "WT") %>% fct_recode("HRD" = "OTHER") %>% relevel("WT")
          modglm = cv.glmnet(x = x,
                          y = y, family = family,
                          lower.limits = rep_len(0, NCOL(x)))
      }
      return(modglm)
    }


    # mylambdas <- function(M = 2, N = 10, step = 0.005) {
    #     return(10^(seq(M, -N, -abs(step))))
    # }

    mylambdas <- function(M = 2, N = 10, step = 0.1) {
        return(10^(seq(M, -N, -abs(step))))
    }


    fitcvglmnet = function(var, lambda = mylambdas(), family = "multinomial", ycol, dat,
                           only_positive = TRUE, nlambda = 100, maxit = 1e5, thresh = 1e-07,
                           type.multinomial = "grouped",
                           standardize = FALSE, standardize.response = FALSE,
                           type.measure = "default", alpha = 1,
                           ...
                           ) {
      if (isFALSE(standardize)) {
        message("standardize is FALSE")
        message("you must standardize the features yourself")
      }
      if (isTRUE(standardize)) message("standardize is TRUE")
      if (missing(ycol)) ycol = "fmut"
      if (missing(dat)) {
        message("attempting to grab variable 'fulldat' from env stack")
        fulldat = dg(fulldat)
      } else {
        fulldat = dat
      }
      ## source("/gpfs/commons/groups/imielinski_lab/home/khadi/git/khscripts/.Rprofile")
      ## require3(glmnet)
      set.seed(10);
      x = select(fulldat, !!var) %>% asm
      ## x.train = model.matrix(~., (select(traind, !!vars)))
      y = fulldat[[ycol]]
      if (only_positive)
          llim = rep_len(0.0, NCOL(x))
      else
          llim = rep_len(-Inf, NCOL(x))
      if (family == "multinomial") {
          ## y = fulldat[[ycol]]
          modglm = cv.glmnet(x = x,
                          y = y, family = family, type.multinomial = type.multinomial,
                          lower.limits = llim,
                          lambda = lambda,
                          nlambda = nlambda,
                          maxit = maxit,
                          thresh = thresh,
                          standardize = standardize,
                          standardize.response = standardize.response,
                          type.measure = type.measure, alpha = alpha, ...
                          )
      } else if (family == "binomial") {
          ## y = refactor(fulldat[[ycol]], "WT") %>% fct_recode("HRD" = "OTHER") %>% relevel("WT")
          modglm = cv.glmnet(x = x,
                          y = y, family = family,
                          lower.limits = llim,
                          lambda = lambda,
                          nlambda = nlambda,
                          maxit = maxit,
                          thresh = thresh,
                          standardize = standardize,
                          standardize.response = standardize.response,
                          type.measure = type.measure, alpha = alpha, ...
                          )
      }
      return(modglm)
    }


    fitglmnet = function(var, lambda = 0, family = "multinomial", ycol, dat,
                         only_positive = TRUE, maxit = 1e5, thresh = 1e-07, standardize = FALSE) {
      if (isFALSE(standardize)) {
        message("standardize is FALSE")
        message("you must standardize the features yourself")
      }
      if (isTRUE(standardize)) message("standardize is TRUE")
      if (missing(ycol)) ycol = "fmut"
      if (missing(dat)) {
        message("attempting to grab variable 'fulldat' from env stack")
        fulldat = dg(fulldat)
      } else {
        fulldat = dat
      }
      set.seed(10);
      x = select(fulldat, !!var) %>% asm
      y = fulldat[[ycol]]
      if (only_positive)
          llim = rep_len(0.0, NCOL(x))
      else
          llim = rep_len(-Inf, NCOL(x))
      if (family == "multinomial")
      {
          modglm = glmnet(
              x = x,
              y = y, family = family,
              lower.limits = llim,
              lambda = lambda, ## thresh = 1e-14, maxit = 10e6,
              trace.it = 10,
              maxit = maxit,
              thresh = thresh,
              standardize = standardize
          )
      } else if (family == "binomial") {
          modglm = glmnet(
              x = x,
              y = y, family = family,
              lower.limits = llim,
              lambda = lambda, ## thresh = 1e-14, maxit = 10e6,
              trace.it = 10,
              maxit = maxit,
              thresh = thresh,
              standardize = standardize
          )
      }
      return(modglm)
    }


    fitrandomforest <- function(var, ycol, dat, ...) {
      if (missing(ycol)) ycol = "fmut"
      if (missing(dat)) {
        message("attempting to grab variable 'fulldat' from env stack")
        fulldat = dg(fulldat)
      } else {
        fulldat = dat
      }
      ## source("/gpfs/commons/groups/imielinski_lab/home/khadi/git/khscripts/.Rprofile")
      ## require3(randomForest)
      set.seed(10);
      x = select(fulldat, !!var) %>% asm
      ## x.train = model.matrix(~., (select(traind, !!vars)))
      y = fulldat[[ycol]]
      rfdat = et(sprintf("cbind(asdf(x), %s = y)", ycol))
      ## rfdat = cbind(asdf(x), mutstat = y)
      form = formula(sprintf("%s ~ %s", ycol, paste(var, collapse = "+")))
      message("formula: ", deparse(form))
      rfmod = randomForest(form, data = rfdat, ntree = 1000, importance = TRUE, ...)
      return(rfmod)
    }



    fit_classifier <- function(mod, data, s = 0, alpha = 1, exact = FALSE) {

      data = copy(data)

      cls.mod = class(mod)[1]
      if (identical(cls.mod, "cv.glmnet")) {
        gmod = mod$glmnet.fit
      } else {
        gmod = mod
      }
      classnames = gmod$classnames

      for (thisc in c(classnames, "NRSUM"))
        data[[thisc]] = NULL


      ## debugonce(predict.multnet)
      if (inherits(gmod$beta, "list"))
        vars = rownames2(gmod$beta[[1]])
      else
        vars = rownames2(gmod$beta)

      scores = get_pred(mod, data, vars, type = "response", s = s, alpha = alpha, exact = exact)


      classes = factor(
        get_pred(mod, data, vars, type = "class", s = s, alpha = alpha, exact = exact),
        classnames
      )

      data$classes = classes
      scores = asdf(scores)

      if (NCOL(scores) == 1)
        colnames(scores) = classnames[2]

      for (i in 1:NCOL(scores))
        data[[colnames(scores)[[i]]]] = scores[[i]]

      data$NRSUM = rowSums(qmat(data,,classnames[-1]))

      return(data)

    }



    fit_rforest <- function (mod, data) {
        data = copy(data)
        classnames = mod$classes
        scores = as.data.frame.matrix(predict(mod, data, type = "prob"))
        out = copy(data)
        for (i in 1:NCOL(scores))
            out[[colnames(scores)[[i]]]] = scores[[i]]
        ## out = out %>% mutate(scores)
        out$classes = predict(mod, data, type = "class")
        out$NRSUM = rowSums(qmat(out,,classnames[-1]))

        ## data$NRSUM = rowSums(qmat(data, , classnames[-1]))
        return(out)
    }



    glmnet_aic = function(fit) {
      tLL <- fit$nulldev - deviance(fit)
      k <- fit$df
      n <- fit$nobs
      AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
      AICc
    }


    get_pred <- function(mod, newdata, feats, type = "response", s = 0, alpha = 1, exact = FALSE) {
      cls.mod = class(mod)
      if (any(cls.mod == "cv.glmnet")) cls.mod = class(mod$glmnet.fit)[1]
      if (any(cls.mod == "multnet")) {
        prd = predict(mod, select(newdata, !!feats) %>% asm, type = type, s = s, alpha = alpha, type.multinomial = "grouped", exact = exact)
        if (length(dim(prd)) == 3) {
          dnm2 = dimnames(prd)[[2]]
          dim(prd) = dim(prd)[1:2]
          dimnames(prd) = list(NULL, dnm2)
        }
      } else if (any(cls.mod == "lognet")) {
        prd = predict(mod, select(newdata, !!feats) %>% asm, type = type, s = s, alpha = alpha, exact = exact) %>%
          setColnames(mod$classnames[2])
      } else if (any(cls.mod == "randomForest")) {
          ## type means different things for random forest
          if (identical(type, "class"))
              type = "response"
          else if (identical(type, "response"))
              type = "prob"
          prd = predict(mod, newdata, type = type)
      }
      return(prd)
    }


    make_roc <- function(dat, lab = "lab", score = "HRD", include_group = FALSE) {
      mrlab = mroclab(dat[[lab]])
      mrpred = mrocpred(et(sprintf("dat[, %s, drop=F]", mkst(paste0("'", score, "'")))))
      colslab = gsub("(\\w+)[ ](\\w+)", "\\1", colnames(mrlab))
      colspred = gsub("(\\w+)[ ](\\w+)", "\\1", colnames(mrpred))
      goodcols = intersect(colslab, colspred)
      mrlab = asdf(et(sprintf("mrlab[, %s,drop=F]", mkst(colslab %in% goodcols))))
      mrpred = asdf(et(sprintf("mrpred[, %s,drop=F]", mkst(colspred %in% goodcols))))

      mrocdat = function (lbl, prd) {
        prd00 = rbind(0, asm(prd))
        mg = setcols(with((melt(prd00)), g2()[order(Var2, value),
          ]), c("Var2", "value"), c("Group", "prd"))[, c("Group",
            "prd"), drop = F]
        cb = cbind(lbl, prd)
        roc_res = multiROC::multi_roc(cb, force_diag = T)
        plot_roc_df <- multiROC::plot_roc_data(roc_res)
        if (!include_group) {
          gdat = asdt(plot_roc_df)[Group %nin% c("Macro", "Micro")]## [,
            ## `:=`(prd, rev(mg$prd))]
        } else {
          gdat = asdt(plot_roc_df)## [Group %nin% c("Macro", "Micro"),
            ## `:=`(prd, rev(mg$prd))]
        }
        gdat[, `:=`(Group, trimws(Group))]
        return(gdat)
      }

      mrdat = suppressWarnings(mrocdat(
        mrlab,
        mrpred
      ))
      return(withAutoprint(mrdat, echo = F)$value)
    }



    score_glmnet = function(data, features, family, labels = "lab", lambda = 0, mod) {
      if (missing(mod))
        mod = fitglmnet(features, family = family, dat = data, ycol = labels, lambda = lambda)
      if (class(mod)[1] == "lognet")
        ## data[[levels(data[[labels]])[-1]]] = get_pred(mod, newdata = data, features)
        data[[mod$classnames[-1]]] = as.vector(get_pred(mod, newdata = data, features))
      else if (class(mod)[1] == "multnet") {
        scores = get_pred(mod, newdata = data, features)
        ## lvs = levels(data[[labels]])
        lvs = mod$classnames[-1]
        for (lv in lvs) {
          data[[lv]] = scores[,lv]
        }
      }
      list(data = data, model = mod, features = features)
    }


    score_randomforest = function(data, features, labels = "lab", lambda = 0, mod) {
      if (missing(mod)) {
        mod = fitrandomforest(features, dat = data, ycol = labels)
      }
      if (class(mod)[1] == "lognet")
        ## data[[levels(data[[labels]])[-1]]] = get_pred(mod, newdata = data, features)
        data[[mod$classnames[-1]]] = as.vector(get_pred(mod, newdata = data, features))
      else if (class(mod)[1] == "multnet") {
        scores = get_pred(mod, newdata = data, features)
        ## lvs = levels(data[[labels]])
        lvs = mod$classnames[-1]
        for (lv in lvs) {
          data[[lv]] = scores[,lv]
        }
      } else if (any(grepl("randomForest", class(mod)))) {
        scores = predict(mod, newdata = asdf(data), type = "prob")
        lab = predict(mod, newdata = asdf(data))
        if (length(mod$classes) == 2)
          lvs = mod$classes[-1]
        else
          lvs = mod$classes
        for (lv in lvs) {
          data[[lv]] = scores[,lv]
        }
        data[["class_lab"]] = lab
      }
      list(data = data, model = mod, features = features)
    }

    feature_importance <- function(ixfold, xfolds, dat, lambda = 0, vars, seed = 10, ycol, debug = FALSE,
                                   rpart = FALSE, permute_varlist = NULL, added_testd = NULL,
                                   caret_preprocess = do.pp(traind, vars = vars, method = "range", prefun = log1p),
                                   only_positive = FALSE, standardize = FALSE, model.type = c("glm", "glmnet", "randomforest"),
                                   s = 0, alpha = 1
                                   ) {
        tryCatch({
            if (debug) browser()
            nm = names(xfolds[ixfold])
            xfold = xfolds[[ixfold]]
            message("testing split: ", ixfold, " ", nm)
            traind = copy(dat[xfold$train,])
            testd = copy(dat[xfold$test,])

            if (!is.null(added_testd)) {
                added_testd = copy(added_testd[pair %in% setdiff(pair, c(traind$pair, testd$pair))])
                if (NROW(added_testd) > 0)
                    testd = rbind(testd[, added := FALSE], added_testd[, added := TRUE], fill = T)
            } else {
                testd$added = FALSE
            }

            set.seed(seed);

            traind = copy3(traind)
            traind$fmut = traind[[ycol]]
            testd = copy3(testd)
            testd$fmut = testd[[ycol]]

            if (!is.null(caret_preprocess)) {
                if (is.character(caret_preprocess) && caret_preprocess[1] %in% c("range", "scale")) {
                    pp = preProcess(select(traind, !!vars), method = caret_preprocess)
                    traind = predict(pp, traind)
                    testd = predict(pp, testd)
                } else if (is.function(caret_preprocess)) {
                    traind = traind %>% mutate_at(vars(!!vars), caret_preprocess)
                    testd = testd %>% mutate_at(vars(!!vars), caret_preprocess)
                } else if (inherits(caret_preprocess, "carprep")) {
                    traind = caret_preprocess$dat
                    prep.res = predict.pp(caret_preprocess, testd)
                    testd = prep.res$newdat
                } else if (inherits(caret_preprocess, c("expression", "call"))) {
                    caret_preprocess = eval(caret_preprocess)
                    traind = caret_preprocess$dat
                    prep.res = predict.pp(caret_preprocess, testd)
                    testd = prep.res$newdat
                }
            }


            if (length(levels(traind$fmut)) > 2) {
                family = "multinomial"
            } else {
                family = "binomial"
            }

            model.type = intersect(model.type, c("glm", "glmnet", "randomforest"))

            if (length(model.type) > 1) {
                model.type = model.type[1]
            }



            if (identical(model.type, "glm")) {
                trainglm = fitcvglmnet(vars, family = family, dat = traind, ycol = ycol, only_positive = only_positive, standardize = standardize, lambda = mylambdas())
                testd = fit_classifier(trainglm, testd, s = s, alpha = alpha)
            } else if (identical(model.type, "glmnet")) {
                trainglm = fitglmnet(vars, family = family, dat = traind, ycol = ycol, only_positive = only_positive, standardize = standardize, lambda = mylambdas())
                testd = fit_classifier(trainglm, testd, s = s, alpha = alpha)
            } else if (identical(model.type, "randomforest")) {
                trainglm = fitrandomforest(var = vars, dat = traind, ycol = ycol)
                testd = fit_rforest(trainglm, testd)
            }

            if (any(class(trainglm) %in% c("cv.glmnet", "glmnet"))) {
                co = as.data.frame(as.matrix(coef(trainglm, s = s))) %>% rn2col("feat") %>% as.data.table
                co$nm = nm
                # co$xfold = xfold
            } else {
                co = NULL
            }


            cls.mod = class(trainglm)[1]
            if (identical(cls.mod, "cv.glmnet")) {
                gmod = trainglm$glmnet.fit
            } else {
                gmod = trainglm
            }
            if (any(class(gmod) == "glmnet"))
                classnames = gmod$classnames
            else if (any(class(gmod) == "randomForest"))
                classnames = gmod$classes

            if (length(classnames) == 2) {
                classnames = classnames[2]
            }
            scores = qmat(testd,,classnames)


            testd$fold_id = nm
            testd$Method = "all"

            off_diag = function(x) {
                diag(x) = NA
                as.vector(x) %>% na.omit
            }

            get_accuracy = function(confus_mat) {
                correct = diag(confus_mat)
                sum(correct) / sum(off_diag(confus_mat), correct)
            }

            indiv.levels = levels(testd[[ycol]])[-1]
            full_mod_indiv_accuracy = c()
            if (length(indiv.levels) > 1) {
                for (i in indiv.levels) {
                    indiv.mat = as.matrix(table(refactor(testd[[ycol]], i),
                      refactor(factor(testd$classes, levels = levels(testd[[ycol]])), i)))
                    full_mod_indiv_accuracy = c(full_mod_indiv_accuracy, setNames(get_accuracy(indiv.mat), i))
                }
            }

            confus_mat = as.matrix(table(testd[[ycol]], factor(testd$classes, levels = levels(testd[[ycol]]))))
            full_mod_accuracy = c(FULL = get_accuracy(confus_mat))

            lvars = as.list(vars)

            if (!is.null(permute_varlist)) {
                ## .NotYetImplemented()
                if (is.character(permute_varlist))
                    pvar = list(permute_varlist)
                else if (is.list(permute_varlist)) {
                    pvar = permute_varlist
                } else if (! is.list(permute_varlist)) {
                    warning("ignoring permute_list")
                    pvar = list()
                }
                lvars = c(lvars, pvar)
            }

            permute_lst = purrr::transpose(lapply(lvars, function(v, seed = 10) {
                testd_permute = copy(testd)
                testd_permute$fold_id = nm
                set.seed(seed)
                for (i in 1:NROW(v))
                    testd_permute[[v[i]]] = sample(testd_permute[[v[i]]])
                testd_permute$Method = paste(v, collapse = ", ")
                permute_classes = get_pred(trainglm, testd_permute, vars, type = "class")
                permut_s = asdf(get_pred(trainglm, testd_permute, vars, type = "response"))
                if (NCOL(permut_s) == 1)
                    colnames(permut_s) = levels(testd[[ycol]])[2]
                for (i in seq_len(NCOL(permut_s)))
                    testd_permute[[names(permut_s)[i]]] = permut_s[[i]]
                permute_roc = make_roc(testd_permute, ycol, colnames(permut_s))[
                   ,Method := testd_permute$Method[1]][
                   ,fold_id := nm]

                indiv.levels = levels(testd[[ycol]])[-1]
                indiv.accuracy = c()
                if (length(indiv.levels) > 1) {
                    for (i in indiv.levels) {
                        indiv.mat = as.matrix(table(refactor(testd[[ycol]], i),
                      refactor(factor(permute_classes, levels = levels(testd[[ycol]])), i)))
                        indiv.accuracy = c(indiv.accuracy, setNames(get_accuracy(indiv.mat), i))
                    }
                }
                indiv.dec.accuracy = full_mod_indiv_accuracy[names(indiv.accuracy)] - indiv.accuracy
                confus_mat_permute = as.matrix(table(testd[[ycol]], factor(permute_classes, levels = levels(testd[[ycol]]))))
                dec_accuracy = data.table(feature = paste(v, collapse = ", "), decrease_accuracy = full_mod_accuracy - get_accuracy(confus_mat_permute))[, which := "FULL"]
                dec_accuracy = dec_accuracy[rep_len(1:NROW(dec_accuracy), NROW(dec_accuracy) + NROW(indiv.accuracy))]
                dec_accuracy[-1, decrease_accuracy := indiv.dec.accuracy]
                dec_accuracy[-1, which := names(indiv.dec.accuracy)]
                list(dec_accuracy = dec_accuracy,
                     roc = permute_roc,
                     testd = testd_permute)
            }))

            dec_accuracy = rbindlist(permute_lst$dec_accuracy)
            dec_accuracy$fold_id = nm

            bigroc = rbind(rbindlist(permute_lst$roc),
                           make_roc(testd, ycol, classnames)[, Method := "all"][, fold_id := nm])


            ## mlab = mroclab(testd[[ycol]])
            ## mpred = mrocpred(get_pred(trainglm, testd, vars))

            ## if (class(trainglm)[1] == "lognet")
            ##     mlab = mlab[,2,drop=F]

            outd = rbind(rbindlist(permute_lst$testd), testd)

            ret = list(
                dec_accuracy = dec_accuracy,
                ## lab = mlab, mpred = mpred,
                roc = bigroc,
                classes = testd$classes,
                scores = scores,
                outd = outd,
                coef = co
            )
            return(ret)
        }, error = function(e) printerr(ixfold))
    };  ifun <- feature_importance


    ## feature_importance <- function(ixfold, xfolds, dat, lambda = 0, vars, seed = 10, ycol, debug = FALSE, rpart = FALSE, permute_varlist = NULL, added_testd = NULL, caret_preprocess = NULL, only_positive = TRUE) {
    ##     tryCatch({
    ##     if (debug) browser()
    ##     nm = names(xfolds[ixfold])
    ##     xfold = xfolds[[ixfold]]
    ##     traind = copy(dat[xfold$train,])
    ##     testd = copy(dat[xfold$test,])

    ##     if (!is.null(added_testd)) {
    ##         added_testd = copy(added_testd[pair %in% setdiff(pair, c(traind$pair, testd$pair))])
    ##         if (NROW(added_testd) > 0)
    ##             testd = rbind(testd[, added := FALSE], added_testd[, added := TRUE])
    ##     } else {
    ##         testd$added = FALSE
    ##     }

    ##     set.seed(seed);

    ##     traind = copy3(traind)
    ##     traind$fmut = traind[[ycol]]
    ##     testd = copy3(testd)
    ##     testd$fmut = testd[[ycol]]

    ##     if (!is.null(caret_preprocess)) {
    ##         if (is.character(caret_preprocess) && caret_preprocess[1] %in% c("range", "scale")) {
    ##             pp = preProcess(select(traind, !!vars), method = caret_preprocess)
    ##             traind = predict(pp, traind)
    ##             testd = predict(pp, testd)
    ##         } else if (is.function(caret_preprocess)) {
    ##             traind = traind %>% mutate_at(vars(!!vars), caret_preprocess)
    ##             testd = testd %>% mutate_at(vars(!!vars), caret_preprocess)
    ##         }
    ##     }


    ##     if (length(levels(traind$fmut)) > 2) {
    ##         family = "multinomial"
    ##     } else {
    ##         family = "binomial"
    ##     }

    ##     trainglm = fitcvglmnet(vars, family = family, dat = traind, ycol = ycol, only_positive = only_positive)

    ##     classes = get_pred(trainglm, testd, vars, type = "class")
    ##     scores = get_pred(trainglm, testd, vars, type = "response")

    ##     testd$fold_id = nm
    ##     testd$Method = "all"



    ##     testd$classes = classes
    ##     if (!inherits(scores, "data.frame"))
    ##         scores = asdf(scores)

    ##     if (NCOL(scores) == 1)
    ##         colnames(scores) = levels(testd[[ycol]])[2]

    ##     for (i in seq_len(NCOL(scores)))
    ##         testd[[names(scores)[i]]] = scores[[i]]

    ##     off_diag = function(x) {
    ##         diag(x) = NA
    ##         as.vector(x) %>% na.omit
    ##     }

    ##     get_accuracy = function(confus_mat) {
    ##         correct = diag(confus_mat)
    ##         sum(correct) / sum(off_diag(confus_mat), correct)
    ##     }

    ##     indiv.levels = levels(testd[[ycol]])[-1]
    ##     full_mod_indiv_accuracy = c()
    ##     if (length(indiv.levels) > 1) {
    ##         for (i in indiv.levels) {
    ##             indiv.mat = as.matrix(table(refactor(testd[[ycol]], i),
    ##               refactor(factor(testd$classes, levels = levels(testd[[ycol]])), i)))
    ##             full_mod_indiv_accuracy = c(full_mod_indiv_accuracy, setNames(get_accuracy(indiv.mat), i))
    ##         }
    ##     }

    ##     confus_mat = as.matrix(table(testd[[ycol]], factor(testd$classes, levels = levels(testd[[ycol]]))))
    ##     full_mod_accuracy = c(FULL = get_accuracy(confus_mat))

    ##     lvars = as.list(vars)

    ##     if (!is.null(permute_varlist)) {
    ##         ## .NotYetImplemented()
    ##         if (is.character(permute_varlist))
    ##             pvar = list(permute_varlist)
    ##         else if (is.list(permute_varlist)) {
    ##             pvar = permute_varlist
    ##         } else if (! is.list(permute_varlist)) {
    ##             warning("ignoring permute_list")
    ##             pvar = list()
    ##         }
    ##         lvars = c(lvars, pvar)
    ##     }

    ##     permute_lst = purrr::transpose(lapply(lvars, function(v, seed = 10) {
    ##         testd_permute = copy(testd)
    ##         testd_permute$fold_id = nm
    ##         set.seed(seed)
    ##         for (i in 1:NROW(v))
    ##             testd_permute[[v[i]]] = sample(testd_permute[[v[i]]])
    ##         testd_permute$Method = paste(v, collapse = ", ")
    ##         permute_classes = get_pred(trainglm, testd_permute, vars, type = "class")
    ##         permut_s = asdf(get_pred(trainglm, testd_permute, vars, type = "response"))
    ##         if (NCOL(permut_s) == 1)
    ##             colnames(permut_s) = levels(testd[[ycol]])[2]
    ##         for (i in seq_len(NCOL(permut_s)))
    ##             testd_permute[[names(permut_s)[i]]] = permut_s[[i]]
    ##         permute_roc = make_roc(testd_permute, ycol, colnames(permut_s))[
    ##            ,Method := testd_permute$Method[1]][
    ##            ,fold_id := nm]

    ##         indiv.levels = levels(testd[[ycol]])[-1]
    ##         indiv.accuracy = c()
    ##         if (length(indiv.levels) > 1) {
    ##             for (i in indiv.levels) {
    ##                 indiv.mat = as.matrix(table(refactor(testd[[ycol]], i),
    ##               refactor(factor(permute_classes, levels = levels(testd[[ycol]])), i)))
    ##                 indiv.accuracy = c(indiv.accuracy, setNames(get_accuracy(indiv.mat), i))
    ##             }
    ##         }
    ##         indiv.dec.accuracy = full_mod_indiv_accuracy[names(indiv.accuracy)] - indiv.accuracy
    ##         confus_mat_permute = as.matrix(table(testd[[ycol]], factor(permute_classes, levels = levels(testd[[ycol]]))))
    ##         dec_accuracy = data.table(feature = paste(v, collapse = ", "), decrease_accuracy = full_mod_accuracy - get_accuracy(confus_mat_permute))[, which := "FULL"]
    ##         dec_accuracy = dec_accuracy[rep_len(1:NROW(dec_accuracy), NROW(dec_accuracy) + NROW(indiv.accuracy))]
    ##         dec_accuracy[-1, decrease_accuracy := indiv.dec.accuracy]
    ##         dec_accuracy[-1, which := names(indiv.dec.accuracy)]
    ##         list(dec_accuracy = dec_accuracy,
    ##              roc = permute_roc,
    ##              testd = testd_permute)
    ##     }))

    ##     dec_accuracy = rbindlist(permute_lst$dec_accuracy)
    ##     dec_accuracy$fold_id = nm

    ##     bigroc = rbind(rbindlist(permute_lst$roc),
    ##                    make_roc(testd, ycol, colnames(scores))[, Method := "all"][, fold_id := nm])


    ##     ## mlab = mroclab(testd[[ycol]])
    ##     ## mpred = mrocpred(get_pred(trainglm, testd, vars))

    ##     ## if (class(trainglm)[1] == "lognet")
    ##     ##     mlab = mlab[,2,drop=F]

    ##     outd = rbind(rbindlist(permute_lst$testd), testd)

    ##     ret = list(dec_accuracy = dec_accuracy,
    ##          ## lab = mlab, mpred = mpred,
    ##          roc = bigroc,
    ##          classes = classes,
    ##          scores = scores,
    ##          outd = outd)
    ##     return(ret)
    ##     }, error = function(e) printerr(ixfold))
    ## };  ifun <- feature_importance




































    ## ifun = function(ixfold, xfolds, dat, lambda = 0, vars, seed = 10, ycol, debug = FALSE, rpart = FALSE, permute_varlist = NULL, added_testd = NULL) {
    ##     if (debug) browser()
    ##     nm = names(xfolds[ixfold])
    ##     xfold = xfolds[[ixfold]]
    ##     traind = copy(dat[xfold$train,])
    ##     testd = copy(dat[xfold$test,])

    ##     if (!is.null(added_testd)) {
    ##         added_testd = copy(added_testd[pair %in% setdiff(pair, c(traind$pair, testd$pair))])
    ##         if (NROW(added_testd) > 0)
    ##             testd = rbind(testd[, added := FALSE], added_testd[, added := TRUE])
    ##     }

    ##     set.seed(seed);

    ##     traind = copy3(traind)
    ##     traind$fmut = traind[[ycol]]
    ##     testd = copy3(testd)
    ##     testd$fmut = testd[[ycol]]



    ##     if (length(levels(traind$fmut)) > 2) {
    ##         family = "multinomial"
    ##     } else {
    ##         family = "binomial"
    ##     }

    ##     trainglm = fitcvglmnet(vars, family = family, dat = traind, ycol = ycol)

    ##     classes = get_pred(trainglm, testd, vars, type = "class")
    ##     scores = get_pred(trainglm, testd, vars, type = "response")

    ##     testd$fold_id = nm
    ##     testd$Method = "all"



    ##     testd$classes = classes
    ##     if (!inherits(scores, "data.frame"))
    ##         scores = asdf(scores)

    ##     if (NCOL(scores) == 1)
    ##         colnames(scores) = levels(testd[[ycol]])[2]

    ##     for (i in seq_len(NCOL(scores)))
    ##         testd[[names(scores)[i]]] = scores[[i]]

    ##     off_diag = function(x) {
    ##         diag(x) = NA
    ##         as.vector(x) %>% na.omit
    ##     }

    ##     get_accuracy = function(confus_mat) {
    ##         correct = diag(confus_mat)
    ##         sum(correct) / sum(off_diag(confus_mat), correct)
    ##     }

    ##     confus_mat = as.matrix(table(testd[[ycol]], factor(testd$classes, levels = levels(testd[[ycol]]))))
    ##     full_mod_accuracy = get_accuracy(confus_mat)

    ##     lvars = as.list(vars)

    ##     if (!is.null(permute_varlist)) {
    ##         ## .NotYetImplemented()
    ##         if (is.character(permute_varlist))
    ##             pvar = list(permute_varlist)
    ##         else if (is.list(permute_varlist)) {
    ##             pvar = permute_varlist
    ##         } else if (! is.list(permute_varlist)) {
    ##             warning("ignoring permute_list")
    ##             pvar = list()
    ##         }
    ##         lvars = c(lvars, pvar)
    ##     }

    ##     permute_lst = purrr::transpose(lapply(lvars, function(v, seed = 10) {
    ##         testd_permute = copy(testd)
    ##         testd_permute$fold_id = nm
    ##         set.seed(seed)
    ##         for (i in 1:NROW(v))
    ##             testd_permute[[v[i]]] = sample(testd_permute[[v[i]]])
    ##         testd_permute$Method = paste(v, collapse = ", ")
    ##         permute_classes = get_pred(trainglm, testd_permute, vars, type = "class")
    ##         permut_s = asdf(get_pred(trainglm, testd_permute, vars, type = "response"))
    ##         if (NCOL(permut_s) == 1)
    ##             colnames(permut_s) = levels(testd[[ycol]])[2]
    ##         for (i in seq_len(NCOL(permut_s)))
    ##             testd_permute[[names(permut_s)[i]]] = permut_s[[i]]
    ##         permute_roc = make_roc(testd_permute, ycol, colnames(permut_s))[
    ##            ,Method := testd_permute$Method[1]][
    ##            ,fold_id := nm]
    ##         confus_mat_permute = as.matrix(table(testd[[ycol]], factor(permute_classes, levels = levels(testd[[ycol]]))))
    ##         list(dec_accuracy = data.table(feature = paste(v, collapse = ", "), decrease_accuracy = full_mod_accuracy - get_accuracy(confus_mat_permute)),
    ##              roc = permute_roc,
    ##              testd = testd_permute)
    ##     }))

    ##     dec_accuracy = rbindlist(permute_lst$dec_accuracy)
    ##     dec_accuracy$fold_id = nm

    ##     bigroc = rbind(rbindlist(permute_lst$roc),
    ##                    make_roc(testd, ycol, colnames(scores))[, Method := "all"][, fold_id := nm])


    ##     ## mlab = mroclab(testd[[ycol]])
    ##     ## mpred = mrocpred(get_pred(trainglm, testd, vars))

    ##     ## if (class(trainglm)[1] == "lognet")
    ##     ##     mlab = mlab[,2,drop=F]

    ##     outd = rbind(rbindlist(permute_lst$testd), testd)

    ##     ret = list(dec_accuracy = dec_accuracy,
    ##          ## lab = mlab, mpred = mpred,
    ##          roc = bigroc,
    ##          classes = classes,
    ##          scores = scores,
    ##          outd = outd)
    ##     return(ret)
    ## }
    ### end glm_utils.R
    pp.res = readRDS(paste0(opt$libdir, "/", "robust.pp.res.rds"))

    pred.res = predict.pp(pp.res, tmp)

    ## newdata = data.matrix(expl_variables[, c("tib","ihdel","qrdup","qrdel","RS3","RS5","del.mh.prop","hrd","SNV8","SNV3"),drop=F]) %>% log1p

    ot_scores = fit_rforest(mod, pred.res$newdat)
    ot_scores$SUM12 = ot_scores$BRCA1 + ot_scores$BRCA2


    message("Predicting Oneness Twoness scores")
    ## ot_scores = predict(mod, newdata, type = "response", s = 0, alpha = 1, type.multinomial = "grouped")
    ## ot_scores = reshape2::melt(ot_scores) %>% dcast(. ~ Var2, value.var = "value")
    ## ot_scores$SUM12 = ot_scores$BRCA1 + ot_scores$BRCA2
    ## ot_scores$OTHER = NULL

    outputs = list(expl_variables = expl_variables,
               ot_scores = ot_scores)

    message("Caching Oneness/Twoness results")
    saveRDS(outputs, "onenesstwoness_results.rds")

    quit("no", 0)
}
