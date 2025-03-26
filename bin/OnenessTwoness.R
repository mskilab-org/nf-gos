{
    options(error={expression({traceback(2); quit('no', 1)})})

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
            make_option(c("--outdir"), type = "character", default = './',
                        help = "output directory")
        )
        message("Found input files")

        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)
        message("Parsed input arguments")


        if (
            is.null(opt$complex)
            || (is.null(opt$homeology))
            || is.null(opt$homeology_stats)
            || is.null(opt$hrdetect_results)
            )  # was opt$multinomial
            stop(print_help(parseobj))

        print(opt)

        print(.libPaths())


        message("Establishing record files")
        ## keep record of run
        writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }

    message("Loading Packages")
    suppressWarnings(expr = {
        suppressPackageStartupMessages(expr = {
                library(gGnome)
                library(onenesstwoness)
            })
    })

    message("Done Loading Packages")

    system(paste('mkdir -p',  opt$outdir))

    ########## Oneness Twoness #########

    predict_hrd(complex = opt$complex,
        homeology = opt$homeology,
        homeology_stats = opt$homeology_stats,
        hrdetect_results = opt$hrdetect_results,
        save = TRUE)

    quit("no", 0)
}
