{
    options(error={expression({traceback(2); quit('no', 1)})})

    library(optparse)
    message("Loading input arguments")

    if (!exists('opt'))
    {
        option_list = list(
            make_option(c("--complex"), type = "character", help = "output of Events.task - path to .rds gGraph object annotated with complex SV event classes"),
            make_option(c("--hrdetect_results"), type = "character", help = "output of hrdetect_og.task - preferably run with Strelka2 snv and SvAbA indels"),
            make_option(c("--genome"), type = "character", help = "path to the reference genome fasta file used for alignment (e.g. hg19, hg38, mm10)"),
            make_option(c("--model"), type = "character", help = "path to the model weights for prediction"),
            make_option(c("--cores"), type = "integer", default = 1, help = "number of cores to use for parallel processing [default 1]"),
            make_option(c("--outdir"), type = "character", default = './', help = "output directory")
        )
        message("Found input files")

        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)
        message("Parsed input arguments")


        if (
            is.null(opt$complex)
            || is.null(opt$hrdetect_results)
            || is.null(opt$genome)
            || is.null(opt$model)
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
                library(GxG)
            })
    })

    message("Done Loading Packages")

    system(paste('mkdir -p',  opt$outdir))

    ########## Oneness Twoness #########

    predict_B1_2(
        complex = opt$complex,
        hrdetect_results = opt$hrdetect_results,
        genome = opt$genome,
        model = opt$model,
        cores = opt$cores
    )

    quit("no", 0)
}
