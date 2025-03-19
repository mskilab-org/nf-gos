{
    library(optparse)

    if (!exists('opt'))
    {
      option_list = list(

        make_option(c("-s", "--somatic_snv"), type = "character", help = "Path to somatic snv file"),
        make_option(c("-g", "--germline_snv"), type = "character", help = "Path to germline snv file"),
        make_option(c("-j", "--jabba"), type = "character", help = "Path to jabba file"),
        make_option(c("--snpeff_path"), type = "character", help = "Path to SnpEff directory containing scripts and modules"),
        make_option(c("--tumor_name"), type = "character", default = as.character(NA), help = "Name of tumor sample"),
        make_option(c("--normal_name"), type = "character", default = as.character(NA), help = "Name of normal sample"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "output directory"),
        make_option(c("-l", "--libdir"), type = "character", default = "~/modules/snv_multiplicity2", help = "Directory containing this R file")
      )

      parseobj = OptionParser(option_list=option_list)
      opt = parse_args(parseobj)

      if ((is.null(opt$somatic_snv) && is.null(opt$germline_snv)) | is.null(opt$jabba))
      	stop(print_help(parseobj))

      print(opt)

      print(.libPaths())
      options(error=function() {traceback(2); quit("no", 1)})

        ## keep record of run
      writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
      saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }

    message("Libraries loading...")
    suppressPackageStartupMessages({
      suppressWarnings({
            library(gGnome)
            library(gUtils)
            library(skitools)
            library(bamUtils)
            library(dplyr)
            library(rtracklayer)
            library(khtools)
            library(skidb)
            library(multiplicity)
        })
    })

    if(opt$somatic_snv == "/dev/null")
      opt$somatic_snv <- NULL

    if(opt$germline_snv == "/dev/null")
      opt$germline_snv <- NULL

    snvplicity.run <- snvplicity(somatic_snv = opt$somatic_snv,
                    germline_snv = opt$germline_snv,
                    jabba_rds = opt$jabba,
                    tumor_name = opt$tumor_name,
                    normal_name = opt$normal_name,
                    snpeff_path = opt$snpeff_path,
                    verbose = T)

    if(!is.null(snvplicity.run[[1]])){
      saveRDS(snvplicity.run[[1]], paste0(opt$outdir, "est_snv_cn_somatic.rds"))
    }

    if(!is.null(snvplicity.run[[2]])){
      saveRDS(snvplicity.run[[2]], paste0(opt$outdir, "est_snv_cn_germline.rds"))
    }

}
