{
    library(optparse)

    if (!exists('opt'))
    {
      option_list = list(
        
        make_option(c("-s", "--somatic_snv"), type = "character", default = NULL, help = "Path to somatic snv file"),
        make_option(c("-g", "--germline_snv"), type = "character", default = NULL, help = "Path to germline snv file"),
        make_option(c("-p", "--het_pileups_wgs"), type = "character", default = NULL, help = "Path to Het SNPs pileups WGS file"),
        make_option(c("-j", "--jabba"), type = "character", help = "Path to jabba file"),
        make_option(c("--filter_pass"), type = "logical", default = FALSE, help = "FILTER == PASS variants?"),
        make_option(c("--tumor_cbs"), type = "character", default = NULL, help = "Path to drycleaned CBS file with seg means contained in field seg.mean"),
        make_option(c("--tumor_dryclean"), type = "character", default = NULL, help = "Path to drycleaned cov file"),
        make_option(c("--dryclean_field"), type = "character", default = "foreground", help = "field that references binned coverage in the dryclean file; usually foreground"),
        make_option(c("--read_size"), type = "integer", default = 151L, help = "expected read size for binned coverage to average base coverage"),
        make_option(c("--tumor_name"), type = "character", default = as.character(NA), help = "Name of tumor sample"),
        make_option(c("--normal_name"), type = "character", default = as.character(NA), help = "Name of normal sample"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "output directory")
      )

      parseobj = OptionParser(option_list=option_list)
      opt = parse_args(parseobj)



      if (is.null(opt$jabba)) {
        stop("Error: The jabba argument is required and must not be NULL.")
      }
      
      if (is.null(opt$tumor_dryclean)) {
        stop("Error: Dryclean argument is required and must not be NULL.")
      }

      if (is.null(opt$somatic_snv) && is.null(opt$germline_snv) && is.null(opt$het_pileups_wgs)) {
        stop("Error: At least one of somatic_snv, germline_snv, or het_pileups_wgs must be provided.")
      }

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

    snvplicity.run <- multiplicity(somatic_snv = opt$somatic_snv,
                                      germline_snv = opt$germline_snv,
                                      het_pileups_wgs = opt$het_pileups_wgs,
                                      jabba_rds = opt$jabba,
                                      tumor_cbs = opt$tumor_cbs,
                                      tumor_dryclean = opt$tumor_dryclean,
                                      dryclean_field = opt$dryclean_field,
                                      read_size = opt$read_size,
                                      tumor_name = opt$tumor_name,
                                      normal_name = opt$normal_name,
                                      filterpass = opt$filter_pass,
                                      tau_in_gamma = F,
                                      verbose = T)

    if(!is.null(snvplicity.run$somatic_variants)){
      saveRDS(snvplicity.run$somatic_variants, paste0(opt$outdir, "est_snv_cn_somatic.rds"))
    }

    if(!is.null(snvplicity.run$germline_variants)){
      saveRDS(snvplicity.run$germline_variants, paste0(opt$outdir, "est_snv_cn_germline.rds"))
    }

    if(!is.null(snvplicity.run$het_pileups)){
      saveRDS(snvplicity.run$het_pileups, paste0(opt$outdir, "est_snv_cn_hets.rds"))
    }

}
