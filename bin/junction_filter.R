withAutoprint(
{
    library(optparse)
    options(bitmapType='cairo')
    options(error = function() {traceback(2); quit("no", 1)})

    if (!exists('opt')) ## if opt already exists allow over-ride of command line arg processing for debugging purposes
    {
        option_list = list(
            make_option(c("-l", "--libdir"), type = "character", help = "Libdir"),
            make_option(c("-p", "--pon"), type = "character", help = "Junction PON file as GRanges .rds file"),
            make_option(c("-i", "--sv"), type = "character", help = "Filtered by PASS tumor only SV vcf from a junction caller"),
            make_option(c("-a", "--padding"), type = "numeric", help = "Padding to provide so that junction that fall within the amount of padding of PON will get removed"),
            make_option(c("-g", "--gnomAD"), type = "character", help = "gnomAD .rds file in breakend grl format."),
            make_option(c("-o", "--outfile_path"), type = "character", default = './', help = "Path to save the output")
        )
        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)
    }
    print(opt)
    writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outfile_path, '/cmd.args', sep = '/'))
    saveRDS(opt, paste(opt$outfile_path, 'cmd.opts', sep = '/'))

    setDTthreads(1)
    library(skitools)
    library(gUtils)
    library(gTrack)
    library(gGnome)

    out.file.somatic.sv            = paste(opt$outfile_path, 'somatic.filtered.sv.rds', sep = '/')
    out.file.somatic.sv.gnoMAD     = paste(opt$outfile_path, 'somatic.filtered.gnoMAD.sv.rds', sep = '/')

    if(file.exists(opt$pon)) {
            ## PON filtering
            message(paste0("Reading in the provided Junction PON from ", opt$pon))
            pon_object                     = readRDS(opt$pon)
            this_path                      = opt$sv
            this_filt                      = gGnome:::read.juncs(this_path, chr.convert=FALSE)
            message("read in the filtered SV file from:", this_path)
            within_pon                     = suppressWarnings(suppressMessages(ra.overlaps(this_filt, pon_object, pad = opt$padding)))
            filter_these                   = unique(within_pon[,"ra1.ix"])
            return_this_filtered_somatic   = this_filt[-filter_these]
            message("Finished filtering of the provided SV vcf with Junction PON, now will filter by gnoMAD...")
            rm(pon_object, filter_these, this_filt)
            ## gnoMAD filtering
            gnoMAD_object                  = readRDS(opt$gnoMAD)
            message("Loaded the gnoMAD object, filtering the provided VCF...")
            within_gnoMAD                  = suppressWarnings(suppressMessages(ra.overlaps(return_this_filtered_somatic, gnoMAD_object, pad = opt$padding)))
            filter_these_gnoMAD            = unique(within_gnoMAD[,"ra1.ix"])
            return_this_filt               = return_this_filtered_somatic[-filter_these_gnoMAD]
            ## Saving output
            message("Saving the results in .rds format. You can plug the rds file in JaBbA as input. :)")
            saveRDS(return_this_filtered_somatic, out.file.somatic.sv)
            saveRDS(return_this_filt, out.file.somatic.sv.gnoMAD)
    } else {
            stop(print_help(parseobj))
    }

    cat('Done filtering SV calls from PON and gnoMAD!\n')

    quit("no", 0)
}, echo = FALSE)
