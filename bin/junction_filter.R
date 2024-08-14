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

    library(data.table)
    library(skitools)
    library(gUtils)
    library(gTrack)
    library(gGnome)

    setDTthreads(1)
    out.file.somatic.sv            = paste(opt$outfile_path, 'somatic.filtered.sv.rds', sep = '/')
    out.file.somatic.sv.gnomAD     = paste(opt$outfile_path, 'somatic.filtered.gnomAD.sv.rds', sep = '/')

    if(file.exists(opt$pon)) {
            ## PON filtering
            message(paste0("Reading in the provided Junction PON from ", opt$pon))
            pon_object                     = readRDS(opt$pon)
            this_path                      = opt$sv
            message("Reading in the filtered SV file from:", this_path)
            this_filt                      = gGnome:::read.juncs(this_path, verbose=TRUE)
            message("Checking overlaps between the provided SV vcf and Junction PON...")
            within_pon                     = suppressWarnings(suppressMessages(ra.overlaps(this_filt, pon_object, pad = opt$padding)))
            message("Unique SVs that overlap with the PON:")
            filter_these                   = unique(within_pon[,"ra1.ix"])
            message("NA values in the filtered SVs:")
            filter_these = filter_these[!is.na(filter_these)]  # Remove NAs
            message("Filtering the provided SV vcf...")
            return_this_filtered_somatic   = this_filt[-filter_these]
            rm(pon_object, filter_these, this_filt)
            ## gnomAD filtering
            gnomAD_object                  = readRDS(opt$gnomAD)
            message("Loaded the gnomAD object, filtering the provided VCF...")
            within_gnomAD                  = suppressWarnings(suppressMessages(ra.overlaps(return_this_filtered_somatic, gnomAD_object, pad = opt$padding)))
            filter_these_gnomAD            = unique(within_gnomAD[,"ra1.ix"])
            filter_these_gnomAD = filter_these_gnomAD[!is.na(filter_these_gnomAD)]  # Remove NAs
            return_this_filt               = return_this_filtered_somatic[-filter_these_gnomAD]
            ## Saving output
            message("Saving the results in .rds format. You can plug the rds file in JaBbA as input. :)")
            saveRDS(return_this_filtered_somatic, out.file.somatic.sv)
            saveRDS(return_this_filt, out.file.somatic.sv.gnomAD)
    } else {
            stop(print_help(parseobj))
    }

    cat('Done filtering SV calls from PON and gnomAD!\n')

    quit("no", 0)
}, echo = FALSE)
