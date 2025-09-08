library(optparse)
options(error = function() {traceback(2); quit("no", 1)})
options(scipen = 9999)

if (!exists('opt')) ## if opt already exists allow over-ride of command line arg processing for debugging purposes
{
    option_list = list(
        make_option(c("--ramatch_file"), type = "character", default = "./ra_match___AA_BB.tsv", help = "Path to file containing RA matches"),
        make_option(c("--input_vcf"), type = "character", help = "Path to SV vcf file"),
        make_option(c("-o", "--outfile_path"), type = "character", default = './', help = "Path to save the output")
    )
    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
}

print(opt)

library(gGnome)


writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outfile_path, '/cmd.args', sep = '/'))
saveRDS(opt, paste(opt$outfile_path, 'cmd.opts', sep = '/'))

jun = gGnome:::read.juncs(opt$input_vcf)


ramatch = fread(opt$ramatch_file, header = FALSE)
junmatch = unique(ramatch[[1]])
nr_junmatch = NROW(junmatch)
any_junmatch = nr_junmatch > 0
ix_in_germline = integer(0)
if (any_junmatch) {
    ix_in_germline = as.integer(stringr::str_split_fixed(junmatch, "___", 3)[,1])
}
# ix_in_germline = unique(as.integer(ramatch$V1))

nr = NROW(ix_in_germline)
any_ix_in_germline = nr > 0
jun_filt = jun
if (any_ix_in_germline) {
    jun_filt = jun[-ix_in_germline]
}

mcols(jun_filt)$is_in_gnomAD_sv_pon = rep_len(FALSE, NROW(jun_filt))

saveRDS(jun_filt, "./somatic.filtered.sv.no.gnomAD.rds")


