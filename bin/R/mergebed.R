library(optparse)
options(error = function() {traceback(2); quit("no", 1)})
options(scipen = 9999)

if (!exists('opt')) ## if opt already exists allow over-ride of command line arg processing for debugging purposes
{
    option_list = list(
        make_option(c("-l", "--libdir"), type = "character", help = "Libdir"),
        make_option(c("--overlap_aa"), type = "character", help = "Overlaps from breakend A <> breakend A"),
        make_option(c("--overlap_bb"), type = "character", help = "Overlaps from breakend B <> breakend B"),
        make_option(c("-o", "--outfile_path"), type = "character", default = './', help = "Path to save the output")
    )
    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
}

print(opt)

writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outfile_path, '/cmd.args', sep = '/'))
saveRDS(opt, paste(opt$outfile_path, 'cmd.opts', sep = '/'))

library(data.table)
setDTthreads(1)

oaa = data.table::fread(
    opt$overlap_aa, 
    header = FALSE
)

obb = data.table::fread(
    opt$overlap_bb, 
    header = FALSE
)

bigmerge = merge.data.table(oaa[, .(V4, V10)], obb[, .(V4, V10)], by = c("V4", "V10"))

fwrite(bigmerge, "./ra_match___AA_BB.tsv", sep = "\t", col.names = FALSE)

quit("no", 0)


