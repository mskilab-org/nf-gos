library(optparse)
options(scipen = 9999)
options(error = function() {traceback(2); quit("no", 1)})

if (!exists('opt')) ## if opt already exists allow over-ride of command line arg processing for debugging purposes
{
    option_list = list(
        make_option(c("-l", "--libdir"), type = "character", help = "Libdir"),
        make_option(c("--sv"), type = "character", help = "SV VCF"),
        make_option(c("--padding"), type = "numeric", help = "Padding to provide so that junction that fall within the amount of padding of PON will get removed"),
        make_option(c("-o", "--outfile_path"), type = "character", default = './', help = "Path to save the output"),
        make_option(c("--cores"), type = "integer", default = 1L, help = "Number of cores")

    )
    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
}

print(opt)

writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outfile_path, '/cmd.args', sep = '/'))
saveRDS(opt, paste(opt$outfile_path, 'cmd.opts', sep = '/'))


library(GenomicRanges)
library(gUtils)
library(gGnome)

setDTthreads(1)

out.file.somatic.sv            = paste(opt$outfile_path, 'somatic.filtered.sv.rds', sep = '/')

jun = gGnome:::read.juncs(opt$sv)
GenomeInfoDb::seqlevelsStyle(jun) = "NCBI"
jun = GenomeInfoDb::sortSeqlevels(jun)

grlunl = as(
(
    setDT(as.data.frame(jun))
    [, c("group_n", "group_iix") := list(.N, seq_len(.N)), by = group]
    [order(group, seqnames, start, end)]
    []
),
    "GRanges"
)
mcunl = mcols(grlunl)
mcols(grlunl) = DataFrame(seq_len(NROW(grlunl)))[,0]
jun_sorted = GenomicRanges::split(grlunl + opt$padding, mcunl$group)

mcols(jun_sorted) = mcols(jun)

pivjun = gUtils::grl.pivot(jun_sorted)
ra_a = pivjun[[1]]
ra_b = pivjun[[2]]

ra_a = setDT(as.data.frame(ra_a))[]
ra_b = setDT(as.data.frame(ra_b))[]

nr = NROW(ra_a)
numdigits = function(x) floor(log10(x)) + 1
spad = numdigits(nr)
## ix_padded = stringr::str_pad(seq_len(nr), width = spad, pad = "0")
ix_padded2 = sprintf(glue::glue("%0{spad}d"), seq_len(nr))

ra_sv_name = paste(
    ix_padded2, "___",
    ra_a$seqnames, "_",
    ra_a$start, "_",
    ra_a$end, "_",
    ra_a$strand, "___",
    ra_b$seqnames, "_",
    ra_b$start, "_",
    ra_b$end, "_",
    ra_b$strand,
    sep = ""
)

ra_a$width = NULL
ra_b$width = NULL

names(ra_a) = c("#chrom", "chromStart", "chromEnd", "strand")
names(ra_b) = c("#chrom", "chromStart", "chromEnd", "strand")

ra_a$chromstart = ra_a$chromStart - 1
ra_b$chromstart = ra_b$chromStart - 1

ra_a$score = 0L
ra_b$score = 0L

ra_a$name = ra_sv_name
ra_b$name = ra_sv_name

ra_a$chromStart = pmax(ra_a$chromStart, 0)
ra_b$chromStart = pmax(ra_b$chromStart, 0)

fwrite(
    base::subset(ra_a, select = c("#chrom", "chromStart", "chromEnd", "name", "score", "strand")),
    "ra_A.bed",
    sep = "\t",
    nThread = opt$cores
)


fwrite(
    base::subset(ra_b, select = c("#chrom", "chromStart", "chromEnd", "name", "score", "strand")),
    "ra_B.bed",
    sep = "\t",
    nThread = opt$cores
)


quit("no", 0)


