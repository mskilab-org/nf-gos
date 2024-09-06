library(optparse)
options(bitmapType='cairo')

options(error = function() {traceback(2); quit("no", 1)})
if (!exists('opt'))
{
    option_list = list(
        make_option(c("-l", "--libdir"), type = "character", help = "libdir"),
        make_option(c("-i", "--id"), type = "character", help = "sample id"),
        make_option(c("-g", "--gGraph"), type = "character", help = "an RDS file contains a gGraph or JaBbA graph with cn annotation on nodes and edges"),
        make_option(c("-r", "--gencode"), type = "character", default = "~/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds", help = "an RDS or GTF file of GENCODE"),
        make_option(c("--do_only_canonical"), type = "logical", default = TRUE, help = "Use only canonical (longest per gene), protein-coding transcripts"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into"),
        make_option(c("--cores"), type = "integer", default = 1L, help = "Number of cores")
    )
    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
    saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
}

suppressPackageStartupMessages({
    library(gGnome)
    library(gUtils)
    library(parallel)  ## needed for mc.cores
    library(rtracklayer)
    library(GenomeInfoDb)
})


## setDTthreads(10)
if (grepl('.rds$', opt$gencode))
{
    gencode = readRDS(as.character(opt$gencode))
} else {
    gencode = rtracklayer::import(opt$gencode)
    GenomeInfoDb::seqlevelsStyle(gencode) = "NCBI" # remove chr and deal with chrM properly
}

if (opt$do_only_canonical) {
    message("Using only canonical transcripts")
    canonical_transcripts = (
        gUtils::gr2dt(
            gencode[
                gencode$type == "transcript"
                & gencode$gene_type == "protein_coding"
                & gencode$transcript_type == "protein_coding"
            ]
        )[, .(transcript_id = transcript_id[which.max(end-start+1)]), by = gene_name]
    )
    gencode = gencode[
        gencode$transcript_id %in% canonical_transcripts$transcript_id
        | (gencode$type == "gene" & gencode$gene_name %in% canonical_transcripts$gene_name)
    ]
}


if (grepl(".rds$", opt$gGraph)) {
    gg_input = readRDS(opt$gGraph)
    if (!identical(c("gGraph", "R6"), class(gg_input)) && is.list(gg_input) && any(c("segstats", "purity", "ploidy") %in% names(gg_input)))
        gg_input = gGnome::gG(jabba = gg_input) # is jabba list object, convert
} else if (grepl(".vcf(.gz|.bz|.bz2|.xz)?$", opt$gGraph)) {
    gg_input = gGnome::gG(junctions = opt$gGraph)
} else {
    stop("Unknown input file type: ", opt$gGraph)
}
fus = fusions(gg_input, gencode, verbose = TRUE, mc.cores = opt$cores)

## update events with sample id
if (length(fus))
{
    fus$set(id = opt$id)
    fus$set(mincn = fus$eval(edge = min(cn, na.rm = TRUE)))
}

saveRDS(fus, paste0(opt$outdir, '/', 'fusions.rds'))

quit("no", 0)
