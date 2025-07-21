library(optparse)
options(bitmapType='cairo')

message("R Libraries:")
print(.libPaths())

options(error = function() {traceback(2); quit("no", 1)})
if (!exists('opt'))
{
    option_list = list(
        make_option(c("-l", "--libdir"), type = "character", help = "libdir"),
        make_option(c("-i", "--id"), type = "character", help = "sample id"),
        make_option(c("-j", "--junctions"), type = "character", help = "an RDS file contains a gGraph or JaBbA graph with cn annotation on nodes and edges, or is a BND formatted SV VCF"),
        make_option(c("-r", "--gencode"), type = "character", default = "~/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds", help = "an RDS or GTF file of GENCODE"),
        make_option(c("--genes"), type = "character", default = "/dev/null", help = "Gene list to filter junction breakends"),
        make_option(c("--do_only_canonical"), type = "logical", default = FALSE, help = "Use only canonical (longest per gene), protein-coding transcripts"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into"),
        make_option(c("--cores"), type = "integer", default = 1L, help = "Number of cores")
    )
    parseobj = OptionParser(option_list=option_list)
    # print(parseobj)
    opt = parse_args(parseobj)
    print(opt)
    saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
}

suppressPackageStartupMessages({
    library(gGnome)
    library(gUtils)
    library(parallel)  ## needed for mc.cores 
    library(rtracklayer)
    library(GenomeInfoDb)
})


message("Loading GENCODE: ", opt$gencode)
if (grepl('.rds$', opt$gencode))
{
    gencode = readRDS(as.character(opt$gencode))
} else {
    gencode = rtracklayer::import(opt$gencode)
    GenomeInfoDb::seqlevelsStyle(gencode) = "NCBI" # remove chr and deal with chrM properly
}

if (opt$do_only_canonical) {
    message("Prefiltering to only 1 canonical transcript per gene is deprecated. gGnome::fusions() internally selects the longest transcript with a break inside")
    # message("Using only canonical transcripts")
    # canonical_transcripts = (
    #     gUtils::gr2dt(
    #         gencode[
    #             gencode$type == "transcript"
    #             & gencode$gene_type == "protein_coding"
    #             & gencode$transcript_type == "protein_coding"
    #         ]
    #     )[, .(transcript_id = transcript_id[which.max(end-start+1)]), by = gene_name]
    # )
    # gencode = gencode[
    #     gencode$transcript_id %in% canonical_transcripts$transcript_id 
    #     | (gencode$type == "gene" & gencode$gene_name %in% canonical_transcripts$gene_name)
    # ]
}

message("Loading SV input: ", opt$junctions)
if (grepl(".rds$", opt$junctions)) {
    gg_input = readRDS(opt$junctions)
    if (!identical(c("gGraph", "R6"), class(gg_input)) && is.list(gg_input) && any(c("segstats", "purity", "ploidy") %in% names(gg_input)))
        gg_input = gGnome::gG(jabba = gg_input) # is jabba list object, convert
    message("Refreshing gGraph input")
    gg_input = gGnome::refresh(gg_input)
} else if (grepl(".vcf(.gz|.bz|.bz2|.xz)?$", opt$junctions)) {
    gg_input = gGnome::gG(junctions = opt$junctions)
} else {
    stop("Unknown input file type: ", opt$junctions)
}

if (!identical(opt$genes, "/dev/null") && file.exists(opt$genes)) {
    message("Loading gene list to filter input SVs:")
    cat(normalizePath(opt$genes))
    gene_list = readLines(opt$genes)
    gencode_by_genes = split(gencode, gencode$gene_name)
    genes_not_present = setdiff(gene_list, names(gencode_by_genes))
    if (length(genes_not_present)) {
        message("The following genes are not present in gencode definition")
        dput(genes_not_present)
    }
    genes_present = intersect(gene_list, names(gencode_by_genes))
    message("Number of unique genes to query in gencode: ", length(genes_present))
    message("\n")
    message("Creating new gGraph from junctions overlapping with gene set")
    alt_juncs = gg_input$edges$junctions[type == "ALT"]
    which_juncs_in_genes = which(alt_juncs %^% unlist(gencode_by_genes[genes_present]))
    message(length(which_juncs_in_genes), " ALT junctions in gene set")
    gg_input = gGnome::gG(junctions=alt_juncs[which_juncs_in_genes])
    rm(gencode_by_genes)
}

fus_lst = fusions(gg_input, gencode, verbose = TRUE, mc.cores = opt$cores, return.all.alt.edges = TRUE)
fus = fus_lst$fus
allaltedges.ann = fus_lst$allaltedges.ann

if (length(fus)) {
    fus$set(id = opt$id)
    is_cn_attribute_present = NROW(fus$edges$dt$cn)
    if (is_cn_attribute_present) fus$set(mincn = fus$eval(edge = min(cn, na.rm = TRUE)))
}

saveRDS(fus, paste0(opt$outdir, '/', 'fusions.rds'))
data.table::fwrite(allaltedges.ann, paste0(opt$outdir, '/', 'altedge.annotations.tsv'), sep = "\t")

quit("no", 0)
