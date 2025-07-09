library(optparse)

options(error = function() {traceback(2); quit("no", 1)})
if (!exists('opt'))
{
    option_list = list(
      make_option(c("--vcf"), type = "character", default = "./", help = "path to small muations VCF"),
      make_option(c("--heme_db"), type = "character", default = "/gpfs/data/imielinskilab/projects/Clinical_NYU/db/master_heme_database_anno.rds", help = "path to Heme DB master list"),
      make_option(c("--genome"), type = "character", default = "hg19", help = "genome build"),
      make_option(c("--samp_id"), type = "character", default = "X", help = "Sample ID"),
      make_option(c("--outdir"), type = "character", default = './', help = "output dir")
    )
    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
    saveRDS(opt, paste(opt$outdir, 'rescue.rds', sep = '/'))
}

suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
})



message("Reading unfiltered VCF")
unfiltered_vcf = readVcf(opt$vcf, genome = opt$genome)
heme_db = readRDS(opt$heme_db)
rescue_vcf = unfiltered_vcf[rownames(unfiltered_vcf) %in% unique(heme_db$KEY_M)]
writeVcf(rescue_vcf, paste0(opt$outdir, "/", opt$samp_id, ".sage.pass_rescued.vcf"))
message("Rescued vcf saved")
