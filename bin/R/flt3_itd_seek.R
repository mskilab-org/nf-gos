#!/usr/bin/env Rscript

library(optparse)

options(error = function() {traceback(2); quit("no", 1)})
if (!exists('opt'))
{
    option_list = list(
      make_option(c("--bam"), type = "character", default = "./", help = "path to the bam file"),
      make_option(c("--fasta"), type = "character", default = "/gpfs/data/imielinskilab/projects/Clinical_NYU/addy_debug/flt3_itd/hg19_nochr.fasta", help = "path to ref fasta"),
      make_option(c("--samtools"), type = "character", default = "samtools",help = "path to samtools"),
      make_option(c("--build"), type = "character", default = "hg19", help = "Ref build"),
      make_option(c("--samp_id"), type = "character", default = "X", help = "Sample ID"),
      make_option(c("--outdir"), type = "character", default = './', help = "output dir")
    )
    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
    saveRDS(opt, paste(opt$outdir, 'flt3_itd_seek.rds', sep = '/'))
}

suppressPackageStartupMessages({
    library(S4Vectors)
    library(data.table)
})



## Setting bam and making sure the index file is in the same dir
this_bam = opt$bam
this_bam_dir = dirname(this_bam)

if (!file.exists(gsub("bam$", "bam\\.bai", this_bam))){
  system(paste0("cp ", this_bam_dir, "fq2bam/*.bam.bai ", this_bam_dir, "."))
  message("moved bai")
}

## ITDseek invocation
system(paste0("/root/git/itdseek/itdseek.sh ", this_bam, " ", opt$fasta, " ", opt$samtools, " ", opt$build, " > ", opt$outdir, "/", opt$samp_id, "_flt3_itd.vcf"))

message("Ran ITDseek")

this_itd = gGnome:::read_vcf(paste0(opt$outdir, "/", opt$samp_id, "_flt3_itd.vcf"))
this_itd$REF = as.character(this_itd$REF)
this_itd$ALT = as.character(unlist(this_itd$ALT))
this_itd = gUtils::gr2dt(this_itd)

if (nrow(this_itd[QUAL >= 2]) > 0){
  ITD_status = "Positive"
} else {
  ITD_status = "Negative"
}

saveRDS(ITD_status,  paste0(opt$outdir, "/", opt$samp_id, "_flt3_itd_status.rds"))

quit("no", 0)

