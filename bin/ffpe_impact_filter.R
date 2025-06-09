library(optparse)
options(bitmapType='cairo')

message("R Libraries:")
print(.libPaths())

options(error = function() {traceback(2); quit("no", 1)})
if (!exists('opt'))
{
    option_list = list(
		make_option(c("--id"), type = "character", help = "id of sample"),
        make_option(c("--somatic_mut_vcf"), type = "character", help = "VCF"),
        make_option(c("--signatures_post_prob_sbs"), type = "character", help = "SBS post-probability signatures"),
		make_option(c("--signatures_post_prob_indel"), type = "character", help = "Indel post-probability signatures"),
		make_option(c("--filter_threshold"), type = "numeric", default = 0.5, help = "Indel post-probability signatures"),
		make_option(c("--cores"), type = "numeric", help = "Number of cores"),
        make_option(c("--outdir"), type = "character", default = './', help = "Directory to dump output into")
    )
    parseobj = OptionParser(option_list=option_list)
    # print(parseobj)
    opt = parse_args(parseobj)
    print(opt)
    saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
}

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(data.table)
    library(parallel)  ## needed for mc.cores 
})





sbs_pp = fread(opt$signatures_post_prob_sbs)
indel_pp = fread(opt$signatures_post_prob_indel)

is_sbs_pp_empty = NROW(sbs_pp) == 0
is_indel_pp_empty = NROW(indel_pp) == 0

library(VariantAnnotation)

vcf_path = opt$somatic_mut_vcf
vcf = VariantAnnotation::readVcf(vcf_path)

vcf_header = VariantAnnotation::header(vcf)
vcf_header_info = VariantAnnotation::info(vcf_header)

info_header_to_add = S4Vectors::DataFrame(
  Number = "1",
  Type = "Float",
  Description = "FFPE Impact (SBS57 + SBSFFPE) posterior activity",
  row.names = "FFPEIMPACT"
)

VariantAnnotation::info(vcf_header) = rbind(
  vcf_header_info,
  info_header_to_add
)

VariantAnnotation::header(vcf) = vcf_header


vcf_info_fields = VariantAnnotation::info(vcf)
vcf_info_fields$FFPEIMPACT = NA_real_

rr = MatrixGenerics::rowRanges(vcf)
rr_ckey = paste(seqnames(rr), start(rr))

ref = as.character(rr$REF)
nrows_alt = elementNROWS(rr$ALT)
ix_nrows_multiallelic = which(nrows_alt > 1)
alt = as.character(S4Vectors::unstrsplit(rr$ALT))

ref_char = nchar(ref)
alt_char = nchar(alt)
is_ref_alt_same = ref_char == alt_char
is_ref_char_one = alt_char == 1
is_alt_char_one = alt_char == 1

is_sbs = is_ref_alt_same & is_ref_char_one
is_indel = !is_ref_alt_same
is_mnv = is_ref_alt_same & !is_ref_char_one

is_no_sbs_in_vcf = !any(is_sbs)
is_no_indel_in_vcf = !any(is_indel)


ref_bs = Biostrings::DNAStringSet(ref)
alt_bs = Biostrings::DNAStringSet(alt)

is_collapse = ref_bs %in% c("A", "G") & is_sbs

ref_bs[is_collapse] = Biostrings::reverseComplement(ref_bs[is_collapse])
alt_bs[is_collapse] = Biostrings::reverseComplement(alt_bs[is_collapse])


rr_cra_key = paste(rr_ckey, paste(ref_bs, ">", alt_bs, sep = ""))

is_any_sbs_present = ! (is_sbs_pp_empty || is_no_sbs_in_vcf)

if (is_any_sbs_present) {
	sbs_pp[, reftoalt_collapsed := gsub(".*\\[([AGCT>]+)\\].*", "\\1", MutationType)]
	sbs_pp$seqnames = sbs_pp$Chr
	sbs_pp$start = sbs_pp$Pos
	sbs_pp[, cra_key := paste(seqnames, start, reftoalt_collapsed)]
	mapped_sbs_pp = sbs_pp[match(rr_cra_key, sbs_pp$cra_key)[is_sbs],]
	vcf_info_fields$FFPEIMPACT[is_sbs] = mapped_sbs_pp[, SBSFFPE + SBS57]
}

is_any_indel_present = ! (is_indel_pp_empty || is_no_indel_in_vcf)
if (is_any_indel_present) {

	indel_pp[, ckey := paste(Chr, Pos)]

	mapped_indel_pp = indel_pp[match(rr_ckey, indel_pp$ckey)[is_indel],]

	vcf_info_fields$FFPEIMPACT[is_indel] = mapped_indel_pp$IDFFPE
}

VariantAnnotation::info(vcf) = vcf_info_fields

is_sbs_annotation_missing = any(is.na(vcf_info_fields[is_sbs,]$FFPEIMPACT))
is_indel_annotation_missing = any(is.na(vcf_info_fields[is_indel,]$FFPEIMPACT))
is_either_missing = is_sbs_annotation_missing || is_indel_annotation_missing

msg = NULL

if (is_sbs_annotation_missing) {
  msg = paste(msg, "ERROR: SBS fields have NA.", sep = "\n")
}
if (is_indel_annotation_missing) {
  msg = paste(msg, "ERROR: Indel fields have NA.", sep = "\n")
}
msg = gsub("^\n{1,}|\n{1,}$", "", msg, perl = TRUE)

if (is_either_missing) cat("\n", msg, "\n\n")

is_ffpe_pass = (!is.na(vcf_info_fields$FFPEIMPACT) & vcf_info_fields$FFPEIMPACT <= 0.5)
is_other = is.na(vcf_info_fields$FFPEIMPACT)
is_filtered = is_ffpe_pass | is_other


VariantAnnotation::writeVcf(vcf, paste(opt$id, "_ffpe_annotated.vcf", sep = ""))
VariantAnnotation::writeVcf(vcf[is_filtered], paste(opt$id, "_ffpe_annotated_filtered.vcf", sep = ""))
