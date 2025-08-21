library(optparse)
withAutoprint({
message("Loading input arguments")
if (!exists('opt'))
{
    option_list = list(
        make_option(c("--snv"), type = "character", help = "Somatic snv vcf"),
        make_option(c("--indel"), type = "character", help = "indel vcf, if not present in snv"),
        make_option(c("--jabba"), type = "character",
                    help = "jabba rds file with $asegstats"),
        make_option(c("--sv"), type = "character", help = "SV input (vcf/rds)"),
        make_option(c("--mask"), type = "character", help = "path to universal mask"),
        make_option(c("--outdir"), type = "character", default = './',
                    help = "output directory"),
        make_option(c("--hets"), type = "character",
                    help = "path to het pileups input"),
        make_option(c("--genome"), type = "character", default = "hg19",
                    help = "string that should be one of 'hg19', 'hg38', 'mm10', 'canFam3'"),
        make_option(c("--ref"), type = "character",
                    help = "character indicating reference location")
    )
    message("Found input files")

    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)
    message("Parsed input arguments")

    if ((is.null(opt$snv)) | is.null(opt$jabba))
        stop(print_help(parseobj))

    print(opt)

    print(.libPaths())
    options(error={expression({traceback(2); quit('no', 1)})})

    message("Establishing record files")
    ## keep record of run
    writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
    saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
}

message("Loading Packages")
suppressWarnings(expr = {
    suppressMessages(expr = {
        suppressPackageStartupMessages(expr = {
            library(plyr)
            library(dplyr)
            library(data.table)
            library(gGnome)
            library(gUtils)
            library(tidyverse)
            library(devtools)
            library(GenomeInfoDb)
            library(skitools)
            library(skidb)
            library(khtools)
            library(tools)
            ## library(StructuralVariantAnnotation)
            library(signature.tools.lib) ## HRDetect package
            library(ffTrack)
        })
    })
})

### kat added because usually this function doesn't recognize the .bgz format of the indels below:
bcfindex <- function (vcf, force = TRUE)
{
    if (!grepl(".[bv]cf(.gz)?$", vcf) & !grepl(".[bv]cf(.bgz)?$", vcf)) {
        stop("check if you have a valid bcf/vcf file")
    }
    if (!file.exists(paste0(vcf, ".tbi")) & !file.exists(paste0(vcf,
        ".csi")) || isTRUE(force)) {
        system(sprintf("bcftools index --tbi %s", vcf))
    }
    vcf
}


## khtools::file.exists2(x) == !khtools::file.not.exists(x)

message("Done Loading Packages")

setDTthreads(1)

system(paste('mkdir -p',  opt$outdir))

#' xtYao #' Thursday, May 06, 2021 01:44:08 PM
## the reference genome should be consistent with the input VCF
v = VariantAnnotation::scanVcfHeader(opt$snv)
sl = seqlengths(v)
# Sys.setenv(DEFAULT_GENOME = opt$ref)
# sl = hg_seqlengths()

## wg = si2gr(sl) ## modified line

wg = si2gr(sl[which(!grepl("GL|MT|M|NC", names(sl)))]) ## modified line

wg = keepStandardChromosomes(wg, pruning.mode = "coarse")
wg = sort(sortSeqlevels(wg), ignore.strand = T)

gr.mask = withv(opt$mask,  {
    if (!is.null(x) && file.exists(x)) {
        message("mask provided")
        if (grepl("bed(.gz)?$", x))
            import(x)
        else if (grepl(".rds$", x))
            readRDS(x)
        else {
            message("mask file in incorrect format")
            message("no region filtering will be done")
            wg[0]
        }
    } else {
        wg[0]
    }
})

if (!inherits(gr.mask, "GRanges")) {
    warning("mask specified incorrectly\nproceeding without filtering")
    gr.mask = wg[0]
}

if (!identical(gr.mask, GRanges())) {
    message("parsing mask file, and saving to regions file for subsetting VCF\n")
} else {
    message("No mask regions provided, using whole SNV / InDel VCFs")
}

good.wg = gr.setdiff(wg, gr.mask) %>% keepStandardChromosomes(pruning.mode = "coarse")

regions.file = "./good_rfile.txt"
regions.bed = "./good_rfile.bed"
mask.file = "./mask_rfile.txt"

fwrite(gr2dt(good.wg)[, 1:3, with = FALSE], regions.file, sep = "\t", col.names = FALSE)
## fwrite(gr2dt(gr.mask)[, 1:3, with = FALSE], mask.file, sep = "\t", col.names = FALSE)

## fwrite(setnames(gr2dt(gr.chr(good.wg))[, 1:3, with = FALSE][,start := start - 1],
##                 c("#chrom", "chromStart", "chromEnd")), regions.bed, sep = "\t", col.names = TRUE)
#' xtYao #' Thursday, May 06, 2021 01:19:45 PM
## trying to see if it breaks HRDetect if we don't do anything to the chr prefix
fwrite(setnames(gr2dt(good.wg)[, 1:3, with = FALSE][,start := start - 1],
                c("#chrom", "chromStart", "chromEnd")), regions.bed, sep = "\t", col.names = TRUE)



########## processing SNV / indel
message("\n\n")
message("processing SNV / indel input")
snv = opt$snv
indel = opt$indel

## snv_filetype = ifelse(grepl("vcf(.gz){0,255}$", snv),
##                       "vcf",
##                ifelse(grepl(".rds$|.RDS$", snv),
##                       "rds",
##                ifelse(grepl(".csv$|.txt$|.tsv$", snv),
##                       "txt", "")))


if (file.exists(snv) && file.info(snv)$size) {
    snv.tmp = "./snv.vcf.gz"
    ## if (file.exists(snv.tmp)) system2("rm", c(paste0(snv.tmp, "*")))

    system2("rm", c(paste0(snv.tmp, "*")), stderr = "/dev/null")

    writeLines(system(sprintf("bcftools query -l %s", snv), intern = T),
               "./excls.txt")

    ## v_query = "-v snps"
    ## cmd = sprintf("(bcftools view %s -S ./excls.txt -i 'FILTER==\"PASS\" | FILTER==\".\"' -T %s %s| bcftools norm -Oz -m-any) > %s", v_query, regions.file, bcfindex(snv), snv.tmp)



    ## PCAWG vcf are malformed
    ## had to use vcftools initially due to strange errors with PCAWG samples... and then
    ## removed sample level information which was breaking reading in the vcf
    ## this works generally now for different vcf types
    ## don't use --gzvcf for non-gzipped vcfs
    message("Processing SNV VCF")
    if (grepl("gz$", snv))
    {
        cmd = sprintf("{ vcftools --gzvcf %s --remove-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/){ gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b %s -header | vcf-sort -c | bcftools view -S ^./excls.txt | bcftools view -v snps | bcftools norm -Oz -m-any; } > %s", snv, normalizePath(regions.bed), snv.tmp)
    }
    else
    {
        cmd = sprintf("{ vcftools --vcf %s --remove-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/) { gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b %s -header | vcf-sort -c | bcftools view -v snps | bcftools norm -Oz -m-any; } > %s", snv, normalizePath(regions.bed), snv.tmp)
    }
    print(cmd)

    out = system(cmd)
    if (out != 0)  {
        message("trying snv processing again")
        if (grepl("gz$", snv))
        {
        cmd = sprintf("(vcftools --gzvcf %s --remove-indels --remove-filtered-all --recode --stdout | vcf-sort -c | bcftools view -S ^./excls.txt | bcftools norm -Oz -m-any) > %s", snv, snv.tmp)
        }
        else
        {
            cmd = sprintf("(vcftools --vcf %s --remove-indels --remove-filtered-all --recode --stdout | vcf-sort -c | bcftools view -S ^./excls.txt | bcftools norm -Oz -m-any) > %s", snv, snv.tmp)
        }
        out = system(cmd)
        if (out != 0) warning("could not be intersected with non-mask regions bed")
    }
    bcfindex(snv.tmp)
    names(snv.tmp) = "sample_1"
} else {
    snv.tmp = NULL
}

if (is.null(indel) || !file.exists(indel) || identical(indel, "/dev/null")) {
    message("no proper indel file provided... assuming indels are in snv vcf file")
    indel = snv
}


## indel_filetype = ifelse(grepl("vcf(.gz){0,255}$", indel),
##                         "vcf",
##                  ifelse(grepl(".rds$|.RDS$", indel),
##                         "rds",
##                  ifelse(grepl(".csv$|.txt$|.tsv$", indel),
##                         "txt",
##                         "")))


## if (file.exists(indel) && identical(indel, "/dev/null")) {
## check that file is non-empty
if (!is.null(indel) && file.exists(indel) && file.info(indel)$size) {
    indel.tmp = "./indel.vcf.bgz"
    system2("rm", c(paste0(indel.tmp, "*")), stderr = "/dev/null")

    writeLines(system(sprintf("bcftools query -l %s", indel), intern = T),
               "./excls.txt")

    ## v_query = "-v indels"
    ## cmd = sprintf("(bcftools view %s -i 'FILTER==\"PASS\" | FILTER==\".\"' -T %s %s | bcftools norm -Oz -m-any) > %s", v_query, regions.file, bcfindex(indel), indel.tmp)

    ## PCAWG vcf are malformed
    ## had to use vcftools initially due to strange errors with PCAWG samples... and then
    ## removed sample level information which was breaking reading in the vcf
    ## this works generally now for different vcf types
    message("Processing indel input")
    if (grepl("gz$", indel))
    {
        ## cmd = sprintf("{ vcftools --gzvcf %s --keep-only-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/) { gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b ./good_rfile.bed -header | vcf-sort -c | bcftools view -S ^./excls.txt | bcftools view -v indels | bcftools norm -Ov -m-any | bcftools norm -f human_g1k_v37_decoy.fasta --check-ref s | bgzip -c; } > ./indel.vcf.bgz", indel)
        cmd = sprintf("{ vcftools --gzvcf %s --keep-only-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/) { gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b %s -header | vcf-sort -c | bcftools view -v indels | bcftools norm -Ov -m-any | bcftools norm -f %s --check-ref s | bgzip -c; } > %s", indel, regions.bed, opt$ref, indel.tmp)
    }
    else
    {
        cmd = sprintf("{ vcftools --vcf %s --keep-only-indels --remove-filtered-all --recode --stdout | awk '{ if ($0 ~ /##contig=<ID=chr/) { gsub(/##contig=<ID=chr/, \"##contig=<ID=\"); print } else if ($0 !~ /^#/) { gsub(/^chr/, \"\"); print } else { print } }' | bedtools intersect -a stdin -b %s -header | vcf-sort -c | bcftools view -v indels | bcftools norm -Ov -m-any | bcftools norm -f %s --check-ref s | bgzip -c; } > %s", indel, regions.bed, opt$ref, indel.tmp)
    }
    print(cmd)

    system(cmd)

    bcfindex(indel.tmp)

    ## fixing up indels - PCAWG formats are strange...
    ## some REF are empty, with ALT field (insertion)
    ## some ALT are empty with REF field (deletion)
    ## need to fix the coordinates and the ALT/REF fields
    ## HRDetect requires a VCF in left-aligned formats to do their microhomology calculations

    indel_fix = readVcf(TabixFile(indel.tmp))
    rr = rowRanges(indel_fix)
    ins = which(width(rr) == 0 & nchar(as.character(rr$REF)) == 0)
    del = which(width(rr) > 0 & nchar(as.character(unlist(rr$ALT))) == 0)

    if (length(c(ins, del))) {
        message("fixing indels")
    }

    if (length(ins)) {
        rowRanges(indel_fix)[ins] = GenomicRanges::shift(
            gr.resize(rr[ins],
                      1,
                      pad = FALSE,
                      fix = "start"), 1L)
        insref = unname(subseq(unlist(VariantAnnotation::alt(indel_fix)[ins]), 1, 1))
        VariantAnnotation::ref(indel_fix)[ins] = insref
    }

    if (length(del)) {
        rowRanges(indel_fix)[del] = gr.resize(rr[del],
                                              width(rr[del]) + 1,
                                              pad = FALSE, fix = "end")
        delref = unname(get_seq(TwoBitFile("~/DB/GATK/human_g1k_v37.fasta.2bit"),
                                rowRanges(indel_fix)[del]))
        delalt = unname(split(subseq(delref, 1, 1), seq_along(delref)))
        VariantAnnotation::alt(indel_fix)[del] = delalt
        VariantAnnotation::ref(indel_fix)[del] = delref
    }

    suppressWarnings(rm("rr"))

    if (length(c(ins, del))) {
        message("indels fixed... writing to vcf")
        writeVcf(indel_fix, file_path_sans_ext(indel.tmp), index = TRUE)
    }

    suppressWarnings(rm("indel_fix"))
    names(indel.tmp) = "sample_1"

} else {
    indel.tmp = NULL
}


########## SV processing
message("\n\n")
message("processing SV")

jabba = opt$jabba
sv = opt$sv
if (!(identical(jabba, "/dev/null") && identical(sv, "/dev/null"))) {
    fe.jabba = !file.not.exists(jabba)
    fe.sv = !file.not.exists(sv)

    ## overwriting mskilab strand convention with theirs
    ## Sanger BRASS convention, our convention
    strand.conv = setnames(rbind(
        data.table(cbind("+", "+"), cbind("-", "+")),
        data.table(cbind("-", "-"), cbind("+", "-")),
        data.table(cbind("+", "-"), cbind("-", "-")),
        data.table(cbind("-", "+"), cbind("+", "+"))
    ), c("new.strand1", "new.strand2", "strand1", "strand2"))

    if (fe.jabba & fe.sv) {
        warning("both jabba and sv inputs provided.")
    }

    if (fe.jabba | fe.sv) {
        gg = if (fe.jabba) {
                 gG(jabba = jabba)
             } else {
                 gG(junctions = gGnome:::read.juncs(sv))
             }
    } else {
        stop("must provide at least one of jabba or sv vcf inputs")
    }

    altgrl = gg$edges[type == "ALT"]$grl
    sv.tmp = "./sv.bedpe"

    if (length(altgrl)) {
        bdpe = grl2bedpe(sort.GRangesList(gr.noval(altgrl), ignore.strand = TRUE))
        these.cols = colnames(bdpe)
        fwrite(select(left_join(bdpe, strand.conv), !!these.cols, everything()) %>%
               mutate(strand1 = new.strand1,
                      strand2 = new.strand2,
                      new.strand1 = NULL,
                      new.strand2 = NULL,
                      sample = "sample_1"),
               sv.tmp, sep = "\t")
    } else {
        empt = data.table(chrom1 = "", start1 = "", end1 = "",
                          chrom2 = "", start2 = "", end2 = "",
                          name = "", score = "", strand1 = "", strand2 = "",
                          sample = "")[0,]
        fwrite(empt, sv.tmp, sep = "\t")

    }
    names(sv.tmp) = "sample_1"
} else {
    sv.tmp = NULL
}



########## CNV processing
message("\n\n")
message("processing CNV input")

if (file.exists(jabba) && !identical(jabba, "/dev/null"))  {

    jabd.simple  = readRDS(jabba)
    ## jabd.simple = readRDS(jabba.simple.rds)
    jabba.has.hets = all(c("asegstats","aadj","agtrack") %in% names(jabd.simple))
    het.file = opt$hets
    ## hets.exist = !is.na(het.file) && !is.null(opt$hets) && file.exists(het.file)
    hets.exist = !file.not.exists(het.file)

    out.file = opt$jabba_rds
    if (isTRUE(hets.exist) && isFALSE(jabba.has.hets)) {
        message("Pulling het_pileups data")
        verbose = TRUE
        if (!file.exists(paste(dirname(opt$outdir), 'hets.gr.rds', sep = '/'))) {
            if (!is.null(het.file))
            {
                if (is.character(het.file))
                {
                    if (grepl(".rds$", het.file)){
                        hets = readRDS(het.file)
                    } else {
                        hets = fread(het.file)
                    }
                }
                else
                    hets = het.file
                if (!is.data.table(hets))
                    hets = as.data.table(hets)
                if (verbose)
                {
                    message('loaded hets')
                }
                if (inherits(hets, "data.frame")){
                    if (!is.null(hets$alt.count.n) & !is.null(hets$ref.count.n)){
                        ## old format, apply het filter ourselves
                        hets$ref.frac.n = hets$alt.count.n / (hets$alt.count.n + hets$ref.count.n)
                        hets.gr = dt2gr(hets[pmin(ref.frac.n, 1-ref.frac.n) > 0.2 & (ref.count.n + alt.count.n)>=2, ])
                        hets.gr$alt = hets.gr$alt.count.t
                        hets.gr$ref = hets.gr$ref.count.t
                    }
                    else ## new, standard format, with $alt and $ref field
                    {
                        hets.gr = dt2gr(hets)
                        if (all(c("alt", "ref") %in% colnames(hets))){
                            message("Valid hets already")
                            ## hets.gr$alt.count.t = hets.gr$alt
                            ## hets.gr$ref.count.t = hets.gr$ref
                        } else if (all(c("alt.count.t", "ref.count.t") %in% colnames(hets))){
                            hets.gr$alt = hets.gr$alt.count.t
                            hets.gr$ref = hets.gr$ref.count.t
                            hets.gr = hets.gr %Q% (!is.na(alt)) %Q% (!is.na(ref))
                        } else {
                            message("hets is not in valid format, ignore")
                            hets.gr = NULL
                        }
                    }
                } else if (inherits(hets, "GRanges")){
                    if (all(c("alt.count.t", "ref.count.t") %in% colnames(values(hets)))){
                        hets.gr = hets
                        hets.gr$alt = hets$alt.count.t
                        hets.gr$ref = hets$ref.count.t
                        hets.gr = hets.gr %Q% (!is.na(alt) & ! is.na(ref))
                    }
                } else {
                    message("hets is neither data.table nor GRanges, ignore.")
                }
                if (!is.null(hets.gr)){
                    ## save hets object for later
                    saveRDS(hets.gr, paste(opt$outdir, '/hets.gr.rds', sep = '/'))
                }
            }
        } else {
            hets.gr = readRDS("hets.gr.rds")
        }
        if (!all(c('asegstats', 'aadj', 'agtrack') %in% names(jabd.simple))) {
            jaba = tryCatch(JaBbA:::jabba.alleles(jabd.simple, hets.gr, verbose = TRUE, uncoupled=TRUE)[c('asegstats', 'aadj', 'agtrack')], error = function(e) NULL)
            if (!is.null(jaba)) {
                jabd.simple = c(jabd.simple, jaba)
            } else {
                if (sum(mcols(jabd.simple$junctions)$cn > 0, na.rm = T) == 0) {
                    message("could not estimate allelic CN due to no aberrant edges in graph")
                } else {
                    message("could not estimate allelic CN... check inputs")
                }
                message("assuming high and low are 1/2 of total cn")
                jabd.simple$asegstats = grbind(within(jabd.simple$segstats, {cn = ceiling(cn / 2); type = "high"}),
                                               within(jabd.simple$segstats, {cn = floor(cn / 2); type = "low"}))
            }
        }
    }
    cnv = as.data.table(keepStandardChromosomes(jabd.simple$asegstats, pruning.mode = "coarse") %Q%
                        (strand == "+" & !grepl("M|MT", seqnames) & width != 1) %>% #JR edited 4/16 removed 'chr' from filter
                        sortSeqlevels() %>%
                        sort(ignore.strand = TRUE) %>%
                        gr.nochr)[!duplicated(paste(seqnames, start, end, width, strand, type))]
    cnv = setDT(rename_at(dcast.data.table(cnv, seqnames + start + end + width + strand ~ type, value.var = "cn"), vars(high, low), ~paste0(., "_cn")))
    cnv = cnv[, .(seg_no = seq_len(.N),
                  Chromosome = seqnames,
                  chromStart = start,
                  chromEnd = end,
                  total.copy.number.inNormal = 2,
                  minor.copy.number.inNormal = 1,
                  total.copy.number.inTumour = high_cn + low_cn,
                  minor.copy.number.inTumour = low_cn)] %>%
        subset2(complete.cases(x))
}

if (!exists('cnv') || is.null(cnv)) {
    stop("no allelic CN available...")
}

if (is.character(cnv) && cnv == "/dev/null") {
    cnv.tmp = NULL
} else {
    cnv.tmp = "./cna.txt"
    fwrite(cnv, cnv.tmp, sep = "\t")
    names(cnv.tmp) = "sample_1"
}



########## run pipeline
res = signature.tools.lib::HRDetect_pipeline(SNV_vcf_files = snv.tmp,
                                             Indels_vcf_files = indel.tmp,
                                             SV_bedpe_files = sv.tmp,
                                             CNV_tab_files = cnv.tmp,
                                             genome.v = opt$genome, SNV_signature_version = 'COSMICv3.2')

saveRDS(res, './hrdetect_results.rds')

if (NROW(res$hrdetect_output) > 0)
    fwrite(res$hrdetect_output, './hrdetect_output.txt')
else {
    message("HRDetect score was not calculated!!")
    message("Inputs may not be available")
    ## fwrite(data.table(), './hrdetect_output.txt')
}

quit("no", 0)
}, echo = FALSE)
