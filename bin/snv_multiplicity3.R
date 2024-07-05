{
    library(optparse)

    if (!exists('opt'))
    {
        option_list = list(
            make_option(c("-s", "--somatic_snv"), type = "character", help = "Path to somatic snv file"),
            make_option(c("-g", "--germline_snv"), type = "character", help = "Path to germline snv file"),
            make_option(c("-j", "--jabba"), type = "character", help = "Path to jabba file"),
            make_option(c("-f", "--fasta"), type = "character", help = "Path to ref fasta file"),
            make_option(c("-k", "--cores"), type = "integer", help = "Number of cores"),
            make_option(c("-o", "--outdir"), type = "character", default = './', help = "output directory"),
        )

        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)

        if (is.null(opt$somatic_snv) | is.null(opt$jabba) | is.null(opt$germline_snv))
            stop(print_help(parseobj))

        print(opt)

        print(.libPaths())
        options(error=function() {traceback(2); quit("no", 1)})

        ## keep record of run
        writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }

    message("Libraries loading...")
    suppressPackageStartupMessages({
      suppressWarnings({
            #library(Flow)
            library(gGnome)
            library(gUtils)
            library(skitools)
            library(bamUtils)
            library(dplyr)
            library(rtracklayer)
            library(khtools)
            library(skidb)
            devtools::load_all("~/git/zitools")
        })
    })

    # debug.R contents
    est_snv_cn_stub =
        function (vcf, jab, tumbam = NULL, germ_subsample = 2e+05, somatic = FALSE,
                  saveme = FALSE)
    {
        oldsaf = options()$stringsAsFactors
        options(stringsAsFactors = FALSE)
        oldscipen = options()$scipen
        options(scipen = 999)
        on.exit({
            options(scipen = oldscipen)
            options(stringsAsFactors = oldsaf)
            unlink(tmpvcf)
            unlink(tmpvcf2)
        })
        tmpvcf = tempfile(fileext = ".vcf")
        tmpvcf2 = tempfile(fileext = ".vcf")
        if (!somatic) {
            message("starting germline processing")
            system2("bcftools", c("view -i 'FILTER==\"PASS\"'", vcf),
                    stdout = tmpvcf)
            system2("java",
                    sprintf("-jar ~/software/jvarkit/dist/downsamplevcf.jar -N 10 -n %s %s",
                            germ_subsample, tmpvcf),
                    stdout = tmpvcf2,
                    env = "module unload java; module load java/1.8;")
            gvcf = parsesnpeff(
                tmpvcf, coding_alt_only = TRUE, keepfile = FALSE,
                altpipe = TRUE)
            gvcf_subsam = parsesnpeff(tmpvcf2, coding_alt_only = FALSE,
                                      keepfile = FALSE, altpipe = TRUE)
            gvcf = unique(dt2gr(rbind(gr2dt(gvcf_subsam), gr2dt(gvcf))))
            input = gvcf
            rm("gvcf", "gvcf_subsam")
            fif = file.info(dir("./"))
            fif = arrange(cbind(path = rownames(fif), fif), desc(mtime))
            tmp.t = grep("reg_.*.tsv", fif$path, value = TRUE)[1]
            tmp.b = grep("reg_.*.bed", fif$path, value = TRUE)[1]
            callout = grep("mpileup_", fif$path, value = TRUE)[1]
            if (!file.exists(tmp.t))
                tmp.t = tempfile(pattern = "reg_", fileext = ".tsv",
                                 tmpdir = ".")
            if (!file.exists(tmp.b))
                tmp.b = tempfile(pattern = "reg_", fileext = ".bed",
                                 tmpdir = ".")
            if (!file.exists(callout))
                callout = tempfile(pattern = "mpileup_", fileext = ".vcf",
                                   tmpdir = ".")
            ## xtYao ## Wednesday, Feb 24, 2021 09:15:58 AM
            ## Remove column explicitly with $ outside within function
            input = within(input, {
                nref = nchar(REF)
                nalt = nchar(ALT)
                vartype = ifelse(nref > 1 & nalt == 1, "DEL", ifelse(nref ==
                                                                     1 & nalt > 1, "INS", ifelse(nref == 1 & nalt ==
                                                                                                 1, "SNV", NA_character_)))
                maxchar = pmax(nref, nalt)
                ## nref = NULL
                ## nalt = NULL
            })
            input$nref = NULL
            input$nalt = NULL
            input2 = GenomicRanges::reduce(
                gr.resize(input, ifelse(input$maxchar > 1, 201, input$maxchar), pad = FALSE) %>% gr.sort)
            fwrite(gr2dt(input2[, c()])[, 1:3, with = F][, `:=`(start,
                                                                pmax(start, 1))][, `:=`(end, pmax(end, 1))], tmp.t,
                   sep = "\t", col.names = FALSE)
            fwrite(gr2dt(input2[, c()])[, 1:3, with = F][, `:=`(start,
                                                                pmax(start - 1, 0))][, `:=`(end, pmax(end, 1))],
                   tmp.b, sep = "\t", col.names = FALSE)
            if (!file.exists(callout)) {
                message("starting germline mpileup to call variants in tumor")
                clock(system(sprintf("(bcftools mpileup -d 8000 -Q 0 -q 0 -B -R %s -f %s %s | bcftools call -m --prior 0 -v) > %s",
                                     tmp.t, opt$fasta, tumbam, callout)))
            }
            excls = tempfile(pattern = "excludesam_", fileext = ".txt",
                             tmpdir = tempdir())
            writeLines(system(sprintf("bcftools query -l %s", callout),
                              intern = T), excls)
            cntmp = tempfile(pattern = "cntmp_", fileext = ".vcf.gz",
                             tmpdir = "./")
            system(sprintf("(bcftools view -S ^%s %s | bcftools norm -c f -f %s | bcftools norm -Ov -m-any | bgzip -c) > %s",
                           excls, callout, opt$fasta, cntmp))
            cnfin = S4Vectors::expand(readVcf(cntmp))
            dp4mat = do.call(rbind, as.list(info(cnfin)$DP4))
            altv = dp4mat[, 3] + dp4mat[, 4]
            idv = info(cnfin)$IDV
            refv = dp4mat[, 1] + dp4mat[, 2]
            idrefv = (altv + refv) - idv
            gr4est = rowRanges(cnfin)
            gr4est$ref = coalesce(idrefv, refv)
            gr4est$alt = coalesce(idv, altv)
            gr4est$ALT = as.character(gr4est$ALT)
            gr4est$REF = as.character(gr4est$REF)
            gr4est = merge(gr2dt(input), gr2dt(gr4est)[, `:=`(pileupfound,
                                                              TRUE)], by = c("seqnames", "start", "end", "REF",
                                                                             "ALT"), suffixes = c("_normal", ""), all = TRUE,
                           allow.cartesian = TRUE)
            gr4est = unique(gr4est)
            gr4est[, `:=`(pileupfound, pileupfound %in% TRUE)]
            if (saveme)
                saveRDS(gr4est, "gr4est_germline.rds")
            hold = gr4est[is.na(ref)]
            hold[, `:=`(pileupnotfound, TRUE)]
            germbin = gr4est[is.na(ref_normal)]
            if (saveme)
                saveRDS(germbin, "germpileupbin.rds")
            gr4est = gr4est[!is.na(ref)][!is.na(ref_normal)]
        }
        else {
            message("reading in somatic variants")
            gr4est = parsesnpeff(vcf, coding_alt_only = FALSE, keepfile = FALSE,
                                 altpipe = TRUE)
            if (saveme)
                saveRDS(gr4est, "gr4est_somatic.rds")
            hold = NULL
        }
        out = est_snv_cn(gr4est, jab, somatic = somatic)
        out = rbind(out, hold, fill = TRUE)
        if (saveme)
            if (somatic)
                saveRDS(out, "est_snv_cn_somatic.rds")
            else saveRDS(out, "est_snv_cn_germline.rds")
        return(out)
    }

    est_snv_cn =
        function (gr, jabba, somatic = FALSE)
    {
        if (length(gr) == 0)
            return(NULL)
        if (is.character(jabba))
            jab = readRDS(jabba)
        gg = gG(jabba = jabba)
        lpp = with(jab, list(purity = purity, ploidy = ploidy))
        if (inherits(gr, "data.table"))
            gr = dt2gr(gr)
        gr = gr %*% within(gg$nodes$gr[, c("snode.id", "cn")], {
            segwid = width
        })
        dt = gr2dt(gr)
        dt = dt[order(seqnames, start, end, -alt)][, `:=`(rtot, ref +
                                                                alt)]
        dt = cbind(dt, with(dt, as.data.table(rleseq(seqnames, start,
                                                     REF, ALT, clump = T))))
        dt[, `:=`(rtot, sum(ref[1], alt[1])), by = idx]
        dt[, `:=`(seg_rtot, {
            u = !duplicated(idx)
            mean(rtot[u]) %>% round
        }), by = snode.id]
        dt[, `:=`(vaf, alt/rtot)]
        dt[, `:=`(vaf_segt, pinch(alt/seg_rtot, 0, 1))]
        if (isFALSE(somatic)) {
            message("calculating cn of somatic variants")
            dt[, `:=`(norm_term, ifelse(grepl("^0[/|]1$", gt), 2,
                                 ifelse(grepl("^1[/|]1$", gt), 1, NA_integer_)))]
            dt[!is.na(cn), `:=`(c("est_cn", "est_cn_rm", "est_cn_ll", "est_cn_llrm"), {
                cn = cn[1]
                rtot = rtot[1]
                seg_rtot = seg_rtot[1]
                vaf_segt = vaf_segt[1]
                alt = alt[1]
                norm_term = norm_term[1]
                estcn = round(
                (cn *
                 (
                     (norm_term * vaf) - (1 - lpp$purity)
                 )
                ) / lpp$purity
                )
                estcnrm = round((cn * ((norm_term * vaf_segt) - ((1 -
                                                                  lpp$purity))))/lpp$purity)
                centers = pinch((lpp$purity * (0:cn)/cn) + ((1 -
                                                             lpp$purity)/2))
                ifun = function(cnid, rtot, vaf, alt) {
                    dbinom(alt, rtot, prob = centers[cnid + 1], log = T)
                }
                estllcn = which.max(
                    withv(sapply((0:cn), ifun, rtot = rtot,
                                 vaf = vaf, alt = alt), x - min(x))) - 1
                estllrcn = which.max(withv(sapply((0:cn), ifun, rtot = seg_rtot,
                                                  vaf = vaf_segt, alt = alt), x - min(x))) - 1
                list(estcn, estcnrm, estllcn, estllrcn)
            }), by = .(snode.id, idx)]
        }
        else {
            message("calculating cn of normal variants")
            dt[!is.na(cn), `:=`(c("est_cn", "est_cn_rm", "est_cn_ll", "est_cn_llrm"), {
                cn = cn[1]
                rtot = rtot[1]
                seg_rtot = seg_rtot[1]
                vaf_segt = vaf_segt[1]
                alt = alt[1]
                estcn = round((cn * (2 * vaf))/lpp$purity)
                estcnrm = round((cn * (2 * vaf_segt))/lpp$purity)
                centers = pinch((lpp$purity * (0:cn)/cn))
                ifun = function(cnid, rtot, vaf, alt) {
                    out = dbinom(alt, rtot, prob = centers[cnid +
                                                           1], log = T)
                    out = replace(out, is.infinite(out) & sign(out) <
                                       0, -2e+09)
                    out = replace(out, is.infinite(out) & sign(out) >
                                       0, 2e+09)
                    return(out)
                }
                estllcn = which.max(withv(sapply((0:cn), ifun, rtot = rtot,
                                                 vaf = vaf, alt = alt), x - min(x))) - 1
                estllrcn = which.max(withv(sapply((0:cn), ifun, rtot = seg_rtot,
                                                  vaf = vaf_segt, alt = alt), x - min(x))) - 1
                list(estcn, estcnrm, estllcn, estllrcn)
            }), by = .(snode.id, idx)]
        }
        return(dt)
    }

    ## copied from /gpfs/commons/groups/imielinski_lab/home/khadi/git/khtools/R/package.R
    #' @name parsesnpeff
    #' @title parse snpeff output into granges
    #'
    #'
    #' @param vcf path to snpeff vcf
    #' @param pad Exposed argument to skitools::ra.overlaps()
    #' @return GRangesList of breakpoint pairs with junctions that overlap removed
    #' @export
    parsesnpeff = function (vcf, id = NULL, filterpass = TRUE, coding_alt_only = TRUE,
                            geno = NULL, gr = NULL, keepfile = FALSE, altpipe = FALSE,
                            debug = FALSE)
    {
        if (debug)
            browser()
        out.name = paste0("tmp_", rand.string(), ".vcf.gz")
        tmp.path = paste0(tempdir(), "/", out.name)
        if (!keepfile)
            on.exit(unlink(tmp.path))
        try2({
            catcmd = if (grepl("(.gz)$", vcf)) "zcat" else "cat"
            onepline = "/gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/scripts/vcfEffOnePerLine.pl"
            if (coding_alt_only) {
                filt = "java -Xmx20m -Xms20m -XX:ParallelGCThreads=1 -jar /gpfs/commons/groups/imielinski_lab/git/mskilab/flows/modules/SnpEff/source/snpEff/SnpSift.jar filter \"( ANN =~ 'chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame' )\""
                if (filterpass)
                    cmd = sprintf(paste(catcmd, "%s | %s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"),
                                  vcf, onepline, filt, tmp.path)
                else cmd = sprintf("cat %s | %s | %s | bcftools norm -Ov -m-any | bgzip -c > %s",
                                   vcf, onepline, filt, tmp.path)
            }
            else {
                filt = ""
                if (filterpass)
                    cmd = sprintf(paste(catcmd, "%s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"),
                                  vcf, onepline, tmp.path)
                else cmd = sprintf(paste(catcmd, "%s | %s | bcftools norm -Ov -m-any | bgzip -c > %s"),
                                   vcf, onepline, tmp.path)
            }
            system(cmd)
        })
        if (!altpipe)
            out = grok_vcf(tmp.path, long = TRUE, geno = geno, gr = gr)
        else {
            vcf = readVcf(tmp.path)
            vcf = S4Vectors::expand(vcf)
            rr = within(rowRanges(vcf), {
                REF = as.character(REF)
                ALT = as.character(ALT)
            })
            ann = as.data.table(tstrsplit(unlist(info(vcf)$ANN),
                                          "\\|"))[, 1:15, with = FALSE, drop = FALSE]
            fn = c("allele", "annotation", "impact", "gene", "gene_id",
                   "feature_type", "feature_id", "transcript_type",
                   "rank", "variant.c", "variant.p", "cdna_pos", "cds_pos",
                   "protein_pos", "distance")
            data.table::setnames(ann, fn)
            if ("AD" %in% names(geno(vcf))) {
                adep = setnames(as.data.table(geno(vcf)$AD[, , 1:2]),
                                c("ref", "alt"))
                gt = geno(vcf)$GT
            }
            else if (all(c("AU", "GU", "CU", "TU", "TAR", "TIR") %in%
                         c(names(geno(vcf))))) {
                this.col = dim(geno(vcf)[["AU"]])[2]
                d.a = geno(vcf)[["AU"]][, , 1, drop = F][, this.col,
                                                         1]
                d.g = geno(vcf)[["GU"]][, , 1, drop = F][, this.col,
                                                         1]
                d.t = geno(vcf)[["TU"]][, , 1, drop = F][, this.col,
                                                         1]
                d.c = geno(vcf)[["CU"]][, , 1, drop = F][, this.col,
                                                         1]
                mat = cbind(A = d.a, G = d.g, T = d.t, C = d.c)
                rm("d.a", "d.g", "d.t", "d.c")
                refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(mat))
                refid = ifelse(!isSNV(vcf), NA_integer_, refid)
                altid = match(as.character(VariantAnnotation::fixed(vcf)$ALT), colnames(mat))
                altid = ifelse(!isSNV(vcf), NA_integer_, altid)
                refsnv = mat[cbind(seq_len(nrow(mat)), refid)]
                altsnv = mat[cbind(seq_len(nrow(mat)), altid)]
                this.icol = dim(geno(vcf)[["TAR"]])[2]
                refindel = d.tar = geno(vcf)[["TAR"]][, , 1, drop = F][,
                                                                       this.icol, 1]
                altindel = d.tir = geno(vcf)[["TIR"]][, , 1, drop = F][,
                                                                       this.icol, 1]
                adep = data.table(ref = coalesce(refsnv, refindel),
                                  alt = coalesce(altsnv, altindel))
                gt = NULL
            }
            else {
                message("ref and alt count columns not recognized")
                adep = NULL
                gt = NULL
            }
            ## mcols(rr) = BiocGenerics::cbind(mcols(rr), ann, adep,
            ##                                 gt = gt[, 1])
            mcols(rr) = cbind(as.data.table(mcols(rr)),
                              ann,
                              adep,
                              gt = gt[, 1])
            out = rr
        }
        this.env = environment()
        return(this.env$out)
    }


    #' @name rand.string
    #' @title make a random string
    #'
    #' @return random string
    #' @author Someone from Stackoverflow
    #' @export rand.string
    rand.string = function(n=1, length=12)
    {
        randomString <- c(1:n)                  # initialize vector
        for (i in 1:n)
        {
            randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                            length, replace=TRUE),
                                     collapse="")
        }
        return(randomString)
    }

    #' @name try2
    #' @title wrapper around tryCatch - robust to parallel:: functions
    #'
    #' A slightly more robust version of try that works within the parallel:: set of functions
    #' that pre-deploy a cluster.
    #'
    #' @export
    try2 = function(expr, ..., finally) {
        tryCatch(expr,
                 error = function(e) {
                     msg = structure(paste(conditionMessage(e), conditionCall(e), sep = "\n"), class = "err")
                     cat("Error: ", msg, "\n\n")
                     return(msg)
                 },
                 finally = finally,
                 ... = ...)
    }

    # end of debug.R

    message("\nLibraries loaded!")

    data.table::setDTthreads(1)
    options(scipen = 999)
    options(stringsAsFactors = FALSE)

    message("Loading snv")

    gg = gG(jabba = opt$jabba)

    somatic.snv = parsesnpeff(opt$somatic_snv,
                              coding_alt_only = FALSE,
                              filterpass = FALSE,
                         keepfile = FALSE,
                         altpipe = TRUE)
    names(somatic.snv) <- NULL
    ## germline.snv = parsesnpeff(opt$germline_snv,
    ##                            coding_alt_only = FALSE,
    ##                            filterpass = FALSE,
    ##                            keepfile = FALSE,
    ##                            altpipe = TRUE)
    ## names(germline.snv) <- NULL

    #normalization
    unique.somatic.snv = somatic.snv[,c('ref', 'alt')] %>%
      gr.nochr %Q%
      (seqnames %in% c(1:22, "X")) %>%
      unique
    unique.somatic.snv$major.count <- pmax(unique.somatic.snv$ref, unique.somatic.snv$alt)
    unique.somatic.snv$minor.count <- pmin(unique.somatic.snv$ref, unique.somatic.snv$alt)
    m = length(unique.somatic.snv)
    sf = sum(unique.somatic.snv$major.count + unique.somatic.snv$minor.count) / m
    unique.somatic.snv$major.count <- unique.somatic.snv$major.count / sf
    unique.somatic.snv$minor.count <- unique.somatic.snv$minor.count / sf

    #constitutional_cn assignment
    #c_subj == 1 for major allele
    #c_subj == 1 for autosomes and X chromosome in females, 0 for X and Y in males
    ncn.x = gg$nodes$dt[(seqnames == "X" | seqnames == "chrX"),
                         weighted.mean(cn,
                                       w = end - start + 1,
                                       na.rm = TRUE)]
    message("mean CN of X: ", ncn.x)
    ncn.vec = rep(2, length(unique.somatic.snv))
    if (ncn.x < 1.4) { # if male
      message("Adjusting ncn for XY")
      ncn.vec[which(as.character(seqnames(unique.somatic.snv)) %in% c("chrX", "chrY", "X", "Y"))] = 1
    }
    #values(unique.somatic.snv)[, "ncn"] = ncn.vec
    values(unique.somatic.snv)[, "minor_constitutional_cn"] = ncn.vec - 1

    jab = gg2jab(gg)
    ss.p = jab$segstats[ as.logical( strand(jab$segstats)=='+' ) ]

    unique.somatic.snv = gr.val(unique.somatic.snv, ss.p, "cn", na.rm = T)

    tau_hat = mean(unique.somatic.snv$cn)
    alpha = jab$purity
    beta = alpha / (alpha * tau_hat + 2*(1 - alpha))
    gamma = 2*(1 - alpha) / (alpha * tau_hat + 2*(1 - alpha))

    mcols(unique.somatic.snv)$major_snv_copies =
    (2 * unique.somatic.snv$major.count - gamma) / (2 * beta)

    mcols(unique.somatic.snv)$minor_snv_copies =
    (2 * unique.somatic.snv$major.count - gamma * unique.somatic.snv$minor_constitutional_cn)/ 2 * beta

    somatic.variants = gr.val(somatic.snv %>% gr.nochr,
                              unique.somatic.snv[,c('major.count',
                                                    'minor.count',
                                                    'major_snv_copies',
                                                    'minor_snv_copies')],
                              c('major.count', 'minor.count',
                                'major_snv_copies', 'minor_snv_copies')) %Q%
      (FILTER == "PASS") %>%
      gr2dt

    somatic.variants[alt >= ref, total_copies := major_snv_copies]
    somatic.variants[alt < ref, total_copies := minor_snv_copies]
    somatic.variants[, VAF := alt / (alt + ref)]
    somatic.variants = somatic.variants %>% dt2gr
    message('done')

    ## save file
    saveRDS(somatic.variants, paste0(opt$outdir, "est_snv_cn_somatic.rds"))

    quit("no", 0)
}
