## withAutoprint(
## {
library(optparse)

if (!exists('opt'))
{
    option_list = list(
        make_option(c("-j", "--junctions"), type = "character", help = "Path to .vcf or bedpe or rds file of junctions or gGraph from which alt edges will be taken "),
        make_option(c("-w", "--width"), type = "integer", default = 50, help = "width around junction breakpoint around which to search for homeology"),
        make_option(c("-p", "--pad"), type = "integer", default = 0, help = "number of bases of padding around each sequence position (bin) to use when computing homeology, i.e. we then will be comparing 1 + 2*pad -mer sequences for edit distance "),
        make_option(c("-t", "--thresh"), type = "integer", default = 0, help = "string distance threshold for calling homeology in a bin"),
        make_option(c("-s", "--stride"), type = "integer", default = 0, help = "distance in bases between consecutive bins in which we will be measuring homeology"),
        make_option(c("-G", "--genome"), type = "character", help = "Path to .2bit or ffTrack .rds containing genome sequence"),
        make_option(c("-c", "--cores"), type = "integer", help = "How many cores to use"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "output directory"),
        make_option("--flip", type = "logical", default = FALSE, help = "if flip = FALSE, homeology search for -/- and +/+ junctions is done between a sequence and its reverse complement"),
        make_option("--bidirectional", type = "logical", default = TRUE, help = "adding padding on both sides of each breakpoint (TRUE) or only in the direction of the fused side (FALSE)"),
        make_option("--annotate", type = "logical", default = TRUE, help = "annotate edges in gGraph object and save it in working directory"),
        make_option("--savegMatrix", type = "logical", default = TRUE, help = "save gMatrix object of edit distances"),
    )

    parseobj = OptionParser(option_list=option_list)
    opt = parse_args(parseobj)

    if ((is.null(opt$junctions)) | is.null(opt$genome))
        stop(print_help(parseobj))

    print(opt)

    print(.libPaths())
    options(error=function() { traceback(2); quit("no", 1) })

    ## keep record of run
    writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
    saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
}

suppressWarnings(expr = {
    suppressMessages(expr = {
        suppressPackageStartupMessages(expr = {
            library(skitools)
            library(skidb)
            library(GxG)
            library(data.table)
            library(gGnome)
            library(gUtils)
            library(khtools)
            library(purrr)
            library(imager)
        })
    })
})

setDTthreads(1)

system(paste('mkdir -p',  opt$outdir))

ref = readDNAStringSet(gsub("^(.*(fasta)|(fa)).*", "\\1", opt$genome))
names(ref) = sapply(names(ref), function(x) {strsplit(x, split = " ")[[1]][1]})

hom = function (event, pad = 100, thresh = 2, stride = 1, pad2 = 5,
                flip = FALSE, mc.cores = 1, anchor = TRUE, deanchor_gm = TRUE,
                mat = FALSE,
                genome = "/gpfs/commons/home/khadi/DB/GATK/human_g1k_v37.fasta.2bit",
                bidirected_search = TRUE,
                save_gm = TRUE)
{
    if (!NROW(event)) {
        return(list(gm = list(), rawres = data.table(), res = data.table()))
    }
    event = copy2(event)
    ## event = data.table(bp1 = gr.string(event$left), bp2 = gr.string(event$right))
    gmfeat = function(gm, thresh = 3, op = "<=") {
        if (is(gm, "gMatrix")) {
            mat = gm$mat
            mat = (mat + t(mat))/2
        }
        else mat = gm
        ## dehom = deanchor(jJ(GRangesList(c(dg(evbp1)[dg(i)], dg(evbp2)[dg(i)]))),  dg(hom))
        ## demat = dehom$mat[dg(ix1), dg(ix2)]
        im = imager::as.cimg(as.matrix(mat))
        cmd = sprintf("im %s %s", op, thresh)
        px = eval(parse(text = cmd))
        if (sum(px) == 0)
            return(data.table())
        sp = imager::split_connected(px)
        if (NROW(sp) == 0)
            return(data.table())
        relmat = (1 - (mat / max(mat, na.rm = T)))

        res = rbindlist(lapply(
            1:length(sp),
            function(k) as.data.table(which(as.matrix(sp[[k]]),
                                            arr.ind = TRUE))[
                           ,.(k = k, i = pmin(row, col),
                              j = pmax(row, col), levd = mat[cbind(row, col)],
                              relsim = relmat[cbind(row, col)],
                              iw = rownames(mat)[row],
                              jw = colnames(mat)[col])]), fill = TRUE)
        if (NROW(res) > 0) {
            griw = gr2dt(with(parse.gr2(res$iw, meta = res), {
                SHIFT = abs(min(pmin(start, end))) + 1;
                GenomicRanges::shift(dynGet("data"), SHIFT) %>% gr.spreduce(k)
            }))[, .(iwid = sum(width)), by = .(k = as.integer(k))]
            grjw = gr2dt(with(parse.gr2(res$jw, meta = res), {
                SHIFT = abs(min(pmin(start, end))) + 1;
                GenomicRanges::shift(dynGet("data"), SHIFT) %>% gr.spreduce(k)
            }))[, .(jwid = sum(width)), by = .(k = as.integer(k))]
            res = merge.data.table(merge.data.table(res, griw, by = 'k', all.x = TRUE),
                                   grjw, by = "k", all.x = TRUE)[, minpx := pmin(iwid, jwid)]
            res[, `:=`(N = .N, r = cor(i, j)), by = k]
        }
        return(res)
    }
    gmstats = function(res) {
        if (NROW(res) > 0) {
            res[!duplicated(k)][
               ,.(numfeat = sum(N > 0), numfeat2 = sum(N > 2), numfeat5 = sum(N > 5),
                  numfeat10 = sum(N > 10), maxfeat = max(N),
                  numlines5 = sum(r > 0.5, na.rm = TRUE),
                  maxlines5 = max(c(0, N[r > 0.5]), na.rm = TRUE),
                  numlines = sum(r > 0.9, na.rm = TRUE),
                  maxlines = max(c(0, N[r > 0.9]), na.rm = TRUE),
                  maxcor = max(r, na.rm = TRUE),
                  numglines = sum(na2false(r > 0.9) &
                                  na2false(N >= 24) &
                                  na2false(floor(N / minpx) <= 4)),
                  numglines_pm5 = sum(na2false(r > 0.9) & na2false(minpx >= 8)),
                  numglines_p10 = sum(na2false(r > 0.9) & na2false(minpx >= 15)),
                  numglines_p20 = sum(na2false(r > 0.9) & na2false(minpx >= 25)),
                  hlen = max(ifelse(na2false(r > 0.9), minpx, 0L)))]
        }
        else data.table(numfeat = 0, maxfeat = 0)
    }
    evbp1 = gr.end(gr.flipstrand(parse.gr(event$bp1)))
    if (flip)
        evbp2 = gr.flipstrand(gr.end(parse.gr(event$bp2)))
    else evbp2 = gr.end(parse.gr(event$bp2))



    if (isTRUE(bidirected_search)) {
        win = c(GRanges("Left:0") + pad, GRanges("Right:0") + pad)
        query.bp1 = evbp1 + pad
        query.bp2 = evbp2 + pad
    } else {
        win = c(gr.resize(GRanges("Left:0"), fix = "end", pad * 2, pad = F),
                gr.resize(GRanges("Right:0"), fix = "start", pad * 2, pad = F))
        query.bp1 = gr.resize(evbp1, pad * 2, pad = FALSE, fix = "end") # fix at end because was flipped
        bp2dir = if (flip) "end" else "start"
        query.bp2 = gr.resize(evbp2, pad * 2, pad = FALSE, fix = bp2dir)
    } # querying by both sides of junction, or just looking in the direction of the fused side of junction
    event$query.bp1 = gr.string(query.bp1)
    event$query.bp2 = gr.string(query.bp2)
    seq1 = ffTrack::get_seq(genome,
                            query.bp1)
    seq2 = ffTrack::get_seq(genome,
                            query.bp2)

    ## win = c(GRanges("Left:0") + pad, GRanges("Right:0") + pad)
    ## evbp1 = gr.end(gr.flipstrand(parse.gr(event$bp1)))
    ## if (flip)
    ##     evbp2 = gr.flipstrand(gr.end(parse.gr(event$bp2)))
    ## else evbp2 = gr.end(parse.gr(event$bp2))
    ## seq1 = ffTrack::get_seq(genome,
    ##                         evbp1 + pad)
    ## seq2 = ffTrack::get_seq(genome,
    ##                         evbp2 + pad)
    ## define reference genome

    ifun = function(i) {
        tryCatch({
            message(i)
            this.env = environment()
            if (!anchor)
                win = c(evbp1[i], evbp2[i]) + pad
            win$seq = c(seq1[i], seq2[i])
            ## hom = homeology(win, stride = stride, pad = pad2)
            ## browser()
            refseq = win$seq
            names(refseq) = seqnames(win)
            homeology.ref = DNAStringSet(x = refseq)
                                         ## start = GenomicRanges::start(win),
                                         ## end = GenomicRanges::end(win),
                                         ## use.names = TRUE)
            hom = GxG::homeology(ref = homeology.ref,##win$seq,
                                 stride = stride,
                                 pad = pad2)
            ## ahh replace the sequence names
            ## browser()
            ix1 = which(seqnames(hom$gr) == seqnames(win[1]))## which(hom$gr %^% win[1])
            ix2 = which(seqnames(hom$gr) == seqnames(win[2]))## which(hom$gr %^% win[2])
            gm = gmfeat(hom$mat[ix1, ix2], thresh = thresh)
            gms = gmstats(gm)
            gm[, seq := this.env$i]
            gms[, seq := this.env$i]
            if (anchor && deanchor_gm)
                hom = deanchor(jJ(GRangesList(grbind(evbp1[i], evbp2[i]))), hom)
            return(list(hom = hom, gm = gm, gms = gms))
        }, error = function(e) printerr(i))
    }
    ## browser()
    lst = mclapply(1:length(seq1), ifun, mc.cores = mc.cores)
    lst = purrr::transpose(lst)

    rawres = merge.data.table(event[, seq := seq_len(.N)], as.data.table(rbindlist(lst[[2]], fill = T)), by = "seq", all.x = TRUE)
    res = cbind(event, rbindlist(lst[[3]], fill = TRUE))

    if (save_gm)
        return(list(gm = lst[[1]], rawres = rawres, res = res))
    else
        return(list(gm = NULL, rawres = rawres, res = res))
}

deanchor = function(jj, gm, flip = FALSE) {
    ## browser()
    if (flip)
        bps = gr2dt(c(gr.flipstrand(jj$left), jj$right))[, .(sn = seqnames, st = start, en = end, str = strand, key = c('Left', 'Right'))]
    else
        bps = gr2dt(c(jj$left, jj$right))[, .(sn = seqnames, st = start, en = end, str = strand, key = c('Left', 'Right'))]
    setkey(bps, "key") ##ah it turns from 1-2 to left/right
    newgr = dt2gr(gr2dt(gm$gr)[, .(seqnames = bps[.(seqnames), sn],
                                   start = ifelse(bps[.(seqnames), str] == '+',
                                                  start + bps[.(seqnames), st],
                                                  bps[.(seqnames), st]-end),
                                   end = ifelse(bps[.(seqnames), str] == '+',
                                                end + bps[.(seqnames), st],
                                                bps[.(seqnames), st]-start))])
    ## newgr = dt2gr(gr2dt(gm$gr)[, .(seqnames = bps[, sn],
    ##                                start = ifelse(bps[, str] == '+',
    ##                                               start + bps[, st],
    ##                                               bps[, st]-end),
    ##                                end = ifelse(bps[, str] == '+',
    ##                                             end + bps[, st],
    ##                                             bps[, st]-start))])
    new.gm = gM(newgr, gm$dat)
    return(new.gm)
}


## junctions = tryCatch(gG(junction = opt$junctions)$edges[type == "ALT"]$junctions, error = function(e) NULL)

## junctions = tryCatch(gGnome:::read.juncs(opt$junctions), error = function(e) NULL)


## if (is.null(junctions))
##     junctions = tryCatch(readRDS(opt$junctions)$edges[type == 'ALT']$junctions, error = function(e) NULL)

## if (is.null(junctions))
##     junctions = tryCatch(jJ(readRDS(opt$junctions)), error = function(e) NULL)

## if (is.null(junctions))
##     junctions = tryCatch(readRDS(opt$junctions), error = function(e) NULL)

gg = tryCatch(gG(junction = opt$junctions), error = function(e) NULL)

if (is.null(gg))
    gg = tryCatch(readRDS(opt$junctions), error = function(e) NULL)

if (!inherits(gg, "gGraph")) gg = gG(junctions = gg)

junctions = tryCatch(gg$edges[type == "ALT"]$junctions, error = function(e) NULL)


if (is.null(junctions) || !inherits(junctions, "Junction"))
    stop('Error reading junctions - input should be either vcf, bedpe, .rds to Junction object, GRangesList, or gGraph')


standardchr = which(as.logical(seqnames(junctions$left) %in% GenomeInfoDb::standardChromosomes(junctions$left) &
                               seqnames(junctions$right) %in% GenomeInfoDb::standardChromosomes(junctions$right)))


message("junctions with either breakends mapped outside of standard chromosomes will be thrown out")
message("TRUE = junctions with breakends mapped to standard chromosomes")
print(NROW(standardchr))

junctions = junctions[standardchr]

si = seqinfo(TwoBitFile(opt$genome))

ll = gr.nochr(junctions$left)
rr = gr.nochr(junctions$right)


conform_si = function(x, si) {
    ans = copy3(x)
    osi = seqinfo(x)
    osn = as.character(seqnames(x))
    newslev = union(seqlevels(si), seqlevels(osi))
    new = rn2col(as.data.frame(si[newslev]), "seqnames")
    old = rn2col(as.data.frame(osi[newslev]), "seqnames")
    newsi = merge.repl(
        new, old, force_y = FALSE, overwrite_x = FALSE,
        by = "seqnames", keep_order = T, all = T
    )
    newsi = as(col2rn(asdf(newsi), "seqnames"), "Seqinfo")
    ans@seqnames = Rle(factor(osn, seqlevels(newsi)))
    ans@seqinfo = newsi
    return(ans)
}


outofbounds = union(
    GenomicRanges:::get_out_of_bound_index(conform_si(ll, si) + (ceiling(opt$width) / 2)),
    GenomicRanges:::get_out_of_bound_index(conform_si(rr, si) + (ceiling(opt$width) / 2))
)

if (NROW(outofbounds) > 0) {
    message(NROW(outofbounds), " out of bound junctions found")
    junctions = junctions[-outofbounds]
}

events = data.table(bp1 = gr.string(gr.nochr(junctions$left)),
                    bp2 = gr.string(gr.nochr(junctions$right)))

print(events)

cmd=sprintf("hom(events, pad = ceiling(%s/2), thresh = %s, stride = %s, pad2 = %s, genome = '%s', mc.cores = %s, bidirected_search = %s, flip = %s, save_gm = %s)", opt$width, opt$thresh, opt$stride, opt$pad, opt$genome, opt$cores, opt$bidirectional, opt$flip, opt$savegMatrix)

if (!file.exists("res.rds")) {
    message("running ", cmd)
    res = et(cmd)
    message("finished querying for homeology, saving as intermediate")
    saveRDS(res, 'res.rds', compress = FALSE)
} else {
    message("results already exist...")
    message("reading in to perform post-processing")
    res = readRDS("res.rds")
}


message("Post-processing")
stat = res[[3]]

print(stat)
stat = as.data.table(cbind(stat, junctions$dt)) %>% dedup.cols

keep = which(sapply(stat, class) %in% c('character', 'integer', 'numeric', 'factor', 'logical'))
stat = stat[, keep, with = FALSE]

##  ,.(numfeat = sum(N > 0), numfeat2 = sum(N > 2), numfeat5 = sum(N > 5),
## numfeat10 = sum(N > 10), maxfeat = max(N),
## numlines5 = sum(r > 0.5, na.rm = TRUE),
## maxlines5 = max(c(0, N[r > 0.5]), na.rm = TRUE),
## numlines = sum(r > 0.9, na.rm = TRUE),
## maxlines = max(c(0, N[r > 0.9]), na.rm = TRUE),
## maxcor = max(r, na.rm = TRUE),
## numglines = sum(na2false(r > 0.9) &
##                 na2false(N >= 24) &
##                 na2false(floor(N / minpx) <= 4)),
## numglines_pm5 = sum(na2false(r > 0.9) & na2false(minpx >= 8)),
## numglines_p10 = sum(na2false(r > 0.9) & na2false(minpx >= 15)),
## numglines_p20 = sum(na2false(r > 0.9) & na2false(minpx >= 25)))]

if (opt$annotate) {
    message("annotating gGraph edges with homeology features")
    added_cols = c("numfeat", "numfeat2", "numfeat5",
                   "maxfeat", "numfeat10", "numlines5",
                   "maxlines5", "numlines", "maxlines",
                   "maxcor", "numglines", "numglines_pm5",
                   "numglines_p10", "numglines_p20", "hlen")
    if (NROW(stat)) {
        ## ed.annotate = dplyr::select(khtools::dedup.cols(stat, remove = T), edge.id, matches("numglines"))
        ed.annotate = dplyr::select(khtools::dedup.cols(stat, remove = T), edge.id, one_of(added_cols))
        ## cols = grep("numglines", colnames(ed.annotate), value = T)
        cols = colnames(ed.annotate)[colnames(ed.annotate) %in% added_cols]
        for (col in cols)
            et(sprintf("gg$edges[as.character(ed.annotate$edge.id)]$mark(%s = ed.annotate$%s)", col, col))
    }
    message("saving annotated gGraph")
    saveRDS(gg, "./marked_gGraph.rds")
}

message("Saving junction-level stats")
fwrite(stat, paste(opt$outdir, 'junctions.txt', sep = '/'), sep = '\t')
saveRDS(stat, paste(opt$outdir, 'junctions.rds', sep = '/'))

rawstat = res[[2]]

if (!inherits(rawstat, "data.table")) setDT(rawstat)

if (NROW(rawstat)) {
    rawstat = rawstat[
      , which(sapply(rawstat, class) %in%
              c('character', 'integer', 'numeric', 'factor', 'logical')),with = FALSE]
    if (is.null(rawstat$k)) rawstat$k = NA_integer_
    rawstat[, featuid := rleseq(seq, k, na.ignore = T, clump = T)$idx]
    if (is.null(rawstat$featuid)) rawstat$featuid = NA_integer_
    rawstat = merge(rawstat, stat[, .(seq, edge.id)], by = "seq")
}

message("Saving junction feature-level stats")
fwrite(rawstat, paste(opt$outdir, 'rawstats.txt', sep = "/"), sep = '\t')

saveRDS(rawstat, paste(opt$outdir, 'rawstats.rds', sep = '/'))

if (NROW(res[[1]])) message("Saving gMatrix as list")
saveRDS(res[[1]], paste(opt$outdir, 'gMatrixList.rds', sep = '/'))

quit("no", 0)
## }, echo = FALSE)

