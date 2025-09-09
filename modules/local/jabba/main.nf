process JABBA {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/jabba:0.0.1':
        'mskilab/jabba:0.0.1' }"

    input:
    tuple val(meta), path(junction), path(cov_rds), val(j_supp), val(het_pileups_wgs), val(purity), val(ploidy), val(cbs_seg_rds), val(cbs_nseg_rds)
    val(blacklist_junctions)    // this is declared as val to allow for "NULL" default value, but is treated like a path
    val(geno)
    val(indel)
    val(tfield)
    val(iter)
    val(rescue_window)
    val(rescue_all)
    val(nudgebalanced)
    val(edgenudge)
    val(strict)
    val(allin)
    val(field)
    val(maxna)
    path(blacklist_coverage)
    val(pp_method)
    val(cnsignif)
    val(slack)
    val(linear)
    val(tilim)
    val(epgap)
    val(fix_thres)
    val(lp)
    val(ism)
    val(filter_loose)
    val(gurobi)
    val(verbose)

    output:
    tuple val(meta), path("*jabba.simple.rds")      , emit: jabba_rds, optional: true
    tuple val(meta), path("*jabba.simple.gg.rds")   , emit: jabba_gg, optional: true
    tuple val(meta), path("*jabba.simple.vcf")      , emit: jabba_vcf, optional: true
    tuple val(meta), path("*jabba.raw.rds")         , emit: jabba_raw_rds, optional: true
    tuple val(meta), path("*opt.report.rds")        , emit: opti, optional: true
    tuple val(meta), path("*jabba.seg")             , emit: jabba_seg, optional: true
    tuple val(meta), path("*karyograph.rds")        , emit: karyograph, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def geno_switch = geno == "TRUE" ? "--geno" : ""
    def strict_switch = strict == "TRUE" ? "--strict" : ""
    def allin_switch = allin == "TRUE" ? "--allin" : ""
    def linear_switch = linear == "TRUE" ? "--linear" : ""
    def verbose_switch = verbose == "TRUE" ? "--verbose" : ""
    def cbs_seg_rds = cbs_seg_rds == '/dev/null' || cbs_seg_rds == 'NA' || cbs_seg_rds == 'NULL' ? "" : "--seg ${cbs_seg_rds}"
    def cbs_nseg_rds = cbs_nseg_rds == '/dev/null' || cbs_nseg_rds == 'NA' || cbs_nseg_rds == 'NULL' ? "" : "--nseg ${cbs_nseg_rds}"
    def het_pileups_wgs = het_pileups_wgs == '/dev/null' || het_pileups_wgs == 'NA' || het_pileups_wgs == 'NULL' ? "" : "--hets ${het_pileups_wgs}"
    def j_supp = j_supp == '/dev/null' || j_supp == 'NA' || j_supp == 'NULL' ? "" : "--j.supp ${j_supp}"
    def trelim_mem = (task.memory.toGiga() * 0.80).toInteger()


    def VERSION    = '1.1'
    """
    #!/bin/bash

    set -o allexport

    # Check if the environment has the module program installed
    if command -v module &> /dev/null
    then
        # Check if the required modules are available
        if module avail R/4.0.2 &> /dev/null && module avail gcc/9.2.0 &> /dev/null
        then
            ## load correct R and gcc versions
            module unload R
            module load R/4.0.2
            module unload gcc
            module load gcc/9.2.0
        fi
    fi

    ## find R installation
    unset R_HOME

    echo "USING LIBRARIES: \$(Rscript -e 'print(.libPaths())')"

    export jabPath=\$(Rscript -e 'cat(suppressWarnings(find.package("JaBbA")))')
    export jba=\${jabPath}/extdata/jba
    echo \$jba
    set +x


    export cmd="Rscript \$jba $junction $cov_rds    \\
    $j_supp                                         \\
    $het_pileups_wgs                                \\
    --purity				$purity                 \\
    --ploidy				$ploidy                 \\
    $cbs_seg_rds                                    \\
    $cbs_nseg_rds                                   \\
    --blacklist.junctions   $blacklist_junctions    \\
    $geno_switch                                    \\
    --indel					$indel                  \\
    --tfield				$tfield                 \\
    --iterate				$iter                   \\
    --rescue.window			$rescue_window          \\
    --rescue.all			$rescue_all             \\
    --nudgebalanced			$nudgebalanced          \\
    --edgenudge				$edgenudge              \\
    $strict_switch                                  \\
    $allin_switch                                   \\
    --field					$field                  \\
    --maxna					$maxna                  \\
    --blacklist.coverage	$blacklist_coverage     \\
    --ppmethod				$pp_method              \\
    --cnsignif				$cnsignif               \\
    --slack					$slack                  \\
    $linear_switch                                  \\
    --tilim					$tilim                  \\
    --epgap					$epgap                  \\
    --name                  ${meta.patient}         \\
    --cores                 $task.cpus              \\
    --fix.thres				$fix_thres              \\
    --lp					$lp                     \\
    --ism					$ism                    \\
    --filter_loose			$filter_loose           \\
    --gurobi				$gurobi                 \\
    --mem                   $trelim_mem             \\
    $verbose_switch                                 \\
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        JaBbA: ${VERSION}
    END_VERSIONS

    { echo "Running:" && echo "\$(echo \$cmd)" && echo && eval \$cmd; }
    cmdsig=\$?
    if [ "\$cmdsig" = 0 ]; then
        echo "Finish!"
    else
        echo "Broke!"
        exit \$cmdsig
    fi

    ## exit 0
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch jabba.simple.rds
    touch jabba.simple.gg.rds
    touch jabba.simple.vcf
    touch jabba.raw.rds
    touch opt.report.rds
    touch jabba.seg
    touch karyograph.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        JaBbA: ${VERSION}
    END_VERSIONS
    """
}

process COERCE_SEQNAMES {

    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/jabba:0.0.1':
        'mskilab/jabba:0.0.1' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("coerced_chr*"), emit: file, optional: true

    script:
    """
    #!/usr/bin/env Rscript

    fn <- "${file}"
    outputfn <- "coerced_chr_${file.name}"

    if(grepl('.rds', "${file.name}")){
        library(GenomicRanges)
        data <- readRDS(fn)
        seqlevels(data, pruning.mode = "coarse") <- gsub("chr","",seqlevels(data))
        saveRDS(data, file = outputfn)
    } else if (grepl('.vcf|.vcf.gz|.vcf.bgz', "${file.name}")) {
        library(VariantAnnotation)
        data <- readVcf(fn)
        ##seqlevelsStyle(data) <- 'NCBI'
        seqlevels(data) <- sub("^chr", "", seqlevels(data))
        header = header(data)
        rownames(header@header\$contig) = sub("^chr", "", rownames(header@header\$contig))
        header(data) <- header
        data@fixed\$ALT <- lapply(data@fixed\$ALT, function(x) gsub("chr", "", x))
        writeVcf(data, file = outputfn)
    } else {
        data <- read.table(fn, header=T)
        data[[1]] <- gsub("chr","",data[[1]])
        write.table(data, file = outputfn, sep = "\\t", row.names = F, quote = F)
    }
    """
}

process RETIER_WHITELIST_JUNCTIONS___OLD {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/jabba:0.0.1':
        'mskilab/jabba:0.0.1' }"

    input:
    tuple val(meta), path(junctions)
    val(tfield)
    path(whitelist_genes)

    output:
    tuple val(meta), path("*___tiered.rds"), emit: retiered_junctions, optional: true

    script:
    """
    #!/usr/bin/env Rscript

	options(error = function() {traceback(2); quit(save = "no", status = 1)})

    library(gUtils)
    library(dplyr)
    library(VariantAnnotation)
    library(gGnome)

    # Load the whitelist genes
    heme_gen = readRDS("${whitelist_genes}")

    # Define the path to the junctions file
    jpath = "${junctions}"
    jpath_tiered = glue::glue('{tools::file_path_sans_ext(jpath)}___tiered.rds')

    # Read the VCF file
	is_character = is.character(jpath)
	is_len_one = NROW(jpath) == 1
	is_na = is_len_one && (is.na(jpath) || jpath %in% c("NA", base::nullfile()))
	is_possible_path = is_character && is_len_one && !is_na
  	is_existent_path = is_possible_path && file.exists(jpath)
  	is_rds = is_possible_path && grepl(".rds\$", jpath)
	is_vcf = is_possible_path && grepl(".vcf(.bgz|.gz){0,}\$", jpath)

	if (is_existent_path && is_rds) {
		ra.all = readRDS(jpath)
	} else if (is_existent_path && is_vcf) {
		ra.all = gGnome:::read.juncs(jpath)
	} else if (!is_existent_path) {
		stop("jpath does not exist or is invalid path: ", jpath)
	}

	is_properly_formatted_grangeslist = (
	    inherits(ra.all, "GRangesList")
		&& (
			all(S4Vectors::elementNROWS(ra.all) == 2)
			|| (NROW(ra.all) == 0)
		)
	)

	if (!is_properly_formatted_grangeslist) {
		stop("Improperly formatted junctions")
	}

    # Important part is below
    mcols_ra.all = mcols(ra.all)
    mcols_ra.all[["${tfield}"]] = rep_len(2, NROW(mcols_ra.all))
    ix = unique(mcols(grl.unlist(ra.all) %&% heme_gen)[["grl.ix"]])

    if (NROW(ix)) {
      cat("Whitelisted junctions overlapped with provided junctions. Retiering...\n")
      mcols_ra.all[["${tfield}"]][ix] = 1
      mcols(ra.all) = mcols_ra.all

      # Output the path of the saved file
      cat("Retiered junctions saved to:", jpath_tiered, "\n")
    } else {
      cat("No whitelisted junctions overlapped with provided junctions.\n")
    }

    saveRDS(ra.all, jpath_tiered)
    """
    stub:

    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ___tiered.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Retier Junctions: ${VERSION}
    END_VERSIONS
    """
}


process RETIER_WHITELIST_JUNCTIONS___DEV_COPY {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.8-rcpp':
        'mskilab/unified:0.0.8-rcpp' }"

    input:
    tuple val(meta), path(junctions), path(junctions_raw)
    val(tfield)
    path(whitelist_genes)

    output:
    tuple val(meta), path("*___tiered.rds"), emit: retiered_junctions, optional: true

    script:
    """
    #!/usr/bin/env Rscript

	options(error = function() {traceback(2); quit(save = "no", status = 1)})

    library(gUtils)
    library(dplyr)
    library(VariantAnnotation)
    library(gGnome)

    # Load the whitelist genes
    heme_gen = readRDS("${whitelist_genes}")

    heme_gen = IRanges::resize(
        heme_gen,
        width = as.integer(width(heme_gen)) + 150000,
        fix = "end",
        ignore.strand = FALSE
    )

    # Define the path to the junctions file
    jpath = "${junctions}"
    jpath_tiered = glue::glue('{tools::file_path_sans_ext(jpath)}___tiered.rds')

    # Read the VCF file
	is_character = is.character(jpath)
	is_len_one = NROW(jpath) == 1
	is_na = is_len_one && (is.na(jpath) || jpath %in% c("NA", base::nullfile()))
	is_possible_path = is_character && is_len_one && !is_na
  	is_existent_path = is_possible_path && file.exists(jpath)
  	is_rds = is_possible_path && grepl(".rds\$", jpath)
	is_vcf = is_possible_path && grepl(".vcf(.bgz|.gz){0,}\$", jpath)

	if (is_existent_path && is_rds) {
		ra.all = readRDS(jpath)
	} else if (is_existent_path && is_vcf) {
		ra.all = gGnome:::read.juncs(jpath)
	} else if (!is_existent_path) {
		stop("jpath does not exist or is invalid path: ", jpath)
	}

	is_properly_formatted_grangeslist = (
	    inherits(ra.all, "GRangesList")
		&& (
			all(S4Vectors::elementNROWS(ra.all) == 2)
			|| (NROW(ra.all) == 0)
		)
	)

	if (!is_properly_formatted_grangeslist) {
		stop("Improperly formatted junctions")
	}

    # Important part is below
    mcols_ra.all = mcols(ra.all)
    default_tier = rep_len(2, NROW(mcols_ra.all))
    default_tier[mcols_ra.all[["FILTER"]] != "PASS"] = 3
    mcols_ra.all[["${tfield}"]] = default_tier

    ix = unique(mcols(grl.unlist(ra.all) %&% heme_gen)[["grl.ix"]])

    if (NROW(ix)) {
      cat("Whitelisted junctions overlapped with provided junctions. Retiering...\n")
      mcols_ra.all[["${tfield}"]][ix] = 1

      # Output the path of the saved file
      cat("Retiered junctions saved to:", jpath_tiered, "\n")
    } else {
      cat("No whitelisted junctions overlapped with provided junctions.\n")
    }
	mcols(ra.all) = mcols_ra.all


    # Define the path to the raw junctions file
    jpath = "${junctions_raw}"

    # Read the VCF / RDS file
    is_character = is.character(jpath)
    is_len_one = NROW(jpath) == 1
    is_na = is_len_one && (is.na(jpath) || jpath %in% c("NA", base::nullfile()))
    is_possible_path = is_character && is_len_one && !is_na
    is_existent_path = is_possible_path && file.exists(jpath)
    is_rds = is_possible_path && grepl(".rds\$", jpath)
    is_vcf = is_possible_path && grepl(".vcf(.bgz|.gz){0,}\$", jpath)

    if (is_existent_path && is_rds) {
        ra.raw = readRDS(jpath)
    } else if (is_existent_path && is_vcf) {
        ra.raw = gGnome:::read.juncs(jpath)
    } else if (!is_existent_path) {
        stop("junctions file does not exist or is invalid path: ", jpath)
    }

    is_properly_formatted_grangeslist = (
        inherits(ra.raw, "GRangesList")
        && (
            all(S4Vectors::elementNROWS(ra.raw) == 2)
            || (NROW(ra.raw) == 0)
        )
    )

    if (!is_properly_formatted_grangeslist) {
        stop("Improperly formatted junctions")
    }

    ra.raw.reciprocals = gGnome::get_reciprocal_pairs(ra.raw)

    num_reciprocals_rescue = NROW(ra.raw.reciprocals)
    is_reciprocal_present = num_reciprocals_rescue > 0

    ra.out = ra.all
    if (is_reciprocal_present) {
        mcols_ra.raw = mcols(ra.raw.reciprocals)
        mcols_ra.raw[["${tfield}"]] = 3
        ix_raw = unique(mcols(grl.unlist(ra.raw.reciprocals) %&% heme_gen)[["grl.ix"]])
        is_reciprocal_in_rescue = NROW(ix_raw) > 0
        if (is_reciprocal_in_rescue) {
            mcols_ra.raw[ix_raw, "${tfield}"] = 1
        }
        mcols(ra.raw.reciprocals) = mcols_ra.raw
        ra.all = gUtils::gr.fix(ra.all, ra.raw.reciprocals)
        ra.raw.reciprocals = gUtils::gr.fix(ra.raw.reciprocals, ra.all)
        ra.out = c(
            ra.raw.reciprocals,
            ra.all
        )
    }

    is_ra_duplicated = gGnome::fra.duplicated(ra.out, pad = 1000)
    ra.out = ra.out[!is_ra_duplicated]

    saveRDS(ra.out, jpath_tiered)
    """
    stub:

    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ___tiered.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Retier Junctions: ${VERSION}
    END_VERSIONS
    """
}


process RETIER_WHITELIST_JUNCTIONS {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.15':
        'mskilab/unified:0.0.15' }"

    input:
    tuple val(meta), path(junctions)
    val(tfield)
    path(whitelist_genes)

    output:
    tuple val(meta), path("*___tiered.rds"), emit: retiered_junctions, optional: true

    script:
    """
    #!/usr/bin/env Rscript

	options(error = function() {traceback(2); quit(save = "no", status = 1)})

    library(gUtils)
    library(dplyr)
    library(VariantAnnotation)
    library(gGnome)

    # Load the whitelist genes
    heme_gen = readRDS("${whitelist_genes}")

    heme_gen = IRanges::resize(
        heme_gen,
        width = as.integer(width(heme_gen)) + 150000,
        fix = "end",
        ignore.strand = FALSE
    )

    # Define the path to the junctions file
    jpath = "${junctions}"
    jpath_tiered = glue::glue('{tools::file_path_sans_ext(jpath)}___tiered.rds')

    # Read the VCF file
	is_character = is.character(jpath)
	is_len_one = NROW(jpath) == 1
	is_na = is_len_one && (is.na(jpath) || jpath %in% c("NA", base::nullfile()))
	is_possible_path = is_character && is_len_one && !is_na
  	is_existent_path = is_possible_path && file.exists(jpath)
  	is_rds = is_possible_path && grepl(".rds\$", jpath)
	is_vcf = is_possible_path && grepl(".vcf(.bgz|.gz){0,}\$", jpath)

	if (is_existent_path && is_rds) {
		ra.all = readRDS(jpath)
	} else if (is_existent_path && is_vcf) {
		ra.all = gGnome:::read.juncs(jpath)
	} else if (!is_existent_path) {
		stop("jpath does not exist or is invalid path: ", jpath)
	}

	is_properly_formatted_grangeslist = (
	    inherits(ra.all, "GRangesList")
		&& (
			all(S4Vectors::elementNROWS(ra.all) == 2)
			|| (NROW(ra.all) == 0)
		)
	)

	if (!is_properly_formatted_grangeslist) {
		stop("Improperly formatted junctions")
	}

    # Important part is below
    mcols_ra.all = mcols(ra.all)
    default_tier = rep_len(2, NROW(mcols_ra.all))
    default_tier[mcols_ra.all[["FILTER"]] != "PASS"] = 3
    mcols_ra.all[["${tfield}"]] = default_tier

    mcols(ra.all) = mcols_ra.all
    ra.all = ra.all[order(mcols_ra.all[["${tfield}"]], decreasing = FALSE)]
    mcols_ra.all = mcols(ra.all)
    grpiv = gUtils::grl.pivot(ra.all)
    bp1 = grpiv[[1]]
    bp2 = grpiv[[2]]

    bp1_genes = gUtils::gr.val(bp1, heme_gen, val = "gene_name", sep = "___SEP___")
    bp2_genes = gUtils::gr.val(bp2, heme_gen, val = "gene_name", sep = "___SEP___")

    bp1_genes\$genes = IRanges::CharacterList(strsplit(bp1_genes\$gene_name, "___SEP___"))
    bp2_genes\$genes = IRanges::CharacterList(strsplit(bp2_genes\$gene_name, "___SEP___"))

    is_bp1bp2_genes_different = any(bp1_genes\$genes != bp2_genes\$genes)
    is_bp1_in_a_gene = elementNROWS(bp1_genes\$genes) > 0
    is_bp2_in_a_gene = elementNROWS(bp2_genes\$genes) > 0

    is_bp1bp2_different = (
    is_bp1bp2_genes_different
    | (is_bp1_in_a_gene & !is_bp2_in_a_gene)
    | (!is_bp1_in_a_gene & is_bp2_in_a_gene)
    )


    ra.raw.reciprocals.lst = gGnome::get_reciprocal_pairs(
        ra.all[is_bp1bp2_different]
        ,
        distance_pad = 1e4
        ,
        return_inputs = TRUE  
    )

    ra.raw.reciprocals = ra.raw.reciprocals.lst\$jun_recip_pairs
    mcols_ra.raw.reciprocals = mcols(ra.raw.reciprocals)
    nr_ra.raw.reciprocals = NROW(ra.raw.reciprocals)
    mcols_ra.raw.reciprocals[["${tfield}"]] = rep_len(1, nr_ra.raw.reciprocals)
    mcols(ra.raw.reciprocals) = mcols_ra.raw.reciprocals

    ra.raw.reciprocals = gr.fix(ra.raw.reciprocals, ra.all)
    ra.all = gr.fix(ra.all, ra.raw.reciprocals)
    ra.out = c(ra.raw.reciprocals, ra.all)

    if (nr_ra.raw.reciprocals > 0) {
        cat("Whitelisted", nr_ra.raw.reciprocals, "reciprocal junctions via retiering '${tfield}' field to 1.\n")
    }

    ra.out = ra.out[!gGnome::fra.duplicated(ra.out, pad = 1000)]


    saveRDS(ra.out, jpath_tiered)
    """
    stub:

    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ___tiered.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Retier Junctions: ${VERSION}
    END_VERSIONS
    """
}

