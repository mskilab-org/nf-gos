process NON_INTEGER_BALANCE {

    tag "$meta.id"
    label 'process_medium'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     '/gpfs/commons/groups/imielinski_lab/home/sdider/Projects/nf-jabba/tests/test_runs/work/singularity/jabba_cplex_latest.sif':
    //     'mskilab/jabba:latest' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/jabba:latest':
        'mskilab/jabba:latest' }"

    input:
    tuple val(meta), path(jabba_gg), path(decomposed_cov), path(het_pileups_wgs)
    val(field)
    val(hets_thresh)
    path(mask)
    val(overwrite)
    val(lambda)
    val(allin)
    val(fix_thresh)
    val(nodebounds)
    val(ism)
    val(build)
    val(epgap)
    val(tilim)
    val(gurobi)
    path(fasta)     // path to decoy fasta
    path(fasta_fai)     // path to decoy fasta
    path(bwa_index)
    val(pad)

    output:
    tuple val(meta), path("non_integer.balanced.gg.rds")    , emit: non_integer_balance_balanced_gg, optional: true
    tuple val(meta), path("hets.gg.rds")                    , emit: non_integer_balance_hets_gg, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def id          = "${meta.sample}"
    def bwa = bwa_index ? "ln -nfs \$(readlink -f ${bwa_index})/* \$(dirname \$(readlink -f $fasta))/" : ""
    def hets_flag = het_pileups_wgs ? "--hets ${het_pileups_wgs}" : ""
    def hets_thres_flag = het_pileups_wgs ? "--hets_thresh ${hets_thresh}" : ""
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    ${bwa}

    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/non_integer_balance.R")

    Rscript \$RSCRIPT_PATH \\
        --id $id \\
        --jab $jabba_gg \\
        --cov $decomposed_cov \\
        --field $field \\
        ${hets_flag} \\
        ${hets_thres_flag} \\
        --mask $mask \\
        --overwrite $overwrite \\
        --lambda $lambda \\
        --allin $allin \\
        --fix_thresh $fix_thresh \\
        --nodebounds $nodebounds \\
        --ism $ism \\
        --build $build \\
        --epgap $epgap \\
        --tilim $tilim \\
        --gurobi $gurobi \\
        --fasta $fasta \\
        --pad $pad

    mv balanced.gg.rds non_integer.balanced.gg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        non_integer_balance: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch balanced.gg.rds hets.gg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        non_integer_balance: ${VERSION}
    END_VERSIONS
    """
}

process LP_PHASED_BALANCE {

    tag "$meta.id"
    label 'process_medium'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     '/gpfs/commons/groups/imielinski_lab/home/sdider/Projects/nf-jabba/tests/test_runs/work/singularity/jabba_cplex_latest.sif':
    //     'mskilab/jabba:latest' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/jabba:latest':
        'mskilab/jabba:latest' }"

    input:
    tuple val(meta), path(hets_gg, stageAs: "non_integer_balanced.gg.rds"), path(hets) // output from non_integer_balance, sites.txt from hetpileups
    val(lambda)
    val(cnloh)
    val(major)
    val(allin)
    val(marginal)
    val(from_maf)
    path(mask)
    val(ism)
    val(epgap)
    val(hets_thresh)
    val(min_bins)
    val(min_width)
    val(trelim)
    val(reward)
    val(nodefileind)
    val(tilim)

    output:
    tuple val(meta), path("lp_phased.balanced.gg.rds")                , emit: lp_phased_balance_balanced_gg, optional: true
    tuple val(meta), path("binstats.gg.rds")                , emit: lp_phased_balance_binstats_gg, optional: true
    tuple val(meta), path("unphased.gg.rds")                , emit: lp_phased_balance_unphased_allelic_gg, optional: true
    path "versions.yml"                                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def id   = "${meta.sample}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/lp_phased_balance.R")

    # Remove 'chr' from chromosome names in sites.txt (for hg38)
    awk 'BEGIN{OFS=" "} {gsub(/^chr/,"",\$1); print}' sites.txt > sites.tmp && mv sites.tmp sites.txt

    Rscript \$RSCRIPT_PATH \\
        --id $id \\
        --jab $hets_gg \\
        --hets $hets \\
        --lambda $lambda \\
        --cnloh $cnloh \\
        --major $major \\
        --allin $allin \\
        --marginal $marginal \\
        --from_maf $from_maf \\
        --mask $mask \\
        --ism $ism \\
        --epgap $epgap \\
        --hets_thresh $hets_thresh \\
        --min_bins $min_bins \\
        --min_width $min_width \\
        --trelim $trelim \\
        --reward $reward \\
        --nodefileind $nodefileind \\
        --tilim $tilim

    mv balanced.gg.rds lp_phased.balanced.gg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lp_phased_balance: ${VERSION}
    END_VERSIONS
    """
}
