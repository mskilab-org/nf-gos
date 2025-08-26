process SNV_MULTIPLICITY {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/snv_multiplicity:0.0.5':
        'mskilab/snv_multiplicity:0.0.5' }"

    input:
    tuple val(meta), path(somatic_snv, stageAs: "somatic_snv.vcf"), path(germline_snv, stageAs: "germline_snv.vcf"), path(jabba_gg), path(hets, stageAs: "hets_sites.txt"), path(dryclean_cov, stageAs: "dryclean_cov.rds")

    output:
    tuple val(meta), path('*est_snv_cn_somatic.rds'),         emit: snv_multiplicity_rds
    tuple val(meta), path('*est_snv_cn_germline.rds'),        emit: snv_multiplicity_germline_rds, optional: true
    tuple val(meta), path('*est_snv_cn_hets.rds'),            emit: snv_multiplicity_hets_rds
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def germline_flag = germline_snv ? "--germline_snv ${germline_snv}" : "--germline_snv /dev/null"
    def normal_flag = meta.normal_id ? "--normal_name ${meta.normal_id}" : ""
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export RSCRIPT_PATH=\$(echo "\${NEXTFLOW_PROJECT_DIR}/bin/snv_multiplicity.R")

    Rscript \$RSCRIPT_PATH \\
        --somatic_snv ${somatic_snv} \\
        ${germline_flag} \\
        --het_pileups_wgs ${hets} \\
        --tumor_dryclean ${dryclean_cov} \\
        --jabba ${jabba_gg} \\
        --tumor_name ${meta.tumor_id} \\
        ${normal_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snv_multiplicity: ${VERSION}
    END_VERSIONS
    """

    stub:
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch est_snv_cn_somatic.rds
    touch est_snv_cn_germline.rds
    touch est_snv_cn_hets.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snv_multiplicity: ${VERSION}
    END_VERSIONS
    """
}
