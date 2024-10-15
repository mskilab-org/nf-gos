process SNV_MULTIPLICITY {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/snv_multiplicity:0.0.3':
        'mskilab/snv_multiplicity:0.0.3' }"

    input:
    tuple val(meta), path(somatic_snv, stageAs: "somatic_snv.vcf"), path(germline_snv, stageAs: "germline_snv.vcf"), path(jabba_gg)

    output:
    tuple val(meta), path('*est_snv_cn_somatic.rds'), emit: snv_multiplicity_rds
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def tumor_name  = meta.tumor_id ? : ""
    def normal_name  = meta.normal_id ? : ""
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def SCRIPTS_DIR = "${baseDir}/bin/"

    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/snv_multiplicity3.R")

    Rscript \$RSCRIPT_PATH \\
        --somatic_snv ${somatic_snv} \\
        --germline_snv ${germline_snv} \\
        --jabba ${jabba_gg} \\
        --snpeff_path ${SCRIPTS_DIR} \\
        --tumor_name ${tumor_name} \\
        --normal_name ${normal_name} \\
        --cores ${task.cpus} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snv_multiplicity: ${VERSION}
    END_VERSIONS
    """

    stub:
    """
    touch est_snv_cn_somatic.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snv_multiplicity: ${VERSION}
    END_VERSIONS
    """
}
