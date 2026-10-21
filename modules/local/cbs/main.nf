process CBS {

    tag "$meta.id"
    label 'process_medium'

    // TODO add fragcounter container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/cbs:latest':
        'mskilab/cbs:latest' }"

    input:
    path cbs_driver
    path cbs_coverage_helper
    tuple val(meta), path(tumor_dryclean_cov, stageAs: "tumor_drycleaned_cov.rds"), path(normal_dryclean_cov, stageAs: "normal_drycleaned_cov.rds")
    val(cnsignif)
    val(field)
    val(name)

    output:
    tuple val(meta), path("cov.rds")         , emit: cbs_cov_rds
    tuple val(meta), path("seg.rds")         , emit: cbs_seg_rds
    tuple val(meta), path("nseg.rds")        , emit: cbs_nseg_rds
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def normal_cov = normal_dryclean_cov ?: "NULL"
    def VERSION     = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    Rscript ${cbs_driver} \\
        -t ${tumor_dryclean_cov} \
        -n ${normal_cov} \
        --cnsignif ${cnsignif} \
        -m ${meta.id} \
        -f ${field}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cbs: ${VERSION}
    END_VERSIONS

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch cov.rds
    touch seg.rds
    touch nseg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cbs: ${VERSION}
    END_VERSIONS
    """

}

process VALIDATE_CBS_BUNDLE {
    tag "${meta.patient}/${meta.sample}"
    label 'process_low'
    publishDir enabled: false

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/jabba:0.0.10':
        'mskilab/jabba:0.0.10' }"

    input:
    path validator
    path cbs_coverage_helper
    tuple val(meta), val(row), path(cov, stageAs: 'cov.rds'), path(seg, stageAs: 'seg.rds'), path(nseg, stageAs: 'nseg.rds'), val(graph_names), path(graphs, stageAs: 'graphs/graph??.rds')
    val(mode)

    output:
    tuple val(meta), path(cov), path(seg), path(nseg), emit: bundle
    val(row), emit: validated_rows

    script:
    def graph_files = graphs instanceof List ? graphs : [graphs]
    def graph_args = graph_files.withIndex().collect { graph, index -> "'${graph_names[index]}=${graph}'" }.join(' ')
    """
    Rscript ${validator} '${mode}' '${cov}' '${seg ?: "NULL"}' '${nseg ?: "NULL"}' ${graph_args}
    """
}
