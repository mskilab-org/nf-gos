process BAMMETRICS {
    tag "$meta.id"
    label 'process_gpu' // Label for GPU processes

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1' :
        'nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bammetrics.txt"), emit: metrics
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metrics_file = "${prefix}.bammetrics.txt"

    """
    pbrun bammetrics \\
        --ref ${fasta} \\
        --bam ${bam} \\
        --out-metrics-file ${metrics_file} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun --version 2>&1 | grep "Clara Parabricks Version" | sed 's/Clara Parabricks Version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metrics_file = "${prefix}.bammetrics.txt"
    """
    touch ${metrics_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(echo "stub_version")
    END_VERSIONS
    """
}
