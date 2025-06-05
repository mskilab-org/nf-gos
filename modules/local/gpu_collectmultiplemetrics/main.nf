process GPU_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_gpu'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1' :
        'nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path(out_dir), emit: metrics_dir
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_gpu_collectmultiplemetrics_qc"

    """
    pbrun collectmultiplemetrics \\
        --ref ${fasta} \\
        --bam ${bam} \\
        --out-qc-metrics-dir ${out_dir} \\
        --gen-all-metrics \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun --version 2>&1 | grep "Clara Parabricks Version" | sed 's/Clara Parabricks Version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}.gpu_collectmultiplemetrics.qc"
    """
    mkdir -p ${out_dir}
    touch ${out_dir}/placeholder_metric.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(echo "stub_version")
    END_VERSIONS
    """
}
