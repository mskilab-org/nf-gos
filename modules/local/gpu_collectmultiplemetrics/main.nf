process GPU_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1' :
        'nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1' }"
	containerOptions "${ workflow.containerEngine == "singularity" ? "--nv --bind ${NEXTFLOW_C_DIR}/memcap.so:/memcap.so": ( workflow.containerEngine == "docker" ? '--gpus all': null ) }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*gpu_collectmultiplemetrics_qc/**"), emit: metrics
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_gpu_collectmultiplemetrics_qc"
    def process_mem_limit_mb = (task.memory.toMega() * 0.97).toInteger() // Leave some headroom for the process to avoid OOM killer. This is the value that will be passed to fq2bam via --memory-limit, which is different from the MEMLIMIT env variable that is used by memcap.so to limit GPU memory usage.

    """

    export LD_PRELOAD=/memcap.so
    export MEMLIMIT_MB=${process_mem_limit_mb}

    pbrun collectmultiplemetrics \\
        --ref ${fasta} \\
        --bam ${bam} \\
        --out-qc-metrics-dir ${out_dir} \\
        --gen-all-metrics \\
		--tmp-dir ./ \\
		--num-gpus 1 \\
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
