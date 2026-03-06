process PARABRICKS_BAMMETRICS {
    tag "$meta.id"
    label 'process_high' // Label for GPU processes

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1' :
        'nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1' }"
	containerOptions "${ workflow.containerEngine == "singularity" ? "--nv --bind ${NEXTFLOW_C_DIR}/memcap.so:/memcap.so": ( workflow.containerEngine == "docker" ? '--gpus all': null ) }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path  intervallist

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metrics_file = "${prefix}.bammetrics.txt"
	def interval  = intervallist ? "--interval ${intervallist}" : ''
    def process_mem_limit_mb = (task.memory.toMega() * 0.97).toInteger() // Leave some headroom for the process to avoid OOM killer. This is the value that will be passed to fq2bam via --memory-limit, which is different from the MEMLIMIT env variable that is used by memcap.so to limit GPU memory usage.

    """
    export TMPDIR=/tmp
    tdir=\$(mktemp -d)

    export LD_PRELOAD=/memcap.so
    export MEMLIMIT_MB=${process_mem_limit_mb}

    pbrun bammetrics \\
        --ref ${fasta} \\
        --bam ${bam} \\
        --out-metrics-file ${prefix}.bammetrics.coverage_metrics \\
		--coverage-cap 10000 \\
		--tmp-dir \$tdir \\
        $args \\
		$interval

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
