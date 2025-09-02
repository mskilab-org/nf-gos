process CHIMERA_FILTER {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.12':
        'mskilab/unified:0.0.12' }"

    input:
	tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.ffpe_filtered.bam"), path("*.ffpe_filtered.bam.bai"), emit: bambai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_bam = bam.getName().replaceFirst(/\.bam$/, '.ffpe_filtered.bam')
    """
    python \${NEXTFLOW_BIN_DIR}/pysam_chimera_filter.py ${bam} ${out_bam}

    samtools index ${out_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ffpe_filtered.bam
    touch ${prefix}.ffpe_filtered.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
