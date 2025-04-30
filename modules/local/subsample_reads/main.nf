process SAMTOOLS_SUBSAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(aligned_input), path(aligned_index)
    tuple val(meta2), path(fasta)
    path(qname)

    output:
    tuple val(meta), path("*.subsample.bam"), path("*.subsample.bam.bai"),  emit: bam_subsampled,     optional: false
    tuple val(meta), path("*.subsample.cram"), path("*.subsample.cram.crai"), emit: cram_subsampled,    optional: true
    tuple val(meta), path("*.subsample.sam"),  val(null), emit: sam_subsampled,     optional: true
    // tuple val(meta), path("*.bai"),  emit: bai_subsampled,     optional: true
    // tuple val(meta), path("*.csi"),  emit: csi_subsampled,     optional: true
    // tuple val(meta), path("*.crai"), emit: crai_subsampled,    optional: true
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def readnames = qname ? "--qname-file ${qname}": ""
    def file_type = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt bam") ? "bam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    aligned_input.getExtension()
    if ("$aligned_input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
	def num_subsampled_reads = params.num_subsampled_reads ?: 800_000_000 * 2 // 800 million read pairs ~ 80X WGS
    """

	FRAC=\$( samtools idxstats \\
		${aligned_input} | \\
		cut -f3 | \\
		awk 'BEGIN {total=0} {total += \$1} END {frac=${num_subsampled_reads}/total; if (frac > 1) {printf "%.1f\\n", 1.0} else {print frac}}' )

    samtools \\
        view \\
		-s \$FRAC \\
        --threads ${task.cpus-1} \\
        ${reference} \\
        `# ${readnames}` \\
        $args \\
        -o ${prefix}.subsample.${file_type} \\
        $aligned_input \\
        $args2
	
	samtools index ${prefix}.subsample.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.subsample.bam
    touch ${prefix}.subsample.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
