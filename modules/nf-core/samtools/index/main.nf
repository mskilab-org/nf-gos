process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
	tuple val(meta), path("*.bam") , optional:true, emit: bam
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    tuple val(meta), path("*.crai"), optional:true, emit: crai
	tuple val(meta), path("*.cram") , optional:true, emit: cram
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
	def filename = input.getName()   // just the filename
	def ext = input.name.tokenize('.')[-1]
	def base = input.name - ".${ext}"
    """

	## ln -ns $input ${base}___indexed.${ext}
	## touch ${base}___indexed.${ext}

    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
		$input
        
		
	# ${base}___indexed.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
	touch ${base}___indexed.${ext}
    touch ${input}.bai
    touch ${input}.crai
    touch ${input}.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
