process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input_files, stageAs: "?/*")
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}.bam") , optional:true, emit: bam
    tuple val(meta), path("${prefix}.cram"), optional:true, emit: cram
    tuple val(meta), path("*.csi")         , optional:true, emit: csi
    path  "versions.yml"                                  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    def reference = fasta ? "--reference ${fasta}" : ""
    """
	
	are_all_inputs_complete=true
	for input_file in $input_files; do
		is_input_complete=\$(samtools quickcheck \$input_file && echo 'true' || echo 'false')
		! \$is_input_complete && printf 'TRUNCATED INPUT: %s\n' \$input_file
		are_all_inputs_complete=\$({ \$are_all_inputs_complete && \$is_input_complete; } && echo 'true' || echo 'false')
	done

	if ! \$are_all_inputs_complete; then
		exit 1
	fi

	is_output_complete=\$(samtools quickcheck ${prefix}.${file_type} && echo 'true' || echo 'false')
	is_output_existent=\$([ -e ${prefix}.${file_type} ] && echo 'true' || echo 'false')
	is_output_done_but_completed=\$( { \$is_output_existent &&   \$is_output_complete; } && echo 'true' || echo 'false')
	is_output_done_but_truncated=\$( { \$is_output_existent && ! \$is_output_complete; } && echo 'true' || echo 'false')

	if \$is_output_done_but_truncated; then
		printf "Merged output truncated!!!\n%s:\n" "${prefix}.${file_type}"
		rm ${prefix}.${file_type}
	fi

	## FIXME: if samtools merge dies, the output should be deleted.
    samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        $args \\
        ${reference} \\
        ${prefix}.${file_type} \\
        $input_files
	
	is_output_complete=\$(samtools quickcheck ${prefix}.${file_type} && echo 'true' || echo 'false')

	if ! \$is_output_complete; then
		printf "Merged output is truncated!!!\n"
		rm ${prefix}.${file_type}
		exit 1
	fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    """
    touch ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
