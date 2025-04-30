process GATK4_ESTIMATELIBRARYCOMPLEXITY {
    tag "$meta.id"
    // label 'process_medium'
	label 'process_max'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.5.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(input)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path('*.metrics'), emit: metrics
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = input.collect(){"--INPUT $it"}.join(" ")
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
	def pixel_dist = params.optical_duplicate_pixel_distance ?: 100

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK EstimateLibraryComplexity] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        EstimateLibraryComplexity \\
        $input_list \\
        --OUTPUT ${prefix}.metrics \\
        $reference \\
        --TMP_DIR . \\
		--OPTICAL_DUPLICATE_PIXEL_DISTANCE $pixel_dist \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
