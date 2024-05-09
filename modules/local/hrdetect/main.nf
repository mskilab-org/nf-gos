process HRDETECT {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/hrdetect:latest':
        'mskilab/hrdetect:latest' }"

    input:
    tuple val(meta), path(junctions), path(hets), path(snv_somatic), path(jabba_rds)
    val(mask) // declared as val to allow for NULL value, but this is actually a path
    path(ref_fasta)
    val(genome_version)

    output:
    tuple val(meta), path("*hrdetect_results.rds")      , emit: hrdetect_rds, optional: true
    tuple val(meta), path("*hrdetect_output.txt")       , emit: hrdetect_txt, optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def hrdetect_mask = mask ?: "--mask '/dev/null'"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/HRDetect.R")

    Rscript \$RSCRIPT_PATH \\
	--sv $junctions \\
	--hets $hets \\
	--snv $snv_somatic \\
	--jabba $jabba_rds \\
    ${hrdetect_mask} \\
    --ref $ref_fasta \\
    --genome $genome_version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hrdetect: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch hrdetect_output.txt
    touch hrdetect_results.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hrdetect: ${VERSION}
    END_VERSIONS
    """
}
