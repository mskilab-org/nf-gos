process ONENESS_TWONESS {

    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/hrdetect:0.0.4':
        'mskilab/hrdetect:0.0.4' }"

    input:
    tuple val(meta), path(events_output), path(hrdetect_results)
    path(genome_fasta)
    path(genome_fai)
    path(model)

    output:
    tuple val(meta), path("*oneness_twoness_results.rds"), emit: oneness_twoness_results, optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/OnenessTwoness.R")

    Rscript \$RSCRIPT_PATH \\
    --complex $events_output \\
    --hrdetect_results $hrdetect_results \\
    --genome $genome_fasta \\
    --model $model \\
    --cores ${task.cpus} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oneness_twoness: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch oneness_twoness_results.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oneness_twoness: ${VERSION}
    END_VERSIONS
    """
}
