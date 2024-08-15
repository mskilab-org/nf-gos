process ONENESS_TWONESS {

    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/hrdetect:0.0.4':
        'mskilab/hrdetect:0.0.4' }"

    input:
    tuple val(meta), path(events_output), path(homeology), path(homeology_stats), path(hrdetect_results)
    path(model)

    output:
    tuple val(meta), path("*oneness_twoness_results.rds")      , emit: oneness_twoness, optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/Oneness_Twoness.R")

    Rscript \$RSCRIPT_PATH \\
    --complex $events_output \\
    --homeology $homeology \\
    --homeology_stats $homeology_stats \\
    --hrdetect_results $hrdetect_results \\
    --model $model
    --libdir ./ \\

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

process HOMEOLOGY {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/hrdetect:0.0.4':
        'mskilab/hrdetect:0.0.4' }"

    input:
    tuple val(meta), path(events_output)
    val(width)
    val(pad)
    val(thresh)
    val(stride)
    path(genome)
    val(flip)
    val(bidirectional)
    val(annotate)
    val(savegMatrix)

    output:
    tuple val(meta), path("junctions.txt")      , emit: homeology, optional: true
    tuple val(meta), path("rawstats.txt")       , emit: homeology_stats, optional: true
    tuple val(meta), path("gMatrixList.rds")     , emit: homeology_gm, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/Homeology.R")

    Rscript \$RSCRIPT_PATH \\
    --junctions $events_output \\
    --width $width \\
    --pad $pad \\
    --thresh $thresh \\
    --stride $stride \\
    --genome $genome \\
    --cores ${task.cpus} \\
    --flip $flip \\
    --bidirectional $bidirectional \\
    --annotate $annotate \\
    --savegMatrix $savegMatrix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homeology: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch junctions.txt
    touch rawstats.txt
    touch gMatrixList.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homeology: ${VERSION}
    END_VERSIONS
    """
}
