process EVENTS {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:initial':
        'mskilab/unified:initial' }"

    input:
    tuple val(meta), path(gGraph)
    path(ref)

    output:
    tuple val(meta), path("*complex.rds") , emit: events_output
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def id          = "${meta.sample}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """

    export HOME=/root && \
    set +u  # Disable unbound variable check
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate pact

    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/Events.R")

    Rscript \$RSCRIPT_PATH \\
	--id $id \\
	--gGraph $gGraph \\
	--ref $ref \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Events: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch dummy_output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        events: ${VERSION}
    END_VERSIONS
    """
}
