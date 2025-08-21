process FUSIONS {

    tag "$meta.id"
    label 'process_medium'

    // using events container since the dependencies are the same
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.10':
        'mskilab/unified:0.0.10' }"
    // 0.0.6 and above do not require conda sourcing

    input:
    tuple val(meta), path(gGraph), path(sv_vcf), path(sv_vcf_tbi)
    path(gencode)

    output:
    tuple val(meta), path("*fusions.rds") , emit: fusions_output
    tuple val(meta), path("*altedge.annotations.tsv") , emit: altedge_annotations
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def id          = "${meta.sample}"
    def junctions_flag = sv_vcf ? "--junctions ${sv_vcf}" : "--junctions ${gGraph}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    set +u  # Disable unbound variable check

    export RSCRIPT_PATH=\$(echo "\${NEXTFLOW_PROJECT_DIR}/bin/Fusions.R")

    Rscript \$RSCRIPT_PATH \\
	--id $id \\
    ${junctions_flag} \\
	--gencode $gencode \\
    --cores ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusions: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch dummy_output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusions: ${VERSION}
    END_VERSIONS
    """
}
