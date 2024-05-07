process JUNCTION_FILTER {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/junction_filter:latest':
        'mskilab/junction_filter:latest' }"

    input:
    tuple val(meta), path(filtered_sv_vcf)
    path(junction_pon)
    path(gnomAD_sv_db)
    val(padding)

    output:
    tuple val(meta), path("somatic.filtered.gnoMAD.sv.rds")                , emit: final_filtered_sv_rds,      optional:true
    tuple val(meta), path("somatic.filtered.sv.rds")                       , emit: pon_filtered_sv_rds,        optional: false
    path "versions.yml"                                                    , emit: versions,                   optional:false

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    #!/bin/bash

    set -o allexport
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/junction_filter.R")

    Rscript \$RSCRIPT_PATH           \\
    --sv        ${filtered_sv_vcf}   \\
    --pon       ${junction_pon}      \\
    --gnomAD    ${gnomAD_sv_db}      \\
    --padding   ${padding}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Junction_Filter: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch somatic.filtered.gnoMAD.sv.rds somatic.filtered.sv.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Junction_Filter: ${VERSION}
    END_VERSIONS
    """
}