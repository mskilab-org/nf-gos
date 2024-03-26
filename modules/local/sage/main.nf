process SAGE {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.4--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sage:3.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(tumor_bam_wgs, stageAs: "tumor.bam"), path(tumor_bai, stageAs: "tumor.bam.bai"), path(normal_bam_wgs, stageAs: "normal.bam"), path(normal_bai, stageAs: "normal.bam.bai")
    path(ref)
    path(ref_fai)
    val(ref_genome_version)
    val(ensembl_data_dir) // set as val but is actually a path to a directory
    path(somatic_hotspots)
    path(panel_bed)
    path(high_confidence_bed)

    output:
    tuple val(meta), path("*sage.vcf")                     , emit: sage_vcf, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def output      = "${meta.id}_sage.vcf"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def reference_bam_arg = normal_bam_wgs ? "-reference_bam ${normal_bam_wgs}" : ''
    """

    mkdir -p somatic/

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        ${reference_bam_arg} \\
        -tumor_bam ${tumor_bam_wgs} \\
        -ref_genome ${ref} \\
        -ref_genome_version ${ref_genome_version} \\
        -hotspots ${somatic_hotspots} \\
        -panel_bed ${panel_bed} \\
        -high_confidence_bed ${high_confidence_bed} \\
        -ensembl_data_dir ${ensembl_data_dir} \\
        -write_bqr_data \\
        -write_bqr_plot \\
        -threads ${task.cpus} \\
        -output_vcf ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hetpileups: ${VERSION}
    END_VERSIONS
    """
}
