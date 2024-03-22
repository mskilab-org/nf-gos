process SAGE {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/sage:latest':
        'mskilab/sage:latest' }"

    input:
    tuple val(meta), path(tumor_bam_wgs, stageAs: "tumor.bam"), path(normal_bam_wgs, stageAs: "normal.bam")
    path(ref)
    val(ref_genome_version)
    path(ensembl_data_dir)
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

    """

    export SAGE_JAR_PATH=\$(echo "${baseDir}/bin/sage_v3.2.2.jar")

    java -jar -Xmx32G ${SAGE_JAR_PATH},
         -tumor tumor \\
         -tumor_bam ${tumor_bam_wgs} \\
         -reference normal  \\
         -reference_bam ${normal_bam_wgs} \\
         -out ${output} \\
         -ref_genome ${ref} \\
         -ref_genome_version ${ref_genome_version} \\
         -ensembl_data_dir ${ensembl_data_dir} \\
         -hotspots ${somatic_hotspots} \\
         -panel_bed ${panel_bed} \\
         -high_confidence_bed ${high_confidence_bed} \\
         -validation_stringency LENIENT \\
         -threads $num.cores
        )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HetPileups: ${VERSION}
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
