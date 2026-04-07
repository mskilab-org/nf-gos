process GRIPSS_SOMATIC_FILTER {
    tag "$meta.id"
    label 'process_medium'

    // using events container since the dependencies are the same
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.19-fusions':
        'mskilab/unified:0.0.19-fusions' }"

    input:
    tuple val(meta), path(gridss_output), path(gridss_output_tbi)
    path(pon_gridss_bedpe_svs)
    path(pon_gridss_bed_breakends)
    path(pon_gridss_known_hotspots_bedpe)
    path(fasta)
    path(fasta_fai)
    val(pon_gridss_ref_genome_version)

    output:
    tuple val(meta), path("*.gripss.filtered.vcf.gz"), path("*.gripss.filtered.vcf.gz.tbi"), emit: somatic_high_vcf,          optional:false
    tuple val(meta), path("*.gripss.vcf.gz"), path("*.gripss.vcf.gz.tbi"), emit: somatic_all_vcf, optional:false
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def VERSION       = '2.3.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def jvmheap_mem = (task.memory.toGiga() * 0.9).toInteger() // 
    """

    java -Xmx${jvmheap_mem}g -jar \${NEXTFLOW_BIN_DIR}/jar/gripss_v2.3.4.jar \\
         -log_debug \\
         -log_level WARN \\
         -pon_sv_file ${pon_gridss_bedpe_svs} \\
         -pon_sgl_file ${pon_gridss_bed_breakends} \\
         -known_hotspot_file ${pon_gridss_known_hotspots_bedpe} \\
         -vcf ${gridss_output} \\
         -ref_genome ${fasta} \\
         -ref_genome_version ${pon_gridss_ref_genome_version} \\
         -output_dir ./ \\
         -sample ${meta.sample} \\
         ${reference_arg} \\
         ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: ${VERSION}
    END_VERSIONS

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.high_confidence_somatic.vcf.bgz
    touch ${prefix}.high_and_low_confidence_somatic.vcf.bgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}

process GRIDSS_SOMATIC_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"


    input:
    tuple val(meta), path(gridss_output), path(gridss_output_tbi)
    path(pondir_gridss)

    output:
    tuple val(meta), path("*.high_confidence_somatic.vcf.bgz"), path("*.high_confidence_somatic.vcf.bgz.tbi")             , emit: somatic_high_vcf,          optional:true
    tuple val(meta), path("*.high_and_low_confidence_somatic.vcf.bgz")      , emit: somatic_all_vcf,           optional:true
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def VERSION       = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def output1       = "${meta.id}.high_confidence_somatic.vcf"
    def output2       = "${meta.id}.high_and_low_confidence_somatic.vcf"
    //def pon           = pondir_gridss ? "cp -s ${pondir_gridss}/* ." : ""
    //def scriptDir     = "dirname \$(which gridss_somatic_filter)".execute().text.trim()
    """

    gridss_somatic_filter \\
    --pondir ${pondir_gridss} \\
    --input ${gridss_output} \\
    --output $output1 \\
    --fulloutput $output2 \\
    -n 1 \\
    -t 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.high_confidence_somatic.vcf.bgz
    touch ${prefix}.high_and_low_confidence_somatic.vcf.bgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}

