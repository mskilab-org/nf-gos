process SNV_MULTIPLICITY {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/snv_multiplicity:0.0.3':
        'mskilab/snv_multiplicity:0.0.3' }"

    input:
    tuple val(meta), path(somatic_snv, stageAs: "somatic_snv.vcf"), path(germline_snv, stageAs: "germline_snv.vcf"), path(jabba_gg)
    path(ref)
    path(ref_fai)

    output:
    tuple val(meta), path('*est_snv_cn_somatic.rds'), emit: snv_multiplicity_rds
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def DOWNSAMPLE_JAR = "${baseDir}/bin/downsamplevcf.jar"
    def SNPSIFT_JAR = "${baseDir}/bin/SnpSift.jar"
    def VCFEFF_PERL = "${baseDir}/bin/vcfEffOnePerLine.pl"

    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/snv_multiplicity3.R")

    Rscript \$RSCRIPT_PATH \\
        --somatic_snv ${somatic_snv} \\
        --germline_snv ${germline_snv} \\
        --fasta ${ref} \\
        --jabba ${jabba_gg} \\
        --downsample_jar ${DOWNSAMPLE_JAR} \\
        --snpsift_jar ${SNPSIFT_JAR} \\
        --vcfeff_perl ${VCFEFF_PERL} \\
        --cores ${task.cpus} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snv_multiplicity: ${VERSION}
    END_VERSIONS
    """

    stub:
    """
    touch est_snv_cn_somatic.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snv_multiplicity: ${VERSION}
    END_VERSIONS
    """
}
