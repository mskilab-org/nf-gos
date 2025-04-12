process HRDETECT {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/hrdetect:0.0.5':
        'mskilab/hrdetect:0.0.5' }"

    input:
    tuple val(meta), path(junctions), path(hets), path(snv_somatic), path(snv_somatic_tbi), path(jabba_rds)
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
    def bgzipped_snv_somatic = "${meta.id}.sage.somatic.vcf.bgz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    # converting the vcf.gz to bgzipped vcf.bgz;
    gunzip -c $snv_somatic | bgzip > $bgzipped_snv_somatic;
    tabix -p vcf $bgzipped_snv_somatic;

    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/HRDetect.R")

    Rscript \$RSCRIPT_PATH \\
    --sv $junctions \\
    --hets $hets \\
    --snv $bgzipped_snv_somatic \\
    --jabba $jabba_rds \\
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
