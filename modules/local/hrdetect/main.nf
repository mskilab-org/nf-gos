process HRDETECT {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/hrdetect:0.0.5':
        'mskilab/hrdetect:0.0.5' }"

    input:
    tuple val(meta), path(junctions), path(hets), path(snv_somatic), path(jabba_rds)
    path(fasta)
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
	def snv_somatic_stripped    = snv_somatic.toString().replaceAll(/.bgz$|.gz$/, "")
	def snv_somatic_bgz         = "${snv_somatic_stripped}.bgz"
    def bgzipped_snv_somatic = "${meta.id}.sage.somatic.vcf.bgz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    # converting the vcf.gz to bgzipped vcf.bgz;
	bcftools view $snv_somatic -Oz -o $snv_somatic_bgz
	tabix -p vcf $snv_somatic_bgz
    # gunzip -c $snv_somatic | bgzip > $bgzipped_snv_somatic;
    # tabix -p vcf $bgzipped_snv_somatic;

    export RSCRIPT_PATH=\$(echo "\${NEXTFLOW_PROJECT_DIR}/bin/HRDetect.R")

    Rscript \$RSCRIPT_PATH \\
    --sv $junctions \\
    --hets $hets \\
    --snv $snv_somatic_bgz \\
    --jabba $jabba_rds \\
    --ref $fasta \\
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
