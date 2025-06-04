process FFPE_IMPACT_FILTER {

    tag "$meta.id"
    label 'process_medium'

    // using events container since the dependencies are the same
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.2':
        'mskilab/unified:0.0.2' }"

    input:
    tuple val(meta), path(somatic_mut_vcf), path(signatures_post_prob_sbs, stageAs: "sbs_post_prob.txt"), path(signatures_post_prob_indel, stageAs: "indel_post_prob.txt")

    output:
    tuple val(meta), path("*ffpe_annotated.vcf.gz"), path("*ffpe_annotated.vcf.gz.tbi"),  emit: ffpe_impact_vcf
	tuple val(meta), path("*ffpe_annotated_filtered.vcf.gz"), path("*ffpe_annotated_filtered.vcf.gz.tbi"), emit: ffpe_impact_filtered_vcf
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


    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/ffpe_impact_filter.R")

    Rscript \$RSCRIPT_PATH \\
        --id ${prefix} \\
        --somatic_mut_vcf ${somatic_mut_vcf} \\
        --signatures_post_prob_sbs ${signatures_post_prob_sbs} \\
        --signatures_post_prob_indel ${signatures_post_prob_indel} \\
        --cores ${task.cpus}
	
	bgzip ${prefix}_ffpe_annotated.vcf
	bgzip ${prefix}_ffpe_annotated_filtered.vcf
	
	bcftools index --tbi ${prefix}_ffpe_annotated.vcf.gz
	bcftools index --tbi ${prefix}_ffpe_annotated_filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ffpe_impact_filter: ${VERSION}
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
