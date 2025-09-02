process SV_CHIMERA_FILTER {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.12':
        'mskilab/unified:0.0.12' }"

    input:
	tuple val(meta), path(vcf), path(vcf_tbi)

    output:
    tuple val(meta), path("*.ffpe_filtered.vcf.gz"), path("*.ffpe_filtered.vcf.gz.tbi"), emit: vcftbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_vcf = vcf.getName().replaceFirst(/\.vcf(\.gz|\.bgz)?$/, '.ffpe_filtered.vcf.gz')
    """
    bcftools view -Oz -i "(FORMAT/SR[1] + FORMAT/RP[1]) > 6 && (INFO/AS) >= 1 && (FORMAT/QUAL[1]) >= 150" ${vcf} > ${out_vcf}

    bcftools index --tbi ${out_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ffpe_filtered.bam
    touch ${prefix}.ffpe_filtered.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}