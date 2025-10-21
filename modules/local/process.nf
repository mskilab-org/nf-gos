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

process ITDSEEK {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.17-itdseek':
        'mskilab/unified:0.0.17-itdseek' }"

    input:
	tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    val(assembly)

    output:
    tuple val(meta), path("*_flt3_itd.vcf"), path("*_flt3_itd_status.rds"), emit: vcfrds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def out_vcf = vcf.getName().replaceFirst(/\.vcf(\.gz|\.bgz)?$/, '.ffpe_filtered.vcf.gz')
    """

    Rscript \${NEXTFLOW_BIN_DIR}/R/flt3_itd_seek.R \\
        --bam ${bam} \\
        --fasta ${fasta} \\
        --samp_id ${meta.sample} \\
        --build ${assembly} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}


process CHECK_REFERENCE_BAM_MATCH {
    label 'process_medium'

    // using events container since the dependencies are the same
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.18-fusions':
        'mskilab/unified:0.0.18-fusions' }"

    input:
    tuple val(meta), path(bam_path)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), env(is_fully_matched), path(bam_path), emit: is_fully_matched
    tuple val(meta), path("*mismatches.tsv"), emit: mismatched_rnames
    tuple val(meta), path("*all.tsv"), emit: all_rnames

    script:
    """
    #!/usr/bin/env bash

    ## ls -1 bam_paths/** | parallel 'samtools view -H {}' | grep '@SQ' | cut -f2-3 | sed -E 's/SN://g; s/LN://g' | sort -k1,1 | uniq > bam_rname.txt
    
    samtools view -H ${bam_path} | grep '@SQ' | cut -f2-3 | sed -E 's/SN://g; s/LN://g' | sort -k1,1 | uniq > bam_rname.txt

    ${fasta_fai} | cut -f1,2 | sort -k1,1 > fasta_rname.txt

    comm -3 <(sort bam_rname.txt) <(sort fasta_reads.txt) > ${meta.sample}.mismatches.tsv
    comm <(sort bam_rname.txt) <(sort fasta_reads.txt) > ${meta.sample}.all.tsv
    num_mismatch=\$(cat mismatches.tsv | wc -l)
    is_fully_matched=\$( [ "\$num_mismatch" -eq 0 ] && echo true || echo false )

    

    """
}