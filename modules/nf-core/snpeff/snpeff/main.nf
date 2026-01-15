process SNPEFF_SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2' :
    //     'biocontainers/snpeff:5.1--hdfd78af_2' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.21-snpeff':
        'mskilab/unified:0.0.21-snpeff' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val(db)                       // e.g hg19 or GRch37.87
    tuple val(meta2), path(cache)
    path(fasta)

    output:
    tuple val(meta), path("*.ann.vcf"),   emit: vcf
    tuple val(meta), path("*.csv"),       emit: report
    tuple val(meta), path("*.html"),      emit: summary_html
    tuple val(meta), path("*.genes.txt"), emit: genes_txt
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cache_command = cache ? "-dataDir \${PWD}/${cache}" : ""
    """
    snpEff \\
        -Xmx${avail_mem}M \\
        $db \\
        $args \\
        -csvStats ${prefix}.csv \\
        $cache_command \\
        $vcf \\
        > ${prefix}.ann.vcf
    
    tmpvcf=\$( export TMPDIR=./ && mktemp -t tmp_XXXXXXXXXX.vcf )
    tmpvcf2=\$( export TMPDIR=./ && mktemp -t tmp2_XXXXXXXXXX.vcf )
    (
        bcftools sort ${prefix}.ann.vcf | \\
        bcftools norm --rm-dup none --multiallelics -any --site-win 1000000 --fasta-ref ${fasta} | \\
        bcftools sort > \${tmpvcf}
    )
    
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\n' "\$tmpvcf" > old_ids.tsv
    bgzip -f old_ids.tsv 
    tabix -f -s 1 -b 2 -p vcf old_ids.tsv.gz

    (
        bcftools annotate \\
        --header-lines <(echo '##INFO=<ID=OLD_ID,Number=1,Type=String,Description="Original ID before normalization">') \\
        -a old_ids.tsv.gz \\
        -c CHROM,POS,REF,ALT,INFO/OLD_ID \\
        --set-id '%CHROM:%POS\\_%REF\\/%FIRST_ALT' \\
        -Ov "\${tmpvcf}" > \${tmpvcf2} ; \\
        mv \${tmpvcf2} ${prefix}.ann.vcf ; \\
        rm -f \${tmpvcf}
    )


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

}

process SNPEFF_VCF_TO_BCF {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/utils:0.0.2':
        'mskilab/utils:0.0.2' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.ann.bcf"),   emit: bcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${task.name}"
    """
    bcftools view ${vcf} -O b -o ${prefix}.ann.unsorted.bcf &&
    bcftools sort ${prefix}.ann.unsorted.bcf -O b -o ${prefix}.ann.bcf &&
    bcftools index ${prefix}.ann.bcf;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 0.0.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 0.0.1
    END_VERSIONS
    """
}
