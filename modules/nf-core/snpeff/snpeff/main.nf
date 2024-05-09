process SNPEFF_SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2' :
        'biocontainers/snpeff:5.1--hdfd78af_2' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val   db                       // e.g hg19 or GRch37.87
    path(cache)

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
	-v \\
        -Xmx${avail_mem}M \\
        $db \\
        $args \\
        -csvStats ${prefix}.csv \\
        $cache_command \\
        $vcf \\
        > ${prefix}.ann.vcf

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
        'docker://mskilab/utils:0.1':
        'mskilab/utils:0.1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.ann.bcf"),   emit: bcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view ${vcf} -O b -o ${prefix}.ann.unsorted.bcf &&
    bcftools sort ${prefix}.ann.unsorted.bcf -O b -o ${prefix}.ann.bcf &&
    bcftools index ${prefix}.ann.bcf; } ||

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
