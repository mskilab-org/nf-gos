process PAVE {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.6--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-pave:1.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_tbi)
    path genome_fasta
    path genome_fai
    val genome_ver
    path sage_pon
    path sage_blocklist_regions
    path sage_blocklist_sites
    path clinvar_annotations
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_resource

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index_only
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def gnomad_args
    if (genome_ver == '37') {
        gnomad_args = "-gnomad_freq_file ${gnomad_resource}"
    } else if (genome_ver == '38') {
        gnomad_args = "-gnomad_freq_dir ${gnomad_resource}"
    } else {
        log.error "got bad genome version: ${genome_ver}"
        System.exit(1)
    }

    // Targeted mode
    def sage_blocklist_regions_arg = sage_blocklist_regions ? "-blacklist_bed ${sage_blocklist_regions}" : ''
    def sage_blocklist_sites_arg = sage_blocklist_sites ? "-blacklist_vcf ${sage_blocklist_sites}" : ''
    def clinvar_annotations = clinvar_annotations ? "-clinvar_vcf ${clinvar_annotations}" : ''

    """
    pave \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -sample ${meta.id} \\
        -vcf_file ${sage_vcf} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -pon_file ${sage_pon} \\
        ${clinvar_annotations} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -mappability_bed ${segment_mappability} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        ${sage_blocklist_regions_arg} \\
        ${sage_blocklist_sites_arg} \\
        ${gnomad_args} \\
        -threads ${task.cpus} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch sage.pave_somatic.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

process PAVE_FILTER_VCF {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/utils:0.0.2':
        'mskilab/utils:0.0.2' }"

    input:
    tuple val(meta), path(sage_pave_vcf), path(sage_pave_tbi)

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi") , emit: filtered_vcf
    tuple val(meta), path("*.filtered.vcf.gz.tbi"), emit: filtered_index_only
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcftools view -i \'FILTER="PASS"\' -o sage.pave.filtered.vcf ${sage_pave_vcf}
    bgzip sage.pave.filtered.vcf
    tabix -p vcf sage.pave.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave_filtered_vcf: 0.1
    END_VERSIONS
    """

    stub:
    """
    touch sage.pave.filtered.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

