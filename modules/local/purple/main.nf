process PURPLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-purple:4.0.2--hdfd78af_0' :
        'biocontainers/hmftools-purple:4.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(amber), path(cobalt), path(sv_tumor_vcf), path(sv_tumor_tbi), path(smlv_tumor_vcf), path(smlv_tumor_tbi), path(smlv_normal_vcf), path(smlv_normal_tbi)
    path genome_fasta
    val genome_ver
    val highly_diploid_percentage
    val min_purity
    val ploidy_penalty_factor
    path genome_fai
    path genome_dict
    path gc_profile
    path sage_known_hotspots_somatic
    path sage_known_hotspots_germline
    path driver_gene_panel
    path ensembl_data_resources
    path germline_del
    path target_region_bed
    path target_region_ratios
    path target_region_msi_indels

    output:
    tuple val(meta), path('purple/*'), emit: purple_dir
    tuple val(meta), path('purple/*.purple.purity.tsv'), emit: purple_purity
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''

    def sv_tumor_vcf_arg = sv_tumor_vcf ? "-somatic_sv_vcf ${sv_tumor_vcf}" : ''
    // def sv_normal_vcf_arg = sv_normal_vcf ? "-germline_sv_vcf ${sv_normal_vcf}" : ''

    // def sv_tumor_recovery_vcf_arg = sv_tumor_unfiltered_vcf ? "-sv_recovery_vcf ${sv_tumor_unfiltered_vcf}" : ''

    def smlv_tumor_vcf_arg =  smlv_tumor_vcf ? "-somatic_vcf ${smlv_tumor_vcf}" : ''
    def smlv_normal_vcf_arg =  smlv_normal_vcf ? "-germline_vcf ${smlv_normal_vcf}" : ''

    def highly_diploid_percentage_arg = highly_diploid_percentage ? "-highly_diploid_percentage ${highly_diploid_percentage}" : ''
    def min_purity_arg = min_purity ? "-min_purity ${min_purity}" : ''
    def ploidy_penalty_factor_arg = ploidy_penalty_factor ? "-ploidy_penalty_factor ${ploidy_penalty_factor}" : ''

    def sage_known_hotspots_germline_arg = sage_known_hotspots_germline ? "-germline_hotspots ${sage_known_hotspots_germline}" : ''
    def germline_del_arg = germline_del ? "-germline_del_freq_file ${germline_del}" : ''

    def target_region_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''
    def target_region_ratios_arg = target_region_ratios ? "-target_regions_ratios ${target_region_ratios}" : ''
    def target_region_msi_indels_arg = target_region_msi_indels ? "-target_regions_msi_indels ${target_region_msi_indels}" : ''

    """
    purple \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.sample} \\
        ${reference_arg} \\
        -amber ${amber} \\
        -cobalt ${cobalt} \\
        ${sv_tumor_vcf_arg} \\
        ${smlv_tumor_vcf_arg} \\
        ${smlv_normal_vcf_arg} \\
        ${highly_diploid_percentage_arg} \\
        ${min_purity_arg} \\
        ${ploidy_penalty_factor_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -somatic_hotspots ${sage_known_hotspots_somatic} \\
        ${sage_known_hotspots_germline_arg} \\
        ${target_region_bed_arg} \\
        ${target_region_ratios_arg} \\
        ${target_region_msi_indels_arg} \\
        ${germline_del_arg} \\
        -gc_profile ${gc_profile} \\
        -circos \$(which circos) \\
        -threads ${task.cpus} \\
        -output_dir purple/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: \$(purple -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_id}.purple.cnv.gene.tsv
    touch ${meta.tumor_id}.purple.cnv.somatic.tsv
    touch ${meta.tumor_id}.purple.driver.catalog.germline.tsv
    touch ${meta.tumor_id}.purple.driver.catalog.somatic.tsv
    touch ${meta.tumor_id}.purple.germline.vcf.gz
    touch ${meta.tumor_id}.purple.germline.vcf.gz
    touch ${meta.tumor_id}.purple.purity.tsv
    touch ${meta.tumor_id}.purple.qc
    touch ${meta.tumor_id}.purple.somatic.vcf.gz
    touch ${meta.tumor_id}.purple.sv.germline.vcf.gz
    touch ${meta.tumor_id}.purple.sv.vcf.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

process EXTRACT_PURITYPLOIDY {
    input:
    tuple val(meta), path(purple_purity)

    output:
    tuple val(meta), env(purity_val), emit: purity_val
    tuple val(meta), env(ploidy_val), emit: ploidy_val

    script:
    """
    export purity_val=\$(awk 'NR==2 {print \$1}' ${purple_purity})
    export ploidy_val=\$(awk 'NR==2 {print \$5}' ${purple_purity})
    """
}
