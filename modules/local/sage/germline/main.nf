process SAGE_GERMLINE {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.4--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sage:3.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(tumor_bam_wgs, stageAs: "tumor.bam"), path(tumor_bai, stageAs: "tumor.bam.bai"), path(normal_bam_wgs, stageAs: "normal.bam"), path(normal_bai, stageAs: "normal.bam.bai")
    path(ref)
    path(ref_fai)
    val(ref_genome_version)
    path(ensembl_data_dir)
    path(somatic_hotspots)
    path(panel_bed)
    path(high_confidence_bed)

    output
    tuple val(meta), path('*.sage.germline.vcf.gz'), path('*.sage.germline.vcf.gz.tbi'), emit: vcf
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args        = task.ext.args ?: ''
    """

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        ${reference_bam_arg} \\
        -tumor ${meta.id} \\
        -tumor_bam ${tumor_bam_wgs} \\
        -ref_genome ${ref} \\
        -ref_genome_version ${ref_genome_version} \\
        -hotspots ${somatic_hotspots} \\
        -panel_bed ${panel_bed} \\
        -high_confidence_bed ${high_confidence_bed} \\
        -ensembl_data_dir ${ensembl_data_dir} \\
        -hotspot_min_tumor_qual 50 \\
        -panel_min_tumor_qual 75 \\
        -hotspot_max_germline_vaf 100 \\
        -hotspot_max_germline_rel_raw_base_qual 100 \\
        -panel_max_germline_vaf 100 \\
        -panel_max_germline_rel_raw_base_qual 100 \\
        -ref_sample_count 0 \\
        -panel_only \\
        -disable_bqr true \\
        -threads ${task.cpus} \\
        -output_vcf ${meta.tumor_id}.sage.germline.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_id}.sage.germline.vcf.gz
    touch ${meta.tumor_id}.sage.germline.vcf.gz.tbi
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
