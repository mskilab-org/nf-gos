process ESVEE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:2.0.1--hdfd78af_0' :
        'biocontainers/hmftools-esvee:2.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(fasta)
    path(fasta_fai)
    path(fasta_dict)
    path(fasta_img)
    path(esvee_pon_sgl)
    path(esvee_pon_sv)
    path(esvee_known_hotspots)
    path(esvee_repeat_mask)
    val(ref_genome_version)

    output:
    tuple val(meta), path("*.esvee.somatic.vcf.gz"), path("*.esvee.somatic.vcf.gz.tbi"), emit: somatic_vcf
    tuple val(meta), path("*.esvee.unfiltered.vcf.gz"), path("*.esvee.unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf
    tuple val(meta), path("*.esvee.germline.vcf.gz"), path("*.esvee.germline.vcf.gz.tbi"), emit: germline_vcf, optional: true
    tuple val(meta), path("*.esvee.raw.vcf.gz"), path("*.esvee.raw.vcf.gz.tbi"), emit: raw_vcf
    tuple val(meta), path("*.esvee.ref_depth.vcf.gz"), path("*.esvee.ref_depth.vcf.gz.tbi"), emit: ref_depth_vcf
    tuple val(meta), path("*.esvee.prep.junction.tsv"), emit: prep_junctions
    tuple val(meta), path("*.esvee.prep.bam"), path("*.esvee.prep.bam.bai"), emit: prep_bams
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference_arg = normal_bam ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''
    def jvmheap_mem = (task.memory.toGiga() * 0.9).toInteger()
    """
    esvee \\
        -Xmx${jvmheap_mem}g \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_bam} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        -ref_genome ${fasta} \\
        -ref_genome_version ${ref_genome_version} \\
        -write_types 'PREP_STANDARD;ASSEMBLY_STANDARD' \\
        -known_hotspot_file ${esvee_known_hotspots} \\
        -pon_sgl_file ${esvee_pon_sgl} \\
        -pon_sv_file ${esvee_pon_sv} \\
        -repeat_mask_file ${esvee_repeat_mask} \\
        -bamtool \$(command -v sambamba) \\
        -output_dir ./ \\
        -threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_id}.esvee.prep.junction.tsv
    touch ${meta.tumor_id}.esvee.prep.bam ${meta.tumor_id}.esvee.prep.bam.bai
    for kind in raw ref_depth unfiltered somatic germline; do
        touch ${meta.tumor_id}.esvee.\${kind}.vcf.gz ${meta.tumor_id}.esvee.\${kind}.vcf.gz.tbi
    done

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
