process AMBER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-amber:4.0.1--hdfd78af_0' :
        'biocontainers/hmftools-amber:4.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    val genome_ver
    path heterozygous_sites
    path target_region_bed

    output:
    tuple val(meta), path('amber/'), emit: amber_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''

    def target_regions_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''

    """
    amber \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_bam} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        ${target_regions_bed_arg} \\
        -ref_genome_version ${genome_ver} \\
        -loci ${heterozygous_sites} \\
        -threads ${task.cpus} \\
        -output_dir amber/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: \$(amber -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p amber/
    touch amber/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

process MAKE_HET_SITES {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(amber_dir)

    output:
    tuple val(meta), path("*sites.txt"), emit: sites
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def baf_tsv = "${amber_dir}/${meta.tumor_id}.amber.baf.tsv.gz"

    """
    echo "seqnames start end alt.count.t ref.count.t alt.count.n ref.count.n" > sites.txt
    zcat ${baf_tsv} | awk 'NR>1 {
        # Calculate alt.count.t using tumorModifiedBAF
        alt_count_t = int(\$5 * \$4)  # \$4 is tumorModifiedBAF

        # Calculate ref.count.t using tumorModifiedBAF
        ref_count_t = int(\$5 * (1 - \$4))  # \$4 is tumorModifiedBAF

        # Calculate alt.count.n using normalBAF
        alt_count_n = int(\$8 * \$6)

        # Calculate ref.count.n using normalBAF
        ref_count_n = int(\$8 * (1 - \$6))

        # Print the results
        print \$1, \$2, \$2, alt_count_t, ref_count_t, alt_count_n, ref_count_n
    }' >> sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: 0.1
    END_VERSIONS
    """

    stub:
    """
    touch sites.txt

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
