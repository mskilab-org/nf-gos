process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_high'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1"
    containerOptions "${ workflow.containerEngine == "singularity" ? '--nv': ( workflow.containerEngine == "docker" ? '--gpus all': null ) }"

    input:
    tuple val(meta), path(reads)
    path fasta
    path fasta_fai
    path interval_file
    path known_sites
    path known_sites_tbi
    val mark_duplicates
    val optical_duplicate_pixel_distance
    val low_memory

    output:
    tuple val(meta), path("*.bam")                , emit: bam
    tuple val(meta), path("*.bai")                , emit: bai
    path "versions.yml"                           , emit: versions
    path "qc_metrics", optional:true              , emit: qc_metrics
    path("*.table"), optional:true                , emit: bqsr_table
    path("duplicate-metrics.txt"), optional:true  , emit: duplicate_metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_fq_command = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"
    def known_sites_command = known_sites ? known_sites.collect{"--knownSites $it"}.join(' ') : ""
    def known_sites_output = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def mark_duplicates_output = mark_duplicates ? "--out-duplicate-metrics duplicate-metrics.txt" : ""
    def optical_duplicate_pixel_distance_command = optical_duplicate_pixel_distance && mark_duplicates ? "--optical-duplicate-pixel-distance $optical_duplicate_pixel_distance" : ""
    def qc_metrics_output = "--out-qc-metrics-dir qc_metrics"
    def mem_limit = (task.memory.toGiga() * 0.40).toInteger()
    def low_memory_command = low_memory ? "--low-memory" : ""

    known_sites.eachWithIndex { site, idx ->
        def tbi_file = known_sites_tbi[idx]
        """
        realpath_dir=\$(dirname \$(readlink -f $site))
        mv $tbi_file \$realpath_dir/
        """
    }

    """
    pbrun \\
        fq2bam \\
        --ref $fasta \\
        $in_fq_command \\
        --read-group-sm $meta.id \\
        --out-bam ${prefix}.bam \\
        $known_sites_command \\
        $known_sites_output \\
        $interval_file_command \\
        $mark_duplicates_output \\
        $optical_duplicate_pixel_distance_command \\
        $qc_metrics_output \\
        --num-gpus $task.accelerator.request \\
        --memory-limit $mem_limit \\
        $low_memory_command \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def metrics_output_command = args = "--out-duplicate-metrics duplicate-metrics.txt" ? "touch duplicate-metrics.txt" : ""
    def known_sites_output_command = known_sites ? "touch ${prefix}.table" : ""
    def qc_metrics_output_command = args = "--out-qc-metrics-dir qc_metrics " ? "mkdir qc_metrics && touch qc_metrics/alignment.txt" : ""
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    $metrics_output_command
    $known_sites_output_command
    $qc_metrics_output_command
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
