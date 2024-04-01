process SIGPROFILER {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/sigprofiler:0.1':
        'mskilab/sigprofiler:0.1' }"

    input:
    tuple val(meta), path(vcf)
    val(max_signatures)
    val(cosmic_version)
    val(use_gpu)


    output:
    tuple val(meta), path("*.txt")                   , emit: sigprofiler_raw_cov, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def vcf_dir     = "sigprofiler_input_vcfs/${meta.id}/"
    def output_dir  = "sigprofiler_output/${meta.id}/"
    def gpu_flag  = use_gpu ? "--gpu" : ""
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir -p ${vcf_dir} && mv ${vcf} ${vcf_dir}

    python sigprofiler.py \\
    --data-directory ${vcf_dir} \\
    --output-directory ${output_dir} \\
    --maximum-signatures ${max_signatures} \\
    --cosmic-version ${cosmic_version} \\
    ${use_gpu}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigprofiler: ${VERSION}
    END_VERSIONS

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch signatuers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigprofiler: ${VERSION}
    END_VERSIONS
    """

}
