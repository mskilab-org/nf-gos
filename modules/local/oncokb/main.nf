process ONCOKB_ANNOTATOR {
    tag "$meta.id"
    label 'process_medium'
    secret 'ONCOKB_API_KEY'
    maxForks 3 // Limit the number of concurrent tasks to avoid overwhelming the ONCOKB API
    errorStrategy 'retry'
    maxRetries 100
    maxErrors 100

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.11':
        'mskilab/unified:0.0.11' }"

    input:
    tuple val(meta), path(vcf), path(fusions), path(cna)
    val assembly // 'GRCh37' or 'GRCh38', default is 'GRCh37'
    val do_vep
    path vep_dir
    path oncokb_genes
    path gencode

    output:
    path "merged_oncokb.maf", emit: merged_oncokb_vcf
    path "merged_oncokb_fusions.tsv", emit: merged_oncokb_fusions
    path "merged_oncokb_cna.tsv", emit: merged_oncokb_cna
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def libdir = "${baseDir}/oncokb"
    def verbose = "TRUE"
    def tumor_type = "unknown"
    def VERSION = '0.1'
    def multiplicity = "/dev/null" // multiplicity is a required argument for process_singularity.sh
	def safe_tmpdir = System.getenv('TMPDIR') ?: '/tmp'

    """
    export ONCOKB_TOKEN=\$ONCOKB_API_KEY
    export HOME=/root
	export TMPDIR=${safe_tmpdir}
    set +u  # Disable unbound variable check
    source /opt/conda/etc/profile.d/conda.sh
    conda activate pact
    \${NEXTFLOW_PROJECT_DIR}/bin/oncokb/process_singularity.sh ${vcf} ${fusions} ${cna} ${multiplicity} ${assembly} ${task.cpus} ${do_vep} ${verbose} ${tumor_type} ${vep_dir} ${oncokb_genes} ${gencode}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncokb: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch merged_oncokb.maf
    touch merged_oncokb_fusions.tsv
    touch merged_oncokb_cna.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncokb: ${VERSION}
    END_VERSIONS
    """
}
