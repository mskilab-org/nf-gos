process CONPAIR {
    tag "$meta.id"
    label 'process_medium' // Label for GPU processes

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/conpair:0.0.2':
        'mskilab/conpair:0.0.2' }"

    input:
    tuple val(meta), path(bam_tumor), path(bai_tumor), path(bam_normal), path(bai_normal)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*contamination*"), path("*concordance*"), emit: conpair_metrics
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metrics_file = "${prefix}.bammetrics.txt"
	def interval  = intervallist ? "--interval ${intervallist}" : ''

    """

    export CONPAIR_DIR=/home/git/Conpair ## Internal to Singularity container
    export PYTHONPATH=\${CONPAIR_DIR}/modules:\${PYTHONPATH}
    export GATK_JAR=/opt/conda/envs/conpair/opt/gatk-3.8/GenomeAnalysisTK.jar
    export concordance_output=./concordance.txt
    export contamination_output=./contamination.txt

    export GATK=\${GATK_JAR}

    \${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py --reference ${fasta} -B ${bam_tumor} -O TUMOR_pileup
    \${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py --reference ${fasta} -B ${bam_normal} -O NORMAL_pileup

    \${CONPAIR_DIR}/scripts/verify_concordance.py -T TUMOR_pileup -N NORMAL_pileup > \${concordance_output}
    \${CONPAIR_DIR}/scripts/estimate_tumor_normal_contamination.py -T TUMOR_pileup -N NORMAL_pileup > \${contamination_output}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(pbrun --version 2>&1 | grep "Clara Parabricks Version" | sed 's/Clara Parabricks Version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metrics_file = "${prefix}.bammetrics.txt"
    """
    touch ${metrics_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: \$(echo "stub_version")
    END_VERSIONS
    """
}
