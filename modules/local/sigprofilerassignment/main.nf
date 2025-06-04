process SIGPROFILERASSIGNMENT {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/sigprofilerassignment:0.0.4':
        'mskilab/sigprofilerassignment:0.0.4' }"

    input:
    tuple val(meta), path(vcf)
    val(genome)
    val(cosmic_version)
	path(sbs_sigs)
	path(id_sigs)

    output:
    tuple val(meta), path("sbs_results/Assignment_Solution/**/*.txt")    , emit: sbs_sigs, optional: true
    tuple val(meta), path("indel_results/Assignment_Solution/**/*.txt")    , emit: indel_sigs, optional: true
    tuple val(meta), path("sig_inputs/output/**/*.all")    , emit: sig_matrix, optional: true
	tuple val(meta), path("sig_inputs/**")    , emit: sig_inputs, optional: true
	tuple val(meta), path("sig_inputs/**/{ID}")    , emit: indels_dir, optional: true
	tuple val(meta), path("sig_inputs/output/SBS/*SBS96.all")    , emit: sbs_96_catalogue, optional: true
	tuple val(meta), path("sig_inputs/output/ID/*ID83.all")    , emit: indels_83_catalogue, optional: true
	tuple val(meta), path("sbs_results/**/*sbs_Assignment_Solution_Activities*.txt")    , emit: sbs_activities, optional: true
	tuple val(meta), path("indel_results/**/*indel_Assignment_Solution_Activities*.txt")    , emit: indel_activities, optional: true
	tuple val(meta), path("sbs_results/**/*Decomposed_Mutation_Probabilities*.txt")    , emit: sbs_posterior_prob, optional: true
	tuple val(meta), path("indel_results/**/*Decomposed_Mutation_Probabilities*.txt")    , emit: indel_posterior_prob, optional: true
	
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
	def vcf_gz      = vcf.endsWith('.gz') ? vcf : "${vcf}.gz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export SIGPROFILER_PATH=\$(echo "${baseDir}/bin/sigprofilerassignment.py")

	if  [ ! -f ${vcf_gz} ]; then
		bgzip -c ${vcf} > ${vcf_gz}
		bcftools index --tbi ${vcf_gz}
	fi

    python \$SIGPROFILER_PATH \\
    --input-vcf ${vcf_gz} \\
    --genome ${genome} \\
    --cosmic-version ${cosmic_version} \\
	--sbs-signatures ${sbs_sigs} \\
	--id-signatures ${id_sigs}

    # append sbs_ and indel_ to the output file names
    mv sbs_results/Assignment_Solution/Activities/Assignment_Solution_Activities.txt sbs_results/Assignment_Solution/Activities/sbs_Assignment_Solution_Activities.txt
    mv indel_results/Assignment_Solution/Activities/Assignment_Solution_Activities.txt indel_results/Assignment_Solution/Activities/indel_Assignment_Solution_Activities.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigprofilerassignment: ${VERSION}
    END_VERSIONS

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch signatures.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigprofiler: ${VERSION}
    END_VERSIONS
    """

}
