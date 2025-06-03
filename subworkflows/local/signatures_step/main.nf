// SigProfilerAssignment
include { VCF_SIGPROFILERASSIGNMENT } from "${workflow.projectDir}/subworkflows/local/vcf_sigprofilerassignment/main"


workflow SIGNATURES_STEP {
	
  // defining inputs
  take:
  inputs
  filtered_somatic_vcf_for_merge
  tools_used

  // Creating empty channels for output
  main:
  signatures_inputs = Channel.empty()
  signatures_inputs_somatic_vcf = Channel.empty()
  sbs_signatures_existing_outputs = Channel.empty()
  indel_signatures_existing_outputs = Channel.empty()
  signatures_matrix_existing_outputs = Channel.empty()
  
  sbs_activities_existing_outputs = Channel.empty()
  indel_activities_existing_outputs = Channel.empty()

  sbs_posterior_prob_existing_outputs = Channel.empty()
  indel_posterior_prob_existing_outputs = Channel.empty()
  

  versions = Channel.empty()

  filtered_somatic_vcf_for_merge = inputs
	.map { it -> [it.meta, it.snv_somatic_vcf, it.snv_somatic_tbi] }
	.filter { !it[1].isEmpty() && !it[2].isEmpty()}
	.map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, filtered somatic snv vcf, tbi


  if (tools_used.contains("all") || tools_used.contains("signatures")) {
	inputs_tumor_status = inputs.branch{ tumor: it.meta.status == 1 }
	signatures_inputs = inputs_tumor_status.tumor
		.filter { it.sbs_signatures.isEmpty() || it.indel_signatures.isEmpty() || it.signatures_matrix.isEmpty()}
		.map { it -> [it.meta.patient, it.meta + [id: it.meta.patient]] }

	signatures_inputs_somatic_vcf = signatures_inputs
		.join(filtered_somatic_vcf_for_merge)
		.map { it -> [ it[1], it[2], it[3] ] } // meta, somatic snv, tbi

	sbs_signatures_existing_outputs = inputs_tumor_status.tumor.map { it -> [it.meta, it.sbs_signatures] }.filter { !it[1].isEmpty() && it[0].status == 1 }
	indel_signatures_existing_outputs = inputs_tumor_status.tumor.map { it -> [it.meta, it.indel_signatures] }.filter { !it[1].isEmpty() && it[0].status == 1 }
	signatures_matrix_existing_outputs = inputs_tumor_status.tumor.map { it -> [it.meta, it.signatures_matrix] }.filter { !it[1].isEmpty() && it[0].status == 1 }

	VCF_SIGPROFILERASSIGNMENT(signatures_inputs_somatic_vcf)

	sbs_signatures = Channel.empty()
		.mix(VCF_SIGPROFILERASSIGNMENT.out.sbs_signatures)
		.mix(sbs_signatures_existing_outputs)
	indel_signatures = Channel.empty()
		.mix(VCF_SIGPROFILERASSIGNMENT.out.indel_signatures)
		.mix(indel_signatures_existing_outputs)
	signatures_matrix = Channel.empty()
		.mix(VCF_SIGPROFILERASSIGNMENT.out.signatures_matrix)
		.mix(signatures_matrix_existing_outputs)

	sbs_activities = Channel.empty()
		.mix(VCF_SIGPROFILERASSIGNMENT.out.sbs_activities)
		.mix(sbs_activities_existing_outputs)

	indel_activities = Channel.empty()
  		.mix(VCF_SIGPROFILERASSIGNMENT.out.indel_activities)
  		.mix(indel_activities_existing_outputs)

	sbs_posterior_prob = Channel.empty()
		.mix(VCF_SIGPROFILERASSIGNMENT.out.sbs_posterior_prob)
		.mix(sbs_posterior_prob_existing_outputs)

	indel_posterior_prob = Channel.empty()
		.mix(VCF_SIGPROFILERASSIGNMENT.out.indel_posterior_prob)
		.mix(indel_posterior_prob_existing_outputs)

  }

  emit:
  sbs_signatures
  indel_signatures
  signatures_matrix
  sbs_activities
  indel_activities
  sbs_posterior_prob
  indel_posterior_prob

}
