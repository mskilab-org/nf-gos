// SigProfilerAssignment
include { VCF_SIGPROFILERASSIGNMENT } from "${workflow.projectDir}/subworkflows/local/vcf_sigprofilerassignment/main"

// FFPE IMPACT Filter
include { FFPE_IMPACT_FILTER } from "${workflow.projectDir}/modules/local/ffpe_impact_filter/main"


workflow SIGNATURES_STEP {
	
  // defining inputs
  take:
  inputs
  snv_somatic_annotations_for_merge
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

  ffpe_filtered_existing_vcf = Channel.empty()

  inputs_tumor_status = inputs.branch{ tumor: it.meta.status == 1 }
  

  versions = Channel.empty()


  filtered_ffpe_impact_somatic_vcf_inputs = inputs_tumor_status.tumor
   	//   .filter { it.ffpe_impact_filtered_vcf.isEmpty() }
	  .map { it -> [ it.meta.patient, it.meta, it.variant_somatic_ann ] } // patient, meta, annotated_somatic_vcf
  
  
  
  vcf_out = snv_somatic_annotations_for_merge


  if (tools_used.contains("all") || tools_used.contains("signatures")) {
	filtered_ffpe_impact_somatic_vcf_for_merge = Channel.empty()
	// filtered_ffpe_impact_somatic_vcf_for_merge = inputs_tumor_status.tumor
   	//   .filter { ! it.somatic_ffpe_impact_annotated_vcf.isEmpty() }
	//   .map { it -> [ it.meta, it.ffpe_impact_filtered_vcf, it.ffpe_impact_filtered_vcf_tbi ] }
	
	signatures_inputs = inputs_tumor_status.tumor
		.filter { it.sbs_signatures.isEmpty() || it.indel_signatures.isEmpty() || it.signatures_matrix.isEmpty()}
		.map { it -> [it.meta.patient, it.meta + [id: it.meta.patient]] }

	signatures_inputs_somatic_vcf = signatures_inputs
		.join(snv_somatic_annotations_for_merge)
		.map { it -> [ it[1], it[2] ] } // meta, annotated somatic snv

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


	sbs_posterior_prob_join = sbs_posterior_prob
		.map { it -> [it[0].patient, it[1]] } // meta.patient, sbs_posterior_prob
	
	indel_posterior_prob_join = indel_posterior_prob
		.map { it -> [it[0].patient, it[1] ] } // meta.patient, indel_posterior_prob


	merged_inputs = filtered_ffpe_impact_somatic_vcf_inputs // patient, meta, annotated somatic_mut_vcf
		.join(sbs_posterior_prob_join) // patient, sbs_posterior_prob
		.join(indel_posterior_prob_join) // patient, indel_posterior_prob
		.map { it -> [it[1], it[2], it[3], it[4] ] } // meta, annotated somatic_mut_vcf, sbs_posterior_prob, indel_posterior_prob


	FFPE_IMPACT_FILTER(merged_inputs)

	vcf_out = filtered_ffpe_impact_somatic_vcf_for_merge
		.mix(FFPE_IMPACT_FILTER.out.ffpe_impact_filtered_vcf) // meta, ffpe_impact_filtered_vcf, ffpe_impact_filtered_vcf_tbi

  }

  

  emit:
  sbs_signatures
  indel_signatures
  signatures_matrix
  sbs_activities
  indel_activities
  sbs_posterior_prob
  indel_posterior_prob
  vcf_out

}
