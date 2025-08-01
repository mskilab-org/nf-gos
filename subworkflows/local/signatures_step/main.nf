// SigProfilerAssignment
include { VCF_SIGPROFILERASSIGNMENT } from "${workflow.projectDir}/subworkflows/local/vcf_sigprofilerassignment/main"

// FFPE IMPACT Filter
include { FFPE_IMPACT_FILTER } from "${workflow.projectDir}/modules/local/ffpe_impact_filter/main"


workflow SIGNATURES_STEP {
	
  // defining inputs
  take:
  inputs_unlaned
  snv_somatic_annotations_for_merge // meta.patient, annotated somatic snv vcf
  tools_used

  // Creating empty channels for output
  main:

  sbs_signatures = Channel.empty()
  indel_signatures = Channel.empty()
  signatures_matrix = Channel.empty()
  sbs_activities = Channel.empty()
  indel_activities = Channel.empty()
  sbs_posterior_prob = Channel.empty()
  indel_posterior_prob = Channel.empty()
  annotated_vcf_ffpe_impact_or_snpeff = Channel.empty()

  is_filter_ffpe_impact = params.filter_ffpe_impact ?: false

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

//   inputs_unlaned = inputs.map { it ->
//     it + [meta: Utils.remove_lanes_from_meta(it.meta)]
//   }

  inputs_tumor_status = inputs_unlaned.branch{ tumor: it.meta.status.toString() == "1" }
  

  versions = Channel.empty()


  filtered_ffpe_impact_somatic_vcf_inputs = inputs_tumor_status.tumor
   	  .filter { it -> it.ffpe_impact_filtered_vcf.isEmpty() }
	  .map { it -> [ it.meta.patient, it.meta ] } // patient, meta
	  .join(snv_somatic_annotations_for_merge) // patient, annotated_somatic_vcf
	  .unique()
  
  filtered_ffpe_impact_existing_outputs = inputs_tumor_status.tumor
   	  .filter { it -> ! it.ffpe_impact_filtered_vcf.isEmpty() }
	  .map { it -> [ it.meta, it.ffpe_impact_filtered_vcf, it.ffpe_impact_filtered_vcf_tbi ] }
	  .unique()

  
  inputs_tumors_meta_for_merge = inputs_tumor_status.tumor.map{ it -> [ it.meta.patient, it.meta + [id: it.meta.sample]] }.unique()
  
  annotated_vcf_ffpe_impact_or_snpeff_for_merge = snv_somatic_annotations_for_merge // meta.patient, annotated somatic snv vcf

  annotated_vcf_ffpe_impact_or_snpeff = inputs_tumors_meta_for_merge
	.join(annotated_vcf_ffpe_impact_or_snpeff_for_merge, failOnDuplicate: false, failOnMismatch: false)
	.map { it -> [ it[1], it[2]] } // meta, annotated somatic snv vcf




  if (tools_used.contains("all") || tools_used.contains("signatures")) {
	// filtered_ffpe_impact_somatic_vcf_for_merge = Channel.empty()
	filtered_ffpe_impact_somatic_vcf_for_merge = inputs_tumor_status.tumor
   	  .filter { ! it.ffpe_impact_filtered_vcf.isEmpty() }
	  .map { it -> [ it.meta, it.ffpe_impact_filtered_vcf, it.ffpe_impact_filtered_vcf_tbi ] }
	
	signatures_inputs = inputs_tumor_status.tumor
		.filter { it.sbs_signatures.isEmpty() || it.indel_signatures.isEmpty() || it.signatures_matrix.isEmpty()}
		.map { it -> [it.meta.patient, it.meta + [id: it.meta.patient]] }

	signatures_inputs_somatic_vcf = signatures_inputs
		.join(snv_somatic_annotations_for_merge)
		.map { it -> [ it[1], it[2] ] } // meta, annotated somatic snv

	sbs_signatures_existing_outputs = inputs_tumor_status.tumor.map { it -> [it.meta, it.sbs_signatures] }.filter { !it[1].isEmpty() && it[0].status.toString() == "1" }
	indel_signatures_existing_outputs = inputs_tumor_status.tumor.map { it -> [it.meta, it.indel_signatures] }.filter { !it[1].isEmpty() && it[0].status.toString() == "1" }
	signatures_matrix_existing_outputs = inputs_tumor_status.tumor.map { it -> [it.meta, it.signatures_matrix] }.filter { !it[1].isEmpty() && it[0].status.toString() == "1" }

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
		.map { patient, meta, annotated_vcf, sbs_posterior_prob, indel_posterior_prob -> 
			[meta, annotated_vcf, sbs_posterior_prob, indel_posterior_prob ] 
		} // meta, annotated somatic_mut_vcf, sbs_posterior_prob, indel_posterior_prob

	if (is_filter_ffpe_impact) {
		FFPE_IMPACT_FILTER(merged_inputs)
		annotated_vcf_ffpe_impact_or_snpeff = FFPE_IMPACT_FILTER.out.ffpe_impact_filtered_vcf
			.mix(filtered_ffpe_impact_existing_outputs) // meta, ffpe_impact_filtered_vcf, ffpe_impact_filtered_vcf_tbi
	} 	

  }

  

  emit:
  sbs_signatures
  indel_signatures
  signatures_matrix
  sbs_activities
  indel_activities
  sbs_posterior_prob
  indel_posterior_prob
  annotated_vcf_ffpe_impact_or_snpeff

}
