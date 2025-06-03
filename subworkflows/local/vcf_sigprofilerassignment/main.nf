//
// VCF SIGPROFILERASSIGNMENT
//

include { SIGPROFILERASSIGNMENT } from '../../../modules/local/sigprofilerassignment/main.nf'

genome          = WorkflowNfcasereports.create_value_channel(params.sigprofilerassignment_genome)
cosmic_version  = WorkflowNfcasereports.create_value_channel(params.sigprofilerassignment_cosmic_version)
sbs_signatures  = WorkflowNfcasereports.create_file_channel(params.sigprofilerassignment_sbs)
id_signatures   = WorkflowNfcasereports.create_file_channel(params.sigprofilerassignment_id)

println "${params.sigprofilerassignment_id}"

workflow VCF_SIGPROFILERASSIGNMENT {
    // defining inputs
    take:
    input                                             // required: [meta, vcf, tbi]

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    signatures   = Channel.empty()

    SIGPROFILERASSIGNMENT(
        input,
        genome,
        cosmic_version,
		sbs_signatures,
		id_signatures
    )

    sbs_signatures   = SIGPROFILERASSIGNMENT.out.sbs_sigs
    indel_signatures   = SIGPROFILERASSIGNMENT.out.indel_sigs
    signatures_matrix   = SIGPROFILERASSIGNMENT.out.sig_matrix
	sbs_activities         = SIGPROFILERASSIGNMENT.out.sbs_activities
	indel_activities       = SIGPROFILERASSIGNMENT.out.indel_activities
	sbs_posterior_prob     = SIGPROFILERASSIGNMENT.out.sbs_posterior_prob
	indel_posterior_prob   = SIGPROFILERASSIGNMENT.out.indel_posterior_prob


    emit:
    sbs_signatures
    indel_signatures
    signatures_matrix
	signatures_matrix
	sbs_activities
	indel_activities
	sbs_posterior_prob
	indel_posterior_prob


    versions
}
