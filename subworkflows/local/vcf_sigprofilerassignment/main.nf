//
// VCF SIGPROFILERASSIGNMENT
//

include { SIGPROFILERASSIGNMENT } from '../../../modules/local/sigprofilerassignment/main.nf'

genome          = WorkflowNfcasereports.create_value_channel(params.sigprofilerassignment_genome)
cosmic_version  = WorkflowNfcasereports.create_value_channel(params.sigprofilerassignment_cosmic_version)
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
        cosmic_version
    )

    sbs_signatures   = SIGPROFILERASSIGNMENT.out.sbs_sigs
    indel_signatures   = SIGPROFILERASSIGNMENT.out.indel_sigs
    signatures_matrix   = SIGPROFILERASSIGNMENT.out.sig_matrix

    emit:
    sbs_signatures
    indel_signatures
    signatures_matrix

    versions
}
