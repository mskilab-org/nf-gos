//
// VCF SNV_MULTIPLICITY
//

include { SNV_MULTIPLICITY } from '../../../modules/local/snv_multiplicity/main.nf'

ref                               = WorkflowNfcasereports.create_file_channel(params.fasta)
ref_fai                           = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
workflow VCF_SNV_MULTIPLICITY {
    // defining inputs
    take:
    input                                       // required: [meta, somatic_vcf, somatic_tbi, germline_vcf, germline_tbi, jabba_rds]

    //Creating empty channels for output
    main:
    versions        = Channel.empty()
    snv_multiplicity_rds	= Channel.empty()

    SNV_MULTIPLICITY(
        input,
        ref,
        ref_fai
    )

    // initializing outputs from snv_multiplicity
    snv_multiplicity_rds              = SNV_MULTIPLICITY.out.snv_multiplicity_rds
    versions                = SNV_MULTIPLICITY.out.versions

    emit:
    snv_multiplicity_rds
    versions
}
