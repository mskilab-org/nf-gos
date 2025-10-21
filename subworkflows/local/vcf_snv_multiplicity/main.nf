//
// VCF SNV_MULTIPLICITY
//

include { SNV_MULTIPLICITY } from '../../../modules/local/snv_multiplicity/main.nf'

dryclean_field		    = WorkflowNfcasereports.create_value_channel(params.field_cbs) // channel: [mandatory] dryclean field to use for SNV multiplicity (same as for CBS)

workflow VCF_SNV_MULTIPLICITY {
    // defining inputs
    take:
    input                                       // required: [meta, somatic_vcf, somatic_tbi, germline_vcf, germline_tbi, jabba_gg, hets, dryclean_cov]

    //Creating empty channels for output
    main:
    versions        = Channel.empty()
    snv_multiplicity_rds	= Channel.empty()
    snv_multiplicity_germline_rds	= Channel.empty()
    snv_multiplicity_hets_rds	= Channel.empty()

    SNV_MULTIPLICITY(
        input,
        dryclean_field
    )

    // initializing outputs from snv_multiplicity
    snv_multiplicity_rds              = SNV_MULTIPLICITY.out.snv_multiplicity_rds
    snv_multiplicity_germline_rds              = SNV_MULTIPLICITY.out.snv_multiplicity_germline_rds
    snv_multiplicity_hets_rds              = SNV_MULTIPLICITY.out.snv_multiplicity_hets_rds
    versions                = SNV_MULTIPLICITY.out.versions

    emit:
    snv_multiplicity_rds
    snv_multiplicity_germline_rds
    snv_multiplicity_hets_rds
    versions
}
