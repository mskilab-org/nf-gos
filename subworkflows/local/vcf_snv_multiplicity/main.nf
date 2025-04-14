//
// VCF SNV_MULTIPLICITY
//

include { SNV_MULTIPLICITY } from '../../../modules/local/snv_multiplicity/main.nf'

workflow VCF_SNV_MULTIPLICITY {
    // defining inputs
    take:
    input                                       // required: [meta, somatic_vcf, somatic_tbi, germline_vcf, germline_tbi, jabba_gg]

    //Creating empty channels for output
    main:
    versions        = Channel.empty()
    snv_multiplicity_rds	= Channel.empty()
    snv_multiplicity_germline_rds	= Channel.empty()
    snv_multiplicity_hets_rds	= Channel.empty()

    SNV_MULTIPLICITY(
        input
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
