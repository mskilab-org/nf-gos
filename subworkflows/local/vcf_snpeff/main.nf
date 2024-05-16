//
// VCF SNPEFF
//

include { SNPEFF_SNPEFF } from '../../../modules/nf-core/snpeff/snpeff/main.nf'
include { SNPEFF_VCF_TO_BCF } from '../../../modules/nf-core/snpeff/snpeff/main.nf'

workflow VCF_SNPEFF {
    // defining inputs
    take:
    input                                       // required: [meta, vcf, tbi]
    snpeff_db
    snpeff_cache

    //Creating empty channels for output
    main:
    versions        = Channel.empty()
    snpeff_vcf		        = Channel.empty()
    snpeff_bcf             = Channel.empty()
    snpeff_report		    = Channel.empty()
    snpeff_summary_html	= Channel.empty()
    snpeff_genes_txt		= Channel.empty()

    SNPEFF_SNPEFF(
        input,
        snpeff_db,
        snpeff_cache
    )

    // initializing outputs from snpeff
    snpeff_vcf              = SNPEFF_SNPEFF.out.vcf
    snpeff_report		    = SNPEFF_SNPEFF.out.report
    snpeff_summary_html     = SNPEFF_SNPEFF.out.summary_html
    snpeff_genes_txt        = SNPEFF_SNPEFF.out.genes_txt
    versions                = SNPEFF_SNPEFF.out.versions

    SNPEFF_VCF_TO_BCF(snpeff_vcf)

    snpeff_bcf   = SNPEFF_VCF_TO_BCF.out.bcf

    emit:
    snpeff_vcf
    snpeff_bcf
    snpeff_report
    snpeff_summary_html
    snpeff_genes_txt
    versions
}
