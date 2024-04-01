//
// VCF PAVE
//

include { PAVE } from '../../../modules/local/pave/main.nf'
include { PAVE_FILTER_VCF } from '../../../modules/local/pave/main.nf'

workflow VCF_PAVE {
    // defining inputs
    take:
    input                                             // required: [meta, vcf, tbi]
    genome_fasta
    genome_fai
    genome_ver
    sage_pon
    sage_blocklist_regions
    sage_blocklist_sites
    clinvar_annotations
    segment_mappability
    driver_gene_panel
    ensembl_data_resources
    gnomad_resource

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    pave_raw_cov   = Channel.empty()
    pave_cov   = Channel.empty()
    corrected_bw      = Channel.empty()
    rebinned_raw_cov  = Channel.empty()

    PAVE(
        input,
        genome_fasta,
        genome_fai,
        genome_ver,
        sage_pon,
        sage_blocklist_regions,
        sage_blocklist_sites,
        clinvar_annotations,
        segment_mappability,
        driver_gene_panel,
        ensembl_data_resources,
        gnomad_resource
    )

    // initializing outputs from pave
    pave_vcf   = PAVE.out.vcf
    pave_tbi   = PAVE.out.index_only
    versions   = PAVE.out.versions

    PAVE_FILTER_VCF(pave_vcf)

    pave_filtered_vcf   = PAVE.out.filtered_vcf
    pave_filtered_tbi   = PAVE.out.filtered_index_only

    //
    emit:
    pave_vcf
    pave_tbi
    pave_filtered_vcf
    pave_filtered_tbi

    versions
}
