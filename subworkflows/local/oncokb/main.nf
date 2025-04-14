//
// OncoKB Annotator
//

include { ONCOKB_ANNOTATOR } from '../../../modules/local/oncokb/main.nf'

// api_token    = WorkflowNfcasereports.create_value_channel(params.oncokb_api_token)
genome_version  = WorkflowNfcasereports.create_value_channel(params.sigprofilerassignment_genome)
do_vep = WorkflowNfcasereports.create_value_channel(params.do_vep_oncokb)
vep_dir = WorkflowNfcasereports.create_file_channel(params.vep_dir_oncokb)
oncokb_genes = WorkflowNfcasereports.create_file_channel(params.oncokb_genes)
gencode = WorkflowNfcasereports.create_file_channel(params.gencode_oncokb)

workflow VCF_FUSIONS_CNA_ONCOKB_ANNOTATOR {

    take:
    inputs  // [ meta, annnotated_vcf, fusions, jabba_rds]

    main:
    versions            = Channel.empty()
    merged_oncokb_vcf   = Channel.empty()
    merged_oncokb_fusions = Channel.empty()
    merged_oncokb_cna   = Channel.empty()

    ONCOKB_ANNOTATOR(
        inputs,
        genome_version,
        do_vep,
        vep_dir,
        oncokb_genes,
        gencode
    )

    merged_oncokb_vcf   = ONCOKB_ANNOTATOR.out.merged_oncokb_vcf
    merged_oncokb_fusions = ONCOKB_ANNOTATOR.out.merged_oncokb_fusions
    merged_oncokb_cna   = ONCOKB_ANNOTATOR.out.merged_oncokb_cna

    versions          = ONCOKB_ANNOTATOR.out.versions

    emit:
    merged_oncokb_vcf
	merged_oncokb_fusions
	merged_oncokb_cna

    versions
}
