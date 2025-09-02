//
// HRDetect
//

include { HRDETECT } from '../../../modules/local/hrdetect/main.nf'


workflow JUNC_SNV_GGRAPH_HRDETECT {

    take:
    inputs  // [ meta, junction, hets, snv_somatic, snv_somatic_tbi, jabba_gg]

    main:
    versions            = Channel.empty()
    hrdetect_rds           = Channel.empty()
    hrdetect_txt            = Channel.empty()

    ref_fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
    genome_version  = WorkflowNfcasereports.create_value_channel(params.ref_hrdetect)


    HRDETECT(
        inputs,
        ref_fasta,
        genome_version
    )

    hrdetect_rds           = HRDETECT.out.hrdetect_rds
    hrdetect_txt            = HRDETECT.out.hrdetect_txt

    versions          = HRDETECT.out.versions

    emit:
    hrdetect_rds
    hrdetect_txt

    versions
}
