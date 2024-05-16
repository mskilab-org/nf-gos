//
// HRDetect
//

include { HRDETECT } from '../../../modules/local/hrdetect/main.nf'

workflow JUNC_SNV_GGRAPH_HRDETECT {

    take:
    inputs  // [ meta, junction, hets, snv_somatic, jabba_rds]
    ref_fasta
    genome_version

    main:
    versions            = Channel.empty()
    hrdetect_rds           = Channel.empty()
    hrdetect_txt            = Channel.empty()

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
