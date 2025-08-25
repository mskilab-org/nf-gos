include { RESCUE_CH_HEME } from '../../../modules/local/sage/somatic/main.nf'



workflow RESCUE_CH_HEME_STEP {
    // defining inputs
    take:
    heme_rescue_input      // [meta, sage_vcf, sage_vcf_tbi, tumoronly_vcf, tumoronly_vcf_tbi]

    //Creating empty channels for output
    main:
    heme_ref = WorkflowNfcasereports.create_file_channel(params.heme_db)
    sage_tumor_only_rescue_ch_vcf = Channel.empty()

    RESCUE_CH_HEME(
        heme_rescue_input,
        heme_ref
    )

    sage_tumor_only_rescue_ch_vcf = RESCUE_CH_HEME.out.sage_tumor_only_rescue_ch_vcf

    emit:
    sage_tumor_only_rescue_ch_vcf
}
