//
// BAM MSISENSORPRO
//

include { MSISENSORPRO_MSISOMATIC } from '../../../modules/nf-core/msisensorpro/msisomatic/main'

fasta     = WorkflowNfcasereports.create_file_channel(params.fasta)

workflow BAM_MSISENSORPRO {
    // defining inputs
    take:
    input   // required: Format should be [meta, nbam, nbai, tbam, tbai, _] : can also provide cram & crai
    msisensor_scan

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    msi_report        = Channel.empty()
    msi_dis           = Channel.empty()
    msi_somatic       = Channel.empty()
    msi_germline      = Channel.empty()

    MSISENSORPRO_MSISOMATIC(input, fasta, msisensor_scan)

    // initializing outputs from fragcounter
    msi_report      = MSISENSORPRO_MSISOMATIC.out.output_report
    msi_dis         = MSISENSORPRO_MSISOMATIC.out.output_dis
    msi_somatic     = MSISENSORPRO_MSISOMATIC.out.output_somatic
    msi_germline    = MSISENSORPRO_MSISOMATIC.out.output_germline
    versions        = MSISENSORPRO_MSISOMATIC.out.versions

    //
    emit:
    msi_report
    msi_dis
    msi_somatic
    msi_germline

    versions
}
