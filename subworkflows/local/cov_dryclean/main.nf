//
// DRYCLEAN
//

include { DRYCLEAN } from '../../../modules/local/dryclean/main.nf'

pon_dryclean                        = WorkflowNfcasereports.create_file_channel(params.pon_dryclean)
center_dryclean            = WorkflowNfcasereports.create_value_channel(params.center_dryclean)
cbs_dryclean                 = WorkflowNfcasereports.create_value_channel(params.cbs_dryclean)
cnsignif_dryclean            = WorkflowNfcasereports.create_value_channel(params.cnsignif_dryclean)
wholeGenome_dryclean         = WorkflowNfcasereports.create_value_channel(params.wholeGenome_dryclean)
blacklist_dryclean           = WorkflowNfcasereports.create_value_channel(params.blacklist_dryclean)
blacklist_path_dryclean             = WorkflowNfcasereports.create_file_channel(params.blacklist_path_dryclean)
germline_file_dryclean              = WorkflowNfcasereports.create_file_channel(params.germline_file_dryclean)
germline_filter_dryclean     = WorkflowNfcasereports.create_value_channel(params.germline_filter_dryclean)
field_dryclean               = WorkflowNfcasereports.create_value_channel(params.field_dryclean)
build_dryclean               = WorkflowNfcasereports.create_value_channel(params.build_dryclean)

workflow COV_DRYCLEAN {

    take:
    input_dryclean   // channel: [mandatory] [ meta, cov(.rds file) ]

    main:
    versions          = Channel.empty()
    dryclean_cov      = Channel.empty()
    //dryclean_obj      = Channel.empty()

    DRYCLEAN(
        input_dryclean,
        pon_dryclean,
        center_dryclean,
        cbs_dryclean,
        cnsignif_dryclean,
        wholeGenome_dryclean,
        blacklist_dryclean,
        blacklist_path_dryclean,
        germline_filter_dryclean,
        germline_file_dryclean,
        field_dryclean,
        build_dryclean
    )

    dryclean_cov      = DRYCLEAN.out.decomposed_cov
    //dryclean_obj      = DRYCLEAN.out.dryclean_object

    versions          = DRYCLEAN.out.versions

    emit:
    dryclean_cov    // only need to emit the coverage for JaBbA

    versions
}

