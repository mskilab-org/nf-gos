//
// SV Junction Filtering (tumor only)
//

include { JUNCTION_FILTER } from '../../../modules/local/junction_filter/main.nf'
include { JUNCTION_FILTER_BEDTOOLS } from '../../../modules/local/junction_filter/main.nf'

workflow SV_JUNCTION_FILTER {
    take:
    input                  //format: [meta, filtered_sv_vcf, tbi]

    main:
    def junction_pon_gridss                 = WorkflowNfcasereports.create_file_channel(params.junction_pon_gridss)
    def junction_pon_gridss_dir = WorkflowNfcasereports.create_file_channel(params.junction_pon_gridss_dir)
    def gnomAD_sv_db                        = WorkflowNfcasereports.create_file_channel(params.gnomAD_sv_db)
    def padding                             = WorkflowNfcasereports.create_value_channel(params.pad_junc_filter)

    versions                = Channel.empty()
    pon_filtered_sv_rds     = Channel.empty()
    final_filtered_sv_rds   = Channel.empty()

    JUNCTION_FILTER(input, junction_pon_gridss, gnomAD_sv_db, padding)

    final_filtered_sv_rds   = JUNCTION_FILTER.out.final_filtered_sv_rds
    pon_filtered_sv_rds     = JUNCTION_FILTER.out.pon_filtered_sv_rds
    versions                = JUNCTION_FILTER.out.versions

    emit:
    final_filtered_sv_rds
    pon_filtered_sv_rds
    versions
}


workflow SV_JUNCTION_FILTER_BEDTOOLS {
    take:
    input //format: [meta, filtered_sv_vcf, tbi]

    main:
    def junction_pon_dir = WorkflowNfcasereports.create_file_channel(params.junction_pon_dir)
    def padding = WorkflowNfcasereports.create_value_channel(params.pad_junc_filter)
    versions = Channel.empty()
    pon_filtered_sv_rds = Channel.empty()
    final_filtered_sv_rds = Channel.empty()

    JUNCTION_FILTER_BEDTOOLS(input, junction_pon_dir, padding)

    final_filtered_sv_rds   = JUNCTION_FILTER_BEDTOOLS.out.final_filtered_sv_rds
    pon_filtered_sv_rds     = JUNCTION_FILTER_BEDTOOLS.out.pon_filtered_sv_rds
    versions                = JUNCTION_FILTER_BEDTOOLS.out.versions

    emit:
    final_filtered_sv_rds
    pon_filtered_sv_rds
    versions
}
