//
// SV Junction Filtering (tumor only)
//

include { JUNCTION_FILTER } from '../../../modules/local/junction_filter/main.nf'

workflow SV_JUNCTION_FILTER {
    take:
    input                  //format: [meta, filtered_sv_vcf]
    junction_pon
    gnomAD_sv_db
    padding

    main:
    versions                = Channel.empty()
    pon_filtered_sv_rds     = Channel.empty()
    final_filtered_sv_rds   = Channel.empty()

    JUNCTION_FILTER(input, junction_pon, gnomAD_sv_db, padding)

    final_filtered_sv_rds   = JUNCTION_FILTER.out.final_filtered_sv_rds
    pon_filtered_sv_rds     = JUNCTION_FILTER.out.pon_filtered_sv_rds
    versions                = JUNCTION_FILTER.out.versions

    emit:
    final_filtered_sv_rds
    pon_filtered_sv_rds
    versions
}
