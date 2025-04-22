//
// COBALT
//

include { COBALT } from '../../../modules/local/cobalt/main'

gc_profile = WorkflowNfcasereports.create_file_channel(params.gc_profile)
diploid_bed = WorkflowNfcasereports.create_file_channel(params.diploid_bed)

workflow BAM_COBALT {
    take:
    input // [meta, tbam, tbai, nbam, nbai]

    main:
    cobalt_dir        = Channel.empty()
    versions          = Channel.empty()

    if (params.tumor_only) {
        COBALT(
            input,
            gc_profile,
            diploid_bed,
            []
        )
    } else {
        COBALT(
            input,
            gc_profile,
            [],
            []
        )
    }

    cobalt_dir        = cobalt_dir.mix(COBALT.out.cobalt_dir)
    versions          = versions.mix(COBALT.out.versions)

    emit:
    cobalt_dir
    versions
}
