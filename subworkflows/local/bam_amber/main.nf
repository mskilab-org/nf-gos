//
// AMBER
//

include { AMBER } from '../../../modules/local/amber/main'
include { MAKE_HET_SITES } from '../../../modules/local/amber/main'

//AMBER
genome_ver = WorkflowNfcasereports.create_value_channel(params.genome_ver_amber)
het_sites = WorkflowNfcasereports.create_file_channel(params.het_sites_amber)

workflow BAM_AMBER {
    take:
    input // [meta, tbam, tbai, nbam, nbai]

    main:
    amber_dir = Channel.empty()
    sites = Channel.empty()
    versions = Channel.empty()

    AMBER(
        input,
        genome_ver,
        het_sites,
        []
    )

    amber_dir = amber_dir.mix(AMBER.out.amber_dir)
    versions = versions.mix(AMBER.out.versions)

    MAKE_HET_SITES(amber_dir)
    sites = sites.mix(MAKE_HET_SITES.out.sites)

    emit:
    amber_dir
    sites
    versions
}
