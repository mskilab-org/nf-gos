import Utils // Not necessary but included for clarity

//
// AMBER
//

include { AMBER } from '../../../modules/local/amber/main'
include { MAKE_HET_SITES } from '../../../modules/local/amber/main'

//AMBER
genome_ver  = WorkflowNfcasereports.create_value_channel(params.genome_ver_amber)
het_sites   = WorkflowNfcasereports.create_file_channel(params.het_sites_amber)
if (params.target_bed_amber != null) {
    target_bed_input = Channel.fromPath(params.target_bed_amber)
                        .map{ it -> [ [id: 'target_bed'], it ] }
} else {
    target_bed_input = Channel.value([id: 'target_bed'])
                        .map{ meta -> [ meta, [] ] }
}


workflow BAM_AMBER {
    take:
    tuple_input_to_amber // [meta, tbam, tbai, nbam, nbai]
    inputs_unlaned // samplesheet input

    main:
    amber_dir = Channel.empty()
    sites = Channel.empty()
    versions = Channel.empty()

    // This is a paired or tumor-only workflow, so outputs should be tied to the tumor sample,
    // or at the patient level.

    amber_existing_outputs_without_hets = inputs_unlaned
        .filter { it ->
            it.hets.isEmpty() &&
            ! it.amber_dir.isEmpty() &&
            it.meta.status.toString() == "1" 
        }
        .map { it ->
            [ it.meta + [tumor_id: it.meta.sample ], it.amber_dir ] 
        }
        .unique()
        

    // amber_existing_outputs_without_hets = inputs_unlaned
    //     .filter { it.hets.isEmpty() && ! it.amber_dir.isEmpty() && it.status.toString() == "1" }
    //     .map { [ it.meta, it.amber_dir ] }
    //     .unique()
    //     .view { "amber_existing_outputs_without_hets: ${it}" }

    AMBER(
        tuple_input_to_amber,
        genome_ver,
        het_sites,
        target_bed_input // [meta, target_bed]
    )

    amber_dir = amber_dir.mix(AMBER.out.amber_dir)
    versions = versions.mix(AMBER.out.versions)

    amber_dir_inputs_for_hets = amber_dir
        .mix(amber_existing_outputs_without_hets)
        .view { "amber_dir_inputs_for_hets: ${it}" }

    MAKE_HET_SITES(amber_dir_inputs_for_hets)
    sites = sites.mix(MAKE_HET_SITES.out.sites)

    emit:
    amber_dir
    sites
    versions
}
