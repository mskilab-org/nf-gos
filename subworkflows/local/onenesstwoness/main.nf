//
// Oneness Twoness
//

include { ONENESS_TWONESS } from '../../../modules/local/onenesstwoness/main.nf'
include { HOMEOLOGY } from '../../../modules/local/onenesstwoness/main.nf'

// Oneness Twoness
// model = WorkflowNfcasereports.create_file_channel(params.onenesstwoness_model)
// Homeology
width = WorkflowNfcasereports.create_value_channel(params.homeology_width)
pad = WorkflowNfcasereports.create_value_channel(params.homeology_pad)
thresh = WorkflowNfcasereports.create_value_channel(params.homeology_thresh)
stride = WorkflowNfcasereports.create_value_channel(params.homeology_stride)
genome = WorkflowNfcasereports.create_file_channel(params.fasta)
flip = WorkflowNfcasereports.create_value_channel(params.homeology_flip)
bidirectional = WorkflowNfcasereports.create_value_channel(params.homeology_bidirectional)
annotate = WorkflowNfcasereports.create_value_channel(params.homeology_annotate)
savegMatrix = WorkflowNfcasereports.create_value_channel(params.homeology_savegMatrix)

workflow HRD_ONENESS_TWONESS {
    take:
    events_output // [meta, events_output]
    hrdetect_results // [meta, hrdetect_results]

    main:
    versions = Channel.empty()

    meta_for_merge = events_output.map { meta, events_output -> [meta.patient, meta] }
    events_output_for_merge = events_output.map { meta, events_output -> [meta.patient, events_output] }
    hrdetect_results_for_merge = hrdetect_results.map { meta, hrdetect_results -> [meta.patient, hrdetect_results] }

    EVENTS_HOMEOLOGY(events_output)
    homeology_for_merge = EVENTS_HOMEOLOGY.out.homeology.map { meta, homeology -> [meta.patient, homeology] }
    homeology_stats_for_merge = EVENTS_HOMEOLOGY.out.homeology_stats.map { meta, homeology_stats -> [meta.patient, homeology_stats] }

    oneness_twoness_input = meta_for_merge
        .join(events_output_for_merge)
        .join(homeology_for_merge)
        .join(homeology_stats_for_merge)
        .join(hrdetect_results_for_merge)
        .map { patient, meta, events_output, homeology, homeology_stats, hrdetect_results -> [meta, events_output, homeology, homeology_stats, hrdetect_results] }

    ONENESS_TWONESS(oneness_twoness_input)
    oneness_twoness_results = ONENESS_TWONESS.out.oneness_twoness_results

    emit:
    oneness_twoness_results
    versions
}

workflow EVENTS_HOMEOLOGY {
    take:
    events_output // [meta, events_output]

    main:
    versions = Channel.empty()

    HOMEOLOGY(events_output, width, pad, thresh, stride, genome, flip, bidirectional, annotate, savegMatrix)

    homeology = HOMEOLOGY.out.homeology
    homeology_stats = HOMEOLOGY.out.homeology_stats
    homeology_gm = HOMEOLOGY.out.homeology_gm

    emit:
    homeology
    homeology_stats
    homeology_gm
    versions
}
