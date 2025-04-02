//
// Oneness Twoness
//

include { ONENESS_TWONESS } from '../../../modules/local/onenesstwoness/main.nf'

// Oneness Twoness
model = WorkflowNfcasereports.create_file_channel(params.model_oneness_twoness)
genome = WorkflowNfcasereports.create_file_channel(params.fasta)
genome_fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)

workflow HRD_ONENESS_TWONESS {
    take:
    events_output // [meta, events_output]
    hrdetect_results // [meta, hrdetect_results]

    main:
    versions = Channel.empty()

    meta_for_merge = events_output.map { meta, events_output -> [meta.patient, meta] }
    events_output_for_merge = events_output.map { meta, events_output -> [meta.patient, events_output] }
    hrdetect_results_for_merge = hrdetect_results.map { meta, hrdetect_results -> [meta.patient, hrdetect_results] }

    oneness_twoness_input = meta_for_merge
        .join(events_output_for_merge)
        .join(hrdetect_results_for_merge)
        .map { patient, meta, events_output, hrdetect_results -> [meta, events_output, hrdetect_results] }

    ONENESS_TWONESS(oneness_twoness_input, genome, genome_fai, model)
    oneness_twoness_results = ONENESS_TWONESS.out.oneness_twoness_results

    emit:
    oneness_twoness_results
    versions
}
