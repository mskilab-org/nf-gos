include { EVENTS } from '../../../modules/local/events/main.nf'

workflow GGRAPH_EVENTS {
    take:
    input           // required: format [val(meta), path(gGraph)]

    main:
    fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
    versions            = Channel.empty()
    events_output       = Channel.empty()

    EVENTS(input, fasta)

    events_output = EVENTS.out.events_output
    versions = EVENTS.out.versions

    emit:
    events_output

    versions
}
