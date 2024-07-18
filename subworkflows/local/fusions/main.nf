include { FUSIONS } from '../../../modules/local/fusions/main.nf'

gencode                     = WorkflowNfcasereports.create_file_channel(params.gencode_fusions)

workflow GGRAPH_FUSIONS {
    take:
    input           // required: format [val(meta), path(gGraph)]

    main:
    versions            = Channel.empty()
    fusions_output       = Channel.empty()

    FUSIONS(input, gencode)

    fusions_output = FUSIONS.out.fusions_output
    versions = FUSIONS.out.versions

    emit:
    fusions_output

    versions
}
