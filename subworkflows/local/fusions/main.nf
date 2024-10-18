include { FUSIONS } from '../../../modules/local/fusions/main.nf'

gencode                     = WorkflowNfcasereports.create_file_channel(params.gencode_fusions)

workflow GGRAPH_FUSIONS {
    take:
    input           // required: format [val(meta), path(gGraph)] or [val(meta), []. path(sv_vcf), path(sv_vcf_tbi)]

    main:

    FUSIONS(input, gencode)

    fusions_output = FUSIONS.out.fusions_output
    altedge_annotations = FUSIONS.out.altedge_annotations
    versions = FUSIONS.out.versions

    emit:
    fusions_output
    altedge_annotations

    versions
}
