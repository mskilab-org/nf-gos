//
// BAM FRAGCOUNTER
//

include { FRAGCOUNTER } from '../../../modules/local/fragcounter/main.nf'
include { REBIN_RAW_FRAGCOUNTER } from '../../../modules/local/fragcounter/main.nf'

gcmapdir       = WorkflowNfcasereports.create_file_channel(params.gcmapdir_frag)
windowsize     = WorkflowNfcasereports.create_value_channel(params.windowsize_frag)
minmapq        = WorkflowNfcasereports.create_value_channel(params.minmapq_frag)
midpoint       = WorkflowNfcasereports.create_value_channel(params.midpoint_frag)
paired         = WorkflowNfcasereports.create_value_channel(params.paired_frag)
exome          = WorkflowNfcasereports.create_value_channel(params.exome_frag)

workflow BAM_FRAGCOUNTER {
    // defining inputs
    take:
    input                                             // required: Format should be [meta, bam, bai] : can also provide cram & crai

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    fragcounter_raw_cov   = Channel.empty()
    fragcounter_cov   = Channel.empty()
    corrected_bw      = Channel.empty()
    rebinned_raw_cov  = Channel.empty()

    FRAGCOUNTER(input, midpoint, windowsize, gcmapdir, minmapq, [], [], paired, exome) // We are keeping cov empty because we don't use any input coverage for fragcounter
    //FRAGCOUNTER(input, midpoint, windowsize, gcmapdir, minmapq, paired, exome)

    // initializing outputs from fragcounter
    fragcounter_raw_cov   = FRAGCOUNTER.out.fragcounter_raw_cov
    fragcounter_cov   = FRAGCOUNTER.out.fragcounter_cov
    versions          = FRAGCOUNTER.out.versions
    corrected_bw      = FRAGCOUNTER.out.corrected_bw

    // REBIN_RAW_FRAGCOUNTER(fragcounter_cov, "reads", 1000)

    // rebinned_raw_cov  = REBIN_RAW_FRAGCOUNTER.out.raw_fragcounter_cov_1kb
    // rebinned_raw_cov  = fragcounter_cov

    //
    emit:
    fragcounter_raw_cov
    fragcounter_cov
    // rebinned_raw_cov
    corrected_bw

    versions
}
