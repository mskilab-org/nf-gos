//
// GPU ACCELERATED MAPPING
//

include { PARABRICKS_FQ2BAM            } from '../../../modules/local/fq2bam/main'

fasta                               = WorkflowNfcasereports.create_file_channel(params.fasta)
fasta_fai                           = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
intervals                          = WorkflowNfcasereports.create_file_channel(params.intervals)
low_memory                          = params.low_memory

workflow FASTQ_PARABRICKS_FQ2BAM {
    take:
    reads // channel: [mandatory] meta, reads, intervals (intervals is optional)
    known_sites
    known_sites_tbi

    main:

    bam = Channel.empty()
    bai = Channel.empty()
    qc = Channel.empty()
    bqsr_table = Channel.empty()
    duplicate_metrics = Channel.empty()
    versions = Channel.empty()

    PARABRICKS_FQ2BAM(
        reads,
        fasta,
        fasta_fai,
        intervals,
        known_sites,
        known_sites_tbi,
        low_memory
    )

    bam = bam.mix(PARABRICKS_FQ2BAM.out.bam)
    bai = bai.mix(PARABRICKS_FQ2BAM.out.bai)
    qc = qc.mix(PARABRICKS_FQ2BAM.out.qc_metrics)
    bqsr_table = bqsr_table.mix(PARABRICKS_FQ2BAM.out.bqsr_table)
    duplicate_metrics = duplicate_metrics.mix(PARABRICKS_FQ2BAM.out.duplicate_metrics)

    versions = versions.mix(PARABRICKS_FQ2BAM.out.versions)

    emit:
    bam      // channel: [ [meta], bam ]
    bai      // channel: [ [meta], bai ]
    qc
    bqsr_table
    duplicate_metrics
    versions
}
