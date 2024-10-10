//
// QC on BAM
//

include { PICARD_COLLECTWGSMETRICS     } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS           } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'

fasta                               = WorkflowNfcasereports.create_file_channel(params.fasta)
fai                           = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
intervals                          = WorkflowNfcasereports.create_file_channel(params.intervals)

workflow BAM_QC_PICARD {
    take:
    bam                         // channel: [mandatory] [ meta, bam, bai ]

    main:
    versions = Channel.empty()
    reports = Channel.empty()

    PICARD_COLLECTWGSMETRICS(
        bam,
        fasta.map{ it -> [ [ id:'fasta' ], it ] },
        fai.map{ it -> [ [ id:'fai' ], it ] },
        []
    )

    PICARD_COLLECTMULTIPLEMETRICS(
        bam,
        fasta.map{ it -> [ [ id:'fasta' ], it ] },
        fai.map{ it -> [ [ id:'fai' ], it ] },
    )

    // Gather all reports generated
    reports = reports.mix(PICARD_COLLECTWGSMETRICS.out.metrics)
    reports = reports.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)

    // Gather versions of all tools used
    versions = versions.mix(PICARD_COLLECTWGSMETRICS.out.versions)
    versions = versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    emit:
    reports

    versions // channel: [ versions.yml ]
}

