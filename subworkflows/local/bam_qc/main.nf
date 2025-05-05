//
// QC on BAM
//

include { PICARD_COLLECTWGSMETRICS     } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS           } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { GATK4_ESTIMATELIBRARYCOMPLEXITY           } from '../../../modules/nf-core/gatk4/estimatelibrarycomplexity/main'
include { SAMTOOLS_SUBSAMPLE           } from '../../../modules/local/subsample_reads/main'

fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
intervals = WorkflowNfcasereports.create_file_channel(params.intervals)

workflow SUBSAMPLE_BAM {
    take:
	bam                 // channel: [mandatory] [ meta, bam, bai ]
	fasta_input
                 
    main:
	versions = Channel.empty()
	bams_subsampled = Channel.empty()

	SAMTOOLS_SUBSAMPLE(
		bam,
		fasta_input,
		[]
	)

	bams_subsampled = bams_subsampled.mix(SAMTOOLS_SUBSAMPLE.out.bam_subsampled)
	versions.mix(SAMTOOLS_SUBSAMPLE.out.versions)

    emit:
	bams_subsampled       // channel: [mandatory] [ meta, bam, bai ]

	versions
}

workflow BAM_QC {
    take:
    bam                         // channel: [mandatory] [ meta, bam, bai ]
    dict

    main:
    versions = Channel.empty()
    reports = Channel.empty()


	SUBSAMPLE_BAM(
		bam,
		fasta.map{ it -> [ [ id:'fasta' ], it ] },
	)
	bams_subsampled = SUBSAMPLE_BAM.out.bams_subsampled

	// TODO: define subsample_interval as a parameter in default config
	// on NYU: "/gpfs/data/imielinskilab/DB/references/hg19/human_g1k_v37_decoy.fasta.subsampled_0.33.interval_list"
	intervals_file = params.subsample_interval ?: []

    PICARD_COLLECTWGSMETRICS(
        bam,
        fasta.map{ it -> [ [ id:'fasta' ], it ] },
        fai.map{ it -> [ [ id:'fai' ], it ] },
        intervals_file
    )

    PICARD_COLLECTMULTIPLEMETRICS(
        bam,
        fasta.map{ it -> [ [ id:'fasta' ], it ] },
        fai.map{ it -> [ [ id:'fai' ], it ] }
    )

    // bam_only = bam.map{ meta, bam, bai -> [ meta, bam ] }
	bam_only = bams_subsampled.map{ meta, bam, bai -> [ meta, bam ] }
    GATK4_ESTIMATELIBRARYCOMPLEXITY(
        bam_only,
        fasta,
        fai,
        dict
    )

    // Gather all reports generated
    reports = reports.mix(PICARD_COLLECTWGSMETRICS.out.metrics)
    reports = reports.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)
    reports = reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics)

    // Gather versions of all tools used
    versions = versions.mix(PICARD_COLLECTWGSMETRICS.out.versions)
    versions = versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
    versions = versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions)

    emit:
    reports

    versions // channel: [ versions.yml ]
}

