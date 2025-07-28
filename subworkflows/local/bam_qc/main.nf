//
// QC on BAM
//

include { PICARD_COLLECTWGSMETRICS } from "${workflow.projectDir}/modules/nf-core/picard/collectwgsmetrics/main"
include { PARABRICKS_BAMMETRICS } from "${workflow.projectDir}/modules/local/bammetrics/main"
include { PICARD_COLLECTMULTIPLEMETRICS } from "${workflow.projectDir}/modules/nf-core/picard/collectmultiplemetrics/main"
include { GPU_COLLECTMULTIPLEMETRICS } from "${workflow.projectDir}/modules/local/gpu_collectmultiplemetrics/main"
include { GATK4_ESTIMATELIBRARYCOMPLEXITY } from "${workflow.projectDir}/modules/nf-core/gatk4/estimatelibrarycomplexity/main"
include { SAMTOOLS_SUBSAMPLE } from "${workflow.projectDir}/modules/local/subsample_reads/main"

fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
intervals = WorkflowNfcasereports.create_file_channel(params.intervals)

// tools_to_run = WorkflowNfcasereports.toolsToRun

workflow BAM_QC {
    take:
    inputs
    bam // channel: [mandatory] [ meta.sample, meta, bam, bai ]
    dict
	tools_used

    main:
    versions = Channel.empty()
    reports = Channel.empty()

	// TODO: define subsample_interval as a parameter in default config
	// on NYU: "/gpfs/data/imielinskilab/DB/references/hg19/human_g1k_v37_decoy.fasta.subsampled_0.33.interval_list"
	intervals_file = params.subsample_interval ?: []

	use_gpu = params.use_gpu

    if (tools_used.contains("all") || tools_used.contains("collect_wgs_metrics")) {
        collect_wgs_metrics_inputs = inputs
            .filter { it.qc_coverage_metrics.isEmpty() }
            .map { it -> [it.meta.sample] }
            .unique()
        collect_wgs_metrics_bam = collect_wgs_metrics_inputs
            .join(bam)
            .map { id, meta, bam, bai -> [ meta, bam, bai ] }
		ucollect_wgs_metrics_bam = collect_wgs_metrics_bam.unique { it -> it[0].sample }
		if (use_gpu) {
			process_wgsmetrics = PARABRICKS_BAMMETRICS(
				ucollect_wgs_metrics_bam,
				fasta.map{ it -> [ [ id:'fasta' ], it ] },
				fai.map{ it -> [ [ id:'fai' ], it ] },
				intervals_file
			)
		} else {
			process_wgsmetrics = PICARD_COLLECTWGSMETRICS(
				ucollect_wgs_metrics_bam,
				fasta.map{ it -> [ [ id:'fasta' ], it ] },
				fai.map{ it -> [ [ id:'fai' ], it ] },
				intervals_file
			)
		}

		reports = reports.mix(process_wgsmetrics.metrics)
		versions = versions.mix(process_wgsmetrics.versions)
    }


	is_aligner_excludes_multiple_metrics = params.aligner != "fq2bam"
	do_gpu_multiple_metrics = use_gpu && is_aligner_excludes_multiple_metrics // fq2bam already includes this.
	
	// do_gpu_multiple_metrics = use_gpu

    if (tools_used.contains("all") || tools_used.contains("collect_multiple_metrics")) {
        collect_multiple_metrics_inputs = inputs
            .filter { it.qc_alignment_summary.isEmpty() || it.qc_insert_size.isEmpty() }
            .map { it -> [it.meta.sample] }
            .unique()
        collect_multiple_metrics_bam = collect_multiple_metrics_inputs
            .join(bam)
            .map { id, meta, bam, bai -> [ meta, bam, bai ] }

		process_mm_qc = [metrics: Channel.empty(), versions: Channel.empty()]
		// FIXME: GPU multiple metrics failing to pick up GPU
		if (do_gpu_multiple_metrics) {
			process_mm_qc = GPU_COLLECTMULTIPLEMETRICS(
				collect_multiple_metrics_bam,
				fasta.map{ it -> [ [ id:'fasta' ], it ] },
				fai.map{ it -> [ [ id:'fai' ], it ] }
			)
		} else if (is_aligner_excludes_multiple_metrics) {
			process_mm_qc = PICARD_COLLECTMULTIPLEMETRICS(
				collect_multiple_metrics_bam,
				fasta.map{ it -> [ [ id:'fasta' ], it ] },
				fai.map{ it -> [ [ id:'fai' ], it ] }
			)
			
		}
		reports = reports.mix(process_mm_qc.metrics)
        versions = versions.mix(process_mm_qc.versions)
    }

	is_run_qc_duplicates = params.is_run_qc_duplicates ?: false // if parameter doesn't exist, set to false
	do_qc_duplicates = (tools_used.contains("all") || tools_used.contains("estimate_library_complexity")) && is_run_qc_duplicates && ! params.aligner == "fq2bam"
    if (do_qc_duplicates) {
        estimate_library_complexity_inputs = inputs
            .filter { it.qc_dup_rate.isEmpty() }
            .map { it -> [it.meta.sample] }
        estimate_library_complexity_bam = estimate_library_complexity_inputs
            .join(bam)
            .map { id, meta, bam, bai -> [ meta, bam, bai ] }
        // Subsample BAMs for faster estimation of library complexity
        // SAMTOOLS_SUBSAMPLE(
        //     estimate_library_complexity_bam,
        //     fasta.map{ it -> [ [ id:'fasta' ], it ] },
        //     []
        // )
        // bams_subsampled = SAMTOOLS_SUBSAMPLE.out.bam_subsampled

        // bam_only = bams_subsampled.map{ meta, bam, bai -> [ meta, bam ] }
		bam_only = estimate_library_complexity_bam.map{ meta, bam, bai -> [ meta, bam ] }
        GATK4_ESTIMATELIBRARYCOMPLEXITY(
            bam_only,
            fasta,
            fai,
            dict
        )
        reports = reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics)
        versions = versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions)
    }

    emit:
    reports

    versions // channel: [ versions.yml ]
}
