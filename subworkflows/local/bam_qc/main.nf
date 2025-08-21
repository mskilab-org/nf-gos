//
// QC on BAM
//

include { CONPAIR } from "${workflow.projectDir}/modules/local/conpair/main"
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

    inputs_unlaned = inputs.map { it ->
        it + [meta: Utils.remove_lanes_from_meta(it.meta)]
    }


    inputs_unlaned_split = inputs_unlaned
        .branch { it -> 
            tumor: it.meta.status.toString() == "1"
            normal: it.meta.status.toString() == "0"
        }

    def normal_ids = inputs_unlaned_split.normal.map { it.meta.patient }.unique().collect().ifEmpty(["NO_NORMALS_PRESENT___MD7cicQBtB"])
    def tumor_ids = inputs_unlaned_split.tumor.map { it.meta.patient }.unique().collect()

    mixed_ids = tumor_ids
        .concat(normal_ids)
        .collect(flat: false)

    tumor_only_ids = mixed_ids
        .map{ tumor, normal ->
            tumor.findAll { !normal.contains(it) }
        }
        .flatten()
    
    tumor_paired_ids = mixed_ids
        .map{ tumor, normal ->
            tumor.findAll { normal.contains(it) }
        }
        .flatten()

    ids_without_conpair = inputs_unlaned.filter { it ->
        it.conpair_contamination.isEmpty() || it.conpair_concordance.isEmpty()
    }.map { it -> 
        [ it.meta.patient ]
    }.unique().dump(tag: "ids_without_conpair", pretty: true)

    tumor_paired_ids = tumor_paired_ids.join(ids_without_conpair)
    
    
    bam_key_patient = bam.map { meta_sample, meta, bam, bai ->
        [ meta.patient, meta, bam, bai ]
    } // Note patient keys will be duplicated

    conpair_inputs_to_combine = tumor_paired_ids // A Should be unique.. no dups.
        .cross(bam_key_patient) // B can have dups, according to nextflow docs. This will inner join but allow for dups. Essentially a filter by merge step
        .map { patient, bam ->
            def meta_bam = bam[1]
            def bam_file = bam[2]
            def bai_file = bam[3]
            [ patient, meta_bam, bam_file, bai_file ]
        }
        .dump(tag: "conpair_inputs_to_combine", pretty: true)
        .branch { it ->
            tumor: it[1].status.toString() == "1"
            normal: it[1].status.toString() == "0"
        }
    
    conpair_inputs = conpair_inputs_to_combine.tumor
        .combine(conpair_inputs_to_combine.normal, by: 0)
        .map { patient_tumor, meta_tumor, bam_tumor, bai_tumor, meta_normal, bam_normal, bai_normal ->
            def meta_tumor_out = meta_tumor + [id: meta_tumor.patient]
            meta_tumor_out = meta_tumor_out - meta_tumor_out.subMap('read_group')
            [ meta_tumor_out, bam_tumor, bai_tumor, bam_normal, bai_normal ]
        }.dump(tag: "conpair_inputs", pretty: true)

    
    CONPAIR(conpair_inputs, fasta, fai, dict)



    emit:
    reports

    versions // channel: [ versions.yml ]
}
