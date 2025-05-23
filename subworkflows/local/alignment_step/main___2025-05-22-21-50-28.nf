// This is the main workflow for the alignment step of the pipeline.

include { test_robust_absence; test_robust_presence } from "${workflow.projectDir}/lib/NfUtils"

// Run FASTQC
include { FASTQC } from "${workflow.projectDir}/modules/nf-core/fastqc/main"

// TRIM/SPLIT FASTQ Files
include { FASTP } from "${workflow.projectDir}/modules/nf-core/fastp/main"

// Loading the MULTIQC module
include { MULTIQC } from "${workflow.projectDir}/modules/nf-core/multiqc/main"

// Loading the module that dumps the versions of software being used
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "${workflow.projectDir}/modules/nf-core/custom/dumpsoftwareversions/main"

// Map input reads to reference genome
include { FASTQ_ALIGN_BWAMEM_MEM2 } from "${workflow.projectDir}/subworkflows/local/fastq_align_bwamem_mem2/main"
include { FASTQ_PARABRICKS_FQ2BAM } from "${workflow.projectDir}/subworkflows/local/fastq_parabricks_fq2bam/main"

// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS } from "${workflow.projectDir}/subworkflows/local/bam_merge_index_samtools/main"

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM } from "${workflow.projectDir}/modules/nf-core/samtools/convert/main"
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING } from "${workflow.projectDir}/modules/nf-core/samtools/convert/main"

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from "${workflow.projectDir}/modules/nf-core/samtools/convert/main"
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL } from "${workflow.projectDir}/modules/nf-core/samtools/convert/main"
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_FINAL } from "${workflow.projectDir}/modules/nf-core/samtools/convert/main"

// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES } from "${workflow.projectDir}/subworkflows/local/bam_markduplicates/main"

// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD } from "${workflow.projectDir}/subworkflows/local/cram_qc_mosdepth_samtools/main"
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL } from "${workflow.projectDir}/subworkflows/local/cram_qc_mosdepth_samtools/main"

// BAM Picard QC
include { BAM_QC } from "${workflow.projectDir}/subworkflows/local/bam_qc/main"
include { BAM_QC_PICARD_COLLECTMULTIPLEMETRICS } from "${workflow.projectDir}/subworkflows/local/bam_qc/main"
include { BAM_QC_PICARD_COLLECTWGSMETRICS } from "${workflow.projectDir}/subworkflows/local/bam_qc/main"
include { BAM_QC_GATK4_ESTIMATELIBRARYCOMPLEXITY } from "${workflow.projectDir}/subworkflows/local/bam_qc/main"

// Create recalibration tables
include { BAM_BASERECALIBRATOR } from "${workflow.projectDir}/subworkflows/local/bam_baserecalibrator/main"

// Create recalibrated cram files to use for variant calling (+QC)
include { BAM_APPLYBQSR } from "${workflow.projectDir}/subworkflows/local/bam_applybqsr/main"

workflow ALIGNMENT_STEP {
	take:
	inputs
	selected_tools_map
	tools_used
	known_sites_indels
	known_sites_indels_tbi
	index_alignment
	fasta
	fasta_fai
	dict
	intervals_for_preprocessing
	intervals_and_num_intervals

	main:
	reports = Channel.empty()
	versions = Channel.empty()

	input_fastq = inputs.filter { it.bam.isEmpty() }.map { it -> [it.meta, it.fastq_1, it.fastq_2, it.meta.read_group] }
	
	
	input_fastq_qc = input_fastq.map { it -> [it[0], [it[1], it[2]]] }


	// Inverse logic is used for QC
	// Use of inner join downstream.
	bam_qc_duplicates_inputs = inputs
		.map { it -> [it.meta, it.qc_dup_rate] }
		.filter { it[1].isEmpty() }
		.map { it -> [it[0].sample] } // meta.sample

	bam_qc_multiple_metrics_inputs = inputs
		.map { it -> [it.meta, it.qc_alignment_summary, it.qc_insert_size] }
		.filter { 
			it[1].isEmpty()
			&& it[2].isEmpty()
		}
		.map { it -> [it[0].sample] } // meta.sample

	bam_qc_coverage_inputs = inputs
		.map { it -> [it.meta, it.qc_coverage_metrics] }
		.filter { it[1].isEmpty() }
		.map { it -> [it[0].sample] } // meta.sample

	is_run_qc_duplicates = params.is_run_qc_duplicates ?: false // if parameter doesn't exist, set to false
	do_qc_coverage = tools_used.contains("qc_coverage")
	do_qc_multiple_metrics = tools_used.contains("qc_multiple_metrics")
	do_qc_duplicates = tools_used.contains("qc_duplicates") && is_run_qc_duplicates
	do_bamqc = do_qc_coverage || do_qc_multiple_metrics || do_qc_duplicates


	if (tools_used.contains("all") || tools_used.contains("aligner")) {
        
        alignment_existing_outputs = inputs.map { it -> [it.meta, it.bam] }.filter { !it[1].isEmpty() }

		// alignment_existing_outputs.view { log.info "Alignment existing outputs: $it" }

        // inputs.view { log.info "Input samples: ${it.meta} is empty? ${it.bam.isEmpty()}" }
        // input_fastq.view { log.info "Input FASTQ files: $it" }
        // alignment_existing_outputs.view { log.info "Alignment existing outputs: $it" }

        // // QC
        FASTQC(input_fastq_qc)

        // reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
        // versions = versions.mix(FASTQC.out.versions.first())

        //skipping the UMI Conscensus calling step for now
        reads_for_fastp = input_fastq

        // Trimming and/or splitting
        if (params.trim_fastq || params.split_fastq > 0) {
            if (params.trim_fastq) {
                log.warn "You have mentioned trim_fastq to `$params.trim_fastq`, will do trimming"
            }
            save_trimmed_fail = false
            save_merged = false
            FASTP(
                reads_for_fastp,
                [], // we are not using any adapter fastas at the moment
                save_trimmed_fail,
                save_merged
            )

            // reports = reports.mix(FASTP.out.json.collect{ meta, json -> json })
            // reports = reports.mix(FASTP.out.html.collect{ meta, html -> html })

            if (params.split_fastq) {
                log.warn "You have mentioned split_fastq to `$params.split_fastq`, will do splitting"
                reads_for_alignment = FASTP.out.reads.map{ meta, reads ->
                    read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                    [ meta + [ size:read_files.size() ], read_files ]
                }.transpose()
            } else reads_for_alignment = FASTP.out.reads

            // versions = versions.mix(FASTP.out.versions)

        } else {
            println "Skipping fastp since trim_fastq and split_fastq are false"
            reads_for_alignment = reads_for_fastp
        }

		// GPU vs non-GPU alignment

		

		// // STEP 1: MAPPING READS TO REFERENCE GENOME
        // // reads will be sorted
        // reads_for_alignment = reads_for_alignment
        //     .map{ meta, reads, rg ->
        //     // Update meta.id to meta.sample if no multiple lanes or splitted fastqs
        //     if (meta.num_lanes == null || meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads, rg ]
        //     else [ meta, reads, rg ]
        // }

		// Group fastqs by sample to provide all lanes per sample for fq2bam
		reads_for_alignment = reads_for_alignment
			.map { meta, fastq_1, fastq_2, rg ->
				// Ensure meta.sample exists and is used as grouping key
				[meta.sample, [meta, fastq_1, fastq_2, rg]]
			}
			.groupTuple()
			.map { sample, items ->
				// items: list of [meta, fastq_1, fastq_2, rg] for each lane
				def metas = items.collect { it[0] }
				def fastq_1s = items.collect { it[1] }.flatten().collect { file(it) }
				def fastq_2s = items.collect { it[2] }.flatten().collect { file(it) }
				def rg_list = items.collect { it[3] }

				// Use first meta as representative, add lanes info
				def meta = metas[0] + [lanes: metas.size()]

				// Return: meta, fastq_1s, fastq_2s, RGs
				[meta, fastq_1s, fastq_2s, rg_list.flatten()]
			}
			.map { meta, fastq_1s, fastq_2s, rg ->
				// Optional: Update meta.id if needed
				// if (meta.num_lanes == null || meta.size * meta.num_lanes == 1) {
				// 	out = [ meta + [ id:meta.sample ], fastq_1s, fastq_2s, rg ]
				// } else {
				// 	out = [ meta + [ id:meta.sample ], fastq_1s, fastq_2s, rg ]
				// }
				// out
				[ meta + [ id:meta.sample ], fastq_1s, fastq_2s, rg ]
			}

        // GPU Alignment
        if (params.aligner == "fq2bam") {
            FASTQ_PARABRICKS_FQ2BAM(
                reads_for_alignment,
                known_sites_indels,
                known_sites_indels_tbi
            )
            // merge existing BAMs with newly mapped ones
            bam_mapped = alignment_existing_outputs.mix(FASTQ_PARABRICKS_FQ2BAM.out.bam)
			// bam_mapped = FASTQ_PARABRICKS_FQ2BAM.out.bam}
        } else {
            FASTQ_ALIGN_BWAMEM_MEM2(
                reads_for_alignment,
                index_alignment
            )
            // merge existing BAMs with newly mapped ones
            bam_mapped = alignment_existing_outputs.mix(FASTQ_ALIGN_BWAMEM_MEM2.out.bam)
        }


        // Grouping the bams from the same samples not to stall the workflow


        bam_mapped = bam_mapped
            .map { meta, bam ->

                // Update meta.id to be meta.sample, ditching sample-lane that is not needed anymore
                // Update meta.data_type
                // Remove no longer necessary fields:
                //   read_group: Now in the BAM header
                //    num_lanes: only needed for mapping
                //         size: only needed for mapping

                // Ensure meta.size and meta.num_lanes are integers and handle null values
                int numLanes = (meta.num_lanes != null && meta.num_lanes > 1 ? meta.num_lanes : 1) as int
                int numSplits = (meta.size ?: 1) as int
                // int numReads = numLanes * numSplits
				int numReads = 1

                // Use groupKey to make sure that the correct group can advance as soon as it is complete
                // and not stall the workflow until all reads from all channels are mapped
				// cleanedMeta = meta.subMap('patient', 'sample', 'sex', 'status', 'id') + [ id: meta.sample]
				cleanedMeta = meta - meta.subMap('num_lanes', 'read_group', 'size', 'tumor_id') + [ id:meta.sample ]
                [ groupKey( cleanedMeta, numReads), bam ]
            }
        .groupTuple()

		

        // bams are merged (when multiple lanes from the same sample) and indexed
        BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

        // Gather used softwares versions
        // versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)

        alignment_bams_final = BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai.map({ meta, bam, bai -> [ meta.id, meta, bam, bai ] })
		alignment_bams_final.view{ log.info "alignment_bams_final: $it" }
    }

    // BAM Postprocessing
    // ##############################
	do_post_processing_bc_aligner_not_fq2bam = (tools_used.contains("all") || tools_used.contains("aligner")) && params.aligner != "fq2bam"
	do_post_processing_bc_of_tool_or_flag = tools_used.contains("all") || tools_used.contains("postprocessing") || params.is_run_post_processing // FIXME: If bam is provided as input, tools_used currently will never contain postprocessing and only controlled by params, but leaving here as a reminder.
    if (do_post_processing_bc_aligner_not_fq2bam || do_post_processing_bc_of_tool_or_flag) { // fq2bam does not need postprocessing
		
		bam_mapped = alignment_bams_final
            .map { id, meta, bam, bai -> [meta + [data_type: "bam"], bam] }
        cram_markduplicates_no_spark = Channel.empty()

        // STEP 2: markduplicates (+QC) + convert to CRAM

        BAM_MARKDUPLICATES(
            bam_mapped,
            fasta,
            fasta_fai,
            intervals_for_preprocessing
        )

        cram_markduplicates_no_spark = BAM_MARKDUPLICATES.out.cram

        // Gather QC reports
        // reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ meta, report -> report })

        // Gather used softwares versions
        // versions = versions.mix(BAM_MARKDUPLICATES.out.versions)

        // STEP 3: BASE RECALIBRATION
        ch_cram_for_bam_baserecalibrator = Channel.empty().mix(cram_markduplicates_no_spark)
            // Make sure correct data types are carried through
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        ch_table_bqsr_tab    = Channel.empty()

        BAM_BASERECALIBRATOR(
            ch_cram_for_bam_baserecalibrator,
            dict,
            fasta,
            fasta_fai,
            intervals_and_num_intervals,
            known_sites_indels,
            known_sites_indels_tbi
            )

        ch_table_bqsr_tab = BAM_BASERECALIBRATOR.out.table_bqsr

        // versions = versions.mix(BAM_BASERECALIBRATOR.out.versions)

        // ch_table_bqsr contains table from baserecalibrator
        ch_table_bqsr = Channel.empty().mix(ch_table_bqsr_tab)

        // reports = reports.mix(ch_table_bqsr.collect{ meta, table -> table })
        cram_applybqsr = ch_cram_for_bam_baserecalibrator.join(ch_table_bqsr, failOnDuplicate: true, failOnMismatch: true)

        // STEP 4: RECALIBRATING

        cram_variant_calling = Channel.empty()

        BAM_APPLYBQSR(
            cram_applybqsr,
            dict,
            fasta,
            fasta_fai,
            intervals_and_num_intervals
        )

        cram_variant_calling = BAM_APPLYBQSR.out.cram

        // Gather used softwares versions
        // versions = versions.mix(BAM_APPLYBQSR.out.versions)

        CRAM_QC_RECAL(
            cram_variant_calling,
            fasta,
            intervals_for_preprocessing
        )

        // Gather QC
        // reports = reports.mix(CRAM_QC_RECAL.out.reports.collect{ meta, report -> report })

        // Gather software versions
        // versions = versions.mix(CRAM_QC_RECAL.out.versions)

        // convert CRAM files to BAM for downstream processes
        CRAM_TO_BAM_RECAL(cram_variant_calling, fasta, fasta_fai)
        // versions = versions.mix(CRAM_TO_BAM_RECAL.out.versions)

        CRAM_TO_BAM_FINAL(cram_variant_calling, fasta, fasta_fai)
        // versions = versions.mix(CRAM_TO_BAM_FINAL.out.versions)

        alignment_bams_final = Channel.empty()
            .mix(CRAM_TO_BAM_FINAL.out.alignment_index)
            .map{ meta, bam, bai -> [ meta.id, meta + [data_type: "bam"], bam, bai ] }
    }

    // // Post-alignment QC
    // if (tools_used.contains("all") || tools_used.contains("bamqc")) {
    //     bam_qc_inputs = inputs.map { it -> [it.meta.sample] }
    //     bam_qc_calling = bam_qc_inputs
    //         .join(alignment_bams_final)
    //         .map { id, meta, bam, bai -> [ meta, bam, bai ] }

    //     // omit meta since it is not used in the BAM_QC
    //     dict_path = dict.map{ meta, dict -> [dict] }
    //     BAM_QC(bam_qc_calling, dict_path)

    //     Gather QC
    //     reports = reports.mix(BAM_QC.out.reports.collect{ meta, report -> report })
    //     versions = versions.mix(BAM_QC.out.versions)
    // }

	// Post-alignment QC Draft
    if (tools_used.contains("all") || do_bamqc) {

		bam_qc_duplicates_calling = bam_qc_duplicates_inputs
            .join(alignment_bams_final)
            .map { id, meta, bam, bai -> [ meta, bam, bai ] }

		bam_qc_multiple_metrics_calling = bam_qc_multiple_metrics_inputs
            .join(alignment_bams_final)
            .map { id, meta, bam, bai -> [ meta, bam, bai ] }
		
		bam_qc_coverage_calling = bam_qc_coverage_inputs
            .join(alignment_bams_final)
            .map { id, meta, bam, bai -> [ meta, bam, bai ] }

        // omit meta since it is not used in the BAM_QC

        dict_path = dict.map{ meta, dict -> [dict] }
        // BAM_QC(bam_qc_calling, dict_path)

		if (do_qc_multiple_metrics) {
			BAM_QC_PICARD_COLLECTMULTIPLEMETRICS(bam_qc_multiple_metrics_calling)
			// reports = reports.mix(BAM_QC_PICARD_COLLECTMULTIPLEMETRICS.out.reports.collect{ meta, report -> report })
			// versions = versions.mix(BAM_QC_PICARD_COLLECTMULTIPLEMETRICS.out.versions)
		}
		if (do_qc_coverage) {
			BAM_QC_PICARD_COLLECTWGSMETRICS(bam_qc_coverage_calling)
			// reports = reports.mix(BAM_QC_PICARD_COLLECTWGSMETRICS.out.reports.collect{ meta, report -> report })
			// versions = versions.mix(BAM_QC_PICARD_COLLECTWGSMETRICS.out.versions)
		}
		if (do_qc_duplicates) {
			BAM_QC_GATK4_ESTIMATELIBRARYCOMPLEXITY(bam_qc_duplicates_calling, dict_path)
			// reports = reports.mix(BAM_QC_GATK4_ESTIMATELIBRARYCOMPLEXITY.out.reports.collect{ meta, report -> report })
			// versions = versions.mix(BAM_QC_GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions)
			null
		}
	}

	emit:
	alignment_bams_final
	reports
	versions

}
