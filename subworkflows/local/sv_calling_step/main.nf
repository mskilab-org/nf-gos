//GRIDSS
include { BAM_SVCALLING_GRIDSS } from "${workflow.projectDir}/subworkflows/local/bam_svcalling_gridss/main"
include { BAM_SVCALLING_GRIDSS_SOMATIC } from "${workflow.projectDir}/subworkflows/local/bam_svcalling_gridss/main"

include { GRIDSS_PREPROCESS as GRIDSS_PREPROCESS_TUMOR   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"
include { GRIDSS_PREPROCESS as GRIDSS_PREPROCESS_NORMAL   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"
include { GRIDSS_ASSEMBLE_SCATTER   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"
include { GRIDSS_ASSEMBLE_GATHER   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"
include { GRIDSS_CALL   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"


fasta                               = WorkflowNfcasereports.create_file_channel(params.fasta)
fasta_fai                           = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
blacklist_gridss                    = WorkflowNfcasereports.create_file_channel(params.blacklist_gridss)


workflow SV_CALLING_STEP {
    take:
    inputs_unlaned
    alignment_bams_final
    bwa_index
    tools_used

    main:

    versions               = Channel.empty()
    vcf                    = Channel.empty()
    vcf_index              = Channel.empty()
    assembly_bam           = Channel.empty()
    vcf_from_gridss_gridss = Channel.empty()


    parallelize_gridss = params.parallelize_gridss ?: false

    inputs_unlaned_split = inputs_unlaned
        .branch { it -> 
            tumor: it.meta.status.toString() == "1"
            normal: it.meta.status.toString() == "0"
        }

    def normal_ids = inputs_unlaned_split.normal.map { it.meta.patient }.distinct().collect().ifEmpty(["NO_NORMALS_PRESENT___MD7cicQBtB"]).view { "Normal IDs: $it" }
    def tumor_ids = inputs_unlaned_split.tumor.map { it.meta.patient }.distinct().collect().view { "Tumor IDs: $it" }

    def total_jobnodes = 4


    // SV Calling
    // ##############################
    if (tools_used.contains("all") || tools_used.contains("gridss") || params.is_run_junction_filter) {

        gridss_existing_outputs = inputs_unlaned.map {
            it -> [it.meta, it.vcf, it.vcf_tbi] }
            .filter { !it[1].isEmpty() && !it[2].isEmpty() }

        // Filter out bams for which SV calling has already been done

        bam_sv_inputs = inputs_unlaned.filter { it.vcf.isEmpty() }.map { it -> [it.meta.sample] }.distinct()
        bam_sv_calling = alignment_bams_final // meta.sample, meta, bam, bai
            .join(bam_sv_inputs) // 
            .map { it -> [ it[1], it[2], it[3] ] } // meta, bam, bai
            .view { "BAM SV calling input: $it" }

        bam_sv_calling_status = bam_sv_calling.branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }
        
        if (parallelize_gridss) {

            gridss_preprocess_tumor = GRIDSS_PREPROCESS_TUMOR(bam_sv_calling_status.tumor, fasta, fasta_fai, bwa_index, blacklist_gridss)
            gridss_preprocess_normal = GRIDSS_PREPROCESS_NORMAL(bam_sv_calling_status.normal, fasta, fasta_fai, bwa_index, blacklist_gridss)
            
            gridss_preprocess_tumor_for_merge = gridss_preprocess_tumor.gridss_preprocess_dir
                .map { meta, gridss_preprocess_dir ->
                    [ meta.patient, gridss_preprocess_dir ]
                }

            gridss_preprocess_normal_for_merge = gridss_preprocess_normal.gridss_preprocess_dir
                .map { meta, gridss_preprocess_dir ->
                    [ meta.patient, gridss_preprocess_dir ]
                }

            mixed_ids = tumor_ids
                .concat(normal_ids)
                .collect(flat: false)

            tumor_only_ids = mixed_ids
                .map{ tumor, normal ->
                    tumor.findAll { !normal.contains(it) }
                }
                .flatten()
                .view { "Tumor IDs without normal: $it" }
            
            tumor_paired_ids = mixed_ids
                .map{ tumor, normal ->
                    tumor.findAll { normal.contains(it) }
                }
                .flatten()
                .view { "Tumor IDs with normal: $it" }

            assembly_preinput_tumor_paired = bam_sv_calling_status.tumor
                .map{ it ->
                    [it[0].patient, it[0], it[1], it[2]]
                }
                .join(tumor_paired_ids)
                .join(gridss_preprocess_tumor_for_merge) // patient, meta, bam, bai, gridss_preprocess_dir
                // .view { "Assembly preinput tumor paired: $it" }

            assembly_preinput_normal = bam_sv_calling_status.normal
                .map{ it ->
                    [it[0].patient, it[0], it[1], it[2]]
                }
                .join(gridss_preprocess_normal_for_merge) // patient, meta, bam, bai, gridss_preprocess_dir
                // .view { "Assembly preinput normal: $it" }

            assembly_paired_input = assembly_preinput_tumor_paired
                .cross(assembly_preinput_normal) { v -> v[0] }
                .flatMap { tumor, normal ->
                    def meta_id = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                    ( 0..< total_jobnodes ).collect { i ->
                        [
                            patient: tumor[0],
                            meta: tumor[1] + [id: meta_id, tumor_id: tumor[1].sample, normal_id: normal[1].sample], 
                            tumor_bam: tumor[2], 
                            tumor_bai: tumor[3], 
                            tumor_gridss_preprocess: tumor[4],
                            normal_bam: normal[2],
                            normal_bai: normal[3],
                            normal_gridss_preprocess: normal[4],
                            jobnodes: total_jobnodes,
                            jobindex: i
                        ]
                    }
                }
                .map { 
                    it.values()
                }
                

            assembly_tumoronly_input = bam_sv_calling_status.tumor
                .map{ it ->
                    [it[0].patient, it[0], it[1], it[2]]
                }
                .view { "bam_sv_calling_status.tumor.map: ${it}" }
                .join(tumor_only_ids)
                .join(gridss_preprocess_tumor_for_merge)
                .flatMap { tumor ->
                    ( 0..<total_jobnodes ).collect { i ->
                        [
                            patient: tumor[1].patient, 
                            meta: tumor[1] + [id: tumor[1].sample, tumor_id: tumor[1].sample, normal_id: ""], 
                            tumor_bam: tumor[2], 
                            tumor_bai: tumor[3], 
                            tumor_gridss_preprocess: tumor[4],
                            normal_bam: [],
                            normal_bai: [],
                            normal_gridss_preprocess: [],
                            jobnodes: total_jobnodes,
                            jobindex: i
                        ]
                    }
                }
                .map { 
                    it.values()
                }
            
            assembly_mixed_input = assembly_paired_input
                .mix(assembly_tumoronly_input)
                .map{ it.toList()[1..-1] }
                // .view{ "assembly_mixed_input $it" }
            
            GRIDSS_ASSEMBLE_SCATTER(assembly_mixed_input, fasta, fasta_fai, bwa_index, blacklist_gridss)

            collected_assembly_dirs = GRIDSS_ASSEMBLE_SCATTER.out.gridss_workdir
                .map { meta, gridss_scatter_assembly_paths_per_jobnode ->
                    [meta.patient, gridss_scatter_assembly_paths_per_jobnode]
                }
                .groupTuple(by: 0, size: total_jobnodes) // [meta.patient, list[workdirs0, workdirs1, ...]]
                .map { patient, list_of_gridss_scatter_assembly_paths ->
                    def gridss_scatter_assembly_paths = list_of_gridss_scatter_assembly_paths.flatten()
                    assembly_dir = gridss_scatter_assembly_paths.collect{ it.getParent().getName().toString() }.unique()[0]
                    gridss_scatter_assembly_paths = gridss_scatter_assembly_paths.findAll { it =~ /.*chunk.*\.(bam|bai)$/ }
                    [patient, assembly_dir, gridss_scatter_assembly_paths]
                }
                // .view { "collected_assembly_dirs: $it" }
            
            gatherassembly_mixed_input = assembly_mixed_input
                .map{ meta, tumor_bam, tumor_bai, tumor_gridss_preprocess, normal_bam, normal_bai, normal_gridss_preprocess, jobnodes, jobindex ->
                    [
                        meta.patient, 
                        meta, 
                        tumor_bam, 
                        tumor_bai,
                        tumor_gridss_preprocess,
                        normal_bam, 
                        normal_bai,
                        normal_gridss_preprocess
                    ]
                }
                .distinct()
                .join(
                    collected_assembly_dirs
                )
                .map{ patient, meta, tumor_bam, tumor_bai, tumor_gridss_preprocess, normal_bam, normal_bai, normal_gridss_preprocess, gridss_assembly_dir, gridss_scatter_assembly_paths ->
                    [
                        meta, 
                        tumor_bam, 
                        tumor_bai,
                        tumor_gridss_preprocess,
                        normal_bam, 
                        normal_bai,
                        normal_gridss_preprocess,
                        gridss_assembly_dir,
                        gridss_scatter_assembly_paths
                        
                    ]
                }
                // .view { "GRIDSS_ASSEMBLE_GATHER input: $it" }
            
            GRIDSS_ASSEMBLE_GATHER(gatherassembly_mixed_input, fasta, fasta_fai, bwa_index, blacklist_gridss)

            call_input = gatherassembly_mixed_input.map{ meta, tumor_bam, tumor_bai, tumor_gridss_preprocess, normal_bam, normal_bai, normal_gridss_preprocess, gridss_assembly_dir, gridss_scatter_assembly_paths ->
                    [
                        meta.patient, 
                        meta, 
                        tumor_bam, 
                        tumor_bai,
                        tumor_gridss_preprocess, 
                        normal_bam, 
                        normal_bai,
                        normal_gridss_preprocess, 
                        gridss_assembly_dir,
                        gridss_scatter_assembly_paths
                    ]
                }
                .distinct()
                .join(
                    GRIDSS_ASSEMBLE_GATHER.out.gridss_final_assembly
                        .map { meta, gridss_final_assembly, gridss_gather_assembly_paths ->
                            [meta.patient, gridss_final_assembly, gridss_gather_assembly_paths]
                        }
                )
                .map{ patient, meta, tumor_bam, tumor_bai, tumor_gridss_preprocess, normal_bam, normal_bai, normal_gridss_preprocess, gridss_assembly_dir, gridss_scatter_assembly_paths, gridss_final_assembly, gridss_gather_assembly_paths ->
                    [
                        meta, 
                        tumor_bam, 
                        tumor_bai, 
                        tumor_gridss_preprocess, 
                        normal_bam, 
                        normal_bai,
                        normal_gridss_preprocess,
                        gridss_assembly_dir,
                        gridss_scatter_assembly_paths + gridss_gather_assembly_paths,
                        gridss_final_assembly
                    ]
                }
                // .view { "GRIDSS_CALL input: $it" }
            
            GRIDSS_CALL(call_input, fasta, fasta_fai, bwa_index, blacklist_gridss)

            vcf_from_gridss_gridss = vcf_from_gridss_gridss
                .mix(GRIDSS_CALL.out.filtered_vcf)
                .mix(gridss_existing_outputs)
        } else {

            // gridss_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.vcf, it.vcf_tbi] }.filter { !it[1].isEmpty() && !it[2].isEmpty() }
            
            if (params.tumor_only) {
                // bam_sv_calling_status = bam_sv_calling.branch{
                //     tumor:  it[0].status == 1
                // }

                // add empty arrays to stand-in for normals
                bam_sv_calling_pair = bam_sv_calling_status.tumor.map{ meta, bam, bai -> [ meta + [tumor_id: meta.sample], [], [], bam, bai ] }
            } else {
                // getting the tumor and normal cram files separated
                // bam_sv_calling_status = bam_sv_calling.branch{
                //     normal: it[0].status == 0
                //     tumor:  it[0].status == 1
                // }

                // All normal samples
                bam_sv_calling_normal_for_crossing = bam_sv_calling_status.normal.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }.view{" Normal samples for crossing: $it" }

                // All tumor samples
                bam_sv_calling_tumor_for_crossing = bam_sv_calling_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }.view{" Tumor samples for crossing: $it" }

                // Crossing the normal and tumor samples to create tumor and normal pairs
                bam_sv_calling_pair = bam_sv_calling_normal_for_crossing.cross(bam_sv_calling_tumor_for_crossing)
                    .map { normal, tumor ->
                        def meta = [:]

                        meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                        meta.normal_id  = normal[1].sample
                        meta.patient    = normal[0]
                        meta.sex        = normal[1].sex
                        meta.tumor_id   = tumor[1].sample

                        [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
                }
            }

            BAM_SVCALLING_GRIDSS(
                bam_sv_calling_pair,
                bwa_index
            )
            vcf_from_gridss_gridss = Channel.empty()
                .mix(BAM_SVCALLING_GRIDSS.out.vcf)
                .mix(gridss_existing_outputs)
        }

        // vcf_from_gridss_gridss = Channel.empty()
        //     .mix(BAM_SVCALLING_GRIDSS.out.vcf)
        //     .mix(gridss_existing_outputs)
        // versions = versions.mix(BAM_SVCALLING_GRIDSS.out.versions)

        // if (params.tumor_only) {
        //     vcf_from_sv_calling_for_merge = vcf_from_gridss_gridss
        //         .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, vcf, tbi

        //     JUNCTION_FILTER(vcf_from_gridss_gridss)

        //     pon_filtered_sv_rds = Channel.empty().mix(JUNCTION_FILTER.out.pon_filtered_sv_rds)
        //     final_filtered_sv_rds = Channel.empty().mix(JUNCTION_FILTER.out.final_filtered_sv_rds)
        //     final_filtered_sv_rds_for_merge = final_filtered_sv_rds
        //         .map { it -> [ it[0].patient, it[1] ] } // meta.patient, rds
        // } else {
        //     //somatic filter for GRIDSS
        //     BAM_SVCALLING_GRIDSS_SOMATIC(vcf_from_gridss_gridss)

        //     versions = versions.mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.versions)
        //     vcf_somatic_high_conf = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_high_confidence)
        //     vcf_from_sv_calling_for_merge = vcf_somatic_high_conf
        //         .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, vcf, tbi

        //     unfiltered_som_sv = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_all)
        //     unfiltered_som_sv_for_merge = unfiltered_som_sv
        //         .map { it -> [ it[0].patient, it[1] ] } // meta.patient, vcf
        // }
    }
    
    emit:
    vcf_from_gridss_gridss

}