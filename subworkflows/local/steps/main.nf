include { BAM_MSISENSORPRO } from '../bam_msisensorpro/main'
workflow MSISENSORPRO_STEP {
    take: 
    inputs_unlaned
    alignment_bams_final  // [meta, bam, bai]
    tools_used
    msisensorpro_scan

    main:
    msisensorpro_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.msi] }.filter { !it[1].isEmpty() }.unique { meta, msiPath -> meta.sample }
    msisensorpro_existing_outputs_germline = inputs_unlaned.map { it -> [it.meta, it.msi_germline] }.filter { !it[1].isEmpty() }.unique { meta, msiPath -> meta.sample }
    msi_from_msisensorpro = msisensorpro_existing_outputs
    germline_msi_from_msisensorpro = msisensorpro_existing_outputs_germline
    versions = Channel.empty()

    if (tools_used.contains("all") || tools_used.contains("msisensorpro") && !params.tumor_only) {

        bam_msi_inputs = inputs_unlaned.filter { it.msi.isEmpty() }.map { it -> [it.meta.sample + "___sep___" + it.meta.patient] }.unique()

        bam_msi = alignment_bams_final // reduped to have all tumor/normal pairs (even if normals are reused/duped)
            .map { sample, meta, bam, bai ->
                [ sample + "___sep___" + meta.patient, meta, bam, bai ]
            }
            .join(bam_msi_inputs)
            .map { it -> [ it[1], it[2], it[3] ] } // meta, bam, bai

        if (params.tumor_only) {
            bam_msi_status = bam_msi.branch{
                tumor:  it[0].status.toString() == "1"
            }

            // add empty arrays to stand-in for normals
            bam_msi_pair = bam_msi_status.tumor.map{ meta, bam, bai -> [ meta + [tumor_id: meta.sample], [], [], bam, bai, [] ] }
        } else {
            // getting the tumor and normal cram files separated
            bam_msi_status = bam_msi.branch{
                normal: it[0].status.toString() == "0"
                tumor:  it[0].status.toString() == "1"
            }

            // All normal samples
            bam_msi_normal_for_crossing = bam_msi_status.normal.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

            // All tumor samples
            bam_msi_tumor_for_crossing = bam_msi_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

            // Crossing the normal and tumor samples to create tumor and normal pairs
            bam_msi_pair = bam_msi_normal_for_crossing.cross(bam_msi_tumor_for_crossing)
                .map { normal, tumor ->
                    def meta = [:]

                    meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                    meta.normal_id  = normal[1].sample
                    meta.patient    = normal[0]
                    meta.sex        = normal[1].sex
                    meta.tumor_id   = tumor[1].sample

                    [ meta, normal[2], normal[3], tumor[2], tumor[3], [] ] // add empty array for intervals
            }
        }

        BAM_MSISENSORPRO(
            bam_msi_pair,
            msisensorpro_scan
        )

        msi_from_msisensorpro = msi_from_msisensorpro
            .mix(BAM_MSISENSORPRO.out.msi_somatic)
            .mix(msisensorpro_existing_outputs)
        versions = versions.mix(BAM_MSISENSORPRO.out.versions)

        if (!params.tumor_only) {
            germline_msi_from_msisensorpro = germline_msi_from_msisensorpro
                .mix(BAM_MSISENSORPRO.out.msi_germline)
                .mix(msisensorpro_existing_outputs_germline)
        }

    }

    emit:
    msi_from_msisensorpro
    germline_msi_from_msisensorpro
    versions
}

include { BAM_AMBER} from '../bam_amber/main'
workflow AMBER_STEP {
    take:
    inputs_unlaned
    alignment_bams_final  // [meta, bam, bai]
    tools_used

    main:
    versions = Channel.empty()
    // AMBER
    // ##############################
    amber_existing_outputs_hets = inputs_unlaned
        .map { it -> [it.meta, it.hets] }
        .filter { !it[1].isEmpty() }
        .unique()
        .branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }.tumor
    amber_existing_outputs_amber_dirs = inputs_unlaned
        .map { it -> [it.meta, it.amber_dir] }
        .filter { !it[1].isEmpty() }
        .unique()
        .branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }.tumor
    amber_dir = amber_existing_outputs_amber_dirs
    sites_from_het_pileups_wgs = amber_existing_outputs_hets

    if (tools_used.contains("all") || tools_used.contains("amber")) {
        bam_amber_inputs = inputs_unlaned.filter { it.hets.isEmpty() && it.amber_dir.isEmpty() }.map { it -> [it.meta.sample] }.unique()

        bam_amber_calling = alignment_bams_final
            .combine(bam_amber_inputs, by: 0) // .join(bam_amber_inputs)
            .map{ it -> [ it[1], it[2], it[3] ] } // meta, bam, bai



        

        // getting the tumor and normal cram files separated
        bam_amber_status = bam_amber_calling.branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }

        // All tumor samples
        bam_amber_tumor_for_crossing = bam_amber_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

        if (params.tumor_only) {
            // add empty arrays to stand-in for normals
            bam_amber_pair = bam_amber_status.tumor.map{ meta, bam, bai -> [ meta + [tumor_id: meta.sample], bam, bai, [], [] ] }
        } else {
            // All normal samples
            bam_amber_normal_for_crossing = bam_amber_status.normal.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }
            // Crossing the normal and tumor samples to create tumor and normal pairs
            bam_amber_pair = bam_amber_normal_for_crossing.cross(bam_amber_tumor_for_crossing)
                .map { normal, tumor ->
                    def meta = [:]
                    meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                    meta.normal_id  = normal[1].sample
                    meta.patient    = normal[0]
                    meta.sex        = normal[1].sex
                    meta.tumor_id   = tumor[1].sample

                    [ meta, tumor[2], tumor[3], normal[2], normal[3] ]
                }
        }


        BAM_AMBER(bam_amber_pair, inputs_unlaned)
        versions = versions.mix(BAM_AMBER.out.versions)

		// FIXME: Something is wrong with this when it's instantiated.
        amber_dir = Channel.empty()
            .mix(BAM_AMBER.out.amber_dir)
            .mix(amber_existing_outputs_amber_dirs)

        // amber_dir_for_merge = amber_dir
        //     .map { it -> [ it[0].patient, it[1] ] } // meta.patient, amber_dir

        sites_from_het_pileups_wgs = Channel.empty()
            .mix(BAM_AMBER.out.sites)
            .mix(amber_existing_outputs_hets)

        // hets_sites_for_merge = sites_from_het_pileups_wgs
        //     .map { it -> [ it[0].patient, it[1] ] } // meta.patient, hets
    }

    emit:
    amber_dir
    sites_from_het_pileups_wgs
}


include { BAM_FRAGCOUNTER as NORMAL_FRAGCOUNTER } from '../bam_fragCounter/main'
include { BAM_FRAGCOUNTER as TUMOR_FRAGCOUNTER } from '../bam_fragCounter/main'
workflow FRAGCOUNTER_STEP {
    take:
    inputs_unlaned
    alignment_bams_final  // [meta, bam, bai]
    tools_used

    main:
    fragcounter_existing_outputs = inputs_unlaned
        .map { it -> [it.meta , it.frag_cov] }
        .filter { !it[1].isEmpty() }
        .unique()
        .branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }

    tumor_frag_cov = fragcounter_existing_outputs.tumor
    normal_frag_cov = fragcounter_existing_outputs.normal
    versions = Channel.empty()

    if (tools_used.contains("all") || tools_used.contains("fragcounter")) {
        bam_fragcounter_inputs = inputs_unlaned.filter { it.frag_cov.isEmpty() }.map { it -> [it.meta.sample] }.unique().dump(tag: "bam_fragcounter_inputs", pretty: true)
        bam_fragcounter_calling = alignment_bams_final
            .dump(tag: "bam_fragcounter_calling_before_join", pretty: true)
            .join(bam_fragcounter_inputs)
            .map{ it -> [ it[1], it[2], it[3] ] } // meta, bam, bai
            .dump(tag: "bam_fragcounter_calling_after_join", pretty: true)

        
        // getting the tumor and normal bam files separated
        bam_fragcounter_status = bam_fragcounter_calling.branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }

        if (!params.tumor_only) {
            NORMAL_FRAGCOUNTER(bam_fragcounter_status.normal)
            normal_frag_cov = Channel.empty()
                .mix(NORMAL_FRAGCOUNTER.out.fragcounter_cov)
                .mix(fragcounter_existing_outputs.normal)
        }

        TUMOR_FRAGCOUNTER(bam_fragcounter_status.tumor.dump(tag: "TUMOR_FRAGCOUNTER", pretty: true))

        tumor_frag_cov = Channel.empty()
            .mix(TUMOR_FRAGCOUNTER.out.fragcounter_cov)
            .mix(fragcounter_existing_outputs.tumor)

        // Only need one versions because its just one program (fragcounter)
        versions = versions.mix(TUMOR_FRAGCOUNTER.out.versions)
    }

    emit:
    tumor_frag_cov
    normal_frag_cov

}

include { COV_DRYCLEAN as TUMOR_DRYCLEAN } from '../cov_dryclean/main'
include { COV_DRYCLEAN as NORMAL_DRYCLEAN } from '../cov_dryclean/main'
workflow DRYCLEAN_STEP {
    take:
    inputs_unlaned
    tumor_frag_cov_for_merge
    normal_frag_cov_for_merge
    tools_used

    main:
    versions = Channel.empty()
    dryclean_existing_outputs = inputs_unlaned
        .map { it -> [it.meta, it.dryclean_cov] }
        .filter { it -> !it[1].isEmpty() }
        .unique()
        .branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }
    dryclean_tumor_cov = dryclean_existing_outputs.tumor
    dryclean_normal_cov = dryclean_existing_outputs.normal

    cov_dryclean_inputs = inputs_unlaned
        .filter { it -> it.dryclean_cov.isEmpty() }
        .map { it -> [it.meta.sample, it.meta] }
        .unique()
        .branch{
            normal: it[1].status.toString() == "0"
            tumor:  it[1].status.toString() == "1"
        }
    cov_dryclean_tumor_input = tumor_frag_cov_for_merge
        .join(cov_dryclean_inputs.tumor)
        .map{ it -> [ it[1], it[2] ] } // meta, frag_cov
    
    if (tools_used.contains("all") || tools_used.contains("dryclean")) {

        // Dryclean for both tumor & normal
        TUMOR_DRYCLEAN(cov_dryclean_tumor_input)

        dryclean_tumor_cov = Channel.empty()
            .mix(TUMOR_DRYCLEAN.out.dryclean_cov)
            .mix(dryclean_existing_outputs.tumor)
        dryclean_tumor_cov_for_merge = dryclean_tumor_cov
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, dryclean_cov
            .unique{ patient, dryclean_cov -> patient }

        // Only need one versions because it's one program (dryclean)
        versions = versions.mix(TUMOR_DRYCLEAN.out.versions)

        if (!params.tumor_only) {
            cov_dryclean_normal_input = normal_frag_cov_for_merge
                .join(cov_dryclean_inputs.normal) // join will get rid of duplicates, so any normals matched to multiple tumors will be removed
                .map{ it -> [ it[1], it[2] ] } // meta, frag_cov

            NORMAL_DRYCLEAN(cov_dryclean_normal_input)

            dryclean_normal_cov = Channel.empty()
                .mix(NORMAL_DRYCLEAN.out.dryclean_cov)
                .mix(dryclean_existing_outputs.normal)


            sample_meta_map = inputs_unlaned.map { it -> [ it.meta.sample, it.meta ]}.unique()

            // dryclean_normal_cov_for_merge = dryclean_normal_cov
            //     .map { it -> [ it[0].patient, it[1] ] } // meta.patient, dryclean_cov

            dryclean_normal_cov_for_merge = dryclean_normal_cov
                .map { it -> [ it[0].sample, it[1] ] } // meta.sample, dryclean_cov
                .cross(sample_meta_map)
                .map { dryclean_normal, sample_meta ->
                    def meta_complete = sample_meta[1]
                    [ meta_complete.patient, dryclean_normal[1] ]
                }
                .unique{ patient, dryclean_cov -> patient }
                .dump(tag: "dryclean_normal_cov_for_merge", pretty: true)
        }
    }

    emit:
    dryclean_tumor_cov
    dryclean_normal_cov
    
}

include { COV_CBS } from '../cov_cbs/main'
workflow CBS_STEP {
    take:
    inputs_unlaned
    dryclean_tumor_cov_for_merge
    dryclean_normal_cov_for_merge
    tools_used

    main:
    versions = Channel.empty()
    cbs_inputs = inputs_unlaned
        .filter { it.seg.isEmpty() || it.nseg.isEmpty() }
        .map { it -> [it.meta.patient, it.meta] }
        .unique{ it -> it[1].patient + "___sep___" + it[1].sample } // unique by patient and sample
        .branch{
            normal: it[1].status.toString() == "0"
            tumor:  it[1].status.toString() == "1"
        }
    cbs_existing_seg = inputs_unlaned
        .map { it -> [it.meta, it.seg] }
        .filter { !it[1].isEmpty() }

    cbs_existing_nseg = inputs_unlaned
        .map { it -> [it.meta, it.nseg] }
        .filter { !it[1].isEmpty() }

    cbs_seg_rds = cbs_existing_seg
    cbs_nseg_rds = cbs_existing_nseg
    
    if (tools_used.contains("all") || tools_used.contains("cbs")) {
        cbs_tumor_input = cbs_inputs.tumor
            .join(dryclean_tumor_cov_for_merge)
            .map{ it -> [ it[0], it[1], it[2] ] } // meta.patient, meta, dryclean tumor cov

        if (params.tumor_only) {
            cov_cbs = cbs_tumor_input.map { patient, meta, tumor_cov -> [ meta, tumor_cov, [] ] }
        } else {
            cbs_normal_input = cbs_inputs.normal
                .join(dryclean_normal_cov_for_merge)
                .map{ it -> [ it[0], it[1], it[2] ] } // meta.patient, meta, dryclean normal cov

            cov_cbs = cbs_tumor_input.cross(cbs_normal_input)
                .map { tumor, normal ->
                    def meta = [:]
                        meta.id             = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                        meta.sample         = "${tumor[1].sample}".toString()
                        meta.normal_id      = normal[1].sample
                        meta.patient        = normal[0]
                        meta.sex            = normal[1].sex
                        meta.tumor_id       = tumor[1].sample

                        [ meta, tumor[2], normal[2] ]
                }.unique { meta, tumor_cov, normal_cov ->
                    meta.patient
                }
        }
        COV_CBS(cov_cbs)

        versions = versions.mix(COV_CBS.out.versions)

        cbs_seg_rds = Channel.empty()
            .mix(COV_CBS.out.cbs_seg_rds)
            .mix(cbs_existing_seg)

        cbs_nseg_rds = Channel.empty()
            .mix(COV_CBS.out.cbs_nseg_rds)
            .mix(cbs_existing_nseg)
    }

    emit:
    cbs_seg_rds
    cbs_nseg_rds

}

include { BAM_SAGE } from '../bam_sage/main'
include { BAM_SAGE_TUMOR_ONLY_FILTER } from '../bam_sage/main'
include { RESCUE_CH_HEME_STEP } from '../rescue_ch_step/main'
workflow VARIANT_CALLING_STEP {
    take:
    inputs_unlaned
    alignment_bams_final
    dbsnp_tbi
    known_indels_tbi
    tools_used

    main:
    versions = Channel.empty()

    dict = params.dict ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect() : Channel.empty()

    snv_somatic_existing_outputs = inputs_unlaned
        .map { it -> [it.meta, it.snv_somatic_vcf, it.snv_somatic_tbi] }
        .filter { !it[1].isEmpty() && !it[2].isEmpty()}
        .unique()
    
    snv_germline_existing_outputs = inputs_unlaned
        .map { it -> [it.meta, it.snv_germline_vcf, it.snv_germline_tbi] }
        .filter { !it[1].isEmpty() && !it[2].isEmpty()}
        .unique()

    filtered_somatic_vcf = snv_somatic_existing_outputs
    germline_vcf = snv_germline_existing_outputs
    
    

    filtered_somatic_vcf_tumoronly_outputs = inputs_unlaned
        .map { it -> [it.meta, it.snv_somatic_vcf_tumoronly_filtered, it.snv_somatic_vcf_tumoronly_filtered_tbi] }
        .filter { !it[1].isEmpty() && !it[2].isEmpty()}
        .unique()
    filtered_somatic_vcf_tumor_only = filtered_somatic_vcf_tumoronly_outputs

    rescue_ch_existing_outputs = inputs_unlaned
        .map { it -> [it.meta, it.snv_somatic_vcf_rescue_ch_heme, it.snv_somatic_vcf_rescue_ch_heme_tbi] }
        .filter { !it[1].isEmpty() && !it[2].isEmpty()}
        .unique()
    
    rescue_ch_inputs = inputs_unlaned
        .map { it -> [it.meta, it.snv_somatic_vcf_rescue_ch_heme, it.snv_somatic_vcf_rescue_ch_heme_tbi] }
        .filter { it[1].isEmpty() || it[1].isEmpty() }
        .map { it -> [it.meta.patient] }
        .unique()

    // Set is_heme based on is_retier_whitelist_junctions
    // params.is_heme = params.is_retier_whitelist_junctions
    is_paired = ! params.tumor_only
    is_tumor_only = ! is_paired
    is_heme = params.is_heme && is_tumor_only
    if (is_paired) {
        filtered_somatic_vcf = snv_somatic_existing_outputs
    } else if (is_tumor_only) {
        filtered_somatic_vcf = filtered_somatic_vcf_tumor_only
    } else if (is_heme) {
        filtered_somatic_vcf = rescue_ch_existing_outputs
    } else {
        error "Problem with parsing somatic variant calling outputs. Could not determine if run is paired, tumor-only or heme." 
    }


    // Filter out bams for which SNV calling has already been done
    if (params.tumor_only) {
        bam_snv_inputs = inputs_unlaned
            .filter { it.snv_somatic_vcf.isEmpty() && it.snv_somatic_vcf_tumoronly_filtered.isEmpty() }
            .map { it -> [it.meta.sample] }.unique()
    } else {
        bam_snv_inputs = inputs_unlaned
            .filter { it.snv_somatic_vcf.isEmpty() || it.snv_germline_vcf.isEmpty() }
            .map { it -> [it.meta.sample] }.unique()
    }

    if (tools_used.contains("all") || tools_used.contains("sage")) {
        
        bam_snv_calling = alignment_bams_final
            .combine(bam_snv_inputs, by: 0) // .join(bam_snv_inputs) // join will remove any duplicated normals (if more than one tumor corresponds to a normal)
            .map { it -> [ it[1], it[2], it[3] ] } // meta, bam, bai

        // getting the tumor and normal bams separated
        bam_snv_calling_status = bam_snv_calling.branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }

        if (params.tumor_only) {
            bam_snv_calling_pair = bam_snv_calling_status.tumor.map{ meta, bam, bai -> [ meta + [tumor_id: meta.sample], [], [], bam, bai ] }
        } else {
            // All normal samples
            bam_snv_calling_normal_for_crossing = bam_snv_calling_status.normal.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

            // All tumor samples
            bam_snv_calling_tumor_for_crossing = bam_snv_calling_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

            // Crossing the normal and tumor samples to create tumor and normal pairs
            bam_snv_calling_pair = bam_snv_calling_normal_for_crossing.cross(bam_snv_calling_tumor_for_crossing)
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

        dict_sage = dict.map{ id, dictPath -> dictPath }

        BAM_SAGE(
            bam_snv_calling_pair,
            dict_sage
        )

        versions = versions.mix(BAM_SAGE.out.versions)

        filtered_somatic_vcf_sage = Channel.empty()
            .mix(BAM_SAGE.out.sage_pass_filtered_somatic_vcf)
            .mix(snv_somatic_existing_outputs)


        tumor_only_filter_input = filtered_somatic_vcf_sage
            .map{ it ->
                [it[0].patient] + it.toList()
            }
            .join(
                filtered_somatic_vcf_tumoronly_outputs
                    .filter { it -> it[1].isEmpty() }
                    .map { it -> [it[0].patient] }
                    .unique()
            )
            .map { it[1..-1] }

        if (is_paired) {
            germline_vcf = Channel.empty()
                .mix(BAM_SAGE.out.sage_germline_vcf)
                .mix(snv_germline_existing_outputs)

            germline_vcf_for_merge = germline_vcf
                .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, germline snv vcf, tbi
        }

        if (is_tumor_only) {
            BAM_SAGE_TUMOR_ONLY_FILTER(
                tumor_only_filter_input,
                dbsnp_tbi,
                known_indels_tbi
            )

            filtered_somatic_vcf_tumor_only = Channel.empty()
                .mix(BAM_SAGE_TUMOR_ONLY_FILTER.out.sage_filtered_vcf)
                .mix(
                    filtered_somatic_vcf_tumoronly_outputs
                    .filter { it -> !it[1].isEmpty() }
                )

            filtered_somatic_vcf = filtered_somatic_vcf_tumor_only

            if (is_heme) {
                heme_rescue_input = filtered_somatic_vcf_sage
                    .map { [it[0].patient, it[0], it[1], it[2]] } // meta.patient, meta, sage vcf, sage tbi
                    .join(
                        filtered_somatic_vcf_tumor_only
                            .map { [it[0].patient, it[1], it[2]]  } // meta.patient, vcf, tbi
                    )
                    .join(rescue_ch_inputs)
                    .map { it[1..-1] } // meta, sage vcf, sage tbi, tumoronly vcf, tumoronly tbi

                RESCUE_CH_HEME_STEP(
                    heme_rescue_input
                )
                filtered_somatic_vcf = RESCUE_CH_HEME_STEP.out.sage_tumor_only_rescue_ch_vcf
                    .mix(rescue_ch_existing_outputs)
            }


        } else {
            filtered_somatic_vcf = filtered_somatic_vcf_sage
        }

    }

    emit:
    filtered_somatic_vcf
    germline_vcf

}

// SNPEFF
include { VCF_SNPEFF as VCF_SNPEFF_SOMATIC } from '../vcf_snpeff/main'
include { VCF_SNPEFF as VCF_SNPEFF_GERMLINE } from '../vcf_snpeff/main'
workflow VARIANT_ANNOTATION_STEP {
    take:
    inputs_unlaned
    filtered_somatic_vcf_for_merge
    germline_vcf_for_merge
    snpeff_db_full
    snpeff_cache
    tools_used

    main:
    variant_somatic_ann_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.variant_somatic_ann] }.filter { !it[1].isEmpty() }
    variant_somatic_bcf_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.variant_somatic_bcf] }.filter { !it[1].isEmpty() }
    variant_germline_ann_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.variant_germline_ann] }.filter { !it[1].isEmpty() }
    variant_germline_bcf_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.variant_germline_bcf] }.filter { !it[1].isEmpty() }
    snv_somatic_annotations = variant_somatic_ann_existing_outputs
    snv_somatic_bcf_annotations = variant_somatic_bcf_existing_outputs

    snv_germline_annotations = variant_germline_ann_existing_outputs
    snv_germline_bcf_annotations = variant_germline_bcf_existing_outputs

    variant_somatic_ann_inputs = inputs_unlaned
            .filter { it.variant_somatic_ann.isEmpty() || it.variant_somatic_bcf.isEmpty() }
            .map { it -> [it.meta.patient, it.meta + [id: it.meta.sample ]] }.unique()

    variant_ann_input_somatic = variant_somatic_ann_inputs
        .join(filtered_somatic_vcf_for_merge)
        .map { it -> [ it[1], it[2], it[3] ] } // meta, filtered (pass and tumor-only filter) snvs, tbi

    if (tools_used.contains("all") || tools_used.contains("snpeff")) {

        VCF_SNPEFF_SOMATIC(
            variant_ann_input_somatic,
            snpeff_db_full,
            snpeff_cache
        )

        snv_somatic_annotations = Channel.empty()
            .mix(VCF_SNPEFF_SOMATIC.out.snpeff_vcf)
            .mix(variant_somatic_ann_existing_outputs)

        snv_somatic_bcf_annotations = Channel.empty()
            .mix(VCF_SNPEFF_SOMATIC.out.snpeff_bcf)
            .mix(variant_somatic_bcf_existing_outputs)

        if (!params.tumor_only) {
            variant_germline_ann_inputs = inputs_unlaned
                .filter { it.variant_germline_ann.isEmpty() || it.variant_germline_bcf.isEmpty() }
                .map { it -> [it.meta.patient, it.meta + [id: it.meta.sample ]] }

            variant_ann_input_germline = variant_germline_ann_inputs
                .join(germline_vcf_for_merge)
                .map { it -> [ it[1], it[2], it[3] ] } // meta, germline snvs, tbi

            variant_germline_ann_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.variant_germline_ann] }.filter { !it[1].isEmpty() }
            variant_germline_bcf_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.variant_germline_bcf] }.filter { !it[1].isEmpty() }

            // FIXME: Duplicated normals will be run multiple times if they correspond to different tumors.
            VCF_SNPEFF_GERMLINE(
                variant_ann_input_germline,
                snpeff_db_full,
                snpeff_cache
            )

            snv_germline_annotations = Channel.empty()
                .mix(VCF_SNPEFF_GERMLINE.out.snpeff_vcf)
                .mix(variant_germline_ann_existing_outputs)

            snv_germline_bcf_annotations = Channel.empty()
                .mix(VCF_SNPEFF_GERMLINE.out.snpeff_bcf)
                .mix(variant_germline_bcf_existing_outputs)
        }
    }

    emit:
    snv_somatic_annotations
    snv_somatic_bcf_annotations
    snv_germline_annotations
    snv_germline_bcf_annotations
}

// COBALT
include { BAM_COBALT } from '../bam_cobalt/main'
workflow COBALT_STEP {
    take:
    inputs_unlaned
    alignment_bams_final  // [meta, bam, bai]
    tools_used

    main:
    versions = Channel.empty()
    // Existing
    cobalt_existing_outputs_cobalt_dirs = inputs_unlaned
        .map { it -> [it.meta, it.cobalt_dir] }
        .filter { !it[1].isEmpty() }
        .unique()
        .branch{
            normal: it[0].status.toString() == "0"
            tumor:  it[0].status.toString() == "1"
        }
    // Emit
    cobalt_dir = cobalt_existing_outputs_cobalt_dirs.tumor.unique{ it -> it[0].patient + "___sep___" + it[0].sample }
    
    // Inputs to run
    bam_cobalt_inputs = inputs_unlaned.filter { it.cobalt_dir.isEmpty() }.map { it -> [it.meta.sample] }.unique()
    bam_cobalt_calling = alignment_bams_final
        .combine(bam_cobalt_inputs, by: 0) // .join(bam_cobalt_inputs) // join will remove any duplicated normals (if more than one tumor corresponds to a normal)
        .map{ it -> [ it[1], it[2], it[3] ] } // meta, bam, bai
    
    // getting the tumor and normal cram files separated
    bam_cobalt_status = bam_cobalt_calling.branch{
        normal: it[0].status.toString() == "0"
        tumor:  it[0].status.toString() == "1"
    }

    // All tumor samples
    bam_cobalt_tumor_for_crossing = bam_cobalt_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

    if (params.tumor_only) {
        // add empty arrays to stand-in for normals
        bam_cobalt_pair = bam_cobalt_status.tumor.map{ meta, bam, bai -> [ meta + [tumor_id: meta.sample], bam, bai, [], [] ] }
    } else {
        // All normal samples
        bam_cobalt_normal_for_crossing = bam_cobalt_status.normal.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }
        // Crossing the normal and tumor samples to create tumor and normal pairs
        bam_cobalt_pair = bam_cobalt_normal_for_crossing.cross(bam_cobalt_tumor_for_crossing)
            .map { normal, tumor ->
                def meta = [:]
                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ meta, tumor[2], tumor[3], normal[2], normal[3] ]
            }
    }

    if (tools_used.contains("all") || tools_used.contains("cobalt")) {
        BAM_COBALT(bam_cobalt_pair)
        versions = versions.mix(BAM_COBALT.out.versions)

        cobalt_dir = Channel.empty()
            .mix(BAM_COBALT.out.cobalt_dir)
            .mix(cobalt_existing_outputs_cobalt_dirs)
    }
    

    emit:
    cobalt_dir
    versions
}

// PURPLE
include { BAM_COV_PURPLE } from '../bam_cov_purple/main'
workflow PURPLE_STEP {
    take:
    inputs_unlaned
    germline_vcf_for_merge
    filtered_somatic_vcf_for_merge
    cobalt_dir_for_merge
    amber_dir_for_merge
    vcf_from_sv_calling_for_merge
    tools_used

    main:
    versions = Channel.empty()
    // Existing
    purple_existing_outputs_ploidy = inputs_unlaned
        .branch{ it -> 
            tumor: it.meta.status.toString() == "1"
        }
        .tumor
        .map { it -> 
            [it.meta, it.ploidy] 
        }
        .filter { it -> 
            !(it[1] instanceof List && it[1].isEmpty())
        }.unique()
    purple_existing_outputs_purity = inputs_unlaned
        .branch{it -> 
            tumor: it.meta.status.toString() == "1"
        }
        .tumor
        .map { it -> 
            [it.meta, it.purity] 
        }
        .filter { it -> 
            !(it[1] instanceof List && it[1].isEmpty())
        }
        .unique()
    purple_existing_outputs = inputs_unlaned
        .branch{it -> 
            tumor: it.meta.status.toString() == "1"
        }
        .tumor
        .map { it -> 
            [it.meta, it.purple_pp_best_fit, it.purple_pp_range] 
        }
        .filter { it -> 
            !it[1].isEmpty() && !it[2].isEmpty()
        }
        .unique()
    // Emit
    purity = purple_existing_outputs_purity
    ploidy = purple_existing_outputs_ploidy
    purple_out = purple_existing_outputs

    // need a channel with patient and meta for merging with rest
    purple_inputs_for_merge = inputs_unlaned
        .filter { it -> it.ploidy.isEmpty() }
        .map { it -> [it.meta.patient, it.meta - it.meta.subMap('tumor_id')] }
        .unique()

    meta_purple = purple_inputs_for_merge
        .branch{
            normal: it[1].status.toString() == "0"
            tumor:  it[1].status.toString() == "1"
        }
        .tumor
        .map {
        patient, meta ->
        [patient, meta + [tumor_id: meta.sample, id: meta.sample] ]
    }

    purple_inputs_snv_germline = Channel.empty()
    if (!params.tumor_only) {
        if (params.purple_use_smlvs) {
            purple_inputs_snv_germline = purple_inputs_for_merge
                .join(germline_vcf_for_merge)
                .map { it -> [ it[0], it[2], it[3] ] } // patient, vcf, tbi
        }
    }

    purple_inputs_cobalt_dir = purple_inputs_for_merge
        .join(cobalt_dir_for_merge)
        .map { it -> [ it[0], it[2] ] } // patient, cobalt_dir

    purple_inputs_amber_dir = purple_inputs_for_merge
        .join(amber_dir_for_merge)
        .map { it -> [ it[0], it[2] ] } // patient, amber_dir

    if (params.purple_use_svs) {
        purple_inputs_sv = purple_inputs_for_merge
            .join(vcf_from_sv_calling_for_merge)
            .map { it -> [ it[0], it[2], it[3] ] } // patient, vcf, tbi
    }

    if (params.purple_use_smlvs) {
        println "Using Purple small variants"
        purple_inputs_snv = purple_inputs_for_merge
            .join(filtered_somatic_vcf_for_merge)
            .map { it -> [ it[0], it[2], it[3] ] } // patient, vcf, tbi
    }

    purple_inputs = meta_purple
        .join(purple_inputs_amber_dir)
        .join(purple_inputs_cobalt_dir)
        .map { patient, meta, amber_dir, cobalt_dir ->
            [meta, amber_dir, cobalt_dir, [], [], [], [], [], []]
        }

    if (params.tumor_only) {
        if (params.purple_use_svs && params.purple_use_smlvs) {
            purple_inputs = meta_purple
            .join(purple_inputs_amber_dir)
            .join(purple_inputs_cobalt_dir)
            .join(purple_inputs_sv)
            .join(purple_inputs_snv)
            .map { patient, meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi ->
                [meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi, [], []]
            }
        }
    } else {
        if (params.purple_use_svs && params.purple_use_smlvs) {
            println "Purple SVS and small variants are being used"
            purple_inputs = meta_purple
                .join(purple_inputs_amber_dir)
                .join(purple_inputs_cobalt_dir)
                .join(purple_inputs_sv)
                .join(purple_inputs_snv)
                .join(purple_inputs_snv_germline)
                .map { patient, meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi, germ_snv_vcf, germ_snv_tbi ->
                    [meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi, germ_snv_vcf, germ_snv_tbi]
                }
        }
    }
    
    if (tools_used.contains("all") || tools_used.contains("purple")) {

        BAM_COV_PURPLE(
            purple_inputs
        )

        versions = versions.mix(BAM_COV_PURPLE.out.versions)
        // // Allow purity to be overwritten by provided values from inputs
        purity = Channel.empty()
            .mix(BAM_COV_PURPLE.out.purity)
            .mix(purple_existing_outputs_purity)
            .unique{ it -> it[0].patient}


        // def use_existing_purity = params.get("use_existing_purity", false)
        // purity = purple_existing_outputs_purity
        //     .map {
        //         [ it[0].patient ] + [ it.toList() ] // [patient, [meta, purity]]
        //     }
        //     .join(
        //         BAM_COV_PURPLE.out.purity
        //             .map {
        //                 [ it[0].patient ] + [ it.toList() ] // [patient, [meta, purity]]
        //             }
        //         ,
        //         remainder: true
        //     )
        //     .map { it -> // patient, meta_existing, purity_existing, meta_purple, purity_purple ->
        //         def (key, existing, purple) = (it + [null, null])[0..2]
        //         def (meta_existing, purity_existing) = existing ?: [null, null]
        //         def (meta_purple, purity_purple) = purple ?: [null, null]
        //         def purity_out = purity_existing ?: purity_purple
        //         def meta_out = meta_existing ?: meta_purple
        //         if (!use_existing_purity) {
        //             purity_out = purity_purple ?: purity_existing
        //         }
        //         [ meta_out, purity_out ]
        //     }
        //     .dump(tag: "purity", pretty: true)


        // // Allow ploidy to be overwritten by provided values from inputs
        ploidy = Channel.empty()
            .mix(BAM_COV_PURPLE.out.ploidy)
            .mix(purple_existing_outputs_ploidy)
            .unique{ it -> it[0].patient}

        // def use_existing_ploidy = params.get("use_existing_ploidy", false)
        // ploidy = purple_existing_outputs_ploidy
        //     .map {
        //         [ it[0].patient ] + [ it.toList() ]
        //     }
        //     .join(
        //         BAM_COV_PURPLE.out.ploidy
        //             .map {
        //                 [ it[0].patient ] + [ it.toList() ]
        //             }
        //         ,
        //         remainder: true
        //     )
        //     .map { it -> // patient, meta_existing, ploidy_existing, meta_purple, ploidy_purple ->
        //         def (key, existing, purple) = (it + [null, null])[0..2]
        //         def (meta_existing, ploidy_existing) = existing ?: [null, null]
        //         def (meta_purple, ploidy_purple) = purple ?: [null, null]
        //         def ploidy_out = ploidy_existing ?: ploidy_purple
        //         def meta_out = meta_existing ?: meta_purple
        //         if (!use_existing_ploidy) {
        //             ploidy_out = ploidy_purple ?: ploidy_existing
        //         }
        //         [ meta_out, ploidy_out ]
        //     }
        //     .dump(tag: "ploidy", pretty: true)

    }

    emit:
    purity
    ploidy
    purple_out

}

// JaBbA
include { COV_JUNC_TUMOR_ONLY_JABBA as JABBA_TUMOR_ONLY } from '../jabba/main'
include { COV_JUNC_JABBA as JABBA } from '../jabba/main'
include { RETIER_JUNCTIONS } from '../jabba/main'
workflow JABBA_STEP {
    take:
    inputs_unlaned
    dryclean_tumor_cov_for_merge
    hets_sites_for_merge
    cbs_seg_for_merge
    cbs_nseg_for_merge
    purity_for_merge
    ploidy_for_merge
    vcf_from_sv_calling_for_merge
    vcf_raw_from_gridss_gridss
    final_filtered_sv_rds_for_merge
    unfiltered_som_sv_for_merge
    tools_used
    


    main:
    versions = Channel.empty()
    // Existing
    jabba_rds_existing_outputs = inputs_unlaned.map { it -> [it.meta + [id: it.meta.sample], it.jabba_rds] }.filter { !it[1].isEmpty() }
    jabba_gg_existing_outputs = inputs_unlaned.map { it -> [it.meta + [id: it.meta.sample], it.jabba_gg] }.filter { !it[1].isEmpty() }
    // Emit
    jabba_rds = jabba_rds_existing_outputs
    jabba_gg = jabba_gg_existing_outputs


    jabba_inputs = inputs_unlaned.filter { (it.jabba_gg.isEmpty() || it.jabba_rds.isEmpty()) && it.meta.status.toString() == "1"}.map { it -> [it.meta.patient, it.meta] }.unique()
    if (tools_used.contains("all") || tools_used.contains("jabba")) {
        

		// Dev block to retier either vcf or filtered retiered junctions
		is_final_filtered_sv_rds_for_merge_retiered = params.is_retier_whitelist_junctions && params.tumor_only
		is_vcf_from_sv_calling_for_merge_retiered = params.is_retier_whitelist_junctions && ! params.tumor_only

		// The variable below will get propagated to jabba if retiering is not done to ! params.tumor_only block
		jabba_vcf_from_sv_calling_for_merge = vcf_from_sv_calling_for_merge

        vcf_raw_from_gridss_gridss_for_merge = vcf_raw_from_gridss_gridss
            .map {
                [ it[0].patient ] + it[1..-1].toList() // patient, vcf_raw, vcf_tbi_raw
            }

		// Variable exists in case retiering is done
		untiered_junctions_for_merge = vcf_from_sv_calling_for_merge
        println "is_final_filtered_sv_rds_for_merge_retiered: ${is_final_filtered_sv_rds_for_merge_retiered}"
		if (is_final_filtered_sv_rds_for_merge_retiered) {
			untiered_junctions_for_merge = final_filtered_sv_rds_for_merge
		}
		if (params.is_retier_whitelist_junctions) {
			untiered_junctions_input = jabba_inputs
				.join(untiered_junctions_for_merge)
                .join(vcf_raw_from_gridss_gridss_for_merge)
				.map { it -> [ it[1], it[2], it[3] ] } // meta, (vcf or rds), vcf_raw

			RETIER_JUNCTIONS(untiered_junctions_input)
			retiered_junctions_output_for_merge = Channel.empty()
				.mix(RETIER_JUNCTIONS.out.retiered_junctions)
				.map { meta, rds -> [ meta.patient, rds ] } // meta.patient, retiered junctions
		}

		if (is_vcf_from_sv_calling_for_merge_retiered) {
			jabba_vcf_from_sv_calling_for_merge = retiered_junctions_output_for_merge
		} else if (is_final_filtered_sv_rds_for_merge_retiered) {
			final_filtered_sv_rds_for_merge = retiered_junctions_output_for_merge
		}

        jabba_inputs_sv = jabba_vcf_from_sv_calling_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, (vcf or rds if retiered)

        if (params.tumor_only) {
            jabba_inputs_junction_filtered_sv = final_filtered_sv_rds_for_merge
                .join(jabba_inputs)
                .map { it -> [ it[0], it[1] ] } // meta.patient, pon and gnomaD filtered sv (for tumor only)
        } else {
            jabba_inputs_unfiltered_sv = unfiltered_som_sv_for_merge
                .join(jabba_inputs)
                .map { it -> [ it[0], it[1] ] } // meta.patient, unfiltered_som_sv
        }

        jabba_inputs_hets = hets_sites_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, hets

        jabba_inputs_cov = dryclean_tumor_cov_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, cov

        jabba_inputs_purity = purity_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, purity

        jabba_inputs_ploidy = ploidy_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, ploidy

        jabba_inputs_cbs_seg = cbs_seg_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, cbs_seg

        jabba_inputs_cbs_nseg = cbs_nseg_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, cbs_nseg

        

        if (params.tumor_only) {
            jabba_input = jabba_inputs
                .join(jabba_inputs_cov)
                .join(jabba_inputs_hets)
                .join(jabba_inputs_purity)
                .join(jabba_inputs_ploidy)
                .join(jabba_inputs_cbs_seg)
                .join(jabba_inputs_cbs_nseg)
                .join(jabba_inputs_junction_filtered_sv)
                .map{ _patient, meta, cov, hets, purity, ploidy, seg, nseg, junction ->
                    [
                        meta,
                        junction,
                        cov,
                        [],
                        hets,
                        purity,
                        ploidy,
                        seg,
                        nseg
                    ]
                }
        } else {
            // join all previous outputs to be used as input for jabba
            jabba_input = jabba_inputs
                .join(jabba_inputs_cov)
                .join(jabba_inputs_hets)
                .join(jabba_inputs_purity)
                .join(jabba_inputs_ploidy)
                .join(jabba_inputs_cbs_seg)
                .join(jabba_inputs_cbs_nseg)
                .join(jabba_inputs_sv)
                .join(jabba_inputs_unfiltered_sv)
                .map{ _patient, meta, cov, hets, purity, ploidy, seg, nseg, junction, j_supp ->
                    [
                        meta,
                        junction,
                        cov,
                        j_supp,
                        hets,
                        purity,
                        ploidy,
                        seg,
                        nseg
                    ]
                }
        }

        jabbaOut = Channel.empty()
        if (params.tumor_only) {
            jabbaOut = JABBA_TUMOR_ONLY(jabba_input)
        } else {
            jabbaOut = JABBA(jabba_input)
        }
        

        jabba_rds = Channel.empty()
            .mix(jabbaOut.jabba_rds)
            .mix(jabba_rds_existing_outputs)
        

        jabba_gg = Channel.empty()
            .mix(jabbaOut.jabba_gg)
            .mix(jabba_gg_existing_outputs)
        versions = versions.mix(jabbaOut.versions)
    }

    emit:
    jabba_rds
    jabba_gg
}


// Allelic CN
include { COV_GGRAPH_NON_INTEGER_BALANCE as NON_INTEGER_BALANCE } from '../allelic_cn/main'
workflow NON_INTEGER_BALANCE_STEP {
    take:
    inputs_unlaned
    jabba_gg_for_merge
    hets_sites_for_merge
    dryclean_tumor_cov_for_merge
    tools_used

    main:
    versions = Channel.empty()
    bwa = WorkflowNfcasereports.create_file_channel(params.bwa)
    index_alignment = bwa

    // Existing
    non_integer_balance_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.ni_balanced_gg] }.filter { !it[1].isEmpty() }

    // Emit
    non_integer_balance_balanced_gg = non_integer_balance_existing_outputs

    // Inputs
    non_integer_balance_inputs = inputs_unlaned.filter { it.ni_balanced_gg.isEmpty() }.map { it -> [it.meta.patient, it.meta + [id: it.meta.sample]] }

    non_integer_balance_inputs_jabba_gg = jabba_gg_for_merge
        .join(non_integer_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, jabba ggraph

    non_integer_balance_inputs_hets = hets_sites_for_merge
        .join(non_integer_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, hets

    non_integer_balance_inputs_cov = dryclean_tumor_cov_for_merge
        .join(non_integer_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, cov

    if (params.tumor_only) {
        non_integer_balance_inputs = non_integer_balance_inputs
            .join(non_integer_balance_inputs_jabba_gg)
            .join(non_integer_balance_inputs_cov)
            .map{ _patient, meta, rds, cov -> [ meta, rds, cov, [] ] }

    } else {
        non_integer_balance_inputs = non_integer_balance_inputs
            .join(non_integer_balance_inputs_jabba_gg)
            .join(non_integer_balance_inputs_hets)
            .join(non_integer_balance_inputs_cov)
            .map{ _patient, meta, rds, hets, cov -> [ meta, rds, cov, hets ] }
    }

     if (tools_used.contains("all") || tools_used.contains("non_integer_balance")) {    
        NON_INTEGER_BALANCE(non_integer_balance_inputs, index_alignment)
        versions = Channel.empty().mix(NON_INTEGER_BALANCE.out.versions)

        non_integer_balance_balanced_gg = Channel.empty()
            .mix(NON_INTEGER_BALANCE.out.non_integer_balance_balanced_gg)
            .mix(non_integer_balance_existing_outputs)
    }
    emit:
    non_integer_balance_balanced_gg

}

include { COV_GGRAPH_LP_PHASED_BALANCE as LP_PHASED_BALANCE } from '../allelic_cn/main'
workflow LP_PHASED_BALANCE_STEP {
    take:
    inputs_unlaned
    non_integer_balance_balanced_gg_for_merge
    hets_sites_for_merge
    tools_used

    main:
    // Existing
    lp_phased_balance_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.lp_balanced_gg] }.filter { !it[1].isEmpty() }
    lp_phased_balance_balanced_gg = lp_phased_balance_existing_outputs

    // Inputs
    lp_phased_balance_inputs = inputs_unlaned.filter { it.lp_balanced_gg.isEmpty() }.map { it -> [it.meta.patient, it.meta + [id: it.meta.sample]] }

    lp_phased_balance_inputs_ni_balanced_gg = non_integer_balance_balanced_gg_for_merge
        .join(lp_phased_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, non integer balanced ggraph
    lp_phased_balance_inputs_hets = hets_sites_for_merge
        .join(lp_phased_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, hets

    lp_phased_balance_inputs = lp_phased_balance_inputs
        .join(lp_phased_balance_inputs_ni_balanced_gg)
        .join(lp_phased_balance_inputs_hets)
        .map{ _patient, meta, balanced_gg, hets -> [ meta, balanced_gg, hets ] }

    if (tools_used.contains("all") || tools_used.contains("lp_phased_balance")) {

        LP_PHASED_BALANCE(lp_phased_balance_inputs)

        lp_phased_balance_balanced_gg = Channel.empty()
            .mix(LP_PHASED_BALANCE.out.lp_phased_balance_balanced_gg)
            .mix(lp_phased_balance_existing_outputs)
    }

    emit:
    lp_phased_balance_balanced_gg

}

// Events
include { GGRAPH_EVENTS as EVENTS } from '../events/main'
workflow EVENTS_STEP {
    take:
    inputs_unlaned
    non_integer_balance_balanced_gg_for_merge
    tools_used

    main:
    // Existing
    events_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.events] }.filter { !it[1].isEmpty() }
    // Emit
    events = events_existing_outputs

    // Input
    events_inputs = inputs_unlaned.filter { it.events.isEmpty() }.map { it -> [it.meta.patient, it.meta + [id: it.meta.sample]] }
    events_input_non_integer_balance = non_integer_balance_balanced_gg_for_merge
            .join(events_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, balanced_gg
    events_input = events_inputs
        .join(events_input_non_integer_balance)
        .map{ patient, meta, balanced_gg -> [ meta, balanced_gg ] }

    if (tools_used.contains("all") || tools_used.contains("events")) {

        EVENTS(events_input)

        // events_versions = Channel.empty().mix(EVENTS.out.versions)
        events = Channel.empty()
            .mix(EVENTS.out.events_output)
            .mix(events_existing_outputs)
    }

    emit:
    events
}


// Fusions
include { GGRAPH_FUSIONS as FUSIONS } from '../fusions/main'
workflow FUSIONS_STEP {
    take:
    inputs_unlaned
    non_integer_balance_balanced_gg_for_merge
    vcf_from_sv_calling_for_merge
    tools_used

    main:
    // Existing
    fusions_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.fusions] }.filter { !it[1].isEmpty() }
    altedge_annotations_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.altedge_annotations] }.filter { !it[1].isEmpty() }

    // Emit
    fusions = fusions_existing_outputs
    altedge_annotations = altedge_annotations_existing_outputs

    fusions_inputs = inputs_unlaned.filter { it.fusions.isEmpty() }.map { it -> [it.meta.patient, it.meta + [id: it.meta.sample]] }
    if (tools_used.contains("non_integer_balance") || tools_used.contains("all")) {
        fusions_input_non_integer_balance = non_integer_balance_balanced_gg_for_merge
            .join(fusions_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, balanced_gg
        fusions_input = fusions_inputs
            .join(fusions_input_non_integer_balance)
            .map{ _patient, meta, balanced_gg -> [ meta, balanced_gg, [], [] ] }
    } else {
        fusions_input_sv = vcf_from_sv_calling_for_merge
            .join(fusions_inputs)
            .map { it -> [ it[0], it[1], it[2] ] } // meta.patient, vcf, vcf_tbi
        fusions_input = fusions_inputs
            .join(fusions_input_sv)
            .map{ _patient, meta, sv, sv_tbi -> [ meta, [], sv, sv_tbi ] }
    }
    if (tools_used.contains("all") || tools_used.contains("fusions")) {

        FUSIONS(fusions_input)
        fusions = Channel.empty()
            .mix(FUSIONS.out.fusions_output)
            .mix(fusions_existing_outputs)
        altedge_annotations = Channel.empty()
            .mix(FUSIONS.out.altedge_annotations)
            .mix(altedge_annotations_existing_outputs)

        // versions = Channel.empty().mix(FUSIONS.out.versions)
    }

    emit:
    fusions
    altedge_annotations
}

// SNV MULTIPLICITY
include { VCF_SNV_MULTIPLICITY } from '../vcf_snv_multiplicity/main'
workflow MULTIPLICITY_STEP {
    take:
    inputs_unlaned
    annotated_vcf_ffpe_impact_or_snpeff_for_merge
    non_integer_balance_balanced_gg_for_merge
    hets_sites_for_merge
    dryclean_tumor_cov_for_merge
    snv_germline_annotations_for_merge
    tools_used

    main:
    // Existing
    snv_multiplicity_fields = inputs_unlaned.map { it -> [it.meta, it.snv_multiplicity, it.snv_multiplicity_hets, it.snv_multiplicity_germline] }
    if (!params.tumor_only) {
        snv_multiplicity_existing_outputs = snv_multiplicity_fields
            .filter { !it[1].isEmpty() && !it[2].isEmpty() && !it[3].isEmpty() }
        snv_multiplicity_to_run = snv_multiplicity_fields
            .filter { it[1].isEmpty() || it[2].isEmpty() || it[3].isEmpty() }
    } else {
        snv_multiplicity_existing_outputs = snv_multiplicity_fields
            .filter { !it[1].isEmpty() && !it[2].isEmpty() }
        snv_multiplicity_to_run = snv_multiplicity_fields
            .filter { it[1].isEmpty() || it[2].isEmpty() }
    }    
    // Emit
    snv_multiplicity = snv_multiplicity_existing_outputs

    // Inputs
    snv_multiplicity_inputs = snv_multiplicity_to_run.map { meta, _snv_multiplicity, _snv_multiplicity_hets, _snv_multiplicity_germline -> [meta.patient, meta - meta.subMap('num_lanes', 'lane', 'read_group', 'id', 'tumor_id') + [id: meta.sample, tumor_id: meta.sample]] }.unique()
    snv_multiplicity_inputs_somatic_vcf = annotated_vcf_ffpe_impact_or_snpeff_for_merge
        .join(snv_multiplicity_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, annotated somatic snv vcf
    snv_multiplicity_inputs_jabba_gg = non_integer_balance_balanced_gg_for_merge
        .join(snv_multiplicity_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, jabba ggraph
    snv_multiplicity_inputs_hets_sites = hets_sites_for_merge
        .join(snv_multiplicity_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, het sites
    snv_multiplicity_inputs_dryclean_tumor_cov = dryclean_tumor_cov_for_merge
        .join(snv_multiplicity_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, dryclean cov
    
    if (params.tumor_only) {
        // tumor/sample id is required for snv multiplicity
        snv_multiplicity_inputs_status = snv_multiplicity_inputs.branch{
            tumor:  it[1].status.toString() == "1"
        }

        // add empty arrays to stand-in for normals
        snv_multiplicity_inputs = snv_multiplicity_inputs_status.tumor.map{ patient, meta -> [ patient, meta + [tumor_id: meta.sample] ] }
        input_snv_multiplicity = snv_multiplicity_inputs
            .join(snv_multiplicity_inputs_somatic_vcf)
            .join(snv_multiplicity_inputs_jabba_gg)
            .join(snv_multiplicity_inputs_hets_sites)
            .join(snv_multiplicity_inputs_dryclean_tumor_cov)
            .map{
                patient, meta, somatic_ann, ggraph, hets, dryclean_cov  ->
                [ meta, somatic_ann, [], ggraph, hets, dryclean_cov ]
            }
    } else {
        // getting the tumor and normal cram files separated
        snv_multiplicity_inputs_status = snv_multiplicity_inputs.branch{
            normal: it[1].status.toString() == "0"
            tumor:  it[1].status.toString() == "1"
        }

        // Crossing the normal and tumor samples to create tumor and normal pairs
        snv_multiplicity_inputs = snv_multiplicity_inputs_status.normal.cross(snv_multiplicity_inputs_status.tumor)
            .map { normal, tumor ->
                def patient = normal[0]
                def meta = [:]

                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ patient, meta ]
        }

        snv_multiplicity_inputs_germline_vcf = snv_germline_annotations_for_merge
            .join(snv_multiplicity_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, annotated germline snv vcf

        input_snv_multiplicity = snv_multiplicity_inputs
            .join(snv_multiplicity_inputs_somatic_vcf)
            .join(snv_multiplicity_inputs_germline_vcf)
            .join(snv_multiplicity_inputs_jabba_gg)
            .join(snv_multiplicity_inputs_hets_sites)
            .join(snv_multiplicity_inputs_dryclean_tumor_cov)
            .map{
                patient, meta, somatic_ann, germline_ann, ggraph, hets, dryclean_cov ->
                [ meta, somatic_ann, germline_ann, ggraph, hets, dryclean_cov ]
            }
    }
    if (tools_used.contains("all") || tools_used.contains("snv_multiplicity")) {

        VCF_SNV_MULTIPLICITY(input_snv_multiplicity)

        // snv_multiplicity = Channel.empty()
        //     .mix(VCF_SNV_MULTIPLICITY.out.snv_multiplicity_rds)
        //     .mix(snv_multiplicity_existing_outputs)

        snv_multiplicity = VCF_SNV_MULTIPLICITY.out.snv_multiplicity_rds.map { meta, multiplicity_rds -> [ meta.patient, meta, multiplicity_rds] }
            .join(
                VCF_SNV_MULTIPLICITY.out.snv_multiplicity_hets_rds.map { meta, het_rds -> [ meta.patient, het_rds] }
            )
            .join(
                VCF_SNV_MULTIPLICITY.out.snv_multiplicity_germline_rds.map { meta, germline_rds -> [ meta.patient, germline_rds] }
            )
        
        if (!params.tumor_only) {
            snv_multiplicity = snv_multiplicity
                .join(
                    VCF_SNV_MULTIPLICITY.out.snv_multiplicity_germline_rds.map { meta, germline_rds -> [ meta.patient, germline_rds] }
                )
        } else {
            snv_multiplicity = snv_multiplicity
                .map{ patient, meta, multiplicity_rds, hets_rds ->
                    [ patient, meta, multiplicity_rds, hets_rds, [] ]
                }
        }
        snv_multiplicity = snv_multiplicity
            .mix(snv_multiplicity_existing_outputs.map { it -> [ it[0].patient ] + it.toList() })
            .map { it -> it[1..-1] }
    }

    emit:
    snv_multiplicity
}

// ONCOKB
include { VCF_FUSIONS_CNA_ONCOKB_ANNOTATOR } from '../oncokb/main'
workflow ONCOKB_STEP {
    take:
    inputs_unlaned
    annotated_vcf_ffpe_impact_or_snpeff_for_merge
    non_integer_balance_balanced_gg_for_merge
    fusions_for_merge
    tools_used

    main:
    versions = Channel.empty()
    // Existing
    oncokb_existing_outputs_maf = inputs_unlaned.map { it -> [it.meta, it.oncokb_maf] }.filter { !it[1].isEmpty() }
    oncokb_existing_outputs_fusions = inputs_unlaned.map { it -> [it.meta, it.oncokb_fusions] }.filter { !it[1].isEmpty() }
    oncokb_existing_outputs_cna = inputs_unlaned.map { it -> [it.meta, it.oncokb_cna] }.filter { !it[1].isEmpty() }

    // Emit
    merged_oncokb_vcf = oncokb_existing_outputs_maf
    merged_oncokb_fusions = oncokb_existing_outputs_fusions
    merged_oncokb_cna = oncokb_existing_outputs_cna

    if ((tools_used.contains ("all") || tools_used.contains("oncokb"))) {
        oncokb_inputs = inputs_unlaned
            .filter { it.oncokb_maf.isEmpty() || it.oncokb_fusions.isEmpty() || it.oncokb_cna.isEmpty() }
            .map { it -> [it.meta.patient, it.meta + [id: it.meta.sample]] }

        oncokb_inputs_annotated_vcf = annotated_vcf_ffpe_impact_or_snpeff_for_merge
            .join(oncokb_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, annotated somatic snv

        oncokb_inputs_fusions = fusions_for_merge
            .join(oncokb_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, fusions

        oncokb_inputs_jabba_gg = non_integer_balance_balanced_gg_for_merge
            .join(oncokb_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, jabba ggraph

        

        oncokb_input = oncokb_inputs
            .join(oncokb_inputs_annotated_vcf)
            .join(oncokb_inputs_fusions)
            .join(oncokb_inputs_jabba_gg)
            .map{
                _patient,
                meta,
                snv_ann,
                fusions,
                jabba ->[ meta, snv_ann, fusions, jabba]

            }

            VCF_FUSIONS_CNA_ONCOKB_ANNOTATOR(oncokb_input)

            versions = versions.mix(VCF_FUSIONS_CNA_ONCOKB_ANNOTATOR.out.versions)
            merged_oncokb_vcf = Channel.empty()
                .mix(VCF_FUSIONS_CNA_ONCOKB_ANNOTATOR.out.merged_oncokb_vcf)
                .mix(oncokb_existing_outputs_maf)

            merged_oncokb_fusions = Channel.empty()
                .mix(VCF_FUSIONS_CNA_ONCOKB_ANNOTATOR.out.merged_oncokb_fusions)
                .mix(oncokb_existing_outputs_fusions)

            merged_oncokb_cna = Channel.empty()
                .mix(VCF_FUSIONS_CNA_ONCOKB_ANNOTATOR.out.merged_oncokb_cna)
                .mix(oncokb_existing_outputs_cna)
    }

    emit:
    merged_oncokb_vcf
    merged_oncokb_fusions
    merged_oncokb_cna
}


include { JUNC_SNV_GGRAPH_HRDETECT } from '../hrdetect/main'
workflow HRDETECT_STEP {
    take:
    inputs_unlaned
    vcf_from_sv_calling_for_merge
    hets_sites_for_merge
    annotated_vcf_ffpe_impact_or_snpeff_for_merge
    non_integer_balance_balanced_gg_for_merge
    tools_used

    main:
    versions = Channel.empty()

    // Existing
    hrdetect_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.hrdetect] }.filter { !it[1].isEmpty() }

    // Emit
    hrdetect_rds = hrdetect_existing_outputs
    
    // Inputs
    hrdetect_inputs = inputs_unlaned
        .filter { it.hrdetect.isEmpty() }
        .map { it -> [it.meta.patient, it.meta + [id: it.meta.sample]] }

    hrdetect_inputs_sv = vcf_from_sv_calling_for_merge
        .join(hrdetect_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, somatic/pon+gnomAD filtered sv

    hrdetect_inputs_hets = hets_sites_for_merge
        .join(hrdetect_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, hets

    hrdetect_inputs_vcf = annotated_vcf_ffpe_impact_or_snpeff_for_merge
        .join(hrdetect_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, somatic snv

    hrdetect_inputs_jabba_rds = non_integer_balance_balanced_gg_for_merge
        .join(hrdetect_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, non integer balanced ggraph

    

    hrdetect_input = hrdetect_inputs
        .join(hrdetect_inputs_sv)
        .join(hrdetect_inputs_hets)
        .join(hrdetect_inputs_vcf)
        .join(hrdetect_inputs_jabba_rds)
        .map{
            _patient,
            meta,
            sv,
            hets,
            snv,
            jabba ->[ meta, sv, hets, snv, jabba ]
        }

    if ((tools_used.contains("all") || tools_used.contains("hrdetect"))) {
            JUNC_SNV_GGRAPH_HRDETECT(hrdetect_input)

            versions = versions.mix(JUNC_SNV_GGRAPH_HRDETECT.out.versions)
            hrdetect_rds = Channel.empty()
                .mix(JUNC_SNV_GGRAPH_HRDETECT.out.hrdetect_rds)
                .mix(hrdetect_existing_outputs)
    }

    emit:
    hrdetect_rds

}

// OnenessTwoness
include { HRD_ONENESS_TWONESS } from '../onenesstwoness/main'
workflow ONENESS_TWONESS_STEP {
    take:
    inputs_unlaned
    events_for_merge
    hrdetect_rds_for_merge
    tools_used

    main:
    versions = Channel.empty()

    // Existing
    onenesstwoness_existing_outputs = inputs_unlaned.map { it -> [it.meta, it.onenesstwoness] }.filter { !it[1].isEmpty() }

    // Emit
    onenesstwoness_rds = onenesstwoness_existing_outputs

    // Inputs
    onenesstwoness_inputs = inputs_unlaned
        .filter { it.onenesstwoness.isEmpty() }
        .map { it -> [it.meta.patient, it.meta] }
        .unique()

    onenesstwoness_inputs_events = events_for_merge
        .join(onenesstwoness_inputs)
        .map { it -> [ it[2], it[1] ] } // meta, events

    onenesstwoness_inputs_hrd = hrdetect_rds_for_merge
        .join(onenesstwoness_inputs)
        .map { it -> [ it[2], it[1] ] } // meta, hrdetect rds


    if ((tools_used.contains("all") || tools_used.contains("onenesstwoness"))) {        

        HRD_ONENESS_TWONESS(onenesstwoness_inputs_events, onenesstwoness_inputs_hrd)

        versions = versions.mix(HRD_ONENESS_TWONESS.out.versions)
        onenesstwoness_rds = Channel.empty()
            .mix(HRD_ONENESS_TWONESS.out.oneness_twoness_results)
            .mix(onenesstwoness_existing_outputs)
    }

    emit:
    onenesstwoness_rds
}