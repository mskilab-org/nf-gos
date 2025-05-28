//
// SV CALLING
//

//GRIDSS
include { BAM_SVCALLING_GRIDSS } from '../bam_svcalling_gridss/main'
include { BAM_SVCALLING_GRIDSS_SOMATIC } from '../bam_svcalling_gridss/main'

// SV Junction Filtering
include { SV_JUNCTION_FILTER as JUNCTION_FILTER } from '../junction_filter/main'

tools_to_run = WorkflowNfcasereports.toolsToRun

workflow SV_CALLING {
    take:
    inputs  // channel: [required] inputs channel
    bam       // channel: [mandatory] [ meta, normalbam, normalbai, tumorbam, tumorbai ]
    index_alignment // channel: [mandatory] index alignment file


    main:
    versions                = Channel.empty()
    gridss_vcf         = Channel.empty()
    final_junction_filtered_sv_rds   = Channel.empty()
    somatic_high_confidence_sv = Channel.empty()
    somatic_unfiltered_sv = Channel.empty()

    if (tools_to_run.contains("all") || tools_to_run.contains("gridss") || params.is_run_junction_filter) {

        // Filter out bams for which SV calling has already been done

		bam_sv_inputs = inputs.filter { it.vcf.isEmpty() }.map { it -> [it.meta.sample] }
        bam_sv_calling = bam
            .join(bam_sv_inputs)
            .map { it -> [ it[1], it[2], it[3] ] } // meta, bam, bai

		gridss_existing_outputs = inputs.map {
			it -> [it.meta, it.vcf, it.vcf_tbi] }
			.filter { !it[1].isEmpty() && !it[2].isEmpty() }

        if (params.tumor_only) {
            bam_sv_calling_status = bam_sv_calling.branch{
                tumor:  it[0].status == 1
            }

            // add empty arrays to stand-in for normals
            bam_sv_calling_pair = bam_sv_calling_status.tumor.map{ meta, bam, bai -> [ meta + [tumor_id: meta.sample], [], [], bam, bai ] }
        } else {
            // getting the tumor and normal cram files separated
            bam_sv_calling_status = bam_sv_calling.branch{
                normal: it[0].status == 0
                tumor:  it[0].status == 1
            }

            // All normal samples
            bam_sv_calling_normal_for_crossing = bam_sv_calling_status.normal.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

            // All tumor samples
            bam_sv_calling_tumor_for_crossing = bam_sv_calling_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

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
            index_alignment
        )

        gridss_vcf = Channel.empty()
            .mix(BAM_SVCALLING_GRIDSS.out.vcf)
            .mix(gridss_existing_outputs)
        versions = versions.mix(BAM_SVCALLING_GRIDSS.out.versions)

        if (params.tumor_only) {
            JUNCTION_FILTER(gridss_vcf)

            pon_filtered_sv_rds = Channel.empty().mix(JUNCTION_FILTER.out.pon_filtered_sv_rds)
            final_junction_filtered_sv_rds = Channel.empty().mix(JUNCTION_FILTER.out.final_filtered_sv_rds)
        } else {
            //somatic filter for GRIDSS
            BAM_SVCALLING_GRIDSS_SOMATIC(gridss_vcf)

            versions = versions.mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.versions)
            somatic_high_confidence_sv = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_high_confidence)

            somatic_unfiltered_sv = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_all)
        }
    }

    emit:
    gridss_vcf                    // channel: [mandatory] [ meta, pass filtered vcf, pass filtered vcf tbi ]
    final_junction_filtered_sv_rds // channel: [optional] [ meta.patient, rds ]
    somatic_high_confidence_sv // channel: [optional] [ meta, high confidence somatic vcf, high confidence somatic vcf tbi ]
    somatic_unfiltered_sv // channel: [optional] [ meta, unfiltered somatic vcf, unfiltered somatic vcf tbi ]

}

