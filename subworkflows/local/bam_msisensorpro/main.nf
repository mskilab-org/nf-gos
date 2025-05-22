//
// BAM MSISENSORPRO
//

include { MSISENSORPRO_MSISOMATIC } from '../../../modules/nf-core/msisensorpro/msisomatic/main'

fasta     = WorkflowNfcasereports.create_file_channel(params.fasta)
tools_to_run = WorkflowNfcasereports.toolsToRun

workflow BAM_MSISENSORPRO {
    // defining inputs
    take:
    inputs
    alignment_bams_final // channel: [mandatory] [ meta.sample, meta, bam, bai ]
    msisensor_scan

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    msi_report        = Channel.empty()
    msi_dis           = Channel.empty()
    msi_somatic       = Channel.empty()
    msi_germline      = Channel.empty()

    if (tools_to_run.contains("msisensorpro") && !params.tumor_only) {

        bam_msi_inputs = inputs.filter { it.msi.isEmpty() }.map { it -> [it.meta.sample] }
        bam_msi = alignment_bams_final
            .join(bam_msi_inputs)
            .map { it -> [ it[1], it[2], it[3] ] } // meta, bam, bai
        msisensorpro_existing_outputs = inputs
            .map { it -> [it.meta, it.msi] }
            .filter { !it[1].isEmpty() }
        msisensorpro_existing_outputs_germline = inputs
            .map { it -> [it.meta, it.msi_germline] }
            .filter { !it[1].isEmpty() }

        // getting the tumor and normal cram files separated
        bam_msi_status = bam_msi.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
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

        MSISENSORPRO_MSISOMATIC(bam_msi_pair, fasta, msisensor_scan)

        // initializing outputs from fragcounter
        msi_report      = MSISENSORPRO_MSISOMATIC.out.output_report
        msi_dis         = MSISENSORPRO_MSISOMATIC.out.output_dis
        msi_somatic     = MSISENSORPRO_MSISOMATIC.out.output_somatic
        msi_germline    = MSISENSORPRO_MSISOMATIC.out.output_germline
        versions        = MSISENSORPRO_MSISOMATIC.out.versions
    }


    //
    emit:
    msi_report
    msi_dis
    msi_somatic
    msi_germline

    versions
}
