//
// GPU ACCELERATED MAPPING
//

include { PARABRICKS_FQ2BAM            } from '../../../modules/local/fq2bam/main'

fasta                               = WorkflowNfcasereports.create_file_channel(params.fasta)
fasta_fai                           = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
intervals                          = WorkflowNfcasereports.create_file_channel(params.intervals)
mark_duplicates                     = params.fq2bam_mark_duplicates
low_memory                          = params.fq2bam_low_memory
optical_duplicate_pixel_distance    = params.optical_duplicate_pixel_distance

workflow FASTQ_PARABRICKS_FQ2BAM {
    take:
    reads // channel: [mandatory] meta, reads
    known_sites
    known_sites_tbi

    main:

    bam = Channel.empty()
    bai = Channel.empty()
    qc = Channel.empty()
    bqsr_table = Channel.empty()
    duplicate_metrics = Channel.empty()
    versions = Channel.empty()

      // if normals are present for multiple patients, these need to be deduplicated
    deduped_reads = reads
        .map { it ->
            def meta = it[0]
            [ meta - meta.subMap("patient") ] + it[1..-1].toList()
        }
        .unique { it -> it[0].sample }
        .view{ "deduped_reads: ${it}" }
    
    patient_sample_map = reads
        .map { it ->
            def meta = it[0]
            [ meta.sample, meta ]
        }
        .unique()
        .view { "patient_sample_map: ${it}" }

    PARABRICKS_FQ2BAM(
        deduped_reads,
        fasta,
        fasta_fai,
        intervals,
        known_sites,
        known_sites_tbi,
        mark_duplicates,
        optical_duplicate_pixel_distance,
        low_memory
    )

    bam = bam.mix(PARABRICKS_FQ2BAM.out.bam)
    bai = bai.mix(PARABRICKS_FQ2BAM.out.bai)
    qc = qc.mix(PARABRICKS_FQ2BAM.out.qc_metrics)
    bqsr_table = bqsr_table.mix(PARABRICKS_FQ2BAM.out.bqsr_table)
    duplicate_metrics = duplicate_metrics.mix(PARABRICKS_FQ2BAM.out.duplicate_metrics)

    bam_to_cross = bam.map { meta, bam ->
        [meta.sample, meta, bam] // meta here has no patient
    }

    bam = bam_to_cross
        .cross(patient_sample_map)
        .map { baminfo, patient_map -> 
            def meta_complete = patient_map[1]
            def bam = baminfo[2]
            [ meta_complete, bam ]
        }
    


    versions = versions.mix(PARABRICKS_FQ2BAM.out.versions)

    emit:
    bam      // channel: [ [meta], bam ]
    bai      // channel: [ [meta], bai ]
    qc
    bqsr_table
    duplicate_metrics
    versions
}
