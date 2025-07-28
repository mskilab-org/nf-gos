//
// MERGE INDEX BAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_MERGE_BAM } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_BAM       } from '../../../modules/nf-core/samtools/merge/main'

workflow BAM_MERGE_INDEX_SAMTOOLS {
    take:
    bam // channel: [mandatory] meta, bam

    main:
    versions = Channel.empty()

    // Figuring out if there is one or more bam(s) from the same sample
    bam_to_merge = bam.branch{ meta, bam ->
        // bam is a list, so use bam.size() to asses number of intervals
        single:   bam.size() <= 1
            return [ meta, bam[0] ]
        multiple: bam.size() > 1
    }

    // Only when using intervals
    MERGE_BAM(bam_to_merge.multiple, [ [ id:'null' ], []], [ [ id:'null' ], []])

    // Mix intervals and no_intervals channels together
    bam_all = MERGE_BAM.out.bam.mix(bam_to_merge.single)

    bam_index_input = bam_all.unique { it -> it[0].sample }

    // Index bam
    INDEX_MERGE_BAM(bam_index_input)

    // Join with the bai file
    // bam_bai = bam_all.join(INDEX_MERGE_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)
    bai_combine = INDEX_MERGE_BAM.out.bai.map{ meta, bai -> [ meta.sample, bai ] }
    bam_bai = bam_all
        .map{ it -> [it[0].sample] + it.toList() }
        .combine(bai_combine, by: 0)
        .map{ it -> it[1..-1] }
        .dump(tag: "bam_bai", pretty: true)

    // Gather versions of all tools used
    versions = versions.mix(INDEX_MERGE_BAM.out.versions)
    versions = versions.mix(MERGE_BAM.out.versions)

    emit:
    bam_bai

    versions
}
