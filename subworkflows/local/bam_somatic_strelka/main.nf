//
// STRELKA2 tumor-normal variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_MERGEVCFS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { STRELKA_SOMATIC } from '../../../modules/nf-core/strelka/somatic/main'

workflow BAM_SOMATIC_STRELKA {
    take:
    bam          // channel: [mandatory] [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]

    main:
    versions = Channel.empty()

    STRELKA_SOMATIC(bam, fasta, fasta_fai)

    vcf_indels = STRELKA_SOMATIC.out.vcf_indels
    vcf_snvs = STRELKA_SOMATIC.out.vcf_snvs

    // create a single channel consist of [meta, indels, snvs] from vcf_indels (which is [meta, vcf]) and vcf_snvs (which is also [meta, vcf])
    vcfs_to_merge = vcf_indels.join(vcf_snvs)

    GATK4_MERGEVCFS(vcfs_to_merge)

    vcf = GATK4_MERGEVCFS.out.vcf
        .map{ meta, vcf -> [ meta + [ variantcaller:'strelka' ], vcf ] }

    versions = versions.mix(GATK4_MERGEVCFS.out.versions)
    versions = versions.mix(STRELKA_SOMATIC.out.versions)

    emit:
    vcf

    versions
}
