//
// STRELKA2 single sample variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { STRELKA_GERMLINE } from '../../../modules/nf-core/strelka/germline/main'

workflow BAM_GERMLINE_STRELKA {
    take:
    bam          // channel: [mandatory] [ meta, bam, bai ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]

    main:
    versions = Channel.empty()

    STRELKA_GERMLINE(bam, fasta, fasta_fai)

    vcf = STRELKA_GERMLINE.out.vcf

    versions = versions.mix(STRELKA_GERMLINE.out.versions)

    emit:
    vcf

    versions
}
