//
// GRIDSS SV CALLING
//
//
//

include { GRIDSS_GRIDSS   } from '../../../modules/local/gridss/gridss/main.nf'
include { 
    GRIDSS_SOMATIC_FILTER; 
    GRIPSS_SOMATIC_FILTER  
} from '../../../modules/local/gridss/somaticFilter/main.nf'

workflow BAM_SVCALLING_GRIDSS {
    take:
    cram                                              // channel: [mandatory] [ meta, normalcram, normalcrai, tumorcram, tumorcrai ]
    bwa_index                                         // channel: [mandatory] bwa index path


    main:
    fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
    fasta_fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
    blacklist_gridss = WorkflowNfcasereports.create_file_channel(params.blacklist_gridss)
    versions               = Channel.empty()
    vcf                    = Channel.empty()
    vcf_index              = Channel.empty()
    assembly_bam           = Channel.empty()

    GRIDSS_GRIDSS(cram, fasta, fasta_fai, bwa_index, blacklist_gridss)

    vcf                    = GRIDSS_GRIDSS.out.filtered_vcf
    vcf_index              = GRIDSS_GRIDSS.out.filtered_vcf_index


    versions = versions.mix(GRIDSS_GRIDSS.out.versions)

    emit:
    vcf     // channel: [mandatory] [ meta, pass filltered vcf, pass filtered vcf tbi ]
    vcf_index


    versions

}



workflow GRIDSS_SOMATIC_FILTER_STEP {
    take:
    vcf

    main:
    // pondir_gridss = WorkflowNfcasereports.create_file_channel(params.pon_gridss)
    pon_gridss_bedpe_svs = WorkflowNfcasereports.create_file_channel(params.pon_gridss_bedpe_svs)
    pon_gridss_bed_breakends = WorkflowNfcasereports.create_file_channel(params.pon_gridss_bed_breakends)
    pon_gridss_known_hotspots_bedpe = WorkflowNfcasereports.create_file_channel(params.pon_gridss_known_hotspots_bedpe)
    pon_gridss_ref_genome_version = WorkflowNfcasereports.create_value_channel(params.pon_gridss_ref_genome_version)
    fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
    fasta_fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)

    versions                = Channel.empty()
    somatic_all             = Channel.empty()
    somatic_high_confidence = Channel.empty()

    GRIPSS_SOMATIC_FILTER(
        vcf, 
        pon_gridss_bedpe_svs,
        pon_gridss_bed_breakends,
        pon_gridss_known_hotspots_bedpe,
        fasta,
        fasta_fai,
        pon_gridss_ref_genome_version
    )

    somatic_high_confidence = GRIPSS_SOMATIC_FILTER.out.somatic_high_vcf
    somatic_all             = GRIPSS_SOMATIC_FILTER.out.somatic_all_vcf
    all_vcf = Channel.empty().mix(somatic_all, somatic_high_confidence)

    versions                = GRIPSS_SOMATIC_FILTER.out.versions

    emit:
    somatic_high_confidence // channel: [mandatory] [ meta, high confidence somatic vcf, high confidence somatic vcf tbi ]
    somatic_all
    all_vcf

    versions

}