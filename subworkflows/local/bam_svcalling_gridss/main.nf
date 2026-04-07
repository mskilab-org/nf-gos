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
    inputs_unlaned

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

    vcf = vcf.map { it -> 
        [ it[0].patient, [ it[0], it[1], it[2] ] ] // meta.patient, [vcf, vcf_index]
    }
    .join(
        inputs_unlaned.filter {it -> it.meta.status.toString() == "0"}.map { it ->  [ it.meta.patient, [ it.meta + [ normal_id: it.meta.sample ] ] ]}.unique(),
        remainder: true
    )
    .dump(tag: "right after join vcf_joined_with_inputs for GRIDSS_SOMATIC_FILTER_STEP", pretty: true)
    .map { it -> // patient, [ meta (vcf), vcf (vcf), vcf_tbi (vcf) ], [ meta with normal id (inputs) ]
        def (_key, existing, meta_input_lst) = (it + [null, null])[0..2]
        def (meta_existing, vcf_existing, tbi_existing) = existing ?: [null, null, null]
        def meta_input = meta_input_lst ? meta_input_lst[0] : null // should only be one entry in the list since we unique by patient, but just in case we take the first
        if (meta_input != null) {
            meta_existing = meta_existing + [ normal_id: meta_input.normal_id ]
        }
        [ meta_existing, vcf_existing, tbi_existing ]
    }
    .dump(tag: "vcf_joined_with_inputs for GRIDSS_SOMATIC_FILTER_STEP", pretty: true)
    .filter { it -> it[0] != null } // filter out any entries that don't have meta (shouldn't be any since we join with remainder: true, but just in case)
    
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