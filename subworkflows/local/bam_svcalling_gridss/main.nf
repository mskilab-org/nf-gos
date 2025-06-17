//
// GRIDSS SV CALLING
//
//
//

include { GRIDSS_GRIDSS   } from '../../../modules/local/gridss/gridss/main.nf'
include { GRIDSS_SOMATIC  } from '../../../modules/local/gridss/somaticFilter/main.nf'

// include { GRIDSS_PREPROCESS as GRIDSS_PREPROCESS_TUMOR   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"
// include { GRIDSS_PREPROCESS as GRIDSS_PREPROCESS_NORMAL   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"
// include { GRIDSS_ASSEMBLE_SCATTER   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"
// include { GRIDSS_ASSEMBLE_GATHER   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"
// include { GRIDSS_CALL   } from "${workflow.projectDir}/modules/local/gridss/gridss/main.nf"

fasta                               = WorkflowNfcasereports.create_file_channel(params.fasta)
fasta_fai                           = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
blacklist_gridss                    = WorkflowNfcasereports.create_file_channel(params.blacklist_gridss)

workflow BAM_SVCALLING_GRIDSS {
    take:
    cram                                              // channel: [mandatory] [ meta, normalcram, normalcrai, tumorcram, tumorcrai ]
    bwa_index                                         // channel: [mandatory] bwa index path


    main:
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

pondir_gridss = WorkflowNfcasereports.create_file_channel(params.pon_gridss)

workflow BAM_SVCALLING_GRIDSS_SOMATIC {
    take:
    vcf

    main:
    versions                = Channel.empty()
    somatic_all             = Channel.empty()
    somatic_high_confidence = Channel.empty()

    GRIDSS_SOMATIC(vcf, pondir_gridss)

    somatic_high_confidence = GRIDSS_SOMATIC.out.somatic_high_vcf
    somatic_all             = GRIDSS_SOMATIC.out.somatic_all_vcf
    all_vcf = Channel.empty().mix(somatic_all, somatic_high_confidence)

    versions                = GRIDSS_SOMATIC.out.versions

    emit:
    somatic_high_confidence // channel: [mandatory] [ meta, high confidence somatic vcf, high confidence somatic vcf tbi ]
    somatic_all
    all_vcf

    versions

}





workflow BAM_SVCALLING_GRIDSS_PARALLEL {
    take:
    inputs_unlaned
    bam_sv_calling // channel: [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
    bwa_index     // channel: [mandatory] bwa index path

    main:
    versions               = Channel.empty()
    vcf                    = Channel.empty()
    vcf_index              = Channel.empty()
    assembly_bam           = Channel.empty()

    bam_sv_calling_status = bam_sv_calling // meta, bam, bai
        .branch{
            tumor:  it[0].status == 1
            normal: it[0].status == 0
        }
    
    

    gridss_preprocess_tumor = GRIDSS_PREPROCESS_TUMOR(bam_sv_calling.tumor, fasta, fasta_fai, bwa_index, blacklist_gridss)
    gridss_preprocess_normal = GRIDSS_PREPROCESS_NORMAL(bam_sv_calling.normal, fasta, fasta_fai, bwa_index, blacklist_gridss)

    gridss_preprocess_tumor_for_merge = gridss_preprocess_tumor.gridss_preprocess_dir.map { meta, gridss_preprocess_dir ->
        def meta_tumor = meta + [id: "${meta.sample}"]
        [ meta_tumor.patient, meta_tumor, gridss_preprocess_dir ]
    }

    gridss_preprocess_normal_for_merge = gridss_preprocess_normal.gridss_preprocess_dir.map { meta, gridss_preprocess_dir ->
        def meta_normal = meta + [id: "${meta.sample}"]
        [ meta_normal.patient, meta_normal, gridss_preprocess_dir ]
    }

    

    BAM_SVCALLING_GRIDSS(bam_sv_calling_pair, index_alignment)

    vcf = vcf.mix(BAM_SVCALLING_GRIDSS.out.vcf)
    vcf_index = vcf_index.mix(BAM_SVCALLING_GRIDSS.out.vcf_index)

    versions = versions.mix(BAM_SVCALLING_GRIDSS.out.versions)

    emit:
    vcf
    vcf_index

    versions
}





