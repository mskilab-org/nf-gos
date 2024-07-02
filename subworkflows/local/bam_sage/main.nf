//
// BAM SAGE
//

include { SAGE_SOMATIC } from '../../../modules/local/sage/somatic/main.nf'
include { SAGE_PASS_FILTER } from '../../../modules/local/sage/somatic/main.nf'
include { SAGE_TUMOR_ONLY_FILTER } from '../../../modules/local/sage/somatic/main.nf'
include { SAGE_GERMLINE } from '../../../modules/local/sage/germline/main.nf'

workflow BAM_SAGE {
    // defining inputs
    take:
    inputs      // [meta, normal bam, normal bai, tumor bam, tumor bai ]
    ref
    ref_fai
    ref_genome_dict
    ref_genome_version
    ensembl_data_dir
    somatic_hotspots
    panel_bed
    high_confidence_bed

    //Creating empty channels for output
    main:
    versions            = Channel.empty()
    sage_pass_filtered_somatic_vcf    = Channel.empty()
    sage_germline_vcf   = Channel.empty()

    SAGE_SOMATIC(
        inputs,
        ref,
        ref_fai,
        ref_genome_dict,
        ref_genome_version,
        ensembl_data_dir,
        somatic_hotspots,
        panel_bed,
        high_confidence_bed
    )

    // initializing outputs from fragcounter
    versions        = SAGE_SOMATIC.out.versions
    sage_somatic_vcf        = SAGE_SOMATIC.out.vcf

    SAGE_PASS_FILTER(sage_somatic_vcf)

    sage_pass_filtered_somatic_vcf = SAGE_PASS_FILTER.out.vcf

    SAGE_GERMLINE(
        inputs,
        ref,
        ref_fai,
        ref_genome_dict,
        ref_genome_version,
        ensembl_data_dir,
        somatic_hotspots,
        panel_bed,
        high_confidence_bed
    )

    sage_germline_vcf        = SAGE_GERMLINE.out.vcf

    emit:
    sage_somatic_vcf
    sage_germline_vcf
    sage_pass_filtered_somatic_vcf

    versions
}


workflow BAM_SAGE_TUMOR_ONLY_FILTER {
    // defining inputs
    take:
    inputs
    dbsnp
    dbsnp_tbi
    gnomAD_snv_db
    gnomeAD_snv_db_tbi
    sage_germ_pon
    sage_germ_pon_tbi
    mills_gold_indel
    mills_gold_indel_tbi

    //Creating empty channels for output
    main:
    versions            = Channel.empty()
    sage_filtered_vcf   = Channel.empty()

    SAGE_TUMOR_ONLY_FILTER(
        inputs,
        dbsnp,
        dbsnp_tbi,
        gnomAD_snv_db,
        gnomeAD_snv_db_tbi,
        sage_germ_pon,
        sage_germ_pon_tbi,
        mills_gold_indel,
        mills_gold_indel_tbi
    )

    // initializing outputs from fragcounter
    versions             = SAGE_TUMOR_ONLY_FILTER.out.versions
    sage_filtered_vcf    = SAGE_TUMOR_ONLY_FILTER.out.sage_filtered_vcf

    emit:
    sage_filtered_vcf
    versions
}
