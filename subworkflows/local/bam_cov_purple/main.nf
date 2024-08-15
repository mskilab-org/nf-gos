//
// PURPLE (+AMBER and COBALT)
//

include { AMBER } from '../../../modules/local/amber/main'
include { COBALT } from '../../../modules/local/cobalt/main'
include { PURPLE } from '../../../modules/local/purple/main'
include { EXTRACT_PURITYPLOIDY } from '../../../modules/local/purple/main'

//AMBER
genome_ver     = WorkflowNfcasereports.create_value_channel(params.genome_ver_amber)
het_sites       = WorkflowNfcasereports.create_file_channel(params.het_sites_amber)
//COBALT
gc_profile = WorkflowNfcasereports.create_file_channel(params.gc_profile)
diploid_bed = WorkflowNfcasereports.create_file_channel(params.diploid_bed)
//PURPLE
genome_fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
genome_fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
genome_dict = WorkflowNfcasereports.create_file_channel(params.dict)
sage_known_hotspots_somatic = WorkflowNfcasereports.create_file_channel(params.somatic_hotspots)
sage_known_hotspots_germline = WorkflowNfcasereports.create_file_channel(params.germline_hotspots)
driver_gene_panel = WorkflowNfcasereports.create_file_channel(params.driver_gene_panel)
ensembl_data_resources = WorkflowNfcasereports.create_file_channel(params.ensembl_data_resources)

workflow BAM_COV_PURPLE {
    // defining inputs
    take:
    bams // required: [meta, tbam, tbai, nbam, nbai]
    tumor_sv // required: [meta, vcf, tbi]
    tumor_snv // required: [meta, vcf, tbi]
    normal_snv // optional [meta, vcf, tbi]

    //Creating empty channels for output
    main:
    versions        = Channel.empty()
    ploidy          = Channel.empty()

    BAM_AMBER(bams)

    amber_dir        = Channel.empty()
        .mix(BAM_AMBER.out.amber_dir)
        .map { meta, amber_dir -> [meta.patient, amber_dir] }

    COV_COBALT(bams)

    cobalt_dir        = Channel.empty()
        .mix(COV_COBALT.out.cobalt_dir)
        .map { meta, cobalt_dir -> [meta.patient, cobalt_dir] }

    tumor_sv        = tumor_sv.map { meta, vcf, tbi -> [meta.patient, vcf, tbi] }
    tumor_snv        = tumor_snv.map { meta, vcf, tbi -> [meta.patient, vcf, tbi] }
    if (!params.tumor_only) {
        normal_snv        = normal_snv.map { meta, vcf, tbi -> [meta.patient, vcf, tbi] }
    }
    meta = bams.map { meta, tbam, tbai, nbam, nbai -> [meta.patient, meta] }

    if (params.tumor_only) {
        purple_inputs = meta
            .join(amber_dir)
            .join(cobalt_dir)
            .join(tumor_sv)
            .join(tumor_snv)
            .map { patient, meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi ->
                [meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi, [], []]
            }
    } else {
        purple_inputs = meta
            .join(amber_dir)
            .join(cobalt_dir)
            .join(tumor_sv)
            .join(tumor_snv)
            .join(normal_snv)
            .map { patient, meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi, germ_snv_vcf, germ_snv_tbi ->
                [meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi, germ_snv_vcf, germ_snv_tbi]
            }

    }

    PURPLE(
        purple_inputs,
        genome_fasta,
        genome_ver,
        genome_fai,
        genome_dict,
        gc_profile,
        sage_known_hotspots_somatic,
        sage_known_hotspots_germline,
        driver_gene_panel,
        ensembl_data_resources,
        [],
        [],
        [],
        []
    )


    // initializing outputs from fragcounter
    purple_dir        = Channel.empty().mix(PURPLE.out.purple_dir)
    versions          = versions.mix(PURPLE.out.versions)

    EXTRACT_PURITYPLOIDY(purple_dir)

    purity = Channel.empty().mix(EXTRACT_PURITYPLOIDY.out.purity_val) // meta, purity
    ploidy = ploidy.mix(EXTRACT_PURITYPLOIDY.out.ploidy_val) // meta, ploidy

    emit:
    ploidy

    versions
}


workflow BAM_AMBER {
    take:
    input // [meta, tbam, tbai, nbam, nbai]

    main:
    amber_dir        = Channel.empty()
    versions          = Channel.empty()

    AMBER(
        input,
        genome_ver,
        het_sites,
        []
    )

    amber_dir        = amber_dir.mix(AMBER.out.amber_dir)
    versions         = versions.mix(AMBER.out.versions)

    emit:
    amber_dir
    versions
}

workflow COV_COBALT {
    take:
    input // [meta, tbam, tbai, nbam, nbai]

    main:
    cobalt_dir        = Channel.empty()
    versions          = Channel.empty()

    if (params.tumor_only) {
        COBALT(
            input,
            gc_profile,
            diploid_bed,
            []
        )
    } else {
        COBALT(
            input,
            gc_profile,
            [],
            []
        )
    }

    cobalt_dir        = cobalt_dir.mix(COBALT.out.cobalt_dir)
    versions          = versions.mix(COBALT.out.versions)

    emit:
    cobalt_dir
    versions
}
