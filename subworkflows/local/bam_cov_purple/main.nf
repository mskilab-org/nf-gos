//
// PURPLE (+AMBER and COBALT)
//

include { AMBER } from '../../../modules/local/amber/main'
include { COBALT } from '../../../modules/local/cobalt/main'
include { PURPLE } from '../../../modules/local/purple/main'

//AMBER
genome_ver     = WorkflowNfcasereports.create_value_channel(params.genome_ver_amber)
het_sites       = WorkflowNfcasereports.create_file_channel(params.het_sites_amber)
//COBALT
gc_profile = WorkflowNfcasereports.create_file_channel(params.gc_profile)
//PURPLE
genome_fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
genome_dict = WorkflowNfcasereports.create_file_channel(params.dict)
sage_known_hotspots_somatic = WorkflowNfcasereports.create_file_channel(params.somatic_hotspots)
driver_gene_panel = WorkflowNfcasereports.create_file_channel(params.driver_gene_panel)
ensembl_data_resources = WorkflowNfcasereports.create_file_channel(params.ensembl_data_resources)

workflow BAM_COV_PURPLE {
    // defining inputs
    take:
    bams // required: [meta, tbam, tbai, nbam, nbai]
    tumor_sv // required: [meta, vcf, tbi]
    tumor_snv // required: [meta, vcf, tbi]

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    purple_dir        = Channel.empty()

    BAM_AMBER(bams)

    amber_dir        = AMBER.out.amber_dir

    COV_COBALT(bams)

    cobalt_dir        = COBALT.out.cobalt_dir

    purple_inputs = amber_dir
        .join(cobalt_dir)
        .join(tumor_sv)
        .join(tumor_snv)

    PURPLE(
        purple_inputs,
        genome_fasta,
        genome_ver,
        genome_fai,
        genome_dict,
        gc_profile,
        sage_known_hotspots_somatic,
        [],
        driver_gene_panel,
        ensembl_data_resources,
        [],
        [],
        [],
        []
    )

    // initializing outputs from fragcounter
    purple_dir        = PURPLE.out.purple_dir
    versions          = PURPLE.out.versions

    emit:
    purple_dir

    versions
}


workflow BAM_AMBER {
    take:
    input // [meta, tbam, tbai, nbam, nbai]

    main:
    versions          = Channel.empty()
    amber_dir        = Channel.empty()

    AMBER(
        input,
        genome_ver,
        het_sites,
        []
    )

    amber_dir        = AMBER.out.amber_dir
    versions         = AMBER.out.versions

    emit:
    amber_dir
    versions
}

workflow COV_COBALT {
    take:
    input // [meta, tbam, tbai, nbam, nbai]

    main:
    versions          = Channel.empty()
    cobalt_dir        = Channel.empty()

    COBALT(
        input,
        gc_profile,
        [],
        []
    )

    cobalt_dir        = COBALT.out.cobalt_dir
    versions          = COBALT.out.versions

    emit:
    cobalt_dir
    versions
}
