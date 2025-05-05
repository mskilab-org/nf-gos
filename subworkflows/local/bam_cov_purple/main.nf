//
// PURPLE (COBALT)
//

include { PURPLE } from '../../../modules/local/purple/main'
include { EXTRACT_PURITYPLOIDY } from '../../../modules/local/purple/main'

//PURPLE
gc_profile = WorkflowNfcasereports.create_file_channel(params.gc_profile)
genome_fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
genome_fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
genome_ver     = WorkflowNfcasereports.create_value_channel(params.genome_ver_amber)
highly_diploid_percentage = WorkflowNfcasereports.create_value_channel(params.purple_highly_diploid_percentage)
min_purity   = WorkflowNfcasereports.create_value_channel(params.purple_min_purity)
ploidy_penalty_factor = WorkflowNfcasereports.create_value_channel(params.purple_ploidy_penalty_factor)
genome_dict = WorkflowNfcasereports.create_file_channel(params.dict)
sage_known_hotspots_somatic = WorkflowNfcasereports.create_file_channel(params.somatic_hotspots)
sage_known_hotspots_germline = WorkflowNfcasereports.create_file_channel(params.germline_hotspots)
driver_gene_panel = WorkflowNfcasereports.create_file_channel(params.driver_gene_panel)
ensembl_data_resources = WorkflowNfcasereports.create_file_channel(params.ensembl_data_resources)

workflow BAM_COV_PURPLE {
    // defining inputs
    take:
    purple_inputs

    //Creating empty channels for output
    main:
    versions        = Channel.empty()
    ploidy          = Channel.empty()

    PURPLE(
        purple_inputs,
        genome_fasta,
        genome_ver,
        highly_diploid_percentage,
        min_purity,
        ploidy_penalty_factor,
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
    purple_purity     = Channel.empty().mix(PURPLE.out.purple_purity)
    versions          = versions.mix(PURPLE.out.versions)

    EXTRACT_PURITYPLOIDY(purple_purity)

    purity = Channel.empty().mix(EXTRACT_PURITYPLOIDY.out.purity_val) // meta, purity
    ploidy = ploidy.mix(EXTRACT_PURITYPLOIDY.out.ploidy_val) // meta, ploidy

    emit:
    purity
    ploidy

    versions
}
