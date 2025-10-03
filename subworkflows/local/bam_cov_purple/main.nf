//
// PURPLE (COBALT)
//

include { PURPLE } from '../../../modules/local/purple/main'
include { EXTRACT_PURITYPLOIDY } from '../../../modules/local/purple/main'
include { REVISE_PURITYPLOIDY } from '../../../modules/local/purple/main'

//PURPLE

workflow BAM_COV_PURPLE {
    // defining inputs
    take:
    purple_inputs
    inputs_unlaned

    //Creating empty channels for output
    main:
    gc_profile = WorkflowNfcasereports.create_file_channel(params.gc_profile)
    genome_fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
    genome_fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
    genome_ver = WorkflowNfcasereports.create_value_channel(params.genome_ver_amber)
    highly_diploid_percentage = WorkflowNfcasereports.create_value_channel(params.purple_highly_diploid_percentage)
    // min_purity   = WorkflowNfcasereports.create_value_channel(params.purple_min_purity)
    ploidy_penalty_factor = WorkflowNfcasereports.create_value_channel(params.purple_ploidy_penalty_factor)
    genome_dict = WorkflowNfcasereports.create_file_channel(params.dict)
    sage_known_hotspots_somatic = WorkflowNfcasereports.create_file_channel(params.somatic_hotspots)
    sage_known_hotspots_germline = WorkflowNfcasereports.create_file_channel(params.germline_hotspots)
    driver_gene_panel = WorkflowNfcasereports.create_file_channel(params.driver_gene_panel)
    ensembl_data_resources = WorkflowNfcasereports.create_file_channel(params.ensembl_data_resources)

    versions = Channel.empty()
    ploidy = Channel.empty()


    minmaxvals = inputs_unlaned.filter{ it -> it.meta.status.toString() == "1" }.map { it -> [it.meta.patient, it.subMap(["min_purity_limit", "max_purity_limit", "min_ploidy_limit", "max_ploidy_limit"]) ] }.unique{ it -> it[0] }
    purple_inputs_w_minmaxpurity = purple_inputs
        .map{it -> [ it[0].patient ] + it.toList() }
        .cross(minmaxvals)
        .map { tuple_purple, tuple_minmaxvals ->
            def minmaxes = tuple_minmaxvals[1]  // min_purity, max_purity, min_ploidy, max_ploidy
            def is_minmaxes_empty = minmaxes instanceof List && minmaxes.isEmpty()
            if (is_minmaxes_empty) {
                minmaxes = [min_purity_limit: [], max_purity_limit: [], min_ploidy_limit: 1.5, max_ploidy_limit: 5.5]
            }
            def rewrite_min_purity = params.is_heme && minmaxes.min_purity_limit instanceof List && minmaxes.min_purity_limit.isEmpty()
            if (rewrite_min_purity) {
                minmaxes = minmaxes + [min_purity_limit: params.purple_min_purity]
            }
            def purple_out_tuple = tuple_purple.toList() + minmaxes.values()
            purple_out_tuple[1..-1] // meta, amber_dir, cobalt_dir, sv_vcf, sv_tbi, snv_vcf, snv_tbi, germ_snv_vcf or [], germ_snv_tbi or [], min_purity, max_purity, min_ploidy, max_ploidy
        }
        .dump(tag: "purple_inputs_w_minmaxpurity", pretty: true)

    PURPLE(
        purple_inputs_w_minmaxpurity.map{ it -> it[0..-5] },
        genome_fasta,
        genome_ver,
        highly_diploid_percentage,
        purple_inputs_w_minmaxpurity.map{ it -> [it[0].patient, it[-4]] }.dump(tag: "min purity val input to PURPLE", pretty: true).map {it -> it[1]},
        // min_purity,
        purple_inputs_w_minmaxpurity.map{ it -> [it[0].patient, it[-3]] }.dump(tag: "max purity val input to PURPLE", pretty: true).map {it -> it[1]},
        purple_inputs_w_minmaxpurity.map{ it -> [it[0].patient, it[-2]] }.dump(tag: "min ploidy val input to PURPLE", pretty: true).map {it -> it[1]},
        purple_inputs_w_minmaxpurity.map{ it -> [it[0].patient, it[-1]] }.dump(tag: "max ploidy val input to PURPLE", pretty: true).map {it -> it[1]},
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
    purple_purity = Channel.empty().mix(PURPLE.out.purple_purity)
    purple_purity_to_extract = purple_purity
    versions          = versions.mix(PURPLE.out.versions)


    if (params.is_heme || params.purple_revise_purity_ploidy) {
        input_purple_revision = (
            purple_purity
                .map { meta, purity_bestfit -> [ meta.patient, meta, purity_bestfit ] }
                .unique { it -> it[0] }
                .cross(
                    purple_dir.map { meta, purple_paths_list ->
                        def paths = purple_paths_list.flatten()
                        def purity_range_path = paths.findAll { it =~ /.*\.purple.purity.range.tsv$/ }
                        [ meta.patient, meta, purity_range_path ]
                        .unique { it -> it[0] }
                    }
                )
        )

        purple_revise_with_range = WorkflowNfcasereports.create_value_channel(params.purple_revise_with_range)
        
        REVISE_PURITYPLOIDY(
            input_purple_revision.map{ it -> it[0][1..-1]}, 
            purple_revise_with_range,
            input_purple_revision.map{ it -> it[1][1..-1]}
        )
        purple_purity_to_extract = Channel.empty().mix(REVISE_PURITYPLOIDY.out.purple_purity_revised)
        // Use revised purity/ploidy values for heme
    }

    
    EXTRACT_PURITYPLOIDY(purple_purity_to_extract)

    purity = Channel.empty().mix(EXTRACT_PURITYPLOIDY.out.purity_val) // meta, purity
    ploidy = ploidy.mix(EXTRACT_PURITYPLOIDY.out.ploidy_val) // meta, ploidy

    emit:
    purity
    ploidy

    versions
}
