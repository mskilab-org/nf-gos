#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mskilab-org/nf-jabba
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://tanubrata/mskilab-org/nf-jabba
    Website: https://nf-co.re/nfcasereports
    Slack  : https://nfcore.slack.com/channels/nfcasereports
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.ascat_alleles         = WorkflowMain.getGenomeAttribute(params, 'ascat_alleles')
params.ascat_genome          = WorkflowMain.getGenomeAttribute(params, 'ascat_genome')
params.ascat_loci            = WorkflowMain.getGenomeAttribute(params, 'ascat_loci')
params.ascat_loci_gc         = WorkflowMain.getGenomeAttribute(params, 'ascat_loci_gc')
params.ascat_loci_rt         = WorkflowMain.getGenomeAttribute(params, 'ascat_loci_rt')
params.bwa                   = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bwamem2               = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.cf_chrom_len          = WorkflowMain.getGenomeAttribute(params, 'cf_chrom_len')
params.chr_dir               = WorkflowMain.getGenomeAttribute(params, 'chr_dir')
params.dbsnp                 = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi             = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.dbsnp_vqsr            = WorkflowMain.getGenomeAttribute(params, 'dbsnp_vqsr')
params.dict                  = WorkflowMain.getGenomeAttribute(params, 'dict')
//params.dragmap               = WorkflowMain.getGenomeAttribute(params, 'dragmap')
params.fasta                 = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai             = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.msisensorpro_list     = WorkflowMain.getGenomeAttribute(params, 'msisensorpro_list')
//params.germline_resource     = WorkflowMain.getGenomeAttribute(params, 'germline_resource')
//params.germline_resource_tbi = WorkflowMain.getGenomeAttribute(params, 'germline_resource_tbi')
params.intervals             = WorkflowMain.getGenomeAttribute(params, 'intervals')
params.known_snps            = WorkflowMain.getGenomeAttribute(params, 'known_snps')
params.known_snps_tbi        = WorkflowMain.getGenomeAttribute(params, 'known_snps_tbi')
params.known_snps_vqsr       = WorkflowMain.getGenomeAttribute(params, 'known_snps_vqsr')
params.known_indels          = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi      = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')
params.known_indels_vqsr     = WorkflowMain.getGenomeAttribute(params, 'known_indels_vqsr')
//params.mappability           = WorkflowMain.getGenomeAttribute(params, 'mappability')
//params.pon                   = WorkflowMain.getGenomeAttribute(params, 'pon')
//params.pon_tbi               = WorkflowMain.getGenomeAttribute(params, 'pon_tbi')
params.snpeff_db             = WorkflowMain.getGenomeAttribute(params, 'snpeff_db')
params.snpeff_genome         = WorkflowMain.getGenomeAttribute(params, 'snpeff_genome')
params.snpeff_cache          = WorkflowMain.getGenomeAttribute(params, 'snpeff_cache')
params.vep_cache_version     = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
params.vep_genome            = WorkflowMain.getGenomeAttribute(params, 'vep_genome')
params.vep_species           = WorkflowMain.getGenomeAttribute(params, 'vep_species')
params.indel_mask            = WorkflowMain.getGenomeAttribute(params, 'indel_mask')
params.germ_sv_db            = WorkflowMain.getGenomeAttribute(params, 'germ_sv_db')
params.simple_seq_db         = WorkflowMain.getGenomeAttribute(params, 'simple_seq_db')
params.blacklist_gridss      = WorkflowMain.getGenomeAttribute(params, 'blacklist_gridss')
params.pon_gridss            = WorkflowMain.getGenomeAttribute(params, 'pon_gridss')
params.gnomAD_sv_db          = WorkflowMain.getGenomeAttribute(params, 'gnomAD_sv_db')
params.junction_pon_gridss   = WorkflowMain.getGenomeAttribute(params, 'junction_pon_gridss')
params.junction_pon_svaba    = WorkflowMain.getGenomeAttribute(params, 'junction_pon_svaba')
params.gcmapdir_frag         = WorkflowMain.getGenomeAttribute(params, 'gcmapdir_frag')
params.build_dryclean        = WorkflowMain.getGenomeAttribute(params, 'build_dryclean')
params.hapmap_sites          = WorkflowMain.getGenomeAttribute(params, 'hapmap_sites')
params.genome_ver_amber          = WorkflowMain.getGenomeAttribute(params, 'genome_ver_amber')
params.het_sites_amber          = WorkflowMain.getGenomeAttribute(params, 'het_sites_amber')
params.gc_profile          = WorkflowMain.getGenomeAttribute(params, 'gc_profile')
params.diploid_bed          = WorkflowMain.getGenomeAttribute(params, 'diploid_bed')
params.pon_dryclean          = WorkflowMain.getGenomeAttribute(params, 'pon_dryclean')
params.ref_genome_version          = WorkflowMain.getGenomeAttribute(params, 'ref_genome_version')
params.ensembl_data_dir          = WorkflowMain.getGenomeAttribute(params, 'ensembl_data_dir')
params.somatic_hotspots          = WorkflowMain.getGenomeAttribute(params, 'somatic_hotspots')
params.germline_hotspots          = WorkflowMain.getGenomeAttribute(params, 'germline_hotspots')
params.panel_bed          = WorkflowMain.getGenomeAttribute(params, 'panel_bed')
params.high_confidence_bed          = WorkflowMain.getGenomeAttribute(params, 'high_confidence_bed')
params.sage_pon          = WorkflowMain.getGenomeAttribute(params, 'sage_pon')
params.sage_blocklist_regions          = WorkflowMain.getGenomeAttribute(params, 'sage_blocklist_regions')
params.sage_blocklist_sites          = WorkflowMain.getGenomeAttribute(params, 'sage_blocklist_sites')
params.sage_clinvar_annotations          = WorkflowMain.getGenomeAttribute(params, 'sage_clinvar_annotations')
params.gnomAD_snv_db                    = WorkflowMain.getGenomeAttribute(params, 'gnomAD_snv_db')
params.gnomAD_snv_db_tbi                = WorkflowMain.getGenomeAttribute(params, 'gnomAD_snv_db_tbi')
params.sage_germline_pon                = WorkflowMain.getGenomeAttribute(params, 'sage_germline_pon')
params.sage_germline_pon_tbi            = WorkflowMain.getGenomeAttribute(params, 'sage_germline_pon_tbi')
params.sigprofilerassignment_genome          = WorkflowMain.getGenomeAttribute(params, 'sigprofilerassignment_genome')
params.sigprofilerassignment_cosmic_version          = WorkflowMain.getGenomeAttribute(params, 'sigprofilerassignment_cosmic_version')
params.segment_mappability          = WorkflowMain.getGenomeAttribute(params, 'segment_mappability')
params.driver_gene_panel          = WorkflowMain.getGenomeAttribute(params, 'driver_gene_panel')
params.ensembl_data_resources          = WorkflowMain.getGenomeAttribute(params, 'ensembl_data_resources')
params.gnomad_resource          = WorkflowMain.getGenomeAttribute(params, 'gnomad_resource')
params.blacklist_coverage_jabba     = WorkflowMain.getGenomeAttribute(params, 'blacklist_coverage_jabba')
params.whitelist_genes_jabba     = WorkflowMain.getGenomeAttribute(params, 'whitelist_genes_jabba')
params.gencode_fusions     = WorkflowMain.getGenomeAttribute(params, 'gencode_fusions')
params.gencode_oncokb     = WorkflowMain.getGenomeAttribute(params, 'gencode_oncokb')
params.oncokb_genes     = WorkflowMain.getGenomeAttribute(params, 'oncokb_genes')
params.build_non_integer_balance     = WorkflowMain.getGenomeAttribute(params, 'build_non_integer_balance')
params.mask_non_integer_balance     = WorkflowMain.getGenomeAttribute(params, 'mask_non_integer_balance')
params.mask_lp_phased_balance     = WorkflowMain.getGenomeAttribute(params, 'mask_lp_phased_balance')
params.ref_hrdetect     = WorkflowMain.getGenomeAttribute(params, 'ref_hrdetect')
//params.blacklist_junctions_jabba     = WorkflowMain.getGenomeAttribute(params, 'blacklist_junctions_jabba')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NFCASEREPORTS } from './workflows/nfcasereports'

//
// WORKFLOW: Run main mskilab-org/nf-jabba analysis pipeline
//
workflow MSKILABORG_NFCASEREPORTS {
    NFCASEREPORTS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    MSKILABORG_NFCASEREPORTS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
