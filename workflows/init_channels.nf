/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
ascat_alleles      = params.ascat_alleles      ? Channel.fromPath(params.ascat_alleles).collect()     : Channel.empty()
ascat_loci         = params.ascat_loci         ? Channel.fromPath(params.ascat_loci).collect()        : Channel.empty()
ascat_loci_gc      = params.ascat_loci_gc      ? Channel.fromPath(params.ascat_loci_gc).collect()     : Channel.value([])
ascat_loci_rt      = params.ascat_loci_rt      ? Channel.fromPath(params.ascat_loci_rt).collect()     : Channel.value([])
cf_chrom_len       = params.cf_chrom_len       ? Channel.fromPath(params.cf_chrom_len).collect()      : []
chr_dir            = params.chr_dir            ? Channel.fromPath(params.chr_dir).collect()           : Channel.value([])
dbsnp              = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()             : Channel.value([])
fasta              = params.fasta              ? Channel.fromPath(params.fasta).first()               : Channel.empty()
fasta_fai          = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()         : Channel.empty()
germline_resource  = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect() : Channel.value([]) // Mutect2 does not require a germline resource, so set to optional input
known_indels       = params.known_indels       ? Channel.fromPath(params.known_indels).collect()      : Channel.value([])
known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()        : Channel.value([])
mappability        = params.mappability        ? Channel.fromPath(params.mappability).collect()       : Channel.value([])
pon                = params.pon                ? Channel.fromPath(params.pon).collect()               : Channel.value([]) // PON is optional for Mutect2 (but highly recommended)
snpeff_cache       = params.snpeff_cache       ? Channel.fromPath(params.snpeff_cache).collect()      : Channel.empty()

// SVABA
indel_mask         = params.indel_mask         ? Channel.fromPath(params.indel_mask).collect()        : Channel.empty()   // This is the indel mask for SVABA
germ_sv_db         = params.germ_sv_db         ? Channel.fromPath(params.germ_sv_db).collect()        : Channel.empty()   // This is the germline SV mask for Svaba
simple_seq_db      = params.simple_seq_db      ? Channel.fromPath(params.simple_seq_db).collect()     : Channel.empty()   // This is the file containing sites of simple DNA that can confuse the contig re-alignment for SVABA

// GRIDSS
blacklist_gridss   = params.blacklist_gridss   ? Channel.fromPath(params.blacklist_gridss).collect()  : Channel.empty()   // This is the mask for gridss SV calls
pon_gridss         = params.pon_gridss         ? Channel.fromPath(params.pon_gridss).collect()        : Channel.empty()   //This is the pon directory for GRIDSS SOMATIC. (MUST CONTAIN .bed and .bedpe files)

//SV Junction Filtering
junction_pon_svaba  = params.junction_pon_svaba   ? Channel.fromPath(params.junction_pon_svaba).collect()   : Channel.empty()
junction_pon_gridss = params.junction_pon_gridss  ? Channel.fromPath(params.junction_pon_gridss).collect()  : Channel.empty()
gnomAD_sv_db        = params.gnomAD_sv_db         ? Channel.fromPath(params.gnomAD_sv_db).collect()         : Channel.empty()

// FragCounter
gcmapdir_frag      = params.gcmapdir_frag      ? Channel.fromPath(params.gcmapdir_frag).collect()     : Channel.empty()   // This is the GC/Mappability directory for fragCounter. (Must contain gc* & map* .rds files)

// HetPileups
hapmap_sites       = params.hapmap_sites       ? Channel.fromPath(params.hapmap_sites).collect()      : Channel.empty()

// Dryclean
pon_dryclean      = params.pon_dryclean      ? Channel.fromPath(params.pon_dryclean).collect()     : Channel.empty()   // This is the path to the PON for Dryclean.
blacklist_path_dryclean      = params.blacklist_path_dryclean      ? Channel.fromPath(params.blacklist_path_dryclean).collect()     : Channel.empty()   // This is the path to the blacklist for Dryclean (optional).
germline_file_dryclean      = params.germline_file_dryclean      ? Channel.fromPath(params.germline_file_dryclean).collect()     : Channel.empty()   // This is the path to the germline mask for dryclean (optional).

// JaBbA
blacklist_coverage_jabba		= params.blacklist_coverage_jabba		  ? Channel.fromPath(params.blacklist_coverage_jabba).collect() : Channel.empty()

// Fusions
gencode_fusions             = params.gencode_fusions        ? Channel.fromPath(params.gencode_fusions).collect() : Channel.empty()

// Allelic CN
mask_non_integer_balance    = params.mask_non_integer_balance   ? Channel.fromPath(params.mask_non_integer_balance).collect() : Channel.empty()
mask_lp_phased_balance    = params.mask_lp_phased_balance   ? Channel.fromPath(params.mask_lp_phased_balance).collect() : Channel.empty()

// Sage
somatic_hotspots_sage        = params.somatic_hotspots       ? Channel.fromPath(params.somatic_hotspots).collect()       : Channel.empty()
panel_bed_sage               = params.panel_bed              ? Channel.fromPath(params.panel_bed).collect()              : Channel.empty()
high_confidence_bed_sage     = params.high_confidence_bed    ? Channel.fromPath(params.high_confidence_bed).collect()    : Channel.empty()
ensembl_data_dir_sage         = params.ensembl_data_dir       ? Channel.fromPath(params.ensembl_data_dir).collect()    : Channel.empty()
gnomAD_snv_db                 = params.gnomAD_snv_db          ? Channel.fromPath(params.gnomAD_snv_db).collect()    : Channel.empty()
gnomAD_snv_db_tbi             = params.gnomAD_snv_db_tbi      ? Channel.fromPath(params.gnomAD_snv_db_tbi).collect()    : Channel.empty()
sage_germline_pon             = params.sage_germline_pon      ? Channel.fromPath(params.sage_germline_pon).collect()    : Channel.empty()
sage_germline_pon_tbi         = params.sage_germline_pon_tbi  ? Channel.fromPath(params.sage_germline_pon_tbi).collect()    : Channel.empty()

// Pave
sage_pon_pave                   = params.sage_pon ? Channel.fromPath(params.sage_pon).collect() : Channel.empty()
sage_blocklist_regions_pave     = params.sage_blocklist_regions ? Channel.fromPath(params.sage_blocklist_regions).collect() : Channel.empty()
sage_blocklist_sites_pave       = params.sage_blocklist_sites ? Channel.fromPath(params.sage_blocklist_sites).collect() : Channel.empty()
clinvar_annotations_pave        = params.clinvar_annotations ? Channel.fromPath(params.clinvar_annotations).collect() : Channel.empty()
segment_mappability_pave        = params.segment_mappability ? Channel.fromPath(params.segment_mappability).collect() : Channel.empty()
driver_gene_panel_pave          = params.driver_gene_panel ? Channel.fromPath(params.driver_gene_panel).collect() : Channel.empty()
ensembl_data_resources_pave     = params.ensembl_data_resources ? Channel.fromPath(params.ensembl_data_resources).collect() : Channel.empty()
gnomad_resource_pave            = params.gnomad_resource ? Channel.fromPath(params.gnomad_resource).collect() : Channel.empty()


// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
ascat_genome       = params.ascat_genome       ?: Channel.empty()
dbsnp_vqsr         = params.dbsnp_vqsr         ? Channel.value(params.dbsnp_vqsr) : Channel.empty()
known_indels_vqsr  = params.known_indels_vqsr  ? Channel.value(params.known_indels_vqsr) : Channel.empty()
known_snps_vqsr    = params.known_snps_vqsr    ? Channel.value(params.known_snps_vqsr) : Channel.empty()
snpeff_genome      = params.snpeff_genome      ? Channel.value(params.snpeff_genome) : Channel.empty()
snpeff_db          = params.snpeff_db          ? Channel.value(params.snpeff_db) : Channel.empty()
snpeff_db_full     = params.snpeff_db && params.snpeff_genome   ? Channel.value("${params.snpeff_genome}.${params.snpeff_db}") : Channel.empty()
vep_cache_version  = params.vep_cache_version  ?: Channel.empty()
vep_genome         = params.vep_genome         ?: Channel.empty()
vep_species        = params.vep_species        ?: Channel.empty()
error_rate         = params.error_rate         ?: Channel.empty()                         // For SVABA

//SV Junction Filter
junction_padding    = params.pad_junc_filter    ?: Channel.empty()                        //For SV Filtering (tumor only)

// Hetpileups
filter_hets         = params.filter_hets       ?: Channel.empty()
max_depth           = params.max_depth         ?: Channel.empty()

// FragCounter
windowsize_frag    = params.windowsize_frag    ?: Channel.empty()                                                         // For fragCounter
minmapq_frag       = params.minmapq_frag       ?: Channel.empty()                                                         // For fragCounter
midpoint_frag      = params.midpoint_frag      ?: Channel.empty()                                                         // For fragCounter
paired_frag        = params.paired_frag        ?: Channel.empty()                                                         // For fragCounter
exome_frag         = params.exome_frag         ?: Channel.empty()                                                         // For fragCounter

// Dryclean
center_dryclean           = params.center_dryclean          ?: Channel.empty()
cbs_dryclean                = params.cbs_dryclean               ?: Channel.empty()
cnsignif_dryclean           = params.cnsignif_dryclean          ?: Channel.empty()
wholeGenome_dryclean        = params.wholeGenome_dryclean       ?: Channel.empty()
blacklist_dryclean          = params.blacklist_dryclean         ?: Channel.empty()
germline_filter_dryclean    = params.germline_filter_dryclean   ?: Channel.empty()
human_dryclean              = params.human_dryclean             ?: Channel.empty()
field_dryclean              = params.field_dryclean             ?: Channel.empty()
build_dryclean              = params.build_dryclean             ?: Channel.empty()

// ASCAT_seg
field_ascat                 = params.field_ascat                ?: Channel.empty()
hets_thresh_ascat           = params.hets_thresh_ascat          ?: Channel.empty()
penalty_ascat               = params.penalty_ascat              ?: Channel.empty()
gc_correct_ascat            = params.gc_correct_ascat           ?: Channel.empty()
rebin_width_ascat           = params.rebin_width_ascat          ?: Channel.empty()
from_maf_ascat              = params.from_maf_ascat             ?: Channel.empty()

// CBS
cnsignif_cbs                    = params.cnsignif_cbs               ?: Channel.empty()
field_cbs                       = params.field_cbs                  ?: Channel.empty()
name_cbs                        = params.name_cbs                   ?: Channel.empty()

// JaBbA
blacklist_junctions_jabba       = params.blacklist_junctions_jabba      ?: Channel.empty()
geno_jabba					    = params.geno_jabba			            ?: Channel.empty()
indel_jabba					    = params.indel_jabba			        ?: Channel.empty()
tfield_jabba					= params.tfield_jabba			        ?: Channel.empty()
iter_jabba					    = params.iter_jabba			            ?: Channel.empty()
rescue_window_jabba				= params.rescue_window_jabba			?: Channel.empty()
rescue_all_jabba				= params.rescue_all_jabba			    ?: Channel.empty()
nudgebalanced_jabba				= params.nudgebalanced_jabba			?: Channel.empty()
edgenudge_jabba					= params.edgenudge_jabba			    ?: Channel.empty()
strict_jabba					= params.strict_jabba			        ?: Channel.empty()
allin_jabba					    = params.allin_jabba			        ?: Channel.empty()
field_jabba					    = params.field_jabba			        ?: Channel.empty()
maxna_jabba					    = params.maxna_jabba			        ?: Channel.empty()
purity_jabba					= params.purity_jabba                   ?: Channel.empty()
ploidy_jab     					= params.ploidy_jabba                   ?: Channel.empty()
pp_method_jabba					= params.pp_method_jabba                ?: Channel.empty()
cnsignif_jabba					= params.cnsignif_jabba                 ?: Channel.empty()
slack_jabba					    = params.slack_jabba                    ?: Channel.empty()
linear_jabba					= params.linear_jabba                   ?: Channel.empty()
tilim_jabba					    = params.tilim_jabba                    ?: Channel.empty()
epgap_jabba					    = params.epgap_jabba                    ?: Channel.empty()
fix_thres_jabba					= params.fix_thres_jabba			    ?: Channel.empty()
lp_jabba					    = params.lp_jabba			            ?: Channel.empty()
ism_jabba					    = params.ism_jabba			            ?: Channel.empty()
filter_loose_jabba				= params.filter_loose_jabba			    ?: Channel.empty()
gurobi_jabba					= params.gurobi_jabba			        ?: Channel.empty()
nonintegral_jabba				= params.nonintegral_jabba			    ?: Channel.empty()
verbose_jabba					= params.verbose_jabba			        ?: Channel.empty()
help_jabba					    = params.help_jabba			            ?: Channel.empty()

//Allelic CN (Non-integer balance)
field_non_integer_balance = params.field_non_integer_balance ?: Channel.empty()
hets_thresh_non_integer_balance = params.hets_thresh_non_integer_balance ?: Channel.empty()
overwrite_non_integer_balance = params.overwrite_non_integer_balance ?: Channel.empty()
lambda_non_integer_balance = params.lambda_non_integer_balance ?: Channel.empty()
allin_non_integer_balance = params.allin_non_integer_balance ?: Channel.empty()
fix_thresh_non_integer_balance = params.fix_thresh_non_integer_balance ?: Channel.empty()
nodebounds_non_integer_balance = params.nodebounds_non_integer_balance ?: Channel.empty()
ism_non_integer_balance = params.ism_non_integer_balance ?: Channel.empty()
build_non_integer_balance = params.build_non_integer_balance ?: Channel.empty()
epgap_non_integer_balance = params.epgap_non_integer_balance ?: Channel.empty()
tilim_non_integer_balance = params.tilim_non_integer_balance ?: Channel.empty()
gurobi_non_integer_balance = params.gurobi_non_integer_balance ?: Channel.empty()
pad_non_integer_balance = params.pad_non_integer_balance ?: Channel.empty()

// ...(LP Phased balance)
lambda_lp_phased_balance = params.lambda_lp_phased_balance ?: Channel.empty()
cnloh_lp_phased_balance = params.cnloh_lp_phased_balance ?: Channel.empty()
major_lp_phased_balance = params.major_lp_phased_balance ?: Channel.empty()
allin_lp_phased_balance = params.allin_lp_phased_balance ?: Channel.empty()
marginal_lp_phased_balance = params.marginal_lp_phased_balance ?: Channel.empty()
from_maf_lp_phased_balance = params.from_maf_lp_phased_balance ?: Channel.empty()
ism_lp_phased_balance = params.ism_lp_phased_balance ?: Channel.empty()
epgap_lp_phased_balance = params.epgap_lp_phased_balance ?: Channel.empty()
hets_thresh_lp_phased_balance = params.hets_thresh_lp_phased_balance ?: Channel.empty()
min_bins_lp_phased_balance = params.min_bins_lp_phased_balance ?: Channel.empty()
min_width_lp_phased_balance = params.min_width_lp_phased_balance || params.min_width_lp_phased_balance == 0 ? params.min_width_lp_phased_balance : Channel.empty()
trelim_lp_phased_balance = params.trelim_lp_phased_balance ?: Channel.empty()
reward_lp_phased_balance = params.reward_lp_phased_balance ?: Channel.empty()
nodefileind_lp_phased_balance = params.nodefileind_lp_phased_balance ?: Channel.empty()
tilim_lp_phased_balance = params.tilim_lp_phased_balance ?: Channel.empty()

// HRDetect
ref_genome_version_hrdetect = params.ref_hrdetect ?: Channel.empty()

// Sage
ref_genome_version_sage       = params.ref_genome_version   ?: Channel.empty()

// Pave
ref_genome_version_pave       = params.ref_genome_version   ?: Channel.empty()

// SigprofilerAssignment
sigprofilerassignment_genome         = params.sigprofilerassignment_genome   ?: Channel.empty()
sigprofilerassignment_cosmic_version = params.sigprofilerassignment_cosmic_version ?: Channel.empty()


// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
if (params.snpeff_cache && params.tools && params.tools.contains("snpeff")) {
    def snpeff_annotation_cache_key = params.use_annotation_cache_keys ? "${params.snpeff_genome}.${params.snpeff_db}/" : ""
    def snpeff_cache_dir =  "${snpeff_annotation_cache_key}${params.snpeff_genome}.${params.snpeff_db}"
    def snpeff_cache_path_full = file("$params.snpeff_cache/$snpeff_cache_dir", type: 'dir')
    if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
        error("Files within --snpeff_cache invalid. Make sure there is a directory named ${snpeff_cache_dir} in ${params.snpeff_cache}.\nhttps://nf-co.re/sarek/dev/usage#how-to-customise-snpeff-and-vep-annotation")
    }
    snpeff_cache = Channel.fromPath(file("${params.snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
        .map{ cache -> [ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], cache ] }
} else snpeff_cache = []

if (params.vep_cache && params.tools && params.tools.contains("vep")) {
    def vep_annotation_cache_key = params.use_annotation_cache_keys ? "${params.vep_cache_version}_${params.vep_genome}/" : ""
    def vep_cache_dir = "${vep_annotation_cache_key}${params.vep_species}/${params.vep_cache_version}_${params.vep_genome}"
    def vep_cache_path_full = file("$params.vep_cache/$vep_cache_dir", type: 'dir')
    if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
        error("Files within --vep_cache invalid. Make sure there is a directory named ${vep_cache_dir} in ${params.vep_cache}.\nhttps://nf-co.re/sarek/dev/usage#how-to-customise-snpeff-and-vep-annotation")
    }
    vep_cache = Channel.fromPath(file("${params.vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
} else vep_cache = []

vep_extra_files = []
