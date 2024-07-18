/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'
def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//Check input path parameters to see if they exist
def checkPathParamList = [
    params.ascat_alleles,
    params.ascat_loci,
    params.ascat_loci_gc,
    params.ascat_loci_rt,
    params.bwa,
    params.bwamem2,
    params.cf_chrom_len,
    params.chr_dir,
    params.cnvkit_reference,
    params.dbsnp,
    params.dbsnp_tbi,
    params.dict,
    params.fasta,
    params.fasta_fai,
    params.germline_resource,
    params.germline_resource_tbi,
    params.input,
    params.intervals,
    params.known_indels,
    params.known_indels_tbi,
    params.known_snps,
    params.known_snps_tbi,
    params.mappability,
    params.multiqc_config,
]

def toolParamMap = [
    "fragcounter": [
        params.gcmapdir_frag
    ],
    "dryclean": [
        params.pon_dryclean,
    ],
    "gridss": [
        params.blacklist_gridss,
        params.pon_gridss
    ],
    "hetpileups" : [
        params.hapmap_sites
    ],
    "fusions"    : [
        params.gencode_fusions
    ],
    "allelic_cn" : [
        params.mask_non_integer_balance,
        params.mask_lp_phased_balance
    ],
    "vep"        : [
        params.vep_cache
    ],
    "snpeff"     : [
        params.snpeff_cache
    ],
    "sage"       : [
        params.ensembl_data_dir,
        params.somatic_hotspots,
        params.panel_bed,
        params.high_confidence_bed
    ],
]

skip_tools = params.skip_tools ? params.skip_tools.split(',') : []

if (!params.dbsnp && !params.known_indels) {
    if (!params.skip_tools || (params.skip_tools && !params.skip_tools.contains('baserecalibrator'))) {
        log.warn "Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command."
    }
    if (params.skip_tools && (!params.skip_tools.contains('haplotypecaller') || !params.skip_tools.contains('sentieon_haplotyper'))) {
        log.warn "If GATK's Haplotypecaller or Sentieon's Haplotyper is specified, without `--dbsnp` or `--known_indels no filtering will be done. For filtering, please provide at least one of `--dbsnp` or `--known_indels`.\nFor more information see FilterVariantTranches (single-sample, default): https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches\nFor more information see VariantRecalibration (--joint_germline): https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator\nFor more information on GATK Best practice germline variant calling: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-"
    }
}

if ((params.download_cache) && (params.snpeff_cache || params.vep_cache)) {
    error("Please specify either `--download_cache` or `--snpeff_cache`, `--vep_cache`.")
}

toolParamMap.each { tool, toolParams ->
    checkPathParamList.addAll(toolParams)
}


// Initialise the workflow
WorkflowNfcasereports.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check if the path parameters exist
@Grab(group='org.codehaus.gpars', module='gpars', version='1.2.1')
import groovyx.gpars.GParsPool

def isAwsCliInstalled() {
    def command = "aws --version"
    def process = command.execute()
    process.waitFor()
    return process.exitValue() == 0
}

def checkS3PathExists(s3Path) {
    def command = "aws s3 ls ${s3Path} --no-sign-request"
    def process = command.execute()
    process.waitFor()
    return process.exitValue() == 0
}

def expandS3Paths(s3Path) {
    def matcher = s3Path =~ /\{([^}]+)\}/
    if (!matcher.find()) {
        return [s3Path]
    }
    def expandedPaths = []
    matcher.each { match ->
        def prefix = s3Path.substring(0, matcher.start())
        def suffix = s3Path.substring(matcher.end())
        def options = match[1].split(',')
        options.each { option ->
            expandedPaths += expandS3Paths(prefix + option + suffix)
        }
    }
    return expandedPaths
}

def awsCliInstalled = isAwsCliInstalled()
if (!awsCliInstalled) {
    println "AWS CLI is not installed/loaded. Will proceed, but S3 paths will not be checked."
}

println "Checking if parameter paths exist..."
int numThreads = 16 // Specify the number of threads to use in the pool
GParsPool.withPool(numThreads) {
    checkPathParamList.eachParallel { param ->
        if (param == null) {
            println "Skipping null path"
            return
        }

        if (param.startsWith("s3://")) {
            if (awsCliInstalled) {
                def expandedPaths = expandS3Paths(param)
                expandedPaths.each { expandedPath ->
                    if (checkS3PathExists(expandedPath)) {
                        println "Path exists: ${expandedPath}"
                    } else {
                        println "Path does not exist: ${expandedPath}"
                    }
                }
            } else {
                println "Skipping S3 path check for: ${param}"
            }
        } else {
            def file = new File(param)
            if (file.exists()) {
                println "Path exists: ${param}"
            } else {
                println "Path does not exist: ${param}"
            }
        }
    }
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}

def inputType = params.input ? "input" : "input_restart"
def ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet(inputType)

inputs = ch_from_samplesheet.map {
    meta,
    fastq_1,
    fastq_2,
    table,
    cram,
    crai,
    bam,
    bai,
    hets,
    frag_cov,
    dryclean_cov,
    ploidy,
    seg,
    nseg,
    vcf,
    vcf2,
    jabba_rds,
    ni_balanced_gg,
    lp_balanced_gg,
    events,
    fusions,
    snv_somatic_vcf,
    snv_somatic_tbi,
    snv_germline_vcf,
    snv_germline_tbi,
    variant_somatic_ann,
    variant_somatic_bcf,
    variant_germline_ann,
    variant_germline_bcf,
    snv_multiplicity,
    sbs_signatures,
    indel_signatures,
    signatures_matrix,
    hrdetect
    -> [
        meta: meta,
        fastq_1: fastq_1,
        fastq_2: fastq_2,
        table: table,
        cram: cram,
        crai: crai,
        bam: bam,
        bai: bai,
        hets: hets,
        frag_cov: frag_cov,
        dryclean_cov: dryclean_cov,
        ploidy: ploidy,
        seg: seg,
        nseg: nseg,
        vcf: vcf,
        vcf2: vcf2,
        jabba_rds: jabba_rds,
        ni_balanced_gg: ni_balanced_gg,
        lp_balanced_gg: lp_balanced_gg,
        events: events,
        fusions: fusions,
        snv_somatic_vcf: snv_somatic_vcf,
        snv_somatic_tbi: snv_somatic_tbi,
        snv_germline_vcf: snv_germline_vcf,
        snv_germline_tbi: snv_germline_tbi,
        variant_somatic_ann: variant_somatic_ann,
        variant_somatic_bcf: variant_somatic_bcf,
        variant_germline_ann: variant_germline_ann,
        variant_germline_bcf: variant_germline_bcf,
        snv_multiplicity: snv_multiplicity,
        sbs_signatures: sbs_signatures,
        indel_signatures: indel_signatures,
        signatures_matrix: signatures_matrix,
        hrdetect: hrdetect
    ]
}

inputs = inputs
    .map {it -> [
            it.meta.patient + it.meta.sample, // create a patient_sample key
            it
        ]
    }
    .tap{ ch_with_patient_sample } // save the channel
    .groupTuple() //group by patient_sample to get all lanes
    .map { patient_sample, ch_items ->
        // get number of lanes per sample
        [ patient_sample, ch_items.size() ]
    }
    .combine(ch_with_patient_sample, by: 0) // for each entry add numLanes
    .map { patient_sample, num_lanes, ch_items ->

        if (ch_items.meta.lane && ch_items.fastq_2) {
            ch_items.meta   = ch_items.meta + [id: "${ch_items.meta.sample}-${ch_items.meta.lane}".toString()]
            def CN         = params.seq_center ? "CN:${params.seq_center}\\t" : ''

            def flowcell   = flowcellLaneFromFastq(fastq_1)
            // Don't use a random element for ID, it breaks resuming
            def read_group = "\"@RG\\tID:${flowcell}.${ch_items.meta.sample}.${ch_items.meta.lane}\\t${CN}PU:${ch_items.meta.lane}\\tSM:${ch_items.meta.patient}_${ch_items.meta.sample}\\tLB:${ch_items.meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

            ch_items.meta           = ch_items.meta - ch_items.meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), size: 1]

        } else if (ch_items.meta.lane && ch_items.bam) {
            ch_items.meta   = ch_items.meta + [id: "${ch_items.meta.sample}-${ch_items.meta.lane}".toString()]
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${ch_items.meta.sample}_${ch_items.meta.lane}\\t${CN}PU:${ch_items.meta.lane}\\tSM:${ch_items.meta.patient}_${ch_items.meta.sample}\\tLB:${ch_items.meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

            ch_items.meta = ch_items.meta - ch_items.meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), size: 1]
        } else {
            ch_items.meta = ch_items.meta + [id: ch_items.meta.sample.toString()]
        }

        ch_items
    }

// Fails when missing sex information for CNV tools
inputs.map{
    if (it.meta.sex == 'NA') {
        log.warn "Please specify sex information for each sample in your samplesheet when using '--tools' with 'ascat' if known for comparison"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope

ascat_alleles                       = WorkflowNfcasereports.create_file_channel(params.ascat_alleles)
ascat_loci                          = WorkflowNfcasereports.create_file_channel(params.ascat_loci)
ascat_loci_gc                       = WorkflowNfcasereports.create_file_channel(params.ascat_loci_gc)
ascat_loci_rt                       = WorkflowNfcasereports.create_file_channel(params.ascat_loci_rt)
cf_chrom_len                        = WorkflowNfcasereports.create_file_channel(params.cf_chrom_len)
chr_dir                             = WorkflowNfcasereports.create_file_channel(params.chr_dir)
dbsnp                               = WorkflowNfcasereports.create_file_channel(params.dbsnp)
fasta                               = WorkflowNfcasereports.create_file_channel(params.fasta)
fasta_fai                           = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
germline_resource                   = WorkflowNfcasereports.create_file_channel(params.germline_resource)
known_indels                        = WorkflowNfcasereports.create_file_channel(params.known_indels)
known_snps                          = WorkflowNfcasereports.create_file_channel(params.known_snps)
mappability                         = WorkflowNfcasereports.create_file_channel(params.mappability)
pon                                 = WorkflowNfcasereports.create_file_channel(params.pon)

// SVABA
indel_mask                          = WorkflowNfcasereports.create_file_channel(params.indel_mask)
germ_sv_db                          = WorkflowNfcasereports.create_file_channel(params.germ_sv_db)
simple_seq_db                       = WorkflowNfcasereports.create_file_channel(params.simple_seq_db)

// GRIDSS
blacklist_gridss                    = WorkflowNfcasereports.create_file_channel(params.blacklist_gridss)
pon_gridss                          = WorkflowNfcasereports.create_file_channel(params.pon_gridss)

//SV Junction Filtering
junction_pon_svaba                  = WorkflowNfcasereports.create_file_channel(params.junction_pon_svaba)
junction_pon_gridss                 = WorkflowNfcasereports.create_file_channel(params.junction_pon_gridss)
gnomAD_sv_db                        = WorkflowNfcasereports.create_file_channel(params.gnomAD_sv_db)

// FragCounter
gcmapdir_frag                       = WorkflowNfcasereports.create_file_channel(params.gcmapdir_frag)

// HetPileups
hapmap_sites                        = WorkflowNfcasereports.create_file_channel(params.hapmap_sites)

// Dryclean
pon_dryclean                        = WorkflowNfcasereports.create_file_channel(params.pon_dryclean)
blacklist_path_dryclean             = WorkflowNfcasereports.create_file_channel(params.blacklist_path_dryclean)
germline_file_dryclean              = WorkflowNfcasereports.create_file_channel(params.germline_file_dryclean)

// JaBbA
blacklist_coverage_jabba		    = WorkflowNfcasereports.create_file_channel(params.blacklist_coverage_jabba)

// Fusions
gencode_fusions                     = WorkflowNfcasereports.create_file_channel(params.gencode_fusions)

// Allelic CN
mask_non_integer_balance            = WorkflowNfcasereports.create_file_channel(params.mask_non_integer_balance)
mask_lp_phased_balance              = WorkflowNfcasereports.create_file_channel(params.mask_lp_phased_balance)

// Sage
somatic_hotspots_sage               = WorkflowNfcasereports.create_file_channel(params.somatic_hotspots)
panel_bed_sage                      = WorkflowNfcasereports.create_file_channel(params.panel_bed)
high_confidence_bed_sage            = WorkflowNfcasereports.create_file_channel(params.high_confidence_bed)
ensembl_data_dir_sage               = WorkflowNfcasereports.create_file_channel(params.ensembl_data_dir)
gnomAD_snv_db                       = WorkflowNfcasereports.create_file_channel(params.gnomAD_snv_db)
gnomAD_snv_db_tbi                   = WorkflowNfcasereports.create_file_channel(params.gnomAD_snv_db_tbi)
sage_germline_pon                   = WorkflowNfcasereports.create_file_channel(params.sage_germline_pon)
sage_germline_pon_tbi               = WorkflowNfcasereports.create_file_channel(params.sage_germline_pon_tbi)

// Pave
sage_pon_pave                       = WorkflowNfcasereports.create_file_channel(params.sage_pon)
sage_blocklist_regions_pave         = WorkflowNfcasereports.create_file_channel(params.sage_blocklist_regions)
sage_blocklist_sites_pave           = WorkflowNfcasereports.create_file_channel(params.sage_blocklist_sites)
clinvar_annotations_pave            = WorkflowNfcasereports.create_file_channel(params.clinvar_annotations)
segment_mappability_pave            = WorkflowNfcasereports.create_file_channel(params.segment_mappability)
driver_gene_panel_pave              = WorkflowNfcasereports.create_file_channel(params.driver_gene_panel)
ensembl_data_resources_pave         = WorkflowNfcasereports.create_file_channel(params.ensembl_data_resources)
gnomad_resource_pave                = WorkflowNfcasereports.create_file_channel(params.gnomad_resource)


// Initialize value channels based on params, defined in the params.genomes[params.genome] scope

ascat_genome        = WorkflowNfcasereports.create_value_channel(params.ascat_genome)
dbsnp_vqsr          = WorkflowNfcasereports.create_value_channel(params.dbsnp_vqsr)
known_indels_vqsr   = WorkflowNfcasereports.create_value_channel(params.known_indels_vqsr)
known_snps_vqsr     = WorkflowNfcasereports.create_value_channel(params.known_snps_vqsr)
snpeff_genome       = WorkflowNfcasereports.create_value_channel(params.snpeff_genome)
snpeff_db           = WorkflowNfcasereports.create_value_channel(params.snpeff_db)
snpeff_db_full     = params.snpeff_db && params.snpeff_genome   ? Channel.value("${params.snpeff_genome}.${params.snpeff_db}") : Channel.empty()
vep_cache_version   = WorkflowNfcasereports.create_value_channel(params.vep_cache_version)
vep_genome          = WorkflowNfcasereports.create_value_channel(params.vep_genome)
vep_species         = WorkflowNfcasereports.create_value_channel(params.vep_species)
error_rate          = WorkflowNfcasereports.create_value_channel(params.error_rate)

//SV Junction Filter
junction_padding     = WorkflowNfcasereports.create_value_channel(params.pad_junc_filter)

// Hetpileups
filter_hets          = WorkflowNfcasereports.create_value_channel(params.filter_hets)
max_depth            = WorkflowNfcasereports.create_value_channel(params.max_depth)

// FragCounter
windowsize_frag     = WorkflowNfcasereports.create_value_channel(params.windowsize_frag)
minmapq_frag        = WorkflowNfcasereports.create_value_channel(params.minmapq_frag)
midpoint_frag       = WorkflowNfcasereports.create_value_channel(params.midpoint_frag)
paired_frag         = WorkflowNfcasereports.create_value_channel(params.paired_frag)
exome_frag          = WorkflowNfcasereports.create_value_channel(params.exome_frag)

// Dryclean
center_dryclean            = WorkflowNfcasereports.create_value_channel(params.center_dryclean)
cbs_dryclean                 = WorkflowNfcasereports.create_value_channel(params.cbs_dryclean)
cnsignif_dryclean            = WorkflowNfcasereports.create_value_channel(params.cnsignif_dryclean)
wholeGenome_dryclean         = WorkflowNfcasereports.create_value_channel(params.wholeGenome_dryclean)
blacklist_dryclean           = WorkflowNfcasereports.create_value_channel(params.blacklist_dryclean)
germline_filter_dryclean     = WorkflowNfcasereports.create_value_channel(params.germline_filter_dryclean)
field_dryclean               = WorkflowNfcasereports.create_value_channel(params.field_dryclean)
build_dryclean               = WorkflowNfcasereports.create_value_channel(params.build_dryclean)

// ASCAT_seg
field_ascat                  = WorkflowNfcasereports.create_value_channel(params.field_ascat)
hets_thresh_ascat            = WorkflowNfcasereports.create_value_channel(params.hets_thresh_ascat)
penalty_ascat                = WorkflowNfcasereports.create_value_channel(params.penalty_ascat)
gc_correct_ascat             = WorkflowNfcasereports.create_value_channel(params.gc_correct_ascat)
rebin_width_ascat            = WorkflowNfcasereports.create_value_channel(params.rebin_width_ascat)
from_maf_ascat               = WorkflowNfcasereports.create_value_channel(params.from_maf_ascat)

// CBS
cnsignif_cbs                     = WorkflowNfcasereports.create_value_channel(params.cnsignif_cbs)
field_cbs                        = WorkflowNfcasereports.create_value_channel(params.field_cbs)
name_cbs                         = WorkflowNfcasereports.create_value_channel(params.name_cbs)

// JaBbA
blacklist_junctions_jabba        = WorkflowNfcasereports.create_value_channel(params.blacklist_junctions_jabba)
geno_jabba					     = WorkflowNfcasereports.create_value_channel(params.geno_jabba)
indel_jabba					     = WorkflowNfcasereports.create_value_channel(params.indel_jabba)
tfield_jabba					 = WorkflowNfcasereports.create_value_channel(params.tfield_jabba)
iter_jabba					     = WorkflowNfcasereports.create_value_channel(params.iter_jabba)
rescue_window_jabba				 = WorkflowNfcasereports.create_value_channel(params.rescue_window_jabba)
rescue_all_jabba				 = WorkflowNfcasereports.create_value_channel(params.rescue_all_jabba)
nudgebalanced_jabba				 = WorkflowNfcasereports.create_value_channel(params.nudgebalanced_jabba)
edgenudge_jabba					 = WorkflowNfcasereports.create_value_channel(params.edgenudge_jabba)
strict_jabba					 = WorkflowNfcasereports.create_value_channel(params.strict_jabba)
allin_jabba					     = WorkflowNfcasereports.create_value_channel(params.allin_jabba)
field_jabba					     = WorkflowNfcasereports.create_value_channel(params.field_jabba)
maxna_jabba					     = WorkflowNfcasereports.create_value_channel(params.maxna_jabba)
purity_jabba					 = WorkflowNfcasereports.create_value_channel(params.purity_jabba)
ploidy_jab     					 = WorkflowNfcasereports.create_value_channel(params.ploidy_jabba)
pp_method_jabba					 = WorkflowNfcasereports.create_value_channel(params.pp_method_jabba)
cnsignif_jabba					 = WorkflowNfcasereports.create_value_channel(params.cnsignif_jabba)
slack_jabba					     = WorkflowNfcasereports.create_value_channel(params.slack_jabba)
linear_jabba					 = WorkflowNfcasereports.create_value_channel(params.linear_jabba)
tilim_jabba					     = WorkflowNfcasereports.create_value_channel(params.tilim_jabba)
epgap_jabba					     = WorkflowNfcasereports.create_value_channel(params.epgap_jabba)
fix_thres_jabba					 = WorkflowNfcasereports.create_value_channel(params.fix_thres_jabba)
lp_jabba					     = WorkflowNfcasereports.create_value_channel(params.lp_jabba)
ism_jabba					     = WorkflowNfcasereports.create_value_channel(params.ism_jabba)
filter_loose_jabba				 = WorkflowNfcasereports.create_value_channel(params.filter_loose_jabba)
gurobi_jabba					 = WorkflowNfcasereports.create_value_channel(params.gurobi_jabba)
nonintegral_jabba				 = WorkflowNfcasereports.create_value_channel(params.nonintegral_jabba)
verbose_jabba					 = WorkflowNfcasereports.create_value_channel(params.verbose_jabba)
help_jabba					     = WorkflowNfcasereports.create_value_channel(params.help_jabba)

//Allelic CN (Non-integer balance)
field_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.field_non_integer_balance)
hets_thresh_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.hets_thresh_non_integer_balance)
overwrite_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.overwrite_non_integer_balance)
lambda_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.lambda_non_integer_balance)
allin_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.allin_non_integer_balance)
fix_thresh_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.fix_thresh_non_integer_balance)
nodebounds_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.nodebounds_non_integer_balance)
ism_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.ism_non_integer_balance)
build_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.build_non_integer_balance)
epgap_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.epgap_non_integer_balance)
tilim_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.tilim_non_integer_balance)
gurobi_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.gurobi_non_integer_balance)
pad_non_integer_balance  = WorkflowNfcasereports.create_value_channel(params.pad_non_integer_balance)

// ...(LP Phased balance)
lambda_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.lambda_lp_phased_balance)
cnloh_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.cnloh_lp_phased_balance)
major_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.major_lp_phased_balance)
allin_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.allin_lp_phased_balance)
marginal_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.marginal_lp_phased_balance)
from_maf_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.from_maf_lp_phased_balance)
ism_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.ism_lp_phased_balance)
epgap_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.epgap_lp_phased_balance)
hets_thresh_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.hets_thresh_lp_phased_balance)
min_bins_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.min_bins_lp_phased_balance)
min_width_lp_phased_balance = params.min_width_lp_phased_balance || params.min_width_lp_phased_balance == 0 ? params.min_width_lp_phased_balance : Channel.empty()
trelim_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.trelim_lp_phased_balance)
reward_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.reward_lp_phased_balance)
nodefileind_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.nodefileind_lp_phased_balance)
tilim_lp_phased_balance  = WorkflowNfcasereports.create_value_channel(params.tilim_lp_phased_balance)

// HRDetect
ref_genome_version_hrdetect  = WorkflowNfcasereports.create_value_channel(params.ref_hrdetect)

// Sage
ref_genome_version_sage        = WorkflowNfcasereports.create_value_channel(params.ref_genome_version)

// Pave
ref_genome_version_pave        = WorkflowNfcasereports.create_value_channel(params.ref_genome_version)

// SigprofilerAssignment
sigprofilerassignment_genome          = WorkflowNfcasereports.create_value_channel(params.sigprofilerassignment_genome)
sigprofilerassignment_cosmic_version  = WorkflowNfcasereports.create_value_channel(params.sigprofilerassignment_cosmic_version)

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
if (params.snpeff_cache) {
    def snpeff_annotation_cache_key = params.use_annotation_cache_keys ? "${params.snpeff_genome}.${params.snpeff_db}/" : ""
    def snpeff_cache_dir =  "${snpeff_annotation_cache_key}${params.snpeff_genome}.${params.snpeff_db}"
    def snpeff_cache_path_full = file("$params.snpeff_cache/$snpeff_cache_dir", type: 'dir')
    if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
        error("Files within --snpeff_cache invalid. Make sure there is a directory named ${snpeff_cache_dir} in ${params.snpeff_cache}.\nhttps://nf-co.re/sarek/dev/usage#how-to-customise-snpeff-and-vep-annotation")
    }
    snpeff_cache = Channel.fromPath(file("${params.snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
        .map{ cache -> [ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], cache ] }
} else snpeff_cache = []

if (params.vep_cache) {
    def vep_annotation_cache_key = params.use_annotation_cache_keys ? "${params.vep_cache_version}_${params.vep_genome}/" : ""
    def vep_cache_dir = "${vep_annotation_cache_key}${params.vep_species}/${params.vep_cache_version}_${params.vep_genome}"
    def vep_cache_path_full = file("$params.vep_cache/$vep_cache_dir", type: 'dir')
    if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
        error("Files within --vep_cache invalid. Make sure there is a directory named ${vep_cache_dir} in ${params.vep_cache}.\nhttps://nf-co.re/sarek/dev/usage#how-to-customise-snpeff-and-vep-annotation")
    }
    vep_cache = Channel.fromPath(file("${params.vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
} else vep_cache = []

vep_extra_files = []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Create samplesheets to restart from different steps
include { CHANNEL_ALIGN_CREATE_CSV                    } from '../subworkflows/local/channel_align_create_csv/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV           } from '../subworkflows/local/channel_markduplicates_create_csv/main'
include { CHANNEL_BASERECALIBRATOR_CREATE_CSV         } from '../subworkflows/local/channel_baserecalibrator_create_csv/main'
include { CHANNEL_APPLYBQSR_CREATE_CSV                } from '../subworkflows/local/channel_applybqsr_create_csv/main'
include { CHANNEL_SVCALLING_CREATE_CSV                } from '../subworkflows/local/channel_svcalling_create_csv/main'
include { CHANNEL_HETPILEUPS_CREATE_CSV               } from '../subworkflows/local/channel_hetpileups_create_csv/main'

// Download annotation cache if needed
include { PREPARE_CACHE                               } from '../subworkflows/local/prepare_cache/main'

// Build indices if needed
include { PREPARE_GENOME                              } from '../subworkflows/local/prepare_genome/main'

// Build intervals if needed
include { PREPARE_INTERVALS                           } from '../subworkflows/local/prepare_intervals/main'

// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT } from '../subworkflows/local/bam_convert_samtools/main'

// Run FASTQC
include { FASTQC                                      } from '../modules/nf-core/fastqc/main'

// TRIM/SPLIT FASTQ Files
include { FASTP                                       } from '../modules/nf-core/fastp/main'

// Loading the MULTIQC module
include { MULTIQC                                     } from '../modules/nf-core/multiqc/main'

// Loading the module that dumps the versions of software being used
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// Map input reads to reference genome
include { FASTQ_ALIGN_BWAMEM_MEM2                     } from '../subworkflows/local/fastq_align_bwamem_mem2/main'
include { FASTQ_PARABRICKS_FQ2BAM                     } from '../subworkflows/local/fastq_parabricks_fq2bam/main'

// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                    } from '../subworkflows/local/bam_merge_index_samtools/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING     } from '../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL       } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_FINAL       } from '../modules/nf-core/samtools/convert/main'

// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES                          } from '../subworkflows/local/bam_markduplicates/main'

// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'

// Create recalibration tables
include { BAM_BASERECALIBRATOR                        } from '../subworkflows/local/bam_baserecalibrator/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { BAM_APPLYBQSR                               } from '../subworkflows/local/bam_applybqsr/main'

// Svaba
include { BAM_SVCALLING_SVABA                         } from '../subworkflows/local/bam_svcalling_svaba/main'

//GRIDSS
include { BAM_SVCALLING_GRIDSS                        } from '../subworkflows/local/bam_svcalling_gridss/main'
include { BAM_SVCALLING_GRIDSS_SOMATIC                } from '../subworkflows/local/bam_svcalling_gridss/main'

// SV Junction Filtering
include { SV_JUNCTION_FILTER as JUNCTION_FILTER       } from '../subworkflows/local/junction_filter/main'

// HETPILEUPS
include { BAM_HETPILEUPS                              } from '../subworkflows/local/bam_hetpileups/main'

// fragCounter
include { BAM_FRAGCOUNTER as TUMOR_FRAGCOUNTER         } from '../subworkflows/local/bam_fragCounter/main'
include { BAM_FRAGCOUNTER as NORMAL_FRAGCOUNTER        } from '../subworkflows/local/bam_fragCounter/main'

// dryclean
include { COV_DRYCLEAN as TUMOR_DRYCLEAN               } from '../subworkflows/local/cov_dryclean/main'
include { COV_DRYCLEAN as NORMAL_DRYCLEAN              } from '../subworkflows/local/cov_dryclean/main'

//ASCAT
include { COV_ASCAT                                   } from '../subworkflows/local/cov_ascat/main'
include { EXTRACT_PURITYPLOIDY                        } from '../modules/local/ascat/main'

// CBS
include { COV_CBS as CBS                              } from '../subworkflows/local/cov_cbs/main'

// JaBbA
if (params.tumor_only) {
    include { COV_JUNC_TUMOR_ONLY_JABBA as JABBA } from '../subworkflows/local/jabba/main'
} else {
    include { COV_JUNC_JABBA as JABBA } from '../subworkflows/local/jabba/main'
}

// Events
include { GGRAPH_EVENTS as EVENTS                            } from '../subworkflows/local/events/main'

// Fusions
include { GGRAPH_FUSIONS as FUSIONS                            } from '../subworkflows/local/fusions/main'

// Allelic CN
include { COV_GGRAPH_NON_INTEGER_BALANCE as NON_INTEGER_BALANCE                            } from '../subworkflows/local/allelic_cn/main'

include { COV_GGRAPH_LP_PHASED_BALANCE as LP_PHASED_BALANCE                            } from '../subworkflows/local/allelic_cn/main'

// HRDetect
include { JUNC_SNV_GGRAPH_HRDETECT } from '../subworkflows/local/hrdetect/main'

//STRELKA2
include { BAM_SOMATIC_STRELKA                        } from '../subworkflows/local/bam_somatic_strelka/main'
include { BAM_GERMLINE_STRELKA                       } from '../subworkflows/local/bam_germline_strelka/main'

// SAGE
include { BAM_SAGE                              } from '../subworkflows/local/bam_sage/main'
include { BAM_SAGE_TUMOR_ONLY_FILTER            } from '../subworkflows/local/bam_sage/main'

// SNPEFF
include { VCF_SNPEFF as VCF_SNPEFF_SOMATIC                              } from '../subworkflows/local/vcf_snpeff/main'
include { VCF_SNPEFF as VCF_SNPEFF_GERMLINE                             } from '../subworkflows/local/vcf_snpeff/main'

// SNV MULTIPLICITY
include { VCF_SNV_MULTIPLICITY           } from '../subworkflows/local/vcf_snv_multiplicity/main'

// PAVE
include { VCF_PAVE as VCF_PAVE_SOMATIC                              } from '../subworkflows/local/vcf_pave/main'
include { VCF_PAVE as VCF_PAVE_GERMLINE                              } from '../subworkflows/local/vcf_pave/main'

// SigProfilerAssignment
include { VCF_SIGPROFILERASSIGNMENT                              } from '../subworkflows/local/vcf_sigprofilerassignment/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCASEREPORTS {

    // MULTIQC
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // To gather all QC reports for MultiQC
    reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    versions = Channel.empty()

        // Download cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache    ? [] : Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    snpeff_info     = params.snpeff_cache ? [] : Channel.of([ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], params.snpeff_genome, params.snpeff_db ])

    if (params.download_cache) {
        PREPARE_CACHE(ensemblvep_info, snpeff_info)
        snpeff_cache = PREPARE_CACHE.out.snpeff_cache
        vep_cache    = PREPARE_CACHE.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

        versions = versions.mix(PREPARE_CACHE.out.versions)
    }

    // Build indices if needed
    // ##############################
    PREPARE_GENOME()

    // Gather built indices or get them from the params
    // Built from the fasta file:
    dict       = params.dict        ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                    : PREPARE_GENOME.out.dict
    fasta_fai  = WorkflowNfcasereports.create_file_channel(params.fasta_fai, PREPARE_GENOME.out.fasta_fai)
    bwa        = WorkflowNfcasereports.create_file_channel(params.bwa, PREPARE_GENOME.out.bwa)
    bwamem2    = WorkflowNfcasereports.create_file_channel(params.bwamem2, PREPARE_GENOME.out.bwamem2)

    // Gather index for mapping given the chosen aligner
    index_alignment = (params.aligner == "bwa-mem") ? bwa :
        params.aligner == "bwa-mem2" ? bwamem2 : null

    // TODO: add a params for msisensorpro_scan
    msisensorpro_scan      = PREPARE_GENOME.out.msisensorpro_scan

    // For ASCAT, extracted from zip or tar.gz files:
    allele_files           = PREPARE_GENOME.out.allele_files
    chr_files              = PREPARE_GENOME.out.chr_files
    gc_file                = PREPARE_GENOME.out.gc_file
    loci_files             = PREPARE_GENOME.out.loci_files
    rt_file                = PREPARE_GENOME.out.rt_file

    dbsnp_tbi              = WorkflowNfcasereports.create_index_channel(params.dbsnp, params.dbsnp_tbi, PREPARE_GENOME.out.dbsnp_tbi)
    //do not change to Channel.value([]), the check for its existence then fails for Getpileupsumamries
    germline_resource_tbi  = params.germline_resource ? params.germline_resource_tbi ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : []
    known_indels_tbi       = WorkflowNfcasereports.create_index_channel(params.known_indels, params.known_indels_tbi, PREPARE_GENOME.out.known_indels_tbi)
    known_snps_tbi         = WorkflowNfcasereports.create_index_channel(params.known_snps, params.known_snps_tbi, PREPARE_GENOME.out.known_snps_tbi)
    pon_tbi                = WorkflowNfcasereports.create_index_channel(params.pon, params.pon_tbi, PREPARE_GENOME.out.pon_tbi)

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    known_sites_snps       = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi   = dbsnp_tbi.concat(known_snps_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai, params.intervals, params.no_intervals)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined         = WorkflowNfcasereports.create_file_channel(params.no_intervals, PREPARE_INTERVALS.out.intervals_bed_combined)
    intervals_bed_gz_tbi_combined  = WorkflowNfcasereports.create_file_channel(params.no_intervals, PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined)

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = params.wes ?
        intervals_bed_combined.map{it -> [ [ id:it.baseName ], it ]}.collect() :
        Channel.value([ [ id:'null' ], [] ])

    intervals            = PREPARE_INTERVALS.out.intervals_bed        // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
        if ( num_intervals < 1 ) [ [], num_intervals ]
        else [ interval, num_intervals ]
    }

    intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map{ intervals, num_intervals ->
        if ( num_intervals < 1 ) [ [], [], num_intervals ]
        else [ intervals[0], intervals[1], num_intervals ]
    }

    // Gather used softwares versions
    versions = versions.mix(PREPARE_GENOME.out.versions)
    versions = versions.mix(PREPARE_INTERVALS.out.versions)

    // Alignment
    // ##############################

    input_fastq = inputs.filter { it.bam.isEmpty() }.map { it -> [it.meta, [it.fastq_1, it.fastq_2]] }
    alignment_existing_outputs = inputs.map { it -> [it.meta, it.bam] }.filter { !it[1].isEmpty() }

    // QC
    FASTQC(input_fastq)

    reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
    versions = versions.mix(FASTQC.out.versions.first())

    //skipping the UMI Conscensus calling step for now
    reads_for_fastp = input_fastq

    // Trimming and/or splitting
    if (params.trim_fastq && params.split_fastq > 0) {
        log.warn "You have mentioned trim_fastq to `$params.trim_fastq`, will do trimming"
        save_trimmed_fail = false
        save_merged = false
        FASTP(
            reads_for_fastp,
            [], // we are not using any adapter fastas at the moment
            save_trimmed_fail,
            save_merged
        )

        reports = reports.mix(FASTP.out.json.collect{ meta, json -> json })
        reports = reports.mix(FASTP.out.html.collect{ meta, html -> html })

        if (params.split_fastq) {
            reads_for_alignment = FASTP.out.reads.map{ meta, reads ->
                read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                [ meta + [ size:read_files.size() ], read_files ]
            }.transpose()
        } else reads_for_alignment = FASTP.out.reads

        versions = versions.mix(FASTP.out.versions)

    } else {
        println "Skipping trimming since trim_fastq is false"
        reads_for_alignment = reads_for_fastp
    }

    // STEP 1: MAPPING READS TO REFERENCE GENOME
    // reads will be sorted
    reads_for_alignment = reads_for_alignment.map{ meta, reads ->
        // Update meta.id to meta.sample no multiple lanes or splitted fastqs
        if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
        else [ meta, reads ]
    }

    // GPU Alignment
    if (params.aligner == "fq2bam") {
        FASTQ_PARABRICKS_FQ2BAM(
            input_fastq,
            known_sites_indels
        )
        // merge existing BAMs with newly mapped ones
        bam_mapped = alignment_existing_outputs.mix(FASTQ_PARABRICKS_FQ2BAM.out.bam)
    } else {
        FASTQ_ALIGN_BWAMEM_MEM2(
            reads_for_alignment,
            index_alignment
        )
        // merge existing BAMs with newly mapped ones
        bam_mapped = alignment_existing_outputs.mix(FASTQ_ALIGN_BWAMEM_MEM2.out.bam)
    }



    // Grouping the bams from the same samples not to stall the workflow
    bam_mapped = bam_mapped.map{ meta, bam ->

        // Update meta.id to be meta.sample, ditching sample-lane that is not needed anymore
        // Update meta.data_type
        // Remove no longer necessary fields:
        //   read_group: Now in the BAM header
        //    num_lanes: only needed for mapping
        //         size: only needed for mapping

        // Use groupKey to make sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        [ groupKey( meta - meta.subMap('num_lanes', 'read_group', 'size') + [ id:meta.sample ], (meta.num_lanes ?: 1) * (meta.size ?: 1)), bam ]
    }.groupTuple()

    // bams are merged (when multiple lanes from the same sample) and indexed
    BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

    // Gather used softwares versions
    versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)

    // Gather used softwares versions
    versions = versions.mix(FASTQ_ALIGN_BWAMEM_MEM2.out.versions)

    alignment_bams_final = BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai.map({ meta, bam, bai -> [ meta.id, meta, bam, bai ] })

    // BAM Postprocessing
    // ##############################

    if (params.postprocess_bams) {
        cram_markduplicates_no_spark = Channel.empty()

        // STEP 2: markduplicates (+QC) + convert to CRAM

        BAM_MARKDUPLICATES(
            bam_mapped,
            fasta,
            fasta_fai,
            intervals_for_preprocessing
        )

        cram_markduplicates_no_spark = BAM_MARKDUPLICATES.out.cram

        // Gather QC reports
        reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ meta, report -> report })

        // Gather used softwares versions
        versions = versions.mix(BAM_MARKDUPLICATES.out.versions)

        // STEP 3: BASE RECALIBRATION
        ch_cram_for_bam_baserecalibrator = Channel.empty().mix(cram_markduplicates_no_spark)
            // Make sure correct data types are carried through
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        ch_table_bqsr_tab    = Channel.empty()

        BAM_BASERECALIBRATOR(
            ch_cram_for_bam_baserecalibrator,
            dict,
            fasta,
            fasta_fai,
            intervals_and_num_intervals,
            known_sites_indels,
            known_sites_indels_tbi
            )

        ch_table_bqsr_tab = BAM_BASERECALIBRATOR.out.table_bqsr

        versions = versions.mix(BAM_BASERECALIBRATOR.out.versions)

        // ch_table_bqsr contains table from baserecalibrator
        ch_table_bqsr = Channel.empty().mix(ch_table_bqsr_tab)

        reports = reports.mix(ch_table_bqsr.collect{ meta, table -> table })
        cram_applybqsr = ch_cram_for_bam_baserecalibrator.join(ch_table_bqsr, failOnDuplicate: true, failOnMismatch: true)

        // STEP 4: RECALIBRATING

        cram_variant_calling = Channel.empty()

        BAM_APPLYBQSR(
            cram_applybqsr,
            dict,
            fasta,
            fasta_fai,
            intervals_and_num_intervals
        )

        cram_variant_calling = BAM_APPLYBQSR.out.cram

        // Gather used softwares versions
        versions = versions.mix(BAM_APPLYBQSR.out.versions)

        CRAM_QC_RECAL(
            cram_variant_calling,
            fasta,
            intervals_for_preprocessing
        )

        // Gather QC
        reports = reports.mix(CRAM_QC_RECAL.out.reports.collect{ meta, report -> report })

        // Gather software versions
        versions = versions.mix(CRAM_QC_RECAL.out.versions)

        // convert CRAM files to BAM for downstream processes
        CRAM_TO_BAM_RECAL(cram_variant_calling, fasta, fasta_fai)
        versions = versions.mix(CRAM_TO_BAM_RECAL.out.versions)

        CRAM_TO_BAM_FINAL(cram_variant_calling, fasta, fasta_fai)
        versions = versions.mix(CRAM_TO_BAM_FINAL.out.versions)

        alignment_bams_final = Channel.empty()
            .mix(CRAM_TO_BAM_FINAL.out.alignment_index)
            .map{ meta, bam, bai -> [ meta.id, meta + [data_type: "bam"], bam, bai ] }
    }

    // SV Calling
    // ##############################

    // Filter out bams for which SV calling has already been done
    bam_sv_inputs = inputs.filter { it.vcf.isEmpty() }.map { it -> [it.meta.id] }
    bam_sv_calling = alignment_bams_final
        .join(bam_sv_inputs)
        .map { it -> [ it[1], it[2], it[3] ] } // meta, bam, bai
    gridss_existing_outputs = inputs.map { it -> [it.meta, it.vcf] }.filter { !it[1].isEmpty() }

    if (params.tumor_only) {
        // add empty arrays to stand-in for normals
        bam_sv_calling_pair = bam_sv_calling_status.tumor.map{ meta, bam, bai -> [ meta, [], [], bam, bai ] }
    } else {
        // getting the tumor and normal cram files separated
        bam_sv_calling_status = bam_sv_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

        // All normal samples
        bam_sv_calling_normal_for_crossing = bam_sv_calling_status.normal.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

        // All tumor samples
        bam_sv_calling_tumor_for_crossing = bam_sv_calling_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

        // Crossing the normal and tumor samples to create tumor and normal pairs
        bam_sv_calling_pair = bam_sv_calling_normal_for_crossing.cross(bam_sv_calling_tumor_for_crossing)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
        }
    }

    BAM_SVCALLING_GRIDSS(
        bam_sv_calling_pair,
        bwa
    )

    vcf_from_gridss_gridss = Channel.empty()
        .mix(BAM_SVCALLING_GRIDSS.out.vcf)
        .mix(gridss_existing_outputs)
    versions = versions.mix(BAM_SVCALLING_GRIDSS.out.versions)

    if (!params.tumor_only) {
        //somatic filter for GRIDSS
        BAM_SVCALLING_GRIDSS_SOMATIC(vcf_from_gridss_gridss, pon_gridss)

        versions = versions.mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.versions)
        vcf_from_sv_calling = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_high_confidence)
        vcf_from_sv_calling_for_merge = vcf_from_sv_calling
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, vcf

        unfiltered_som_sv = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_all)
        unfiltered_som_sv_for_merge = unfiltered_som_sv
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, vcf
    }

    // if tumor only, filter the junctions
    if (params.tumor_only) {
        JUNCTION_FILTER(
            vcf_from_gridss_gridss,
            junction_padding
        )

        final_filtered_sv_rds = Channel.empty().mix(JUNCTION_FILTER.out.final_filtered_sv_rds)
        final_filtered_sv_rds_for_merge = final_filtered_sv_rds
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, vcf
        pon_filtered_sv_rds = Channel.empty().mix(JUNCTION_FILTER.out.pon_filtered_sv_rds)
    }

    // HETPILEUPS
    // ##############################
    if (!params.tumor_only) {
        bam_hetpileups_inputs = inputs.filter { it.hets.isEmpty() }.map { it -> [it.meta.id] }
        bam_hetpileups_calling = alignment_bams_final
            .join(bam_hetpileups_inputs)
            .map{ it -> [ it[1], it[2], it[3] ] } // meta, bam, bai

        hetpileups_existing_outputs = inputs
            .map { it -> [it.meta, it.hets] }
            .filter { !it[1].isEmpty() }
            .branch{
                normal: it[0].status == 0
                tumor:  it[0].status == 1
            }

        // getting the tumor and normal cram files separated
        bam_hetpileups_status = bam_hetpileups_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

        // All normal samples
        bam_hetpileups_normal_for_crossing = bam_hetpileups_status.normal.map{ meta, bam, bai -> [ meta.patient, meta + [id: meta.sample], bam, bai ] }

        // All tumor samples
        bam_hetpileups_tumor_for_crossing = bam_hetpileups_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta + [id: meta.sample], bam, bai ] }

        // Crossing the normal and tumor samples to create tumor and normal pairs
        bam_hetpileups_pair = bam_hetpileups_normal_for_crossing.cross(bam_hetpileups_tumor_for_crossing)
            .map { normal, tumor ->
                def meta = [:]
                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
            }


        BAM_HETPILEUPS(bam_hetpileups_pair)
        versions = versions.mix(BAM_HETPILEUPS.out.versions)

        sites_from_het_pileups_wgs = Channel.empty()
            .mix(BAM_HETPILEUPS.out.het_pileups_wgs)
            .mix(hetpileups_existing_outputs)

        hets_sites_for_merge = sites_from_het_pileups_wgs
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, hets
    }
    // FRAGCOUNTER
    // ##############################

    bam_fragcounter_inputs = inputs.filter { it.frag_cov.isEmpty() }.map { it -> [it.meta.id] }
    bam_fragcounter_calling = alignment_bams_final
        .join(bam_fragcounter_inputs)
        .map{ it -> [ it[1], it[2], it[3] ] } // meta, bam, bai

    fragcounter_existing_outputs = inputs
        .map { it -> [it.meta, it.frag_cov] }
        .filter { !it[1].isEmpty() }
        .branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

    if (params.tumor_only) {
        bam_fragcounter_status = bam_fragcounter_calling.branch{
            tumor:  it[0].status == 1
        }
    } else {
        // getting the tumor and normal bam files separated
        bam_fragcounter_status = bam_fragcounter_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }
    }

    if (!params.tumor_only) {
        NORMAL_FRAGCOUNTER(bam_fragcounter_status.normal)
        normal_frag_cov = Channel.empty()
            .mix(NORMAL_FRAGCOUNTER.out.fragcounter_cov)
            .mix(fragcounter_existing_outputs.normal)
    }

    TUMOR_FRAGCOUNTER(bam_fragcounter_status.tumor)

    tumor_frag_cov = Channel.empty()
        .mix(TUMOR_FRAGCOUNTER.out.fragcounter_cov)
        .mix(fragcounter_existing_outputs.tumor)

    // Only need one versions because its just one program (fragcounter)
    versions = versions.mix(TUMOR_FRAGCOUNTER.out.versions)

    // Dryclean
    // ##############################
    cov_dryclean_inputs = inputs
        .filter { it.dryclean_cov.isEmpty() }
        .map { it -> [it.meta] }
        .branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

    cov_dryclean_tumor_input = tumor_frag_cov
        .join(cov_dryclean_inputs.tumor)
        .map{ it -> [ it[0], it[1] ] } // meta, frag_cov

    cov_dryclean_normal_input = normal_frag_cov
        .join(cov_dryclean_inputs.normal)
        .map{ it -> [ it[0], it[1] ] } // meta, frag_cov

    dryclean_existing_outputs = inputs
        .map { it -> [it.meta, it.dryclean_cov] }
        .filter { !it[1].isEmpty() }
        .branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

    // Dryclean for both tumor & normal
    TUMOR_DRYCLEAN(cov_dryclean_tumor_input)

    dryclean_tumor_cov = Channel.empty()
        .mix(TUMOR_DRYCLEAN.out.dryclean_cov)
        .mix(dryclean_existing_outputs.tumor)
    dryclean_tumor_cov_for_merge = dryclean_tumor_cov
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, dryclean_cov

    // Only need one versions because it's one program (dryclean)
    versions = versions.mix(TUMOR_DRYCLEAN.out.versions)

    if (!params.tumor_only) {
        NORMAL_DRYCLEAN(cov_dryclean_normal_input)

        dryclean_normal_cov = Channel.empty()
            .mix(NORMAL_DRYCLEAN.out.dryclean_cov)
            .mix(dryclean_existing_outputs.normal)
    }

    // CBS
    // ##############################

    cbs_inputs = inputs
        .filter { it.seg.isEmpty() || it.nseg.isEmpty() }
        .map { it -> [it.meta] }
        .branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

    cbs_tumor_input = dryclean_tumor_cov
        .join(cbs_inputs.tumor)
        .map{ it -> [ it[0].patient, it[0], it[1] ] } // meta.patient, meta, dryclean tumor cov

    cbs_normal_input = dryclean_normal_cov
        .join(cbs_inputs.normal)
        .map{ it -> [ it[0].patient, it[0], it[1] ] } // meta.patient, meta, dryclean normal cov

    cbs_existing_seg = inputs
        .map { it -> [it.meta, it.seg] }
        .filter { !it[1].isEmpty() }

    cbs_existing_nseg = inputs
        .map { it -> [it.meta, it.nseg] }
        .filter { !it[1].isEmpty() }

    cov_cbs = cbs_tumor_input.cross(cbs_normal_input)
        .map { tumor, normal ->
            def meta = [:]
                meta.id             = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.sample         = "${tumor[1].sample}".toString()
                meta.normal_id      = normal[1].sample
                meta.patient        = normal[0]
                meta.sex            = normal[1].sex
                meta.tumor_id       = tumor[1].sample

                [ meta, tumor[2], normal[2] ]
        }

    CBS(cov_cbs)

    versions = versions.mix(CBS.out.versions)

    cbs_seg_rds = Channel.empty()
        .mix(CBS.out.cbs_seg_rds)
        .mix(cbs_existing_seg)
    cbs_seg_for_merge = cbs_seg_rds
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, cbs_seg

    cbs_nseg_rds = Channel.empty()
        .mix(CBS.out.cbs_nseg_rds)
        .mix(cbs_existing_nseg)
    cbs_nseg_for_merge = cbs_nseg_rds
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, cbs_nseg

    // ASCAT
    // ##############################

    if (!params.tumor_only) {
        ascat_inputs = inputs.filter { it.ploidy.isEmpty() }.map { it -> [it.meta.patient, it.meta] }
        ascat_inputs_hets = hets_sites_for_merge
            .join(ascat_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, hets
        ascat_inputs_cov = dryclean_tumor_cov_for_merge
            .join(ascat_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, cov

        ascat_existing_outputs = inputs.map { it -> [it.meta, it.ploidy] }.filter { !it[1].isEmpty() }

        ascat_input = ascat_inputs
            .join(ascat_inputs_hets)
            .join(ascat_inputs_cov)
            .map { patient, meta, hets, cov ->
                [
                    meta,
                    hets,
                    cov,
                ]
            }

        COV_ASCAT(ascat_input)

        versions = versions.mix(COV_ASCAT.out.versions)
        ploidy = Channel.empty()
            .mix(COV_ASCAT.out.ploidy)
            .mix(ascat_existing_outputs)

        ploidy_for_merge = ploidy
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, ploidy
    }

    // JaBbA
    // ##############################
    jabba_inputs = inputs.filter { it.jabba_rds.isEmpty() }.map { it -> [it.meta.patient, it.meta] }
    jabba_inputs_sv = vcf_from_sv_calling_for_merge
        .join(jabba_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, vcf

    jabba_inputs_unfiltered_sv = unfiltered_som_sv_for_merge
        .join(jabba_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, unfiltered_som_sv

    if (params.tumor_only) {
        jabba_inputs_junction_filtered_sv = final_filtered_sv_rds_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, pon and gnomaD filtered sv (for tumor only)
    }

    jabba_inputs_hets = hets_sites_for_merge
        .join(jabba_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, hets

    jabba_inputs_cov = dryclean_tumor_cov_for_merge
        .join(jabba_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, cov

    if (!params.tumor_only) {
        jabba_inputs_ploidy = ploidy_for_merge
            .join(jabba_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, ploidy
    } else {
        jabba_inputs_ploidy = inputs.map { it -> [it.meta.patient, it.ploidy] }.filter { !it[1].isEmpty() }
    }

    jabba_inputs_cbs_seg = cbs_seg_for_merge
        .join(jabba_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, cbs_seg

    jabba_inputs_cbs_nseg = cbs_nseg_for_merge
        .join(jabba_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, cbs_nseg

    jabba_existing_outputs = inputs.map { it -> [it.meta, it.jabba_rds] }.filter { !it[1].isEmpty() }

    if (params.tumor_only) {
        jabba_inputs = jabba_inputs
            .join(jabba_inputs_cov)
            .join(jabba_inputs_ploidy)
            .join(jabba_inputs_junction_filtered_sv)
            .map{ patient, meta, cov, ploidy, junction ->
                [
                    meta,
                    junction,
                    cov,
                    [],
                    [],
                    ploidy,
                    [],
                    []
                ]
            }
    } else {
        // join all previous outputs to be used as input for jabba
        jabba_inputs = jabba_inputs
            .join(jabba_inputs_cov)
            .join(jabba_inputs_hets)
            .join(jabba_inputs_ploidy)
            .join(jabba_inputs_cbs_seg)
            .join(jabba_inputs_cbs_nseg)
            .join(jabba_inputs_sv)
            .join(jabba_inputs_unfiltered_sv)
            .map{ patient, meta, cov, hets, ploidy, seg, nseg, junction, j_supp ->
                [
                    meta,
                    junction,
                    cov,
                    j_supp,
                    hets,
                    ploidy,
                    seg,
                    nseg
                ]
            }
    }

    JABBA(jabba_inputs)

    jabba_rds = Channel.empty()
        .mix(JABBA.out.jabba_rds)
        .mix(jabba_existing_outputs)
    jabba_rds_for_merge = jabba_rds
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, jabba rds
    versions = versions.mix(JABBA.out.versions)

    // Non-integer balance
    // ##############################
    non_integer_balance_inputs = inputs.filter { it.ni_balanced_gg.isEmpty() }.map { it -> [it.meta.patient, it.meta] }

    non_integer_balance_inputs_jabba_rds = jabba_rds_for_merge
        .join(non_integer_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, jabba ggraph

    non_integer_balance_inputs_hets = hets_sites_for_merge
        .join(non_integer_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, hets

    non_integer_balance_inputs_cov = dryclean_tumor_cov_for_merge
        .join(non_integer_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, cov

    non_integer_balance_existing_outputs = inputs.map { it -> [it.meta, it.ni_balanced_gg] }.filter { !it[1].isEmpty() }

    non_integer_balance_inputs = non_integer_balance_inputs
        .join(non_integer_balance_inputs_jabba_rds)
        .join(non_integer_balance_inputs_hets)
        .join(non_integer_balance_inputs_cov)
        .map{ patient, meta, rds, hets, cov -> [ meta, rds, cov, hets ] }

    NON_INTEGER_BALANCE(non_integer_balance_inputs, bwa)
    versions = Channel.empty().mix(NON_INTEGER_BALANCE.out.versions)

    non_integer_balance_balanced_gg = Channel.empty()
        .mix(NON_INTEGER_BALANCE.out.non_integer_balance_balanced_gg)
        .mix(non_integer_balance_existing_outputs)
    non_integer_balance_balanced_gg_for_merge = non_integer_balance_balanced_gg
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, non integer balanced ggraph

    // LP Phased Balance
    // ##############################
    lp_phased_balance_inputs = inputs.filter { it.lp_balanced_gg.isEmpty() }.map { it -> [it.meta.patient, it.meta] }
    lp_phased_balance_inputs_ni_balanced_gg = non_integer_balance_balanced_gg_for_merge
        .join(lp_phased_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, non integer balanced ggraph
    lp_phased_balance_inputs_hets = hets_sites_for_merge
        .join(lp_phased_balance_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, hets

    lp_phased_balance_inputs = lp_phased_balance_inputs
        .join(lp_phased_balance_inputs_ni_balanced_gg)
        .join(lp_phased_balance_inputs_hets)
        .map{ patient, meta, balanced_gg, hets -> [ meta, balanced_gg, hets ] }

    lp_phased_balance_existing_outputs = inputs.map { it -> [it.meta, it.lp_balanced_gg] }.filter { !it[1].isEmpty() }

    LP_PHASED_BALANCE(lp_phased_balance_inputs)

    lp_phased_balance_balanced_gg = Channel.empty()
        .mix(LP_PHASED_BALANCE.out.lp_phased_balance_balanced_gg)
        .mix(lp_phased_balance_existing_outputs)

    // Events
    // ##############################
    events_inputs = inputs.filter { it.events.isEmpty() }.map { it -> [it.meta.patient, it.meta] }
    events_input_jabba_rds = jabba_rds_for_merge
        .join(events_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, jabba ggraph
    events_input_non_integer_balance = non_integer_balance_balanced_gg_for_merge
        .join(events_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, balanced_gg

    events_existing_outputs = inputs.map { it -> [it.meta, it.events] }.filter { !it[1].isEmpty() }

    events_input = events_inputs
        .join(events_input_jabba_rds)
        .join(events_input_non_integer_balance)
        .map{ patient, meta, rds, balanced_gg -> [ meta, rds, balanced_gg ] }

    EVENTS(events_input)

    events_versions = Channel.empty().mix(EVENTS.out.versions)
    events = Channel.empty()
        .mix(EVENTS.out.events_output)
        .mix(events_existing_outputs)

    // Fusions
    // ##############################
    fusions_inputs = inputs.filter { it.fusions.isEmpty() }.map { it -> [it.meta.patient, it.meta] }
    fusions_input_jabba_rds = jabba_rds_for_merge
        .join(fusions_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, jabba ggraph
    fusions_input_non_integer_balance = non_integer_balance_balanced_gg_for_merge
        .join(fusions_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, balanced_gg

    fusions_existing_outputs = inputs.map { it -> [it.meta, it.fusions] }.filter { !it[1].isEmpty() }

    fusions_input = fusions_inputs
        .join(fusions_input_jabba_rds)
        .join(fusions_input_non_integer_balance)
        .map{ patient, meta, rds, balanced_gg -> [ meta, rds, balanced_gg ] }

    FUSIONS(fusions_input)
    fusions = Channel.empty()
        .mix(FUSIONS.out.fusions_output)
        .mix(fusions_existing_outputs)
    versions = Channel.empty().mix(FUSIONS.out.versions)

    // SNV Calling
    // ##############################

    // Filter out bams for which SNV calling has already been done
    if (params.tumor_only) {
        bam_snv_inputs = inputs
            .filter { it.snv_somatic_vcf.isEmpty() }
            .map { it -> [it.meta.id] }
    } else {
        bam_snv_inputs = inputs
            .filter { it.snv_somatic_vcf.isEmpty() || it.snv_germline_vcf.isEmpty() }
            .map { it -> [it.meta.id] }
    }

    bam_snv_calling = alignment_bams_final
        .join(bam_snv_inputs)
        .map { it -> [ it[1].patient, it[1], it[2], it[3] ] } // meta.patient, meta, bam, bai

    snv_somatic_existing_outputs = inputs
        .map { it -> [it.meta, it.snv_somatic_vcf, it.snv_somatic_tbi] }
        .filter { !it[1].isEmpty() && !it[2].isEmpty()}

    if (!params.tumor_only) {
        snv_germline_existing_outputs = inputs
            .map { it -> [it.meta, it.snv_germline_vcf, it.snv_germline_tbi] }
            .filter { !it[1].isEmpty() && !it[2].isEmpty()}
    }

    if (params.tumor_only) {
        bam_snv_calling_status = bam_snv_calling.branch{
            tumor:  it[1].status == 1
        }

        bam_snv_calling_pair = bam_snv_calling_status.tumor.map{ meta, bam, bai -> [ meta, [], [], bam, bai ] }
    } else {
        // getting the tumor and normal cram files separated
        bam_snv_calling_status = bam_snv_calling.branch{
            normal: it[1].status == 0
            tumor:  it[1].status == 1
        }


        // Crossing the normal and tumor samples to create tumor and normal pairs
        bam_snv_calling_pair = bam_snv_calling_status.normal.cross(bam_snv_calling_status.tumor)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
        }
    }

    dict_sage = dict.map{ id, dict -> dict }

    BAM_SAGE(
        bam_snv_calling_pair,
        dict_sage
    )

    versions = versions.mix(BAM_SAGE.out.versions)

    filtered_somatic_vcf = Channel.empty()
        .mix(BAM_SAGE.out.sage_pass_filtered_somatic_vcf)
        .mix(snv_somatic_existing_outputs)
    filtered_somatic_vcf_for_merge = filtered_somatic_vcf
        .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, filtered somatic snv vcf, tbi

    if (!params.tumor_only) {
        germline_vcf = Channel.empty()
            .mix(BAM_SAGE.out.sage_germline_vcf)
            .mix(snv_germline_existing_outputs)
        germline_vcf_for_merge = germline_vcf
            .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, germline snv vcf, tbi
    }

    if (params.tumor_only) {
        BAM_SAGE_TUMOR_ONLY_FILTER(
            filtered_somatic_vcf,
            dbsnp_tbi,
            known_indels_tbi
        )

        filtered_somatic_vcf = Channel.empty()
            .mix(BAM_SAGE_TUMOR_ONLY_FILTER.out.sage_filtered_vcf)
        filtered_somatic_vcf_for_merge = filtered_somatic_vcf
            .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, filtered somatic snv vcf, tbi
    }

    // Variant Annotation
    // ##############################
    variant_somatic_ann_inputs = inputs
        .filter { it.variant_somatic_ann.isEmpty() || it.variant_somatic_bcf.isEmpty() }
        .map { it -> [it.meta.patient, it.meta] }

    variant_ann_input_somatic = variant_somatic_ann_inputs
        .join(filtered_somatic_vcf_for_merge)
        .map { it -> [ it[1], it[2], it[3] ] } // meta, filtered (pass and tumor-only filter) snvs, tbi


    variant_somatic_ann_existing_outputs = inputs.map { it -> [it.meta, it.variant_somatic_ann] }.filter { !it[1].isEmpty() }
    variant_somatic_bcf_existing_outputs = inputs.map { it -> [it.meta, it.variant_somatic_bcf] }.filter { !it[1].isEmpty() }

    VCF_SNPEFF_SOMATIC(
        variant_ann_input_somatic,
        snpeff_db_full,
        snpeff_cache
    )

    snv_somatic_annotations = Channel.empty()
        .mix(VCF_SNPEFF_SOMATIC.out.snpeff_vcf)
        .mix(variant_somatic_ann_existing_outputs)
    snv_somatic_annotations_for_merge = snv_somatic_annotations
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, annotated somatic snv vcf

    snv_somatic_bcf_annotations = Channel.empty()
        .mix(VCF_SNPEFF_SOMATIC.out.snpeff_bcf)
        .mix(variant_somatic_bcf_existing_outputs)

    if (!params.tumor_only) {
        variant_germline_ann_inputs = inputs
            .filter { it.variant_germline_ann.isEmpty() || it.variant_germline_bcf.isEmpty() }
            .map { it -> [it.meta.patient, it.meta] }

        variant_ann_input_germline = variant_germline_ann_inputs
            .join(germline_vcf_for_merge)
            .map { it -> [ it[1], it[2], it[3] ] } // meta, germline snvs, tbi

        variant_germline_ann_existing_outputs = inputs.map { it -> [it.meta, it.variant_germline_ann] }.filter { !it[1].isEmpty() }
        variant_germline_bcf_existing_outputs = inputs.map { it -> [it.meta, it.variant_germline_bcf] }.filter { !it[1].isEmpty() }

        VCF_SNPEFF_GERMLINE(
            variant_ann_input_germline,
            snpeff_db_full,
            snpeff_cache
        )

        snv_germline_annotations = Channel.empty()
            .mix(VCF_SNPEFF_GERMLINE.out.snpeff_vcf)
            .mix(variant_germline_ann_existing_outputs)
        snv_germline_annotations_for_merge = snv_germline_annotations
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, annotated germline snv vcf

        snv_germline_bcf_annotations = Channel.empty()
            .mix(VCF_SNPEFF_GERMLINE.out.snpeff_bcf)
            .mix(variant_germline_bcf_existing_outputs)
    }
    // SNV Multiplicity
    // ##############################
    if (!params.tumor_only) {
        snv_multiplicity_inputs = inputs.filter { it.snv_multiplicity.isEmpty() }.map { it -> [it.meta.patient, it.meta] }

        snv_multiplicity_inputs_somatic_vcf = snv_somatic_annotations_for_merge
            .join(snv_multiplicity_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, annotated somatic snv vcf
        snv_multiplicity_inputs_germline_vcf = snv_germline_annotations_for_merge
            .join(snv_multiplicity_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, annotated germline snv vcf
        snv_multiplicity_inputs_jabba_rds = jabba_rds_for_merge
            .join(snv_multiplicity_inputs)
            .map { it -> [ it[0], it[1] ] } // meta.patient, jabba ggraph

        snv_multiplicity_existing_outputs = inputs.map { it -> [it.meta, it.snv_multiplicity] }.filter { !it[1].isEmpty() }

        input_snv_multiplicity = snv_multiplicity_inputs
            .join(snv_multiplicity_inputs_somatic_vcf)
            .join(snv_multiplicity_inputs_germline_vcf)
            .join(snv_multiplicity_inputs_jabba_rds)
            .map{
                patient, meta, somatic_ann, germline_ann, ggraph ->
                [ meta, somatic_ann, germline_ann, ggraph ]
            }

        VCF_SNV_MULTIPLICITY(input_snv_multiplicity)

        snv_multiplicity = Channel.empty()
            .mix(VCF_SNV_MULTIPLICITY.out.snv_multiplicity_rds)
            .mix(snv_multiplicity_existing_outputs)
    }

    // Signatures
    // ##############################
    signatures_inputs = inputs
        .filter { it.sbs_signatures.isEmpty() || it.indel_signatures.isEmpty() || it.signatures_matrix.isEmpty()}
        .map { it -> [it.meta.patient, it.meta] }

    signatures_inputs_somatic_vcf = signatures_inputs
        .join(filtered_somatic_vcf_for_merge)
        .map { it -> [ it[1], it[2], it[3] ] } // meta, somatic snv, tbi

    sbs_signatures_existing_outputs = inputs.map { it -> [it.meta, it.sbs_signatures] }.filter { !it[1].isEmpty() }
    indel_signatures_existing_outputs = inputs.map { it -> [it.meta, it.indel_signatures] }.filter { !it[1].isEmpty() }
    signatures_matrix_existing_outputs = inputs.map { it -> [it.meta, it.signatures_matrix] }.filter { !it[1].isEmpty() }

    VCF_SIGPROFILERASSIGNMENT(signatures_inputs_somatic_vcf)

    sbs_signatures = Channel.empty()
        .mix(VCF_SIGPROFILERASSIGNMENT.out.sbs_signatures)
        .mix(sbs_signatures_existing_outputs)
    indel_signatures = Channel.empty()
        .mix(VCF_SIGPROFILERASSIGNMENT.out.indel_signatures)
        .mix(indel_signatures_existing_outputs)
    signatures_matrix = Channel.empty()
        .mix(VCF_SIGPROFILERASSIGNMENT.out.signatures_matrix)
        .mix(signatures_matrix_existing_outputs)

    // HRDetect
    // ##############################
    hrdetect_inputs = inputs
        .filter { it.hrdetect.isEmpty() }
        .map { it -> [it.meta.patient, it.meta] }

    hrdetect_inputs_sv = vcf_from_sv_calling_for_merge
        .join(hrdetect_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, somatic/pon+gnomAD filtered sv

    hrdetect_inputs_hets = hets_sites_for_merge
        .join(hrdetect_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, hets

    hrdetect_inputs_vcf = filtered_somatic_vcf_for_merge
        .join(hrdetect_inputs)
        .map { it -> [ it[0], it[1], it[2] ] } // meta.patient, somatic snv, tbi

    hrdetect_inputs_jabba_rds = jabba_rds_for_merge
        .join(hrdetect_inputs)
        .map { it -> [ it[0], it[1] ] } // meta.patient, jabba ggraph

    hrdetect_existing_outputs = inputs.map { it -> [it.meta, it.hrdetect] }.filter { !it[1].isEmpty() }

    hrdetect_input = hrdetect_inputs
        .join(hrdetect_inputs_sv)
        .join(hrdetect_inputs_hets)
        .join(hrdetect_inputs_vcf)
        .join(hrdetect_inputs_jabba_rds)
        .map{
            patient,
            meta,
            sv,
            hets,
            snv,
            snv_tbi,
            jabba ->[ meta, sv, hets, snv, snv_tbi, jabba ]
        }

        JUNC_SNV_GGRAPH_HRDETECT(hrdetect_input)

        versions = versions.mix(JUNC_SNV_GGRAPH_HRDETECT.out.versions)
        hrdetect_rds = Channel.empty()
            .mix(JUNC_SNV_GGRAPH_HRDETECT.out.hrdetect_rds)
            .mix(hrdetect_existing_outputs)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
}
