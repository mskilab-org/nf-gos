include { paramsSummaryLog; paramsSummaryMap; samplesheetToList } from 'plugin/nf-schema'

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

// Function to expand brace-enclosed parts of a path
def expandBraces(String path) {
    def regex = /\{([^}]+)\}/
    def matcher = path =~ regex
    if (!matcher) {
        return [path] // No braces to expand
    }

    def expandedPaths = [path]
    matcher.each { match ->
        def options = match[1].split(',')
        def newPaths = []
        expandedPaths.each { expandedPath ->
            options.each { option ->
                newPaths << expandedPath.replaceFirst(regex, option)
            }
        }
        expandedPaths = newPaths
    }
    return expandedPaths
}


// def expandBraces(String path) {
//     def regex = /\{([^}]+)\}/
//     def matcher = path =~ regex
//     if (!matcher) {
//         return [path] // No braces to expand
//     }

//     def expandedPaths = [path]
//     matcher.each { match ->
//         def options = match[1].split(',')
//         def newPaths = []
//         expandedPaths.each { expandedPath ->
//             options.each { option ->
//                 newPaths << expandedPath.replaceFirst(regex, option)
//             }
//         }
//         expandedPaths = newPaths
//     }
//     return expandedPaths
// }

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

include { validateParameters; paramsHelp } from 'plugin/nf-schema'

// Download annotation cache if needed
include { PREPARE_CACHE } from '../subworkflows/local/prepare_cache/main'

// Build indices if needed
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome/main'

// Build intervals if needed
include { PREPARE_INTERVALS } from '../subworkflows/local/prepare_intervals/main'

// BAM Picard QC
// include { BAM_QC } from '../subworkflows/local/bam_qc/main'

//GRIDSS
include { 
    BAM_SVCALLING_GRIDSS; 
    BAM_SVCALLING_GRIDSS_SOMATIC
} from '../subworkflows/local/bam_svcalling_gridss/main'

// SV Junction Filtering
include { 
    SV_JUNCTION_FILTER as JUNCTION_FILTER;
    SV_JUNCTION_FILTER_BEDTOOLS as JUNCTION_FILTER_BEDTOOLS
} from '../subworkflows/local/junction_filter/main'

include { 
    ALIGNMENT_STEP; 
    BAM_QC;
    MSISENSORPRO_STEP;
    SV_CALLING_STEP;
    FRAGCOUNTER_STEP;
    DRYCLEAN_STEP;
    AMBER_STEP;
    COBALT_STEP;
    PURPLE_STEP;
    CBS_STEP;
    VARIANT_CALLING_STEP;
    VARIANT_ANNOTATION_STEP;
    JABBA_STEP;
    NON_INTEGER_BALANCE_STEP;
    LP_PHASED_BALANCE_STEP;
    EVENTS_STEP;
    FUSIONS_STEP;
    MULTIPLICITY_STEP;
    ONCOKB_STEP;
    SIGNATURES_STEP;
    HRDETECT_STEP;
    ONENESS_TWONESS_STEP
} from '../subworkflows/local/steps.nf'

include {
    SV_CHIMERA_FILTER as SV_CHIMERA_FILTER_RAWVCF;
    SV_CHIMERA_FILTER as SV_CHIMERA_FILTER_VCF
} from '../modules/local/process.nf'


workflow SETUP {

    main:
    logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    summary_params = paramsSummaryMap(workflow)

    // Print parameter summary log to screen
    log.info logo + paramsSummaryLog(workflow) + citation

    //Check input path parameters to see if they exist
    checkPathParamList = [
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

    toolParamMap = [
        "msisensorpro": [
            params.msisensorpro_list
        ],
        "gridss": [
            params.blacklist_gridss,
            params.pon_gridss
        ],
        "hetpileups" : [
            params.hapmap_sites
        ],
        "fragcounter": [
            params.gcmapdir_frag
        ],
        "dryclean": [
            params.pon_dryclean,
        ],
        "fusions"    : [
            params.gencode_fusions
        ],
        "non_integer_balance" : [
            params.mask_non_integer_balance
        ],
        "lp_phased_balance" : [
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
        "cobalt"    : [
            params.gc_profile,
            params.diploid_bed
        ],
        "purple"    : [
            params.het_sites_amber,
            params.gc_profile,
        ],
    ]

    toolParamMap.each { tool, params ->
        params.each { param ->
            if (param) {
                checkPathParamList.add(param)
            }
        }
    }


    if ((params.download_cache) && (params.snpeff_cache || params.vep_cache)) {
        error("Please specify either `--download_cache` or `--snpeff_cache`, `--vep_cache`.")
    }

    def awsCliInstalled = isAwsCliInstalled()
    if (!awsCliInstalled) {
        println "AWS CLI is not installed/loaded. Will proceed, but S3 paths will not be checked."
    }

    println "Checking if parameter paths exist..."

    println "Checking if parameter paths exist..."
    int numThreads = 16 // Specify the number of threads to use in the pool

    def GParsPool = Class.forName('groovyx.gpars.GParsPool')
    GParsPool.withPool(numThreads) {
        checkPathParamList.eachParallel { param ->
            if (param == null) {
                println "Skipping null path"
                return
            }

            def expandedPaths = expandBraces(param)
            expandedPaths.each { expandedPath ->
                if (expandedPath.startsWith("s3://")) {
                    if (awsCliInstalled) {
                        def s3ExpandedPaths = expandS3Paths(expandedPath)
                        s3ExpandedPaths.each { s3Path ->
                            if (checkS3PathExists(s3Path)) {
                                println "Path exists: ${s3Path}"
                            } else {
                                println "Path does not exist: ${s3Path}"
                            }
                        }
                    } else {
                        println "Skipping S3 path check for: ${expandedPath}"
                    }
                } else {
                    def file = new File(expandedPath)
                    if (file.exists()) {
                        println "Path exists: ${expandedPath}"
                    } else {
                        println "Path does not exist: ${expandedPath}"
                    }
                }
            }
        }
    }

    inputType = params.input ? "input" : "input_restart"

    if (params.build_only_index) {
        ch_from_samplesheet = Channel.empty()
        samplesheetList = []
    } else {
        samplesheetList = samplesheetToList(params.get(inputType), "gos-assets/nf-gos/assets/schema_input.json")
        ch_from_samplesheet = Channel.fromList(samplesheetList)
    }
    // def ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet(inputType)

    rowsAsMaps = WorkflowNfcasereports.samplesheetTuplesToMaps(
        samplesheetList, 
        file("${projectDir}/gos-assets/nf-gos/assets/schema_input.json").toFile()
    )

    rowsAsMaps = rowsAsMaps.collect { it ->
        it + [
            crai: it.cram ? it.cram + '.crai' : [],
            bai: it.bam ? it.bam + '.bai': [],
            vcf_tbi: it.vcf ? it.vcf + '.tbi' : [],
            vcf_raw_tbi: it.vcf_raw ? it.vcf_raw + '.tbi' : [],
            snv_somatic_vcf_tumoronly_filtered_tbi: it.snv_somatic_vcf_tumoronly_filtered ? it.snv_somatic_vcf_tumoronly_filtered + ".tbi" : [],
            snv_somatic_vcf_rescue_ch_heme_tbi: it.snv_somatic_vcf_rescue_ch_heme ? it.snv_somatic_vcf_rescue_ch_heme + '.tbi' : [],
            snv_somatic_tbi: it.snv_somatic_vcf ? it.snv_somatic_vcf + '.tbi' : [],
            snv_germline_tbi: it.snv_germline_vcf ? it.snv_germline_vcf + '.tbi' : [],
            ffpe_impact_vcf_tbi: it.ffpe_impact_vcf ? it.ffpe_impact_vcf + '.tbi' : [],
            ffpe_impact_filtered_vcf_tbi: it.ffpe_impact_filtered_vcf ? it.ffpe_impact_filtered_vcf + '.tbi' : []
        ]
    }
    
    println "Settings Globals.rowsAsMaps to:"
    rowsAsMaps.eachWithIndex { m, i -> log.info "Row ${i}: ${m}" }
    Globals.rowsAsMaps = rowsAsMaps

    inputs = Channel.fromList(rowsAsMaps)

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
                ch_items.meta =  ch_items.meta + [id: "${ch_items.meta.sample}-${ch_items.meta.lane}".toString()]
                def CN =  params.seq_center ? "CN:${params.seq_center}\\t" : ''

                def flowcell =  flowcellLaneFromFastq(ch_items.fastq_1)
                // Don't use a random element for ID, it breaks resuming
                def read_group = "\"@RG\\tID:${flowcell}.${ch_items.meta.sample}.${ch_items.meta.lane}\\t${CN}PU:${ch_items.meta.lane}\\tSM:${ch_items.meta.patient}_${ch_items.meta.sample}\\tLB:${ch_items.meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

                ch_items.meta =  ch_items.meta - ch_items.meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), size: 1]

            } else if (ch_items.fastq_2) {
                ch_items.meta =  ch_items.meta + [id: ch_items.meta.sample.toString()]
                def CN =  params.seq_center ? "CN:${params.seq_center}\\t" : ''

                def flowcell =  flowcellLaneFromFastq(ch_items.fastq_1)
                def read_group = "\"@RG\\tID:${flowcell}.${ch_items.meta.sample}\\t${CN}PU:${ch_items.meta.sample}\\tSM:${ch_items.meta.patient}_${ch_items.meta.sample}\\tLB:${ch_items.meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

                ch_items.meta = ch_items.meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), size: 1]
            } else if (ch_items.meta.lane && ch_items.bam) {
                ch_items.meta =  ch_items.meta + [id: "${ch_items.meta.sample}-${ch_items.meta.lane}".toString()]
                def CN =  params.seq_center ? "CN:${params.seq_center}\\t" : ''
                def read_group  = "\"@RG\\tID:${ch_items.meta.sample}_${ch_items.meta.lane}\\t${CN}PU:${ch_items.meta.lane}\\tSM:${ch_items.meta.patient}_${ch_items.meta.sample}\\tLB:${ch_items.meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

                ch_items.meta = ch_items.meta - ch_items.meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), size: 1]
            } else {
                ch_items.meta = ch_items.meta + [id: ch_items.meta.sample.toString()]
            }

            ch_items
        }

    inputs_unlaned = inputs.map { it ->
        it + [meta: Utils.remove_lanes_from_meta(it.meta)]
    }

    emit:
    inputs
    inputs_unlaned

}

workflow TOOLS {

    main:
    tool_input_output_map = Globals.tool_input_output_map
    
    // see lib/Globals.groovy
    sampleList = Globals.rowsAsMaps
    available_inputs = new HashSet()
    present_outputs = new HashSet()
    

    sampleList.each { input_map ->
        input_map.each { key, value ->
            def is_value_present = value && !(value instanceof Collection && value.empty)
            if (is_value_present) {
                available_inputs.add(key)
            }
        }
    }

    println "Provided inputs: ${available_inputs}"

    schemaFile = file("$projectDir/gos-assets/nf-gos/assets/schema_input.json")
    schema = new groovy.json.JsonSlurper().parse(schemaFile)


    props = schema.items.properties
    requiredFields = props.findAll { !it.value.containsKey('meta') }.keySet()
    println "requiredFields: $requiredFields"

    missing_outputs = requiredFields.findAll { field ->
        // Check if this field is missing (null or empty collection) in any sample
        sampleList.any { sample ->
            def value = sample[field]
            // !value || (value instanceof Collection && value.empty) || value.toString() == ""
            def truthy_val = (
                (!value) ||
                (value == null) ||
                (value instanceof Collection && value.empty) ||
                (value.toString().trim() == "") ||
                (value instanceof String && value.replaceAll(/["']/, "").trim() == "")
            )
            truthy_val
        }
    }
    println "Outputs MISSING from at least one sample: $missing_outputs"

    // Iteratively select tools based on available inputs
    skip_tools = params.skip_tools ? params.skip_tools.split(',').collect { it.trim() } : []
    println "Skipping tools: ${skip_tools}"
    // TODO: if GRIDSS - skip if vcf is found, but not if vcf_raw is present.
    selected_tools = []
    tools_qc = ["collect_wgs_metrics", "collect_multiple_metrics", "estimate_library_complexity"]
    selected_tools_map = [:]
    tool_input_output_map.each { tool, io ->
        if (!selected_tools.contains(tool) && !skip_tools.contains(tool)) {

            def inputsRequired = io.inputs
            def inputsPresent = inputsRequired.every { available_inputs.contains(it) }
            def outputsNeeded = io.outputs.any { missing_outputs.contains(it) }

            // special cases
            def is_sage_tumor_only = tool == "sage" && params.tumor_only
            def is_sage_heme = is_sage_tumor_only && params.is_heme
            def is_current_tool_qc = tools_qc.contains(tool)
            def is_current_tool_qc_multiple_metrics = tool == "collect_multiple_metrics" // nested
            // TODO: for later
            def is_current_tool_jabba = tool == "jabba" // separate first input, (vcf or sv_retier) and remaining required inputs
            def is_output_nested_list = io.outputs instanceof List && io.outputs.every { it instanceof List }
            def is_output_generic_case = !is_sage_tumor_only && !is_current_tool_qc_multiple_metrics
            def is_input_generic_case = !is_current_tool_jabba

            // Treat special cases
            if (is_sage_tumor_only) {
                outputsNeeded = ["snv_somatic_vcf", "snv_somatic_vcf_tumoronly_filtered"].any {
                    missing_outputs.contains(it)
                }
            }

            if (is_sage_heme) {
                outputsNeeded = ["snv_somatic_vcf", "snv_somatic_vcf_tumoronly_filtered", "snv_somatic_vcf_rescue_ch_heme"].any {
                    missing_outputs.contains(it)
                }
            }

            if (is_current_tool_qc_multiple_metrics) {
                def is_any_alignment_summary_absent = io.outputs[0].any {
                    missing_outputs.contains(it)
                }
                def is_any_insert_size_absent = io.outputs[1].any {
                    missing_outputs.contains(it)
                }
                outputsNeeded = is_any_alignment_summary_absent || is_any_insert_size_absent
            }

            if (inputsPresent && outputsNeeded) {
                selected_tools.add(tool)
                available_inputs.addAll(io.outputs)
            }




            // selected_tools_map[tool] = inputs.filter { sample ->
            // 	def is_all_input_col_present = false
            // 	def is_any_output_col_empty = false

            // 	is_all_input_col_present = io.inputs.every { field ->
            // 		test_robust_presence(sample[field], test_file = false) // Tests if file exists and is nonzero file size
            // 	}

            // 	// if (is_input_generic_case) {
            // 	// 	is_all_input_col_present = io.inputs.every { field ->
            // 	// 		! sample[field].isEmpty() // Tests if file exists and is nonzero file size
            // 	// 	}
            // 	// }
            // 	// if (is_current_tool_jabba) {
            // 	// 	is_any_sv_input_col_present = io.inputs[0].any { field ->
            // 	// 		! sample[field].isEmpty() // Tests if file exists and is nonzero file size
            // 	// 	}
            // 	// 	is_all_remaining_input_col_present = io.inputs[1].every { field ->
            // 	// 		! sample[field].isEmpty() // Tests if file exists and is nonzero file size
            // 	// 	}
            // 	// 	is_all_input_col_present = is_any_sv_input_col_present && is_all_remaining_input_col_present
            // 	// }



            // 	// Generic output case
            // 	if (is_output_generic_case) {
            // 		is_any_output_col_empty = io.outputs.any { field ->
            // 			test_robust_absence(sample[field], test_file = false)
            // 		}
            // 	}


            // 	def is_any_output_empty_test = io.outputs.any { field ->
            // 		value = sample[field]
            // 		def is_absent = test_robust_absence(sample[field], test_file = false)
            // 		is_absent
            // 	}

            // 	def is_all_input_present_test = io.inputs.every { field ->
            // 		value = sample[field]
            // 		def is_present = test_robust_presence(sample[field], test_file = false)
            // 		is_present
            // 	}

            // 	// Treat special cases
            // 	if (is_sage_tumor_only) {
            // 		is_any_output_col_empty = test_robust_absence(sample["snv_somatic_vcf"], test_file = false)
            // 	}

            // 	if (is_current_tool_qc_multiple_metrics) {
            // 		is_any_alignment_summary_absent = io.outputs[0].any { field ->
            // 			test_robust_absence(sample[field], test_file = false) // Tests if file exists and is nonzero file size
            // 		}
            // 		def is_any_insert_size_absent = io.outputs[1].any { field ->
            // 			test_robust_absence(sample[field], test_file = false) // Tests if file exists and is nonzero file size
            // 		}
            // 		is_any_output_col_empty = is_any_alignment_summary_absent || is_any_insert_size_absent
            // 	}
            // 	// End Treat special cases

            // 	return is_any_output_col_empty && is_all_input_col_present
            // 	// return [tool, is_any_output_col_empty && is_all_input_col_present, sample]
            // }


        }
    }

    // see lib/Globals.groovy
    println "Setting Globals.selected_tools_map to: ${selected_tools_map}"
    Globals.selected_tools_map = selected_tools_map


    tools_used = selected_tools

    tools_used.removeAll(skip_tools)
    // see lib/Globals.groovy
    println "Setting Globals.tools_used to: ${tools_used}"
    Globals.tools_used = tools_used

    println "Tools that will be run based on your inputs: ${tools_used}"

    if (!params.dbsnp && !params.known_indels) {
        if (!params.skip_tools || (params.skip_tools && !params.skip_tools.contains('baserecalibrator'))) {
            log.warn "Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command."
        }
        if (params.skip_tools && (!params.skip_tools.contains('haplotypecaller') || !params.skip_tools.contains('sentieon_haplotyper'))) {
            log.warn "If GATK's Haplotypecaller or Sentieon's Haplotyper is specified, without `--dbsnp` or `--known_indels no filtering will be done. For filtering, please provide at least one of `--dbsnp` or `--known_indels`.\nFor more information see FilterVariantTranches (single-sample, default): https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches\nFor more information see VariantRecalibration (--joint_germline): https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator\nFor more information on GATK Best practice germline variant calling: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-"
        }
    }

    emit:
    tools_used
    selected_tools_map


}


workflow NFGOS {

    main:
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

    SETUP()
    inputs = SETUP.out.inputs
    inputs_unlaned = SETUP.out.inputs_unlaned
    rowsAsMaps = Globals.rowsAsMaps // Written inside SETUP()
    TOOLS()
    // See lib/Globals.groovy
    tools_used = Globals.tools_used // Written inside TOOLS()
    // See lib/Globals.groovy
    selected_tools_map = Globals.selected_tools_map

    // TOOLS(inputs, inputs_unlaned, rowsAsMaps)
    // tools_used = TOOLS.out.tools_used
    // selected_tools_map = TOOLS.out.selected_tools_map

    dbsnp = WorkflowNfcasereports.create_file_channel(params.dbsnp)
    fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
    fasta_fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
    germline_resource = WorkflowNfcasereports.create_file_channel(params.germline_resource)
    known_indels = WorkflowNfcasereports.create_file_channel(params.known_indels)
    known_snps = WorkflowNfcasereports.create_file_channel(params.known_snps)
    pon = WorkflowNfcasereports.create_file_channel(params.pon)

    // Initialize value channels based on params, defined in the params.genomes[params.genome] scope

    snpeff_genome = WorkflowNfcasereports.create_value_channel(params.snpeff_genome)
    snpeff_db = WorkflowNfcasereports.create_value_channel(params.snpeff_db)
    snpeff_db_full = params.snpeff_db && params.snpeff_genome   ? Channel.value("${params.snpeff_genome}.${params.snpeff_db}") : Channel.empty()
    vep_cache_version = WorkflowNfcasereports.create_value_channel(params.vep_cache_version)
    vep_genome = WorkflowNfcasereports.create_value_channel(params.vep_genome)
    vep_species = WorkflowNfcasereports.create_value_channel(params.vep_species)

    ch_multiqc_config = Channel.fromPath("$projectDir/gos-assets/nf-gos/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/gos-assets/nf-gos/assets/methods_description_template.yml", checkIfExists: true)

    // Set is_heme based on is_retier_whitelist_junctions
    params.is_heme = params.is_retier_whitelist_junctions

    // To gather all QC reports for MultiQC
    reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    versions = Channel.empty()

        // Download cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache    ? [] : Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    snpeff_info = params.snpeff_cache ? [] : Channel.of([ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], params.snpeff_genome, params.snpeff_db ])

    if (params.snpeff_cache) {
        def snpeff_annotation_cache_key = params.use_annotation_cache_keys ? "${params.snpeff_genome}.${params.snpeff_db}/" : ""
        def snpeff_cache_dir =  "${snpeff_annotation_cache_key}${params.snpeff_genome}.${params.snpeff_db}"
        def snpeff_cache_path_full = file("$params.snpeff_cache/$snpeff_cache_dir", type: 'dir')
        if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
            error("Files within --snpeff_cache invalid. Make sure there is a directory named ${snpeff_cache_dir} in ${params.snpeff_cache}.\nhttps://nf-co.re/sarek/dev/usage#how-to-customise-snpeff-and-vep-annotation")
        }
        snpeff_cache = Channel.fromPath(file("${params.snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
            .map{ cache -> [ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], cache ] }
    } else {
        snpeff_cache = []
    }



    if (params.download_cache) {
        PREPARE_CACHE(ensemblvep_info, snpeff_info)
        snpeff_cache = PREPARE_CACHE.out.snpeff_cache
        vep_cache = PREPARE_CACHE.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

        versions = versions.mix(PREPARE_CACHE.out.versions)
    }

    // Always build indices
    // ##############################
    PREPARE_GENOME()

    // Gather built indices or get them from the params
    // Built from the fasta file:
    dict = params.dict ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                    : PREPARE_GENOME.out.dict
    fasta_fai = WorkflowNfcasereports.create_file_channel(params.fasta_fai, PREPARE_GENOME.out.fasta_fai)
    bwa = WorkflowNfcasereports.create_file_channel(params.bwa)

    // Gather index for mapping given the chosen aligner
    index_alignment = bwa

    // TODO: add a params for msisensorpro_scan
    msisensorpro_scan = PREPARE_GENOME.out.msisensorpro_scan

    dbsnp_tbi =  WorkflowNfcasereports.create_index_channel(params.dbsnp, params.dbsnp_tbi, PREPARE_GENOME.out.dbsnp_tbi)
    //do not change to Channel.value([]), the check for its existence then fails for Getpileupsumamries
    germline_resource_tbi = params.germline_resource ? params.germline_resource_tbi ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : []
    known_indels_tbi = WorkflowNfcasereports.create_index_channel(params.known_indels, params.known_indels_tbi, PREPARE_GENOME.out.known_indels_tbi)
    known_snps_tbi = WorkflowNfcasereports.create_index_channel(params.known_snps, params.known_snps_tbi, PREPARE_GENOME.out.known_snps_tbi)
    pon_tbi = WorkflowNfcasereports.create_index_channel(params.pon, params.pon_tbi, PREPARE_GENOME.out.pon_tbi)

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    known_sites_snps = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi = dbsnp_tbi.concat(known_snps_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai, params.intervals, params.no_intervals)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined = WorkflowNfcasereports.create_file_channel(params.no_intervals, PREPARE_INTERVALS.out.intervals_bed_combined)
    intervals_bed_gz_tbi_combined = WorkflowNfcasereports.create_file_channel(params.no_intervals, PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined)

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = params.wes ?
        intervals_bed_combined.map{it -> [ [ id:it.baseName ], it ]}.collect() :
        Channel.value([ [ id:'null' ], [] ])

    intervals =  PREPARE_INTERVALS.out.intervals_bed // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
        if ( num_intervals < 1 ) [ [], num_intervals ]
        else [ interval, num_intervals ]
    }

    // Gather used softwares versions
    versions = versions.mix(PREPARE_GENOME.out.versions)
    versions = versions.mix(PREPARE_INTERVALS.out.versions)

    // inputs.view{ "inputs: $it"}

	ALIGNMENT_STEP(
		inputs,
		selected_tools_map,
		tools_used,
		known_sites_indels,
		known_sites_indels_tbi,
		index_alignment,
		fasta,
		fasta_fai,
		dict,
		intervals_for_preprocessing,
		intervals_and_num_intervals
	)

	alignment_bams_final = ALIGNMENT_STEP.out.alignment_bams_final.dump(tag: 'alignment_bams_final', pretty: true)


    dict_path = dict.map{ _meta, dictOut -> [dictOut] }
    BAM_QC(inputs, alignment_bams_final, dict_path, tools_used)

    // MSISensorPro
    // ##############################
    MSISENSORPRO_STEP(
        inputs_unlaned,
        alignment_bams_final,
        tools_used,
        msisensorpro_scan
    )
    msi_from_msisensorpro = MSISENSORPRO_STEP.out.msi_from_msisensorpro
    germline_msi_from_msisensorpro = MSISENSORPRO_STEP.out.germline_msi_from_msisensorpro

    // SV Calling
    // ##############################
    SV_CALLING_STEP(
        inputs_unlaned,
        alignment_bams_final,
        index_alignment,
        tools_used
    )

    vcf_from_gridss_gridss = SV_CALLING_STEP.out.vcf_from_gridss_gridss
    vcf_raw_from_gridss_gridss = SV_CALLING_STEP.out.vcf_raw_from_gridss_gridss

    do_filter_ffpe_chimera = params.filter_ffpe_chimera ?: false
    if (do_filter_ffpe_chimera) {
        SV_CHIMERA_FILTER_VCF(vcf_from_gridss_gridss)
        SV_CHIMERA_FILTER_RAWVCF(vcf_raw_from_gridss_gridss)
        vcf_from_gridss_gridss = SV_CHIMERA_FILTER_VCF.out.vcftbi
        vcf_raw_from_gridss_gridss = SV_CHIMERA_FILTER_RAWVCF.out.vcftbi
    }

    /* FIXME: Junction Filtering step
        Placed here for now as the new JUNCTION_FILTER_BEDTOOLS needs to be implemented.
    */
    final_filtered_sv_rds_for_merge = Channel.empty()
    unfiltered_som_sv_for_merge = Channel.empty()
    if (params.tumor_only) {
        vcf_from_sv_calling_for_merge = vcf_from_gridss_gridss
            .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, vcf, tbi

        vcf_from_gridss_gridss.dump(tag: "vcf_from_gridss_gridss", pretty: true)

        // JUNCTION_FILTER(vcf_from_gridss_gridss)
        JUNCTION_FILTER(vcf_raw_from_gridss_gridss)
        // JUNCTION_FILTER_BEDTOOLS(vcf_raw_from_gridss_gridss)

        pon_filtered_sv_rds = Channel.empty().mix(JUNCTION_FILTER.out.pon_filtered_sv_rds)
        final_filtered_sv_rds = Channel.empty().mix(JUNCTION_FILTER.out.final_filtered_sv_rds)
        final_filtered_sv_rds_for_merge = final_filtered_sv_rds
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, rds
            .mix(
                inputs_unlaned
                    .map { it -> [it.meta, it.vcf, it.vcf_tbi] }
                    .filter {
                        def vcf_or_rds = it[1]
                        def is_rds = vcf_or_rds =~ /\.rds$/
                        def is_vcf_or_rds_filesize_zero_or_nonexistent = vcf_or_rds.isEmpty()
                        ! is_vcf_or_rds_filesize_zero_or_nonexistent && is_rds }
                    .map {
                        it -> [ it[0].patient, it[1] ] } // meta.patient, vcf_or_rds
                    .unique()
            )

    } else {
        //somatic filter for GRIDSS
        BAM_SVCALLING_GRIDSS_SOMATIC(vcf_raw_from_gridss_gridss)

        versions = versions.mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.versions)
        vcf_somatic_high_conf = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_high_confidence)
        vcf_from_sv_calling_for_merge = vcf_somatic_high_conf
            .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, vcf, tbi

        unfiltered_som_sv = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_all)
        unfiltered_som_sv_for_merge = unfiltered_som_sv
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, vcf
    }


    AMBER_STEP(
        inputs_unlaned,
        alignment_bams_final,
        tools_used
    )

    amber_dir_for_merge = AMBER_STEP.out.amber_dir
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, amber_dir

    hets_sites_for_merge = AMBER_STEP.out.sites_from_het_pileups_wgs
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, hets
    
    COBALT_STEP(
        inputs_unlaned,
        alignment_bams_final,
        tools_used
    )

    cobalt_dir_for_merge = COBALT_STEP.out.cobalt_dir
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, cobalt_dir

    FRAGCOUNTER_STEP(
        inputs_unlaned,
        alignment_bams_final, // [meta, bam, bai]
        tools_used
    )
    tumor_frag_cov_for_merge = FRAGCOUNTER_STEP.out.tumor_frag_cov.map { meta, frag_cov -> [ meta.sample, meta, frag_cov ] }
    normal_frag_cov_for_merge = FRAGCOUNTER_STEP.out.normal_frag_cov.map { meta, frag_cov -> [ meta.sample, meta, frag_cov ] }

    DRYCLEAN_STEP(
        inputs_unlaned,
        tumor_frag_cov_for_merge,
        normal_frag_cov_for_merge,
        tools_used
    )

     sample_meta_map = inputs_unlaned.map { it -> [ it.meta.sample, it.meta ]}.unique()

     dryclean_tumor_cov_for_merge = DRYCLEAN_STEP.out.dryclean_tumor_cov
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, dryclean_cov
            .unique{ patient, _dryclean_cov -> patient }

    // Need to recombine normals here with full metadata if a normal is matched across multiple tumor/normal pairs.
     dryclean_normal_cov_for_merge = DRYCLEAN_STEP.out.dryclean_normal_cov
            .map { it -> [ it[0].sample, it[1] ] } // meta.sample, dryclean_cov
            .cross(sample_meta_map)
            .map { dryclean_normal, sample_meta ->
                def meta_complete = sample_meta[1]
                [ meta_complete.patient, dryclean_normal[1] ]
            }
            .unique{ patient, _dryclean_cov -> patient }
            .dump(tag: "dryclean_normal_cov_for_merge", pretty: true)


    CBS_STEP(
        inputs_unlaned,
        dryclean_tumor_cov_for_merge,
        dryclean_normal_cov_for_merge,
        tools_used
    )

    cbs_seg_for_merge = CBS_STEP.out.cbs_seg_rds
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, cbs_seg
    cbs_nseg_for_merge = CBS_STEP.out.cbs_nseg_rds
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, cbs_nseg

    VARIANT_CALLING_STEP(
        inputs_unlaned,
        alignment_bams_final,
        dbsnp_tbi,
        known_indels_tbi,
        tools_used
    )


    filtered_somatic_vcf_for_merge = VARIANT_CALLING_STEP.out.filtered_somatic_vcf
        .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, filtered somatic snv vcf, tbi
    
    germline_vcf_for_merge = VARIANT_CALLING_STEP.out.germline_vcf
                .map { it -> [ it[0].patient, it[1], it[2] ] } // meta.patient, germline snv vcf, tbi

    VARIANT_ANNOTATION_STEP(
        inputs_unlaned,
        filtered_somatic_vcf_for_merge,
        germline_vcf_for_merge,
        snpeff_db_full,
        snpeff_cache,
        tools_used
    )

    snv_somatic_annotations_for_merge = VARIANT_ANNOTATION_STEP.out.snv_somatic_annotations
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, annotated somatic snv vcf
    
    snv_germline_annotations_for_merge = VARIANT_ANNOTATION_STEP.out.snv_germline_annotations
                .map { it -> [ it[0].patient, it[1] ] } // meta.patient, annotated germline snv vcf


    PURPLE_STEP(
        inputs_unlaned,
        germline_vcf_for_merge,
        filtered_somatic_vcf_for_merge,
        cobalt_dir_for_merge,
        amber_dir_for_merge,
        vcf_from_sv_calling_for_merge,
        tools_used
    )
    
    purity_for_merge = PURPLE_STEP.out.purity
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, purity
    ploidy_for_merge = PURPLE_STEP.out.ploidy
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, ploidy

    JABBA_STEP(
        inputs_unlaned,
        dryclean_tumor_cov_for_merge,
        hets_sites_for_merge,
        cbs_seg_for_merge,
        cbs_nseg_for_merge,
        purity_for_merge,
        ploidy_for_merge,
        vcf_from_sv_calling_for_merge,
        vcf_raw_from_gridss_gridss,
        final_filtered_sv_rds_for_merge,
        unfiltered_som_sv_for_merge,
        tools_used
    )

    jabba_gg_for_merge = JABBA_STEP.out.jabba_gg
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, jabba.gg.rds
    jabba_rds_for_merge = JABBA_STEP.out.jabba_rds
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, jabba rds
    

    NON_INTEGER_BALANCE_STEP(
        inputs_unlaned,
        jabba_gg_for_merge,
        hets_sites_for_merge,
        dryclean_tumor_cov_for_merge,
        tools_used
    )
    non_integer_balance_balanced_gg_for_merge = NON_INTEGER_BALANCE_STEP.out.non_integer_balance_balanced_gg
        .map { it -> [ it[0].patient, it[1] ] } // meta.patient, non integer balanced ggraph
    
    LP_PHASED_BALANCE_STEP(
        inputs_unlaned,
        non_integer_balance_balanced_gg_for_merge,
        hets_sites_for_merge,
        tools_used
    )

    EVENTS_STEP(
        inputs_unlaned,
        non_integer_balance_balanced_gg_for_merge,
        tools_used
    )

    events_for_merge = EVENTS_STEP.out.events.map { it -> [ it[0].patient, it[1] ] } // meta.patient, events

    FUSIONS_STEP(
        inputs_unlaned,
        non_integer_balance_balanced_gg_for_merge,
        vcf_from_sv_calling_for_merge,
        tools_used
    )
    fusions_for_merge = FUSIONS_STEP.out.fusions
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, fusions
    altedge_annotations_for_merge = FUSIONS_STEP.out.altedge_annotations
            .map { it -> [ it[0].patient, it[1] ] } // meta.patient, altedge_annotations

    // FFPE Impact Filtering encoded inside SIGNATURES_STEP() workflow
    SIGNATURES_STEP(
		inputs_unlaned,
		snv_somatic_annotations_for_merge, // meta.patient, annotated somatic snv vcf
		tools_used
	)

	sbs_signatures = SIGNATURES_STEP.out.sbs_signatures
	indel_signatures = SIGNATURES_STEP.out.indel_signatures
	signatures_matrix = SIGNATURES_STEP.out.signatures_matrix
	sbs_activities = SIGNATURES_STEP.out.sbs_activities
	indel_activities = SIGNATURES_STEP.out.indel_activities
	sbs_posterior_prob = SIGNATURES_STEP.out.sbs_posterior_prob
	indel_posterior_prob = SIGNATURES_STEP.out.indel_posterior_prob
	annotated_vcf_ffpe_impact_or_snpeff = SIGNATURES_STEP.out.annotated_vcf_ffpe_impact_or_snpeff // meta.patient, annotated vcf, tbi

	annotated_vcf_ffpe_impact_or_snpeff_for_merge = annotated_vcf_ffpe_impact_or_snpeff
		.map { it -> [ it[0].patient, it[1] ] } // meta.patient, annotated vcf
    

    MULTIPLICITY_STEP(
        inputs_unlaned,
        annotated_vcf_ffpe_impact_or_snpeff_for_merge,
        non_integer_balance_balanced_gg_for_merge,
        hets_sites_for_merge,
        dryclean_tumor_cov_for_merge,
        snv_germline_annotations_for_merge,
        tools_used
    )

    snv_multiplicity_for_merge = MULTIPLICITY_STEP.out.snv_multiplicity
            .map { it -> [ it[0].patient, it[1], it[2], it[3] ] } // meta.patient, snv multiplicity rds, snv_multiplicity_germline_rds, snv_multiplicity_hets_rds

    if (params.tumor_only) {
        snv_multiplicity_for_merge.map{
            patient, snv_multiplicity_rds, snv_multiplicity_germline_rds, snv_multiplicity_hets_rds ->
            [ patient, snv_multiplicity_rds, [], snv_multiplicity_hets_rds ]
        }
    }

    ONCOKB_STEP(
        inputs_unlaned,
        annotated_vcf_ffpe_impact_or_snpeff_for_merge,
        non_integer_balance_balanced_gg_for_merge,
        fusions_for_merge,
        tools_used
    )

    HRDETECT_STEP(
        inputs_unlaned,
        vcf_from_sv_calling_for_merge,
        hets_sites_for_merge,
        annotated_vcf_ffpe_impact_or_snpeff_for_merge,
        // non_integer_balance_balanced_gg_for_merge,
        jabba_rds_for_merge,
        tools_used
    )
    hrdetect_rds_for_merge = HRDETECT_STEP.out.hrdetect_rds
                .map { it -> [ it[0].patient, it[1] ] } // meta.patient, hrdetect rds

    ONENESS_TWONESS_STEP(
        inputs_unlaned,
        events_for_merge,
        hrdetect_rds_for_merge,
        tools_used
    )

    onenesstwoness_rds_for_merge = ONENESS_TWONESS_STEP.out.onenesstwoness_rds
                .map { it -> [ it[0].patient, it[1] ] } // meta.patient, onenesstwoness rds
}