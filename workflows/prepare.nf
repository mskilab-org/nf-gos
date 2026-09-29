include { paramsSummaryLog; paramsSummaryMap; samplesheetToList } from 'plugin/nf-schema'
include { validateParameters; paramsHelp } from 'plugin/nf-schema'

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
        "echtvar"  : [
            params.echtvar_dbnsfp
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

    // rowsAsMaps = rowsAsMaps.collect { it ->
    //     it + [
    //         crai: it.cram ? it.cram + '.crai' : [],
    //         bai: it.bam ? it.bam + '.bai': [],
    //         vcf_tbi: it.vcf ? it.vcf + '.tbi' : [],
    //         vcf_raw_tbi: it.vcf_raw ? it.vcf_raw + '.tbi' : [],
    //         structural_variants_chimera_filtered_tbi: it.structural_variants_chimera_filtered ? it.structural_variants_chimera_filtered + '.tbi' : [],
    //         structural_variants_raw_chimera_filtered_tbi: it.structural_variants_raw_chimera_filtered ? it.structural_variants_raw_chimera_filtered + '.tbi' : [],
    //         snv_somatic_vcf_tumoronly_filtered_tbi: it.snv_somatic_vcf_tumoronly_filtered ? it.snv_somatic_vcf_tumoronly_filtered + ".tbi" : [],
    //         snv_somatic_vcf_rescue_ch_heme_tbi: it.snv_somatic_vcf_rescue_ch_heme ? it.snv_somatic_vcf_rescue_ch_heme + '.tbi' : [],
    //         snv_somatic_tbi: it.snv_somatic_vcf ? it.snv_somatic_vcf + '.tbi' : [],
    //         snv_germline_tbi: it.snv_germline_vcf ? it.snv_germline_vcf + '.tbi' : [],
    //         ffpe_impact_vcf_tbi: it.ffpe_impact_vcf ? it.ffpe_impact_vcf + '.tbi' : [],
    //         ffpe_impact_filtered_vcf_tbi: it.ffpe_impact_filtered_vcf ? it.ffpe_impact_filtered_vcf + '.tbi' : []
    //     ]
    // }
    
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
        .map { _patient_sample, num_lanes, ch_items ->

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

    Globals.inputs = inputs
    Globals.inputs_unlaned = inputs_unlaned

    emit:
    inputs
    inputs_unlaned

}

workflow TOOLS {

    main:
    tool_input_output_map = Globals.tool_input_output_map
    
    // see lib/Globals.groovy
    sampleList = Globals.rowsAsMaps
    inputs = Globals.inputs
    inputs_unlaned = Globals.inputs_unlaned

    // ── Parse tool-control params early ───────────────────────────────────
    // Parsed here (before available_inputs / missing_outputs) so that
    // --overwrite_subsequent can use run_tools as DAG seeds to clear
    // sampleList before the main tool-selection scan.
    skip_tools = params.skip_tools ? params.skip_tools.split(',').collect { it.trim() } : []
    log.info "Skipping tools: ${skip_tools}"
    run_tools = params.only_tools ? params.only_tools.split(',').collect { it -> it.trim() } : []
    force_tools = params.force_tools ? params.force_tools.split(',').collect { it -> it.trim() } : []
    run_tools = (run_tools + force_tools).unique()
    is_run_tools_populated = ! run_tools.isEmpty()
    if (is_run_tools_populated) {
        log.info "Running tools: ${run_tools}" 
    }
    
    is_overlapping = run_tools.any { it ->
        skip_tools.contains(it)
    }
    if (is_overlapping) {
        log.info "Overlapping tool sets specified in skip and only tools parameters.. defaulting to running the tool specified"
    }
    if (! ( force_tools.isEmpty() )) {
        log.info "Forcing tools: ${force_tools}" 
    }

    // ── overwrite_subsequent ──────────────────────────────────────────────
    // When --overwrite_subsequent is set, blank out every output column
    // produced by any explicitly requested tool (run_tools = only_tools ∪
    // force_tools) **and every downstream tool in the DAG** across all
    // samplesheet representations.  Clearing happens BEFORE available_inputs
    // and missing_outputs are computed so the tool-selection loop naturally
    // picks up downstream tools whose outputs are now absent.
    //
    // The tool_input_output_map encodes a DAG: tool A's outputs may be tool
    // B's inputs, whose outputs may in turn be tool C's inputs, etc.
    // A fixpoint forward-propagation walks the DAG to collect the full
    // transitive closure of fields that must be cleared, and (in parallel)
    // the set of tools downstream of the seeds — these are the
    // `subsequent_tools` used by is_scenario3 in the selection loop below.
    //
    // Always defined so the tool-selection loop can reference it
    // unconditionally; only populated when overwrite_subsequent triggers.
    subsequent_tools = new LinkedHashSet()

    if (params.overwrite_subsequent && is_run_tools_populated) {
        // Step 1 – seed with direct outputs of the requested tools.
        // .flatten() handles nested output lists (e.g. collect_multiple_metrics).
        def fields_to_clear = new HashSet()
        run_tools.each { tool ->
            def io = tool_input_output_map[tool]
            if (io) { fields_to_clear.addAll(io.outputs.flatten()) }
        }

        // Step 2 – propagate downstream through the DAG until fixpoint.
        // A tool is "subsequent" iff any of its inputs is in fields_to_clear;
        // we record it AND add its outputs to the set, then iterate.
        def changed = true
        while (changed) {
            changed = false
            tool_input_output_map.each { tool, io ->
                def tool_outputs = io.outputs.flatten() as Set
                if (io.inputs.any { fields_to_clear.contains(it) }
                        && !fields_to_clear.containsAll(tool_outputs)) {
                    fields_to_clear.addAll(tool_outputs)
                    subsequent_tools.add(tool)
                    changed = true
                }
            }
        }

        log.info "overwrite_subsequent: clearing output columns ${fields_to_clear} from samplesheet rows"
        log.info "overwrite_subsequent: subsequent tools to re-run: ${subsequent_tools}"

        // Step 3 – build the overlay once, apply to all representations.
        def cleared = fields_to_clear.collectEntries { field -> [(field): []] }

        sampleList = sampleList.collect { row -> row + cleared }
        Globals.rowsAsMaps = sampleList

        inputs         = inputs.map         { row -> row + cleared }
        inputs_unlaned = inputs_unlaned.map { row -> row + cleared }
    }

    available_inputs = new HashSet()
    present_outputs = new HashSet()
    

    // Step 1 – seed available_inputs from the samplesheet: anything with a
    // non-empty value in at least one row is provided.
    sampleList.each { input_map ->
        input_map.each { key, value ->
            def is_value_present = value && !(value instanceof Collection && value.empty)
            if (is_value_present) {
                available_inputs.add(key)
            }
        }
    }

    log.info "Provided inputs (samplesheet only): ${available_inputs}"

    // Step 2 – DAG forward-propagation (fixpoint).
    // A tool's outputs are also "available" if (a) the tool is eligible to
    // run under the same is_scenarioN logic used by the main selection
    // loop below AND (b) all of its inputs are already in available_inputs.
    // This makes the inputsPresent check in the main loop independent of
    // tool_input_output_map iteration order: as long as the tool's
    // dependencies can transitively be produced from the samplesheet by
    // selectable tools, the inputs count as "are or will be present".
    def changed = true
    while (changed) {
        changed = false
        tool_input_output_map.each { tool, io ->
            def is_eligible_s1 = ! is_run_tools_populated && !skip_tools.contains(tool)
            def is_eligible_s2 = run_tools.contains(tool)
            def is_eligible_s3 = params.overwrite_subsequent && subsequent_tools.contains(tool)
            if (!(is_eligible_s1 || is_eligible_s2 || is_eligible_s3)) return

            def tool_outputs = io.outputs.flatten() as Set
            if (io.inputs.every { available_inputs.contains(it) }
                    && !available_inputs.containsAll(tool_outputs)) {
                available_inputs.addAll(tool_outputs)
                changed = true
            }
        }
    }

    log.info "Available inputs (samplesheet + transitively producible via DAG): ${available_inputs}"

    schemaFile = file("$projectDir/gos-assets/nf-gos/assets/schema_input.json")
    schema = new groovy.json.JsonSlurper().parse(schemaFile)


    props = schema.items.properties
    requiredFields = props.findAll { !it.value.containsKey('meta') }.keySet()
    log.info "requiredFields: $requiredFields"

    // Direct per-field missingness probe — does not depend on the schema or on
    // missing_outputs being populated. Utils.robustly_test_if_empty handles null
    // (missing map key), empty collections, blank strings, and missing/empty
    // files/paths uniformly.
    is_field_missing_in_any_sample = { field ->
        sampleList.any { sample -> Utils.robustly_test_if_empty(sample[field]) }
    }

    // Union schema-declared output fields with fields declared as outputs anywhere in
    // tool_input_output_map. This lets us mark a tool as "needed" based on its declared
    // outputs even when those outputs are not yet enumerated in schema_input.json —
    // avoiding the need to update the schema every time a tool is added or renamed.
    // .flatten() handles both flat lists (most tools) and nested lists
    // (e.g. collect_multiple_metrics: [['qc_alignment_summary'], ['qc_insert_size']]).
    all_tool_output_fields = tool_input_output_map.values()
        .collectMany { io -> io.outputs.flatten() as List } as Set
    candidate_output_fields = ((requiredFields as Set) + all_tool_output_fields) as Set

    missing_outputs = candidate_output_fields.findAll(is_field_missing_in_any_sample)
    log.info "Outputs MISSING from at least one sample: $missing_outputs"

    // Iteratively select tools based on available inputs
    // (skip_tools / run_tools / force_tools were parsed above.)
    // TODO: if GRIDSS - skip if vcf is found, but not if vcf_raw is present.
    selected_tools = []
    tools_qc = ["collect_wgs_metrics", "collect_multiple_metrics", "estimate_library_complexity"]
    selected_tools_map = [:]
    // is_scenario1 = run_tools not specified, so try to add to selected tools from the menu (tool_input_output_map)
    // is_scenario2 = run_tools is specified, so try to see if the tool matches the menu
    tool_input_output_map.each { tool, io ->
        def is_scenario1 = ! is_run_tools_populated && !selected_tools.contains(tool) && !skip_tools.contains(tool)
        def is_scenario2 = ! is_scenario1 && run_tools.contains(tool)
        // is_scenario3 = overwrite_subsequent is set and the tool is downstream
        // (in the tool_input_output_map DAG) of any tool in run_tools.
        // Enables re-running of cascading downstream tools even when run_tools
        // is populated (which would otherwise short-circuit is_scenario1).
        def is_scenario3 = params.overwrite_subsequent && subsequent_tools.contains(tool)
        if (is_scenario1 || is_scenario2 || is_scenario3) {

            def inputsRequired = io.inputs
            def inputsPresent = inputsRequired.every { available_inputs.contains(it) }
            
            // Probe each output field directly against sample rows rather than looking it
            // up in missing_outputs — this way fields that are not (yet) in the schema
            // still count as "needed" if absent from any sample.
            def outputsNeeded = io.outputs.flatten().any(is_field_missing_in_any_sample)

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
                outputsNeeded = ["snv_somatic_vcf", "snv_somatic_vcf_tumoronly_filtered"].any(is_field_missing_in_any_sample)
            }

            if (is_sage_heme) {
                outputsNeeded = ["snv_somatic_vcf", "snv_somatic_vcf_tumoronly_filtered", "snv_somatic_vcf_rescue_ch_heme"].any(is_field_missing_in_any_sample)
            }

            if (is_current_tool_qc_multiple_metrics) {
                def is_any_alignment_summary_absent = io.outputs[0].any(is_field_missing_in_any_sample)
                def is_any_insert_size_absent = io.outputs[1].any(is_field_missing_in_any_sample)
                outputsNeeded = is_any_alignment_summary_absent || is_any_insert_size_absent
            }

            if (force_tools && force_tools.contains(tool)) {
                log.info "Tool ${tool} is being forced to run by user request, so it will be added to the selected tools list even if its outputs are not needed or its inputs are not present."
                outputsNeeded = true
                
            }
            
            if (inputsPresent && outputsNeeded) {
                selected_tools.add(tool)
                // available_inputs is no longer mutated here — it was
                // pre-computed via fixpoint propagation above so
                // inputsPresent is order-independent.
            }

            log.info "tool: ${tool} \n inputsPresent: ${inputsPresent} \n outputsNeeded: ${outputsNeeded}"




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

    // Publish (potentially overwritten) channels to Globals so they are
    // accessible from workflow NFTAPS via TOOLS.out.
    Globals.inputs         = inputs
    Globals.inputs_unlaned = inputs_unlaned

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
    inputs           // Channel<Map> — samplesheet rows (lane-aware), output columns cleared if overwrite_subsequent
    inputs_unlaned   // Channel<Map> — same but lane meta stripped


}
