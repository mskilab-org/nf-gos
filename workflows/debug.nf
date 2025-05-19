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

tool_input_output_map = [
    "aligner": [ inputs: ['fastq_1', 'fastq_2'], outputs: ['bam'] ],
    "bamqc": [ inputs: ['bam'], outputs: [] ],
    "msisensorpro": [ inputs: ['bam'], outputs: ['msi', 'msi_germline'] ],
    "gridss": [ inputs: ['bam'], outputs: ['vcf'] ],
    "amber": [ inputs: ['bam'], outputs: ['hets', 'amber_dir'] ],
    "fragcounter": [ inputs: ['bam'], outputs: ['frag_cov'] ],
    "dryclean": [ inputs: ['frag_cov'], outputs: ['dryclean_cov'] ],
    "cbs": [ inputs: ['dryclean_cov'], outputs: ['seg', 'nseg'] ],
    "sage": [ inputs: ['bam'], outputs: ['snv_somatic_vcf', 'snv_germline_vcf'] ],
    "purple": [ inputs: ['bam', 'amber_dir'], outputs: ['purity', 'ploidy'] ],
    "jabba": [ inputs: ['vcf', 'hets', 'dryclean_cov', 'ploidy', 'seg', 'nseg'], outputs: ['jabba_rds', 'jabba_gg'] ],
    "non_integer_balance": [ inputs: ['jabba_gg'], outputs: ['ni_balanced_gg'] ],
    "lp_phased_balance": [ inputs: ['ni_balanced_gg'], outputs: ['lp_balanced_gg'] ],
    "events": [ inputs: ['ni_balanced_gg'], outputs: ['events'] ],
    "fusions": [ inputs: ['ni_balanced_gg'], outputs: ['fusions'] ],
    "snpeff": [ inputs: ['snv_somatic_vcf'], outputs: ['variant_somatic_ann', 'variant_somatic_bcf'] ],
    "snv_multiplicity": [ inputs: ['jabba_gg', 'variant_somatic_ann'], outputs: ['snv_multiplicity'] ],
    "oncokb": [ inputs: ['variant_somatic_ann', 'snv_multiplicity', 'jabba_gg', 'fusions'], outputs: ['oncokb_maf', 'oncokb_fusions', 'oncokb_cna'] ],
    "signatures": [ inputs: ['snv_somatic_vcf'], outputs: ['sbs_signatures', 'indel_signatures', 'signatures_matrix'] ],
    "hrdetect": [ inputs: ['hets', 'vcf', 'jabba_gg', 'snv_somatic_vcf'], outputs: ['hrdetect'] ],
    "onenesstwoness": [ inputs: ['events', 'hrdetect'], outputs: ['onenesstwoness'] ]
]

def samplesheetToList(String filePath) {
    def sampleList = []
    def lines = new File(filePath).readLines()

    if (lines.isEmpty()) {
        return sampleList // Return an empty list if the file is empty
    }

    // Assume the first line contains the headers
    def headers = lines[0].split(',')

    // Process each subsequent line as a data row
    lines.drop(1).each { line ->
        def values = line.split(',')
        def rowMap = [:]

        headers.eachWithIndex { header, index ->
            if (index < values.size()) {
                rowMap[header] = values[index]
            } else {
                rowMap[header] = null // Handle missing values
            }
        }

        sampleList.add(rowMap)
    }

    return sampleList
}

def sampleList = samplesheetToList(params.input)
def available_inputs = new HashSet()
sampleList.each { input_map ->
    input_map.each { key, value ->
        if (value && !(value instanceof Collection && value.empty)) {
            available_inputs.add(key)
        }
    }
}

println "Provided inputs: ${available_inputs}"

// Iteratively select tools based on available inputs
def skip_tools = params.skip_tools ? params.skip_tools.split(',').collect { it.trim() } : []
println "Skipping tools: ${skip_tools}"
def selected_tools = []
boolean changed
do {
    changed = false
    tool_input_output_map.each { tool, io ->
        if (!selected_tools.contains(tool) && !skip_tools.contains(tool)) {
            def inputsRequired = io.inputs
            def inputsPresent = inputsRequired.every { available_inputs.contains(it) }
            def outputsNeeded = io.outputs.any { !available_inputs.contains(it) }
            if (inputsPresent && outputsNeeded) {
                selected_tools.add(tool)
                available_inputs.addAll(io.outputs)
                changed = true
            }
        }
    }
} while (changed)

tools_used = selected_tools

println "Tools that will be run based on your inputs: ${tools_used}"


final_filtered_sv_rds_for_merge = inputs
    .map { it -> [it.meta, it.vcf, it.vcf_tbi] }
	.map { 
		def vcf = it[1]
		// def is_vcf = vcf =~ /\.vcf(\.gz|\.bgz)?$/
		def is_rds = vcf =~ /\.rds$/
		// println it[1]
		// println it[1].getClass()
		// println it[1].isEmpty()
		// println it[1].exists()
		// println is_vcf
		// println is_rds
		// println "tbi"
		// println it[2]
		// println it[2].getClass()
		// println it[2].isEmpty()
		// println it[2].exists()
		return it 
	}
    .filter { 
		def vcf it[1]
		!it[1].isEmpty() && !it[2].isEmpty() }
	.map {
		println "after filter"
		println it[0].patient
		println it[1]
		println it[2]
		return it 
	}
    .map { 
		it -> [ it[0].patient, it[1] ] 
	} // meta.patient, vcf

println "${final_filtered_sv_rds_for_merge.view()}"