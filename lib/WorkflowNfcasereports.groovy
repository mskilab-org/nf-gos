//
// This file holds several functions specific to the workflow/nfcasereports.nf in the mskilab-org/nf-casereports pipeline
//

import nextflow.Nextflow
import nextflow.Channel
import groovy.text.SimpleTemplateEngine

class WorkflowNfcasereports {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        genomeExistsError(params, log)


        if (!params.fasta) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        }

    }

    public static create_file_channel(parameter, else_part = Channel.empty()) {
        return parameter ? Channel.fromPath(parameter).collect() : else_part
    }

    public static create_value_channel(parameter) {
        return parameter ? Channel.value(parameter) : Channel.empty()
    }

    public static create_index_channel(param, param_tbi, prepare_genome_out) {
        if (param) {
            if (param_tbi) {
                return Channel.fromPath(param_tbi).collect()
            } else {
                return prepare_genome_out
            }
        } else {
            return Channel.value([])
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        // TODO Optionally add in-text citation tools to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def citation_text = [
                "Tools used in the workflow included:",
                "FastQC (Andrews 2010),",
                "MultiQC (Ewels et al. 2016)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        // TODO Optionally add bibliographic entries to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        // TODO Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
        //meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        //meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    private static def samplesheetToList(String filePath) {
        def sampleList = []
        if (!filePath || filePath.trim().isEmpty()) {
            return sampleList
        }
        def lines = new File(filePath).readLines()

        if (lines.isEmpty()) {
            return sampleList
        }

        def headers = lines[0].split(',')
        lines.drop(1).each { line ->
            def values = line.split(',')
            def rowMap = [:]
            headers.eachWithIndex { header, index ->
                if (index < values.size()) {
                    rowMap[header] = values[index]
                } else {
                    rowMap[header] = null
                }
            }
            sampleList.add(rowMap)
        }
        return sampleList
    }

    public static List determineToolsToRun(params) {
        def tool_input_output_map = [
            "aligner": [ inputs: ['fastq_1', 'fastq_2'], outputs: ['bam'] ],
            "bamqc": [ inputs: ['bam'], outputs: ['wgs_metrics', 'alignment_metrics', 'insert_size_metrics'] ],
            "postprocessing": [ inputs: ['bam'], outputs: [] ],
            "msisensorpro": [ inputs: ['bam'], outputs: ['msi', 'msi_germline'] ],
            "gridss": [ inputs: ['bam'], outputs: ['vcf'] ],
            "amber": [ inputs: ['bam'], outputs: ['hets', 'amber_dir'] ],
            "fragcounter": [ inputs: ['bam'], outputs: ['frag_cov'] ],
            "dryclean": [ inputs: ['frag_cov'], outputs: ['dryclean_cov'] ],
            "cbs": [ inputs: ['dryclean_cov'], outputs: ['seg', 'nseg'] ],
            "sage": [ inputs: ['bam'], outputs: ['snv_somatic_vcf', 'snv_germline_vcf'] ],
            "cobalt": [ inputs: ['bam'], outputs: ['cobalt_dir'] ],
            "purple": [ inputs: ['cobalt_dir', 'amber_dir'], outputs: ['purity', 'ploidy'] ],
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

        def sampleList = samplesheetToList(params.input)
        def available_inputs = new HashSet()
        sampleList.each { input_map ->
            input_map.each { key, value ->
                available_inputs.add(key)
                if (!value || (value instanceof Collection && value.empty)) {
                    available_inputs.remove(key)
                }
            }
        }

        println "Provided inputs: ${available_inputs}"

        def skip_tools = params.skip_tools ? params.skip_tools.split(',')
            .collect {
            def tool = it.trim()
            if (!tool_input_output_map.containsKey(tool)) {
                println "Tool ${tool} not found in tool_input_output_map. Skipping."
            } else { return tool }
        } : []

        // remove null entries from skip_tools
        skip_tools = skip_tools.findAll { it != null }

        println "Skipping tools: ${skip_tools}"
        // TODO: if GRIDSS - skip if vcf is found, but not if vcf_unfiltered is present.
        def selected_tools = []
        boolean changed
        do {
            changed = false
            tool_input_output_map.each { tool, io ->
                if (!selected_tools.contains(tool) && !skip_tools.contains(tool)) {
                    def inputsRequired = io.inputs
                    def inputsPresent = inputsRequired.every { available_inputs.contains(it) }
                    def outputsNeeded = io.outputs.any { !available_inputs.contains(it) }
                    if (tool == "sage" && params.tumor_only) {
                        outputsNeeded = !available_inputs.contains("snv_somatic_vcf")
                    }
                    if (inputsPresent && outputsNeeded) {
                        selected_tools.add(tool)
                        available_inputs.addAll(io.outputs)
                        changed = true
                    }
                }
            }
        } while (changed)
        return selected_tools
    }
}
