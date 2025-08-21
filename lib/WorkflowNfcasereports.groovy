//
// This file holds several functions specific to the workflow/nfcasereports.nf in the mskilab-org/nf-casereports pipeline
//

import nextflow.Nextflow
import nextflow.Channel
import groovy.json.JsonSlurper
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

	public static map_join(channel_a, channel_b, key) {
        channel_a
			.map{ it -> [it[key], it] }
			.cross(channel_b.map{it -> [it[key], it]})
			.map { it[0][1] + it[1][1] }
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

     public static String retrieveInput(params, log){
        def input = null
            if (!params.input && !params.build_only_index) {
                switch (params.step) {
                    case 'alignment':                 Nextflow.error("Can't start with step $params.step without samplesheet")
                                                    break
                    case 'markduplicates':          log.warn("Using file ${params.outdir}/csv/mapped.csv");
                                                    input = params.outdir + "/csv/mapped.csv"
                                                    break
                    case 'prepare_recalibration':   log.warn("Using file ${params.outdir}/csv/markduplicates_no_table.csv");
                                                    input = params.outdir + "/csv/markduplicates_no_table.csv"
                                                    break
                    case 'recalibrate':             log.warn("Using file ${params.outdir}/csv/markduplicates.csv");
                                                    input = params.outdir + "/csv/markduplicates.csv"
                                                    break
                    case 'variant_calling':         log.warn("Using file ${params.outdir}/csv/recalibrated.csv");
                                                    input = params.outdir + "/csv/recalibrated.csv"
                                                    break
                    // case 'controlfreec':         csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
                    case 'annotate':                log.warn("Using file ${params.outdir}/csv/variantcalled.csv");
                                                    input = params.outdir + "/csv/variantcalled.csv"
                                                    break
                    default:                        log.warn("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
                                                    Nextflow.error("Unknown step $params.step")
                }
            }
            return input

    }

    /**
    * Convert nf-schema samplesheet tuples to List<Map>.
    * - row shape assumed: [ metaMap, v1, v2, ... ]
    * - meta fields are kept under key "meta"
    * - non-meta fields mapped in schema-declared order
    */
    public static List<Map> samplesheetTuplesToMaps(ArrayList tuples, File schemaFile, boolean coerceEmptyToNull = false) {
        def schema = new JsonSlurper().parse(schemaFile)
        assert schema?.items?.properties instanceof Map : "Schema must define items.properties"

        // Preserve schema-declared order
        def propEntries = schema.items.properties.entrySet().asList()

        // Identify meta vs non-meta fields by presence of "meta" tag
        def metaKeys     = propEntries.findAll { it.value?.meta instanceof List && it.value.meta }*.key
        def nonMetaKeys  = propEntries.findAll { !(it.value?.meta instanceof List && it.value.meta) }*.key

        // If there are meta columns, nf-schema puts a combined meta Map at index 0
        def expectsMetaMap = metaKeys && tuples && tuples[0] instanceof List && tuples[0][0] instanceof Map

        tuples.collect { row ->
            assert row instanceof List : "Each row must be a List/tuple"

            Map meta = [:]
            List vals = row
            if (expectsMetaMap) {
                meta = (row[0] as Map) ?: [:]
                vals = row.drop(1)
            }

            // Map non-meta keys to positional values
            def paired = [nonMetaKeys, vals].transpose().collectEntries { k, v ->
                def val = (coerceEmptyToNull && (v instanceof List) && v.isEmpty()) ? null : v
                [(k): val]
            }

            // Attach meta (nested)
            return [meta: meta] + paired
        }
    }
}
