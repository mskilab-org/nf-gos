//
// This file holds several functions specific to the main.nf workflow in the mskilab-org/nf-jabba pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The nf-gos pipeline\n" +
            "  https://github.com/mskilab-org/nf-gos/\n\n"
    }


    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }
    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }

    //
    // Promote every key from params.genomes[params.genome] onto the top-level params map.
    // Existing values on params are preserved, so explicit user overrides (CLI or config)
    // beat the genome preset. Keys present on params but set to null are treated as unset
    // and get filled from the genome block.
    //
    public static void loadGenomeParams(params) {
        if (!(params.genomes && params.genome && params.genomes.containsKey(params.genome))) {
            return
        }
        params.genomes[ params.genome ].each { key, value ->
            if (!params.containsKey(key) || params[key] == null) {
                params[key] = value
            }
        }
    }
}
