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
    // Return a map of every key from params.genomes[params.genome] that is not already
    // set (or is null) on the top-level params map.  Intended to be used with
    // params.putAll() at the Nextflow script top-level so genome attributes are visible
    // before any include / workflow declarations are parsed:
    //
    //   params.putAll( WorkflowMain.genomeParams(params) )
    //
    // Existing non-null values on params are preserved, so explicit user overrides
    // (CLI flags or config entries) always beat the genome preset.
    //
    public static Map genomeParams(params) {
        Map resolved = [:]
        if (!(params.genomes && params.genome && params.genomes.containsKey(params.genome))) {
            return resolved
        }
        params.genomes[ params.genome ].each { key, value ->
            if (!params.containsKey(key) || params[key] == null) {
                resolved[key] = value
            }
        }
        return resolved
    }

    //
    // Promote every key from params.genomes[params.genome] onto the top-level params map.
    // Kept for backwards-compatibility (used inside a workflow block where putAll is not
    // needed, but the in-place mutation still works at that scope).
    //
    public static void loadGenomeParams(params) {
        params.putAll( genomeParams(params) )
    }
}
