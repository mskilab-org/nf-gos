package mylib
class Globals {
  static Map selected_tools_map = [:]
  
  static List tools_used = []

  static List rowsAsMaps = []

  static Map tool_input_output_map = [
        "aligner": [ inputs: ['fastq_1', 'fastq_2'], outputs: ['bam'] ],
        "collect_wgs_metrics": [
            inputs: ['bam'],
            outputs: ['qc_coverage_metrics']
        ],
        "collect_multiple_metrics": [
            inputs: ['bam'],
            outputs: [
                ['qc_alignment_summary'],
                ['qc_insert_size']
            ]
        ],
        "estimate_library_complexity": [ inputs: ['bam'], outputs: ['qc_dup_rate'] ],
        // "bamqc": [ inputs: ['bam'], outputs: ['wgs_metrics', 'alignment_metrics', 'insert_size_metrics', "estimate_library_complexity"] ],
        "postprocessing": [ inputs: ['bam'], outputs: [] ], // FIXME: Postprocessing will never be selected as a tool given the current set of inputs/outputs, empty output means tool will not be selected. postprocessing tool must be controlled by params.is_run_post_processing.
        "msisensorpro": [ inputs: ['bam'], outputs: [
            'msi'
            // 'msi_germline'
            ]
        ],
        // TODO: figure out what the best way to
        // "gridss": [ inputs: ['bam'], outputs: ['vcf_raw'] ],
        // "junctionfilter": [ inputs: ['vcf_raw'], outputs: ['vcf'] ],
        // "retiered_filtered_junctions": [ inputs: ['vcf'], outputs: ['sv_retier'] ], // ?
        "gridss": [ inputs: ['bam'], outputs: ['vcf'] ],
        "amber": [ inputs: ['bam'], outputs: ['hets', 'amber_dir'] ],
        "fragcounter": [ inputs: ['bam'], outputs: ['frag_cov'] ],
        "dryclean": [ inputs: ['frag_cov'], outputs: ['dryclean_cov'] ],
        "cbs": [ inputs: ['dryclean_cov'], outputs: ['seg', 'nseg'] ],
        "sage": [ inputs: ['bam'], outputs: ['snv_somatic_vcf', 'snv_germline_vcf'] ],
        "cobalt": [ inputs: ['bam'], outputs: ['cobalt_dir'] ],
        "purple": [ inputs: ['cobalt_dir', 'amber_dir'], outputs: ['purity', 'ploidy'] ],
        // "jabba": [
        // 	inputs: [
        // 		['vcf', 'sv_retier'],
        // 		['hets', 'dryclean_cov', 'ploidy', 'seg', 'nseg']
        // 	],
        // 	outputs: ['jabba_rds', 'jabba_gg']
        // ],
        "jabba": [ inputs: [ 'vcf', 'hets', 'dryclean_cov', 'ploidy', 'seg', 'nseg'], outputs: ['jabba_rds', 'jabba_gg'] ],
        "non_integer_balance": [ inputs: ['jabba_gg'], outputs: ['ni_balanced_gg'] ],
        "lp_phased_balance": [ inputs: ['ni_balanced_gg'], outputs: ['lp_balanced_gg'] ],
        "events": [ inputs: ['ni_balanced_gg'], outputs: ['events'] ],
        "fusions": [ inputs: ['ni_balanced_gg'], outputs: ['fusions'] ],
        "snpeff": [ inputs: ['snv_somatic_vcf'], outputs: ['variant_somatic_ann', 'variant_somatic_bcf'] ],
        "snv_multiplicity": [ inputs: ['jabba_gg', 'variant_somatic_ann'], outputs: ['snv_multiplicity'] ],
        "oncokb": [ inputs: ['variant_somatic_ann', 'snv_multiplicity', 'jabba_gg', 'fusions'], outputs: ['oncokb_maf', 'oncokb_fusions', 'oncokb_cna'] ],
        "signatures": [ inputs: ['variant_somatic_ann'], outputs: ['sbs_signatures', 'indel_signatures', 'signatures_matrix'] ],
        "hrdetect": [ inputs: ['hets', 'vcf', 'jabba_gg', 'variant_somatic_ann'], outputs: ['hrdetect'] ],
        "onenesstwoness": [ inputs: ['events', 'hrdetect'], outputs: ['onenesstwoness'] ]
    ]
}