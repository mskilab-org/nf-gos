//
// BAM SAGE
//

include { SAGE } from '../../../modules/local/sage/main.nf'

workflow BAM_SAGE {
    // defining inputs
    take:
    inputs      // [meta, tumor bam, tumor bai, normal bam, normal bai]
    ref
    ref_fai
    ref_genome_version
    ensembl_data_dir
    somatic_hotspots
    panel_bed
    high_confidence_bed

    //Creating empty channels for output
    main:
    versions            = Channel.empty()
    sage_vcf            = Channel.empty()

    SAGE(
        inputs,
        ref,
        ref_fai,
        ref_genome_version,
        ensembl_data_dir,
        somatic_hotspots,
        panel_bed,
        high_confidence_bed
    )

    // initializing outputs from fragcounter
    versions        = SAGE.out.versions
    sage_vcf        = SAGE.out.sage_vcf

    //
    emit:
    sage_vcf

    versions
}
