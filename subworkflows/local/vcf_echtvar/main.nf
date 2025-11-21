//
// VCF ECHTVAR
//

include { ECHTVAR_ANNO } from '../../../modules/local/echtvar/main.nf'
include { DECOMPOSE_VCF } from '../../../modules/local/echtvar/main.nf'




workflow VCF_ECHTVAR {
    // defining inputs
    take:
    echtvar_input     // required: [meta, vcf, tbi]

    //Creating empty channels for output
    main:
    fasta = WorkflowNfcasereports.create_file_channel(params.fasta)
    echtvar_dbnsfp = WorkflowNfcasereports.create_file_channel(params.echtvar_dbnsfp) // use dbnsfp database
    echtvar_clinvar = WorkflowNfcasereports.create_file_channel(params.echtvar_clinvar) // use dbnsfp database
    echtvar_civic = WorkflowNfcasereports.create_file_channel(params.echtvar_civic) // use dbnsfp database
    echtvar_bcf  = Channel.empty()

    versions = Channel.empty()

    // decompose VCF first (generally not needed)
    // DECOMPOSE_VCF(
    //     echtvar_input,
    //     fasta
    // )

    // echtvar_input = input.map { meta, vcf, _tbi -> [ meta, vcf ] }

    echtvar_database = echtvar_dbnsfp
        .mix(echtvar_clinvar)
        .mix(echtvar_civic)
        .collect()

    ECHTVAR_ANNO(
        echtvar_input
        .dump(tag: "ECHTVAR_ANNO input", pretty: true),
        echtvar_database
        .dump(tag: "ECHTVAR_ANNO database", pretty: true)
    )

    // initializing outputs from echtvar
    echtvar_bcf         = ECHTVAR_ANNO.out.bcf

    versions                = ECHTVAR_ANNO.out.versions

    emit:
    echtvar_bcf

    versions
}
