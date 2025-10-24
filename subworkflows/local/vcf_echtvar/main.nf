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
    echtvar_database = WorkflowNfcasereports.create_file_channel(params.echtvar_dbnsfp) // use dbnsfp database
    echtvar_bcf  = Channel.empty()

    versions        = Channel.empty()

    // decompose VCF first (generally not needed)
    // DECOMPOSE_VCF(
    //     echtvar_input,
    //     fasta
    // )

    // echtvar_input = input.map { meta, vcf, _tbi -> [ meta, vcf ] }

    ECHTVAR_ANNO(
        echtvar_input,
        echtvar_database
    )

    // initializing outputs from echtvar
    echtvar_bcf         = ECHTVAR_ANNO.out.bcf

    versions                = ECHTVAR_ANNO.out.versions

    emit:
    echtvar_bcf

    versions
}
