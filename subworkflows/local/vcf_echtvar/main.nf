//
// VCF ECHTVAR
//

include { ECHTVAR_ANNO } from '../../../modules/local/echtvar/main.nf'
include { DECOMPOSE_VCF } from '../../../modules/local/echtvar/main.nf'

echtvar_database = WorkflowNfcasereports.create_file_channel(params.echtvar_database)
fasta = WorkflowNfcasereports.create_file_channel(params.fasta)

workflow VCF_ECHTVAR {
    // defining inputs
    take:
    input                                       // required: [meta, vcf, tbi]

    //Creating empty channels for output
    main:
    echtvar_bcf             = Channel.empty()

    versions        = Channel.empty()

    // decompose VCF first
    DECOMPOSE_VCF(
        input,
        fasta
    )

    echtvar_input = DECOMPOSE_VCF.out.vcf // [meta, vcf]

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
