//
// VCF SIGPROFILERASSIGNMENT
//

include { SIGPROFILERASSIGNMENT } from '../../../modules/local/sigprofilerassignment/main.nf'

workflow VCF_SIGPROFILERASSIGNMENT {
    // defining inputs
    take:
    input                                             // required: [meta, vcf, tbi]
    genome
    cosmic_version

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    signatures   = Channel.empty()

    SIGPROFILERASSIGNMENT(
        input,
        genome,
        cosmic_version
    )

    signatures   = SIGPROFILERASSIGNMENT.out.sigs

    emit:
    signatures

    versions
}
