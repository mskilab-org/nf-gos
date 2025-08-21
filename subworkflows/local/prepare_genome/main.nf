//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths


include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { MSISENSORPRO_SCAN                      } from '../../../modules/nf-core/msisensorpro/scan/main'
include { SAMTOOLS_FAIDX                         } from '../../../modules/nf-core/samtools/faidx/main'
include { TABIX_TABIX as TABIX_DBSNP             } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS      } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_SNPS        } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PON               } from '../../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_CHR_DIR                 } from '../../../modules/nf-core/untar/main'
include { UNZIP as UNZIP_ALLELES                 } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GC                      } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_LOCI                    } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_RT                      } from '../../../modules/nf-core/unzip/main'

ascat_alleles                       = WorkflowNfcasereports.create_file_channel(params.ascat_alleles)
ascat_loci                          = WorkflowNfcasereports.create_file_channel(params.ascat_loci)
ascat_loci_gc                       = WorkflowNfcasereports.create_file_channel(params.ascat_loci_gc)
ascat_loci_rt                       = WorkflowNfcasereports.create_file_channel(params.ascat_loci_rt)
chr_dir                             = WorkflowNfcasereports.create_file_channel(params.chr_dir)
dbsnp                               = WorkflowNfcasereports.create_file_channel(params.dbsnp)
fasta                               = WorkflowNfcasereports.create_file_channel(params.fasta)
fasta_fai                           = WorkflowNfcasereports.create_file_channel(params.fasta_fai)
bwa                                 = WorkflowNfcasereports.create_file_channel(params.bwa)
germline_resource                   = WorkflowNfcasereports.create_file_channel(params.germline_resource)
known_indels                        = WorkflowNfcasereports.create_file_channel(params.known_indels)
known_snps                          = WorkflowNfcasereports.create_file_channel(params.known_snps)
pon                                 = WorkflowNfcasereports.create_file_channel(params.pon)
msisensorpro_list                   = WorkflowNfcasereports.create_file_channel(params.msisensorpro_list)

workflow PREPARE_GENOME {

    main:
    fasta = fasta.map{ fasta -> [ [ id:fasta.baseName[0] ], fasta ] }
    fasta.dump(tag: "fasta in PREPARE_GENOME", pretty: true)

    versions = Channel.empty()

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    // only run msisensorpro_scan if the msisensorpro_list is an empty channel
    println "params.msisensorpro_list: ${params.msisensorpro_list}"
    def is_msisensorpro_list_present = params.msisensorpro_list ? Files.exists(Paths.get(params.msisensorpro_list)) : false
    println "is_msisensorpro_list_present: ${is_msisensorpro_list_present}"
    if (! is_msisensorpro_list_present) {
        MSISENSORPRO_SCAN(fasta)
        msisensorpro_list = MSISENSORPRO_SCAN.out.list.map{ meta, list -> list }                // path: genome_msi.list
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }

    SAMTOOLS_FAIDX(fasta, [['id':null], []])


    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [ file1, file2 ] becomes [ [ meta1, file1 ], [ meta2, file2 ] ]
    // outputs are collected to maintain a single channel for relevant TBI files
    TABIX_DBSNP(dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_KNOWN_SNPS(known_snps.flatten().map{ it -> [ [ id:it.baseName ], it ] } )
    TABIX_KNOWN_INDELS(known_indels.flatten().map{ it -> [ [ id:it.baseName ], it ] } )
    TABIX_PON(pon.flatten().map{ it -> [ [ id:it.baseName ], it ] })


    // prepare ascat reference files
    allele_files = ascat_alleles
    if (params.ascat_alleles && params.ascat_alleles.endsWith('.zip')) {
        UNZIP_ALLELES(ascat_alleles.map{ it -> [[id:it[0].baseName], it]})
        allele_files = UNZIP_ALLELES.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_ALLELES.out.versions)
    }

    loci_files = ascat_loci
    if (params.ascat_loci && params.ascat_loci.endsWith('.zip')) {
        UNZIP_LOCI(ascat_loci.map{ it -> [[id:it[0].baseName], it]})
        loci_files = UNZIP_LOCI.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_LOCI.out.versions)
    }
    gc_file = ascat_loci_gc
    if (params.ascat_loci_gc && params.ascat_loci_gc.endsWith('.zip')) {
        UNZIP_GC(ascat_loci_gc.map{ it -> [[id:it[0].baseName], it]})
        gc_file = UNZIP_GC.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_GC.out.versions)
    }
    rt_file = ascat_loci_rt
    if (params.ascat_loci_rt && params.ascat_loci_rt.endsWith('.zip')) {
        UNZIP_RT(ascat_loci_rt.map{ it -> [[id:it[0].baseName], it]})
        rt_file = UNZIP_RT.out.unzipped_archive.map{ it[1] }
        versions = versions.mix(UNZIP_RT.out.versions)
    }


    chr_files = chr_dir
    if (params.chr_dir && params.chr_dir.endsWith('tar.gz')) {
        UNTAR_CHR_DIR(chr_dir.map{ it -> [ [ id:'chr_dir' ], it ] })
        chr_files = UNTAR_CHR_DIR.out.untar.map{ it[1] }
        versions = versions.mix(UNTAR_CHR_DIR.out.versions)
    }


    // Gather versions of all tools used
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    // versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    versions = versions.mix(TABIX_DBSNP.out.versions)
    versions = versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    versions = versions.mix(TABIX_KNOWN_SNPS.out.versions)
    versions = versions.mix(TABIX_KNOWN_INDELS.out.versions)
    versions = versions.mix(TABIX_PON.out.versions)



    emit:
    bwa                  = bwa
    dbsnp_tbi             = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()               // path: dbsnb.vcf.gz.tbi
    dict                  = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
    fasta_fai             = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                      // path: genome.fasta.fai
    germline_resource_tbi = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }.collect()   // path: germline_resource.vcf.gz.tbi
    known_snps_tbi        = TABIX_KNOWN_SNPS.out.tbi.map{ meta, tbi -> [tbi] }.collect()          // path: {known_indels*}.vcf.gz.tbi
    known_indels_tbi      = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
    msisensorpro_scan     = msisensorpro_list
    pon_tbi               = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }.collect()                 // path: pon.vcf.gz.tbi
    allele_files
    chr_files
    gc_file
    loci_files
    rt_file

    versions // channel: [ versions.yml ]
}
