process SAGE_SOMATIC {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.4--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sage:3.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(tumor_bam_wgs, stageAs: "${meta.id}.tumor.bam"), path(tumor_bai, stageAs: "${meta.id}.tumor.bam.bai"), path(normal_bam_wgs, stageAs: "${meta.id}.normal.bam"), path(normal_bai, stageAs: "${meta.id}.normal.bam.bai")
    path(ref)
    path(ref_fai)
    path(ref_genome_dict)
    val(ref_genome_version)
    path(ensembl_data_dir)
    path(somatic_hotspots)
    path(panel_bed)
    path(high_confidence_bed)

    output:
    tuple val(meta), path('*.sage.somatic.vcf.gz'), path('*.sage.somatic.vcf.gz.tbi'), emit: vcf
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def output      = "${meta.id}.sage.vcf.gz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam_wgs ? "-reference_bam ${normal_bam_wgs}" : ''

    """
    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        -tumor ${meta.id} \\
        -tumor_bam ${tumor_bam_wgs} \\
        -ref_genome ${ref} \\
        -ref_genome_version ${ref_genome_version} \\
        -hotspots ${somatic_hotspots} \\
        -panel_bed ${panel_bed} \\
        -high_confidence_bed ${high_confidence_bed} \\
        -ensembl_data_dir ${ensembl_data_dir} \\
        -disable_bqr \\
        -threads ${task.cpus} \\
        -output_vcf ${meta.tumor_id}.sage.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sage.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: ${VERSION}
    END_VERSIONS
    """
}

process SAGE_PASS_FILTER {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/utils:0.0.2':
        'mskilab/utils:0.0.2' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_vcf_tbi)

    output:
    tuple val(meta), path('*.sage.pass_filtered.vcf.gz'), path('*.sage.pass_filtered.vcf.gz.tbi'), emit: vcf

    script:

    """
    bcftools view -i 'FILTER="PASS"' -O z -o ${meta.id}.sage.pass_filtered.vcf.gz ${output} && tabix -p vcf ${meta.id}.sage.pass_filtered.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.sage.pass_filtered.vcf.gz
    touch ${meta.id}.sage.pass_filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools -v | head -n 1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}

process SAGE_FILTER {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.4--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sage:3.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_vcf_tbi)
    path(dbsnp)
    path(gnomAD_snv_db)
    path(sage_germ_pon)
    path(mills_gold_indel)

    output:
    tuple val(meta), path("*.sage.filtered.tumoronly.vcf.gz"), path("*.sage.filtered.tumoronly.vcf.gz.tbi"),                  emit: sage_filtered_vcf
    path "versions.yml"                                                                                                     emit: versions,                 optional:true

    when:
    task.ext.when == null || task.ext.when
    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def output      = "${meta.id}.sage.filtered.tumoronly.vcf.gz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir -p filter_tumoronly/

    if [[ -e ${meta.id}.sage.filtered.vcf.gz ]]; then
        echo "Filtering the tumor-only vcf by dbSNP..."
        bcftools isec -C -O z -o ${meta.id}.sage.filtered.nodbsnp.vcf.gz -p ./ ${meta.id}.sage.filtered.vcf.gz ${dbsnp} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by dbSNP. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.vcf.gz ]]; then
        echo "Now filtering by gnomAD SNV DB ..."
        bcftools isec -C -O z -o ${meta.id}.sage.filtered.nodbsnp.nognoMAD.vcf.gz -p ./ filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.vcf.gz ${gnomAD_snv_db} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by gnomAD SNV database. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.vcf.gz ]]; then
        echo "Now filtering by SAGE Germline PON ..."
        bcftools isec -C -O z -o ${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.vcf.gz -p ./ filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.vcf.gz \\
        ${sage_germ_pon} && mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by SAGE Germline PON. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.vcf.gz ]]; then
        echo "Now filtering by Mills and Gold Standard Indels ..."
        bcftools isec -C -O z -o ${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.noMGIndel.vcf.gz -p ./ \\
        filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.vcf.gz ${mills_gold_indel} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.noMGIndel.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.noMGIndel.vcf.gz.tbi

        cp filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.noMGIndel.vcf.gz.tbi ./${meta.id}.sage.filtered.tumoronly.vcf.gz && \\
        cp filter_tumoronly/${meta.id}.sage.filtered.nodbsnp.nognoMAD.noPON.noMGIndel.vcf.gz.tbi ./${meta.id}.sage.filtered.tumoronly.vcf.gz.tbi

        rm -rf filter_tumoronly/
    else
        echo "Cannot perform SNV filtering by Mills and Gold Standard Indels. Exiting..."
        exit 1;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools -v | head -n 1 | sed 's/^bcftools //')
    END_VERSIONS

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${meta.id}.sage.filtered.tumoronly.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: ${VERSION}
    END_VERSIONS
    """
}
