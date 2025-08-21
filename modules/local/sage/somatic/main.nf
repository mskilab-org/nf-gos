process SAGE_SOMATIC {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.4--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sage:3.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(normal_bam_wgs), path(normal_bai), path(tumor_bam_wgs), path(tumor_bai)
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
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_bam_wgs} \\
        -ref_genome ${ref} \\
        -ref_genome_version ${ref_genome_version} \\
        -hotspots ${somatic_hotspots} \\
        -panel_bed ${panel_bed} \\
        -high_confidence_bed ${high_confidence_bed} \\
        -ensembl_data_dir ${ensembl_data_dir} \\
        -disable_bqr \\
        -threads ${task.cpus} \\
        -output_vcf ${meta.id}.sage.somatic.vcf.gz

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
    bcftools view -i 'FILTER="PASS"' -O z -o ${meta.id}.sage.pass_filtered.vcf.gz ${sage_vcf} && tabix -p vcf ${meta.id}.sage.pass_filtered.vcf.gz
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

process SAGE_TUMOR_ONLY_FILTER {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/utils:0.0.2':
        'mskilab/utils:0.0.2' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_vcf_tbi)
    path(dbsnp)
    path(dbsnp_tbi)
    path(gnomAD_snv_db)
    path(gnomAD_snv_db_tbi)
    path(sage_germ_pon)
    path(sage_germ_pon_tbi)
    path(mills_gold_indel)
    path(mills_gold_indel_tbi)

    output:
    tuple val(meta), path("*.sage.pass_filtered.tumoronly.vcf.gz"), path("*.sage.pass_filtered.tumoronly.vcf.gz.tbi"),    emit: sage_filtered_vcf
    path "versions.yml",                                                                                        emit: versions, optional:true

    when:
    task.ext.when == null || task.ext.when
    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def output      = "${meta.id}.sage.pass_filtered.tumoronly.vcf.gz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir -p filter_tumoronly/

    if [[ -e ${meta.id}.sage.pass_filtered.vcf.gz ]]; then
        echo "Filtering the tumor-only vcf by gnomAD..."
        bcftools isec -C -O z -o ${meta.id}.sage.pass_filtered.nognomAD.vcf.gz -p ./ ${meta.id}.sage.pass_filtered.vcf.gz ${gnomAD_snv_db} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by gnomAD. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.vcf.gz ]]; then
        echo "Now filtering by SAGE Germline PON ..."
        bcftools isec -C -O z -o ${meta.id}.sage.pass_filtered.nognomAD.noPON.vcf.gz -p ./ filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.vcf.gz \\
        ${sage_germ_pon} && mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.noPON.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.noPON.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by SAGE Germline PON. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.noPON.vcf.gz ]]; then
        echo "Now filtering by Mills and Gold Standard Indels ..."
        bcftools isec -C -O z -o ${meta.id}.sage.pass_filtered.nognomAD.noPON.noMGIndel.vcf.gz -p ./ \\
        filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.noPON.vcf.gz ${mills_gold_indel} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.noPON.noMGIndel.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.noPON.noMGIndel.vcf.gz.tbi

        cp filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.noPON.noMGIndel.vcf.gz ./${meta.id}.sage.pass_filtered.tumoronly.vcf.gz && \\
        cp filter_tumoronly/${meta.id}.sage.pass_filtered.nognomAD.noPON.noMGIndel.vcf.gz.tbi ./${meta.id}.sage.pass_filtered.tumoronly.vcf.gz.tbi

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

process SAGE_TUMOR_ONLY_FILTER___DEPRECATED {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/utils:0.0.2':
        'mskilab/utils:0.0.2' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_vcf_tbi)
    path(dbsnp)
    path(dbsnp_tbi)
    path(gnomAD_snv_db)
    path(gnomAD_snv_db_tbi)
    path(sage_germ_pon)
    path(sage_germ_pon_tbi)
    path(mills_gold_indel)
    path(mills_gold_indel_tbi)

    output:
    tuple val(meta), path("*.sage.pass_filtered.tumoronly.vcf.gz"), path("*.sage.pass_filtered.tumoronly.vcf.gz.tbi"),    emit: sage_filtered_vcf
    path "versions.yml",                                                                                        emit: versions, optional:true

    when:
    task.ext.when == null || task.ext.when
    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def output      = "${meta.id}.sage.pass_filtered.tumoronly.vcf.gz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir -p filter_tumoronly/

    if [[ -e ${meta.id}.sage.pass_filtered.vcf.gz ]]; then
        echo "Filtering the tumor-only vcf by dbSNP..."
        bcftools isec -C -O z -o ${meta.id}.sage.pass_filtered.nodbsnp.vcf.gz -p ./ ${meta.id}.sage.pass_filtered.vcf.gz ${dbsnp} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by dbSNP. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.vcf.gz ]]; then
        echo "Now filtering by gnomAD SNV DB ..."
        bcftools isec -C -O z -o ${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.vcf.gz -p ./ filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.vcf.gz ${gnomAD_snv_db} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by gnomAD SNV database. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.vcf.gz ]]; then
        echo "Now filtering by SAGE Germline PON ..."
        bcftools isec -C -O z -o ${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.vcf.gz -p ./ filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.vcf.gz \\
        ${sage_germ_pon} && mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by SAGE Germline PON. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.vcf.gz ]]; then
        echo "Now filtering by Mills and Gold Standard Indels ..."
        bcftools isec -C -O z -o ${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.noMGIndel.vcf.gz -p ./ \\
        filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.vcf.gz ${mills_gold_indel} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.noMGIndel.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.noMGIndel.vcf.gz.tbi

        cp filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.noMGIndel.vcf.gz ./${meta.id}.sage.pass_filtered.tumoronly.vcf.gz && \\
        cp filter_tumoronly/${meta.id}.sage.pass_filtered.nodbsnp.nognomAD.noPON.noMGIndel.vcf.gz.tbi ./${meta.id}.sage.pass_filtered.tumoronly.vcf.gz.tbi

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


process RESCUE_CH_HEME {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.5':
        'mskilab/unified:0.0.5' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_vcf_tbi), path(sage_tumor_only_vcf), path(sage_tumor_only_vcf_tbi)
    path(heme_ref)

    output:
    tuple val(meta), path("*.sage.pass_rescued_concatenated.vcf.gz"), path("*.sage.pass_rescued_concatenated.vcf.gz.tbi"),    emit: sage_tumor_only_rescue_ch_vcf
    path "versions.yml",                                                                                        emit: versions, optional:true

    when:
    task.ext.when == null || task.ext.when
    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def output      = "${meta.id}.sage.pass_filtered.tumoronly.vcf.gz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def safe_tmpdir = System.getenv('TMPDIR') ?: '/tmp'

    """
    export HOME=/root
	export TMPDIR=${safe_tmpdir}
    set +u  # Disable unbound variable check
    source /opt/conda/etc/profile.d/conda.sh
    conda activate pact

    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/rescue_ch.R")

    Rscript \$RSCRIPT_PATH \\
        --vcf ${sage_vcf} \\
        --heme_db ${heme_ref} \\
        --samp_id ${meta.id} \\
        --outdir ./
    
    bgzip ${meta.id}.sage.pass_rescued.vcf

    bcftools index --tbi ${meta.id}.sage.pass_rescued.vcf.gz

    bcftools concat \\
        --allow-overlaps \\
        --remove-duplicates \\
        ${sage_tumor_only_vcf} \\
        ${meta.id}.sage.pass_rescued.vcf.gz \\
        -o ${meta.id}.sage.pass_rescued_concatenated.unsorted.vcf.gz \\
        -Oz
    
    
    bcftools index --tbi ${meta.id}.sage.pass_rescued_concatenated.unsorted.vcf.gz

    bcftools sort \\
        -o ${meta.id}.sage.pass_rescued_concatenated.vcf.gz \\
        ${meta.id}.sage.pass_rescued_concatenated.unsorted.vcf.gz

    rm ${meta.id}.sage.pass_rescued_concatenated.unsorted.vcf.gz*

    bcftools index --tbi ${meta.id}.sage.pass_rescued_concatenated.vcf.gz

    
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

