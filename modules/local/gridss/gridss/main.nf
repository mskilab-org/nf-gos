process GRIDSS_GRIDSS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"


    input:
    tuple val(meta), path(normalbam), path(normalbai), path(tumorbam), path(tumorbai)
    path(fasta)
    path(fasta_fai)
    path(bwa_index)
    path(blacklist_gridss)                                                                    // optional: gridss blacklist bed file based on genome


    output:
    tuple val(meta), path("*gridss.vcf.gz"), path("*gridss.vcf.gz.tbi"), emit: vcf, optional:true
    tuple val(meta), path("*gridss.vcf.gz.tbi"), emit: vcf_index, optional:true
    tuple val(meta), path("*.assembly.bam"), emit: assembly, optional:true
    tuple val(meta), path("*gridss.filtered.vcf.gz"), path("*gridss.filtered.vcf.gz.tbi"), emit: filtered_vcf
    tuple val(meta), path("*gridss.filtered.vcf.gz.tbi"), emit: filtered_vcf_index, optional:true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def assembly_bam = "--assembly ${meta.id}.assembly.bam"
    def bwa = bwa_index ? "cp -s ${bwa_index}/* . || true" : ""
    def blacklist = blacklist_gridss ? "--blacklist ${blacklist_gridss}" : ""

    """
    ${bwa}

    gridss \\
        --labels ${meta.normal_id},${meta.tumor_id} \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        $assembly_bam \\
        $blacklist \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT \\
        --jvmheap 32g \\
        --otherjvmheap 32g \\
        ${normalbam} \\
        ${tumorbam}

    bcftools view -f PASS ${prefix}.vcf.gz -Oz -o ${prefix}.filtered.vcf.gz
    tabix -p vcf ${prefix}.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def steps = args.contains("-s ") ? args.split('-s ')[-1].split(" ")[0] :
                args.contains("--steps ") ? args.split('--steps ')[-1].split(" ")[0] :
                "all"
    def vcf = steps.contains("call") || steps.contains("all") ? "touch ${prefix}.vcf*" : ""
    def assembly_bam = steps.contains("assembly") || steps.contains("all") ? "touch ${meta.id}.assembly.bam" : ""
    """
    ${vcf}
    ${assembly_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}




