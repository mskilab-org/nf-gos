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
    def jvmheap_mem = (task.memory.toGiga() * 0.9).toInteger() // 
    def otherjvmheap_mem = (task.memory.toGiga() * 0.75).toInteger() // 

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
        --jvmheap ${jvmheap_mem}g \\
        --otherjvmheap ${otherjvmheap_mem}g \\
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



process GRIDSS_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"


    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fasta_fai)
    path(bwa_index)
    path(blacklist_gridss)                                                                    // optional: gridss blacklist bed file based on genome


    output:
    tuple val(meta), path("*gridss.working"), emit: gridss_preprocess_dir, optional:false
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
    def jvmheap_mem = (task.memory.toGiga() * 0.9).toInteger() // 
    def otherjvmheap_mem = (task.memory.toGiga() * 0.75).toInteger() // 

    """
    ${bwa}

    rm -rf *gridss.working/*

    gridss \\
        --labels ${meta.sample} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        --steps preprocess \\
        $blacklist \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT \\
        --jvmheap ${jvmheap_mem}g \\
        --otherjvmheap ${otherjvmheap_mem}g \\
        ${bam} \\

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


process GRIDSS_ASSEMBLE_SCATTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"


    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(tumor_gridss_preprocess_dir), path(normalbam), path(normalbai), path(normal_gridss_preprocess_dir), val(jobnodes), val(jobindex)
    path(fasta)
    path(fasta_fai)
    path(bwa_index)
    path(blacklist_gridss)                                                                    // optional: gridss blacklist bed file based on genome


    output:
    tuple val(meta), path("*assembly.bam.gridss.working/**"), emit: gridss_workdir, optional:false
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def sample_labels = normalbam.toString() ? "${meta.normal_id},${meta.tumor_id}" : "${meta.tumor_id}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def assembly_bam = "--assembly ${meta.id}.assembly.bam"
    def bwa = bwa_index ? "cp -s ${bwa_index}/* . || true" : ""
    def blacklist = blacklist_gridss ? "--blacklist ${blacklist_gridss}" : ""
    def jvmheap_mem = (task.memory.toGiga() * 0.9).toInteger() // 
    def otherjvmheap_mem = (task.memory.toGiga() * 0.75).toInteger() // 

    """
    ${bwa}

    gridss \\
        --labels ${sample_labels} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        $assembly_bam \\
        --steps assemble \\
        $blacklist \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT \\
        --jvmheap ${jvmheap_mem}g \\
        --otherjvmheap ${otherjvmheap_mem}g \\
        --jobindex ${jobindex} \\
        --jobnodes ${jobnodes} \\
        ${normalbam} \\
        ${tumorbam}

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


process GRIDSS_ASSEMBLE_GATHER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"


    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(tumor_gridss_preprocess_dir), path(normalbam), path(normalbai), path(normal_gridss_preprocess_dir), val(gridss_assembly_dir), path(gridss_assembly_paths, stageAs: "gridss_assembly_dir___staged/*")
    path(fasta)
    path(fasta_fai)
    path(bwa_index)
    path(blacklist_gridss)                                                                    // optional: gridss blacklist bed file based on genome


    output:
    tuple val(meta), path("*.assembly.bam"), path("*.assembly.bam.gridss.working/**"), emit: gridss_final_assembly, optional:false
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def sample_labels = normalbam.toString() ? "${meta.normal_id},${meta.tumor_id}" : "${meta.tumor_id}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def assembly_bam = "--assembly ${meta.id}.assembly.bam"
    def bwa = bwa_index ? "cp -s ${bwa_index}/* . || true" : ""
    def blacklist = blacklist_gridss ? "--blacklist ${blacklist_gridss}" : ""
    def jvmheap_mem = (task.memory.toGiga() * 0.9).toInteger() // 
    def otherjvmheap_mem = (task.memory.toGiga() * 0.75).toInteger() // 

    """
    mkdir -p ${gridss_assembly_dir}

    echo "gridss_assembly_dir: ${gridss_assembly_dir}"

    cp -P gridss_assembly_dir___staged/* ./${gridss_assembly_dir} || { echo "No assembly files found in gridss_assembly_dir___staged" && exit 1; }

    ls ./${gridss_assembly_dir}

    ${bwa}

    gridss \\
        --labels ${sample_labels} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        $assembly_bam \\
        --steps assemble \\
        $blacklist \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT \\
        --jvmheap ${jvmheap_mem}g \\
        --otherjvmheap ${otherjvmheap_mem}g \\
        ${normalbam} \\
        ${tumorbam}

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


process GRIDSS_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"


    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(tumor_gridss_preprocess_dir), path(normalbam), path(normalbai), path(normal_gridss_preprocess_dir), val(gridss_assembly_dir), path(gridss_assembly_paths, stageAs: "gridss_assembly_dir___staged/*"), path(gridss_final_assembly)
    path(fasta)
    path(fasta_fai)
    path(bwa_index)
    path(blacklist_gridss)                                                                    // optional: gridss blacklist bed file based on genome


    output:
    tuple val(meta), path("*vcf.gz"), path("*vcf.gz.tbi"), emit: vcf, optional:false
    tuple val(meta), path("*vcf.gz.tbi"), emit: vcf_index, optional:false
    tuple val(meta), path("*filtered.vcf.gz"), path("*filtered.vcf.gz.tbi"), emit: filtered_vcf, optional: false
    tuple val(meta), path("*filtered.vcf.gz.tbi"), emit: filtered_vcf_index, optional:false
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def sample_labels = normalbam.toString() ? "${meta.normal_id},${meta.tumor_id}" : "${meta.tumor_id}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def assembly_bam = "--assembly ${meta.id}.assembly.bam"
    def bwa = bwa_index ? "cp -s ${bwa_index}/* . || true" : ""
    def blacklist = blacklist_gridss ? "--blacklist ${blacklist_gridss}" : ""
    def jvmheap_mem = (task.memory.toGiga() * 0.9).toInteger() // 
    def otherjvmheap_mem = (task.memory.toGiga() * 0.75).toInteger() // 

    """
    mkdir -p ${gridss_assembly_dir}

    echo "gridss_assembly_dir: ${gridss_assembly_dir}"

    cp -P gridss_assembly_dir___staged/* ./${gridss_assembly_dir} || { echo "No assembly files found in gridss_assembly_dir___staged" && exit 1; }


    ${bwa}

    gridss \\
        --labels ${sample_labels} \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        $assembly_bam \\
        $blacklist \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT \\
        --jvmheap ${jvmheap_mem}g \\
        --otherjvmheap ${otherjvmheap_mem}g \\
        --steps call \\
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



