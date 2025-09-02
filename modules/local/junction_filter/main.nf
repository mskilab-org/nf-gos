process JUNCTION_FILTER_BEDTOOLS {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.13': 
        'mskilab/unified:0.0.13' }"

    input:
    tuple val(meta), path(sv_vcf), path(sv_vcf_tbi)
    path(junction_pon_dir)
    val(padding)

    output:
    tuple val(meta), path("somatic.filtered.gnoMAD.sv.rds"), emit: pon_filtered_sv_rds, optional:true
    tuple val(meta), path("somatic.filtered.sv.no.gnomAD.rds"), emit: final_filtered_sv_rds, optional: false
    path "versions.yml", emit: versions, optional:false

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def low_memory = task.ext.low_memory ?: false
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_limit = (task.memory.toGiga() / 2).toInteger()
    def ncpus = (task.cpus / 2).toInteger()
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    #!/bin/bash

    export LIBDIR=\${NEXTFLOW_BIN_DIR}/R
    
    if [ "${low_memory}"" = "true" ]; then
        export LOW_MEMORY=true
    else
        export LOW_MEMORY=false
    fi

    export PADDING=${padding}
    printf "Converting SV VCF to BED format with padding %s\n" "\$PADDING"
	Rscript --vanilla \${LIBDIR}/sv2bed.R --sv ${sv_vcf} --padding \${PADDING} --cores ${ncpus * 2}
	

	printf "Sorting SV bed files\n"
	## bed-ops
	sort-bed ra_A.bed > ra_A.sorted.bed
	sort-bed ra_B.bed > ra_B.sorted.bed

	bgzip -f ra_A.sorted.bed
	bgzip -f ra_B.sorted.bed

	tabix -f -0 -p bed ra_A.sorted.bed.gz
	tabix -f -0 -p bed ra_B.sorted.bed.gz
    
    printf "Intersecting SV BED files with PON files\n"
    # pondir_a=\$(find ${junction_pon_dir}/*A.sorted.bed.gz)
	# pondir_b=\$(find ${junction_pon_dir}/*B.sorted.bed.gz)
	bedtools intersect -sorted -s -wa -wb -loj -a ra_A.sorted.bed.gz -b ${junction_pon_dir}/*A.sorted.bed.gz > overlap_AA.tsv
	bedtools intersect -sorted -s -wa -wb -loj -a ra_B.sorted.bed.gz -b ${junction_pon_dir}/*B.sorted.bed.gz > overlap_BB.tsv

	use_data_table=\$( ( ! \${LOW_MEMORY} ) && echo true || echo false )

	exit_merge=0
	if \$use_data_table; then
		printf "Using R data.table to merge SV hits\n"
		Rscript --vanilla \${LIBDIR}/mergebed.R --overlap_aa overlap_AA.tsv --overlap_bb overlap_BB.tsv
		exit_merge=\$?

	else
		printf "Using bash to merge SV hits\n"
		( join -j1 <(<./overlap_AA.tsv awk '{print \$4"______"\$10}' | sort -S ${mem_limit}G --parallel=${ncpus} -k1,1) <(<./overlap_BB.tsv awk '{print \$4"______"\$10}' | sort -S ${mem_limit}G --parallel=${ncpus} -k1,1) | awk -F'______'  '{print \$1"\t"\$2}'  > ra_match___AA_BB.tsv )
		exit_merge=\$?
		# | awk '{print \$1}' | uniq
	fi

	rm -f overlap_*.tsv

	is_failure=\$(test ! \$exit_merge -eq 0 && echo true || echo false)
	if \$is_failure; then
		echo "Failed to merge SV hits, exit code: \$exit_merge"
		exit \$exit_merge
	fi

    Rscript --vanilla \${LIBDIR}/filtersv.R --input_vcf ${sv_vcf} --ramatch_file ./ra_match___AA_BB.tsv
    exit_filter=\$?

    is_failure=\$( (test ! -e ./somatic.filtered.sv.no.gnomAD.rds || test ! \$exit_filter -eq 0) && echo true || echo false)
    if \$is_failure; then
		echo "Failed to produce filtered SV file"
		exit 1
	fi
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Junction_Filter: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch somatic.filtered.gnoMAD.sv.rds somatic.filtered.sv.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Junction_Filter: ${VERSION}
    END_VERSIONS
    """
}


process JUNCTION_FILTER {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/hrdetect:0.0.5':  // using hrdetect container since it has all dependencies
        'mskilab/hrdetect:0.0.5' }"

    input:
    tuple val(meta), path(filtered_sv_vcf), path(filtered_sv_vcf_tbi)
    path(junction_pon)
    path(gnomAD_sv_db)
    val(padding)

    output:
    tuple val(meta), path("somatic.filtered.gnoMAD.sv.rds")                , emit: pon_filtered_sv_rds,      optional:true
    tuple val(meta), path("somatic.filtered.sv.rds")                       , emit: final_filtered_sv_rds,        optional: false
    path "versions.yml"                                                    , emit: versions,                   optional:false

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    #!/bin/bash

    set -o allexport
    export RSCRIPT_PATH=\$(echo "\${NEXTFLOW_PROJECT_DIR}/bin/junction_filter.R")

    Rscript \$RSCRIPT_PATH           \\
    --sv        ${filtered_sv_vcf}   \\
    --pon       ${junction_pon}      \\
    --gnomAD    ${gnomAD_sv_db}      \\
    --padding   ${padding}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Junction_Filter: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch somatic.filtered.gnoMAD.sv.rds somatic.filtered.sv.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Junction_Filter: ${VERSION}
    END_VERSIONS
    """
}
