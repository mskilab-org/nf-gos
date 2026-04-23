process SV_CHIMERA_FILTER {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.12':
        'mskilab/unified:0.0.12' }"

    input:
	tuple val(meta), path(vcf), path(vcf_tbi)

    output:
    tuple val(meta), path("*.ffpe_filtered.vcf.gz"), path("*.ffpe_filtered.vcf.gz.tbi"), emit: vcftbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_vcf = vcf.getName().replaceFirst(/\.vcf(\.gz|\.bgz)?$/, '.ffpe_filtered.vcf.gz')
    """
    bcftools view -Oz -i "(FORMAT/SR[1] + FORMAT/RP[1]) > 6 && (INFO/AS) >= 1 && (FORMAT/QUAL[1]) >= 150" ${vcf} > ${out_vcf}

    bcftools index --tbi ${out_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ffpe_filtered.bam
    touch ${prefix}.ffpe_filtered.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}

process ITDSEEK {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.17-itdseek':
        'mskilab/unified:0.0.17-itdseek' }"

    input:
	tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    val(assembly)

    output:
    tuple val(meta), path("*_flt3_itd.vcf"), path("*_flt3_itd_status.rds"), emit: vcfrds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def out_vcf = vcf.getName().replaceFirst(/\.vcf(\.gz|\.bgz)?$/, '.ffpe_filtered.vcf.gz')
    """

    Rscript \${NEXTFLOW_BIN_DIR}/R/flt3_itd_seek.R \\
        --bam ${bam} \\
        --fasta ${fasta} \\
        --samp_id ${meta.sample} \\
        --build ${assembly} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            samtools: \$(echo \$(samtools version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}


process CHECK_REFERENCE_BAM_MATCH {
    label 'process_medium'

    // using events container since the dependencies are the same
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/unified:0.0.18-fusions':
        'mskilab/unified:0.0.18-fusions' }"

    input:
    tuple val(meta), path(bam_path)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), env(is_fully_matched), path(bam_path), emit: is_fully_matched
    tuple val(meta), path("*mismatches.tsv"), emit: mismatched_rnames
    tuple val(meta), path("*all.tsv"), emit: all_rnames

    script:
    """
    #!/usr/bin/env bash

    ## ls -1 bam_paths/** | parallel 'samtools view -H {}' | grep '@SQ' | cut -f2-3 | sed -E 's/SN://g; s/LN://g' | sort -k1,1 | uniq > bam_rname.txt
    
    samtools view -H ${bam_path} | grep '@SQ' | cut -f2-3 | sed -E 's/SN://g; s/LN://g' | sort -k1,1 | uniq > bam_rname.txt

    ${fasta_fai} | cut -f1,2 | sort -k1,1 > fasta_rname.txt

    comm -3 <(sort bam_rname.txt) <(sort fasta_reads.txt) > ${meta.sample}.mismatches.tsv
    comm <(sort bam_rname.txt) <(sort fasta_reads.txt) > ${meta.sample}.all.tsv
    num_mismatch=\$(cat mismatches.tsv | wc -l)
    is_fully_matched=\$( [ "\$num_mismatch" -eq 0 ] && echo true || echo false )

    

    """
}


process RASTAIR {
    tag "$meta.id"

    // conda '/opt/conda/envs/my-existing-env'
    beforeScript """
        set -a
        . /gpfs/home/hadik01/load_miniforge
        conda activate rna --stack
        conda activate vcftools --stack
        conda activate rastair --stack
        set +a
    """

    input:
	tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(gnomAD_snv_db)
    path(gnomAD_snv_db_tbi)
    path(sage_germline_pon)
    path(sage_germline_pon_tbi)

    // sample_name=${meta.sample}
    // bam_file="$2"
    // reference_fasta="$3"
    // gnomad_vcf="$4"
    // pon_vcf="$5"
    // output_dir="${6:-.}"

    output:
    tuple val(meta), path("*.rastair.final.vcf.gz"), path("*.rastair.final.vcf.gz.tbi"), path("**"), emit: rastair_paths

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus ?: 32
    """

    usage() {
        cat <<'EOF'
        Usage:
        bash run_rastair_taps.sh <sample_name> <bam> <reference_fasta> <gnomad_vcf> <pon_vcf> [output_dir]

        Outputs:
        <sample>.rastair.raw.vcf.gz
        <sample>.rastair.ml_not_m5mc.vcf.gz
        <sample>.rastair.ml_not_m5mc.no_gnomad.vcf.gz
        <sample>.rastair.final.vcf.gz
        <sample>.rastair.cpg_only.bed.gz
EOF
    }

    sample_name=${meta.sample}
    bam_file=${bam}
    reference_fasta=${fasta}
    gnomad_vcf=${gnomAD_snv_db}
    pon_vcf=${sage_germline_pon}
    output_dir=./

    sample_prefix="\$(basename "\$sample_name")"

    for required_path in "\$bam_file" "\$reference_fasta" "\$gnomad_vcf" "\$pon_vcf"; do
        if [[ ! -f "\$required_path" ]]; then
            echo "Missing required file: \$required_path" >&2
            exit 1
        fi
    done

    require_cmd() {
        local cmd="\$1"
        if ! command -v "\$cmd" >/dev/null 2>&1; then
            echo "Required command not found on PATH: \$cmd" >&2
            exit 1
        fi
    }

    for cmd in rastair bcftools awk bgzip tabix; do
        require_cmd "\$cmd"
    done

    mkdir -p "\$output_dir"

    index_vcf() {
        local vcf_path="\$1"
        tabix -f -p vcf "\$vcf_path"
    }

    rastair_raw_vcf="\${output_dir}/\${sample_prefix}.rastair.raw.vcf.gz"
    rastair_ml_vcf="\${output_dir}/\${sample_prefix}.rastair.ml_not_m5mc.vcf.gz"
    rastair_no_gnomad_vcf="\${output_dir}/\${sample_prefix}.rastair.ml_not_m5mc.no_gnomad.vcf.gz"
    rastair_final_vcf="\${output_dir}/\${sample_prefix}.rastair.final.vcf.gz"
    rastair_cpg_only_bed="\${output_dir}/\${sample_prefix}.rastair.cpg_only.bed.gz"

    echo "Step 1/5: Running Rastair in VCF mode..."
    rastair call \\
        -@ ${threads} `# 32` \\
        -q 50 \\
        -Q 20 \\
        --ml 0.8 \\
        --vcf-format-fields GT,DP,M5mC,ML \\
        --vcf "\$rastair_raw_vcf" \\
        -r "\$reference_fasta" \\
        "\$bam_file"
    index_vcf "\$rastair_raw_vcf"

    echo "Step 2/5: Keeping Rastair records with ML but not M5mC..."
    {
    bcftools view -h "\$rastair_raw_vcf"
    bcftools view -H "\$rastair_raw_vcf" | awk -F'\\t' '\$9 ~ /(^|:)ML(:|\$)/ && \$9 !~ /(^|:)M5mC(:|\$)/'
    } | bgzip > "\$rastair_ml_vcf"
    index_vcf "\$rastair_ml_vcf"

    echo "Step 3/5: Removing Rastair calls found in gnomAD..."
    bcftools isec \\
        -C \\
        -c none \\
        -w1 \\
        -Oz \\
        -o "\$rastair_no_gnomad_vcf" \\
        "\$rastair_ml_vcf" \\
        "\$gnomad_vcf"
    index_vcf "\$rastair_no_gnomad_vcf"

    echo "Step 4/5: Removing Rastair calls found in the panel of normals..."
    bcftools isec \\
        -C \\
        -c none \\
        -w1 \\
        -Oz \\
        -o "\$rastair_final_vcf" \\
        "\$rastair_no_gnomad_vcf" \\
        "\$pon_vcf"
    index_vcf "\$rastair_final_vcf"

    echo "Step 5/5: Running Rastair in --cpg-only BED mode..."
    rastair call \\
        --cpgs-only \\
        --threads ${threads} \\
        --bed "\$rastair_cpg_only_bed" \\
        -r "\$reference_fasta" \\
        "\$bam_file"

    echo "Done."
    echo "Rastair raw VCF: \${rastair_raw_vcf}"
    echo "Rastair ML-not-M5mC VCF: \${rastair_ml_vcf}"
    echo "Rastair no-gnomAD VCF: \${rastair_no_gnomad_vcf}"
    echo "Rastair final VCF: \${rastair_final_vcf}"
    echo "Rastair CpG-only BED: \${rastair_cpg_only_bed}"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            rastair: \$(echo \$(rastair --version 2>&1) | sed 's/^rastair //' )
END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            rastair: \$(echo \$(rastair --version 2>&1) | sed 's/^rastair //' )
END_VERSIONS
    """

}

process MUTECT2_TAPS {
    tag "$meta.id"

    // conda '/opt/conda/envs/my-existing-env'
    beforeScript """
        set -a
        . /gpfs/home/hadik01/load_miniforge
        conda activate rna --stack
        conda activate vcftools --stack
        set +a
    """

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(dict)
    path(gnomAD_snv_db)
    path(gnomAD_snv_db_tbi)
    path(sage_germline_pon)
    path(sage_germline_pon_tbi)
    path(skeletal_bed)

    output:
    tuple val(meta), path("*.mutect2.final.vcf.gz"), path("*.mutect2.final.vcf.gz.tbi"), path("**"), emit: mutect2_paths

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus ?: 32
    """

    usage() {
        cat <<'EOF'
        Usage:
        bash run_mutect2_taps.sh <sample_name> <bam> <reference_fasta> <gnomad_vcf> <pon_vcf> <skeletal_bed> [output_dir]

        Outputs:
        <sample>.mutect2.raw.vcf.gz
        <sample>.mutect2.filtered.annotated.vcf.gz
        <sample>.mutect2.filtered.annotated.no_gnomad.vcf.gz
        <sample>.mutect2.filtered.annotated.no_gnomad.no_pon.vcf.gz
        <sample>.mutect2.filtered.annotated.no_gnomad.no_pon.bidir_alt_support.vcf.gz
        <sample>.mutect2.final.vcf.gz
EOF
    }

    require_cmd() {
        local cmd="\$1"
        if ! command -v "\$cmd" >/dev/null 2>&1; then
            echo "Required command not found on PATH: \$cmd" >&2
            exit 1
        fi
    }


    index_vcf() {
        local vcf_path="\$1"
        tabix -f -p vcf "\$vcf_path"
    }

    tag_mutect_bidirectional_support() {
        local input_vcf="\$1"
        local output_vcf="\$2"
        local workdir="\$3"

        mkdir -p "\$workdir"

        bcftools view -h "\$input_vcf" > "\${workdir}/header.txt"
        bcftools view -H "\$input_vcf" > "\${workdir}/body.txt"

        awk '
        BEGIN {
            has_bidir_info = 0
        }

        /^##INFO=<ID=BIDIR_ALT_SUPPORT,/ { has_bidir_info = 1 }
        /^#CHROM/ {
            if (!has_bidir_info) {
                print "##INFO=<ID=BIDIR_ALT_SUPPORT,Number=0,Type=Flag,Description=\\"ALT allele has non-zero support in both F1R2 and F2R1 for the first ALT allele\\">"
                print "##INFO=<ID=F1R2_ALT_COUNT,Number=1,Type=Integer,Description=\\"ALT count for the first ALT allele extracted from FORMAT/F1R2\\">"
                print "##INFO=<ID=F2R1_ALT_COUNT,Number=1,Type=Integer,Description=\\"ALT count for the first ALT allele extracted from FORMAT/F2R1\\">"
            }
            print
            next
        }
        { print }
        ' "\${workdir}/header.txt" > "\${workdir}/header.with_tags.txt"

        awk '
        BEGIN {
            FS = OFS = "\\t"
        }

        function get_format_index(format_string, target,    n, fields, i) {
            n = split(format_string, fields, ":")
            for (i = 1; i <= n; i++) {
                if (fields[i] == target) {
                    return i
                }
            }
            return 0
        }

        function get_first_alt_count(sample_string, idx,    n, sample_fields, allele_counts, allele_n) {
            if (idx == 0) {
                return "."
            }

            n = split(sample_string, sample_fields, ":")
            if (idx > n) {
                return "."
            }

            allele_n = split(sample_fields[idx], allele_counts, ",")
            if (allele_n < 2) {
                return "."
            }

            return allele_counts[2]
        }

        {
            f1r2_idx = get_format_index(\$9, "F1R2")
            f2r1_idx = get_format_index(\$9, "F2R1")

            f1r2_alt = get_first_alt_count(\$10, f1r2_idx)
            f2r1_alt = get_first_alt_count(\$10, f2r1_idx)

            new_info = \$8
            if (new_info == "." || new_info == "") {
                new_info = ""
            } else {
                new_info = new_info ";"
            }

            new_info = new_info "F1R2_ALT_COUNT=" f1r2_alt ";F2R1_ALT_COUNT=" f2r1_alt

            if (f1r2_alt != "." && f2r1_alt != "." && (f1r2_alt + 0) > 0 && (f2r1_alt + 0) > 0) {
                new_info = new_info ";BIDIR_ALT_SUPPORT"
            }

            \$8 = new_info
            print
        }
        ' "\${workdir}/body.txt" > "\${workdir}/body.with_tags.txt"

        cat "\${workdir}/header.with_tags.txt" "\${workdir}/body.with_tags.txt" | bgzip > "\$output_vcf"
        index_vcf "\$output_vcf"
    }

    filter_ct_ga_without_bidir() {
        local input_vcf="\$1"
        local output_vcf="\$2"
        local workdir="\$3"

        mkdir -p "\$workdir"

        bcftools view -h "\$input_vcf" > "\${workdir}/header.txt"
        bcftools view -H "\$input_vcf" > "\${workdir}/body.txt"

        awk '
        BEGIN {
            FS = OFS = "\\t"
        }

        function is_target_conversion(ref, alt) {
            return length(ref) == 1 && length(alt) == 1 && ((ref == "C" && alt == "T") || (ref == "G" && alt == "A"))
        }

        function has_bidir_flag(info) {
            return info ~ /(^|;)BIDIR_ALT_SUPPORT(;|\$)/
        }

        {
            if (is_target_conversion(\$4, \$5) && !has_bidir_flag(\$8)) {
                next
            }
            print
        }
        ' "\${workdir}/body.txt" > "\${workdir}/body.filtered.txt"

        cat "\${workdir}/header.txt" "\${workdir}/body.filtered.txt" | bgzip > "\$output_vcf"
        index_vcf "\$output_vcf"
    }


    sample_name=${meta.sample}
    bam_file=${bam}
    reference_fasta=${fasta}
    gnomad_vcf=${gnomAD_snv_db}
    pon_vcf=${sage_germline_pon}
    skeletal_bed_file=${skeletal_bed}
    output_dir=./

    sample_prefix="\$(basename "\$sample_name")"

    for required_path in "\$bam_file" "\$reference_fasta" "\$gnomad_vcf" "\$pon_vcf" "\$skeletal_bed_file"; do
        if [[ ! -f "\$required_path" ]]; then
            echo "Missing required file: \$required_path" >&2
            exit 1
        fi
    done

    for cmd in gatk bcftools awk bgzip tabix mktemp; do
        require_cmd "\$cmd"
    done

    mkdir -p "\$output_dir"

    mutect2_raw_vcf="\${output_dir}/\${sample_prefix}.mutect2.raw.vcf.gz"
    mutect2_filtered_vcf="\${output_dir}/\${sample_prefix}.mutect2.filtered.annotated.vcf.gz"
    mutect2_no_gnomad_vcf="\${output_dir}/\${sample_prefix}.mutect2.filtered.annotated.no_gnomad.vcf.gz"
    mutect2_no_gnomad_no_pon_vcf="\${output_dir}/\${sample_prefix}.mutect2.filtered.annotated.no_gnomad.no_pon.vcf.gz"
    mutect2_bidir_vcf="\${output_dir}/\${sample_prefix}.mutect2.filtered.annotated.no_gnomad.no_pon.bidir_alt_support.vcf.gz"
    mutect2_final_vcf="\${output_dir}/\${sample_prefix}.mutect2.final.vcf.gz"

    tmp_root="\$(mktemp -d "\${TMPDIR:-/tmp}/mutect2-taps.XXXXXX")"
    cleanup() {
        rm -rf "\$tmp_root"
    }
    trap cleanup EXIT

    echo "Step 1/6: Running Mutect2..."
    gatk Mutect2 \\
        -R "\$reference_fasta" \\
        -I "\$bam_file" \\
        -tumor "\$sample_name" \\
        -L "\$skeletal_bed_file" \\
        --panel-of-normals "\$pon_vcf" \\
        -O "\$mutect2_raw_vcf"
    index_vcf "\$mutect2_raw_vcf"

    echo "Step 2/6: Running FilterMutectCalls..."
    gatk FilterMutectCalls \\
        -R "\$reference_fasta" \\
        -V "\$mutect2_raw_vcf" \\
        --stats "\${mutect2_raw_vcf}.stats" \\
        -O "\$mutect2_filtered_vcf"
    index_vcf "\$mutect2_filtered_vcf"

    echo "Step 3/6: Removing Mutect2 calls found in gnomAD..."
    bcftools isec \\
        -C \\
        -c none \\
        -w1 \\
        -Oz \\
        -o "\$mutect2_no_gnomad_vcf" \\
        "\$mutect2_filtered_vcf" \\
        "\$gnomad_vcf"
    index_vcf "\$mutect2_no_gnomad_vcf"

    echo "Step 4/6: Removing Mutect2 calls found in the panel of normals..."
    bcftools isec \\
        -C \\
        -c none \\
        -w1 \\
        -Oz \\
        -o "\$mutect2_no_gnomad_no_pon_vcf" \\
        "\$mutect2_no_gnomad_vcf" \\
        "\$pon_vcf"
    index_vcf "\$mutect2_no_gnomad_no_pon_vcf"

    echo "Step 5/6: Tagging Mutect2 calls with bidirectional ALT support..."
    tag_mutect_bidirectional_support \\
        "\$mutect2_no_gnomad_no_pon_vcf" \\
        "\$mutect2_bidir_vcf" \\
        "\${tmp_root}/tag_mutect"

    echo "Step 6/6: Filtering Mutect2 C>T and G>A SNVs that lack bidirectional ALT support..."
    filter_ct_ga_without_bidir \\
        "\$mutect2_bidir_vcf" \\
        "\$mutect2_final_vcf" \\
        "\${tmp_root}/filter_mutect"
    
    index_vcf "\$mutect2_final_vcf"

    echo "Done."
    echo "Mutect2 raw VCF: \${mutect2_raw_vcf}"
    echo "Mutect2 filtered VCF: \${mutect2_filtered_vcf}"
    echo "Mutect2 no-gnomAD VCF: \${mutect2_no_gnomad_vcf}"
    echo "Mutect2 no-gnomAD/no-PoN VCF: \${mutect2_no_gnomad_no_pon_vcf}"
    echo "Mutect2 bidirectional-tagged VCF: \${mutect2_bidir_vcf}"
    echo "Mutect2 final VCF: \${mutect2_final_vcf}"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*GATK) v//' )
            bcftools: \$(echo \$(bcftools --version 2>&1) | head -1 | sed 's/^bcftools //' )
END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*GATK) v//' )
            bcftools: \$(echo \$(bcftools --version 2>&1) | head -1 | sed 's/^bcftools //' )
END_VERSIONS
    """

}

process COMBINE_TAPS_VARIANT_CALLS {
    tag "$meta.id"

    // conda '/opt/conda/envs/my-existing-env'
    beforeScript """
        set -a
        . /gpfs/home/hadik01/load_miniforge
        conda activate rna --stack
        conda activate vcftools --stack
        set +a
    """

    input:
    tuple val(meta), path(mutect2_final_vcf), path(mutect2_final_vcf_tbi), path(rastair_final_vcf), path(rastair_final_vcf_tbi)

    output:
    tuple val(meta), path("*.combined.final.vcf.gz"), path("*.combined.final.vcf.gz.tbi"), path("**"), emit: combined_paths

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus ?: 32
    """

    usage() {
        cat <<'EOF'
        Usage:
        bash combine_mutect2_rastair_taps.sh <sample_name> <mutect2_final_vcf> <rastair_final_vcf> [output_dir]

        Outputs:
        <sample>.combined.final.vcf.gz
EOF
    }

    require_cmd() {
        local cmd="\$1"
        if ! command -v "\$cmd" >/dev/null 2>&1; then
            echo "Required command not found on PATH: \$cmd" >&2
            exit 1
        fi
    }

    index_vcf() {
        local vcf_path="\$1"
        tabix -f -p vcf "\$vcf_path"
    }

    combine_union_vcfs() {
        local vcf1="\$1"
        local vcf2="\$2"
        local output_vcf="\$3"
        local merged_sample_name="\$4"
        local workdir="\$5"

        mkdir -p "\$workdir"

        local source1
        local source2
        local sample1_old
        local sample2_old

        source1="\$(basename "\$vcf1" .vcf.gz)"
        source2="\$(basename "\$vcf2" .vcf.gz)"

        sample1_old="\$(bcftools query -l "\$vcf1" | awk 'NR == 1 { print; exit }')"
        sample2_old="\$(bcftools query -l "\$vcf2" | awk 'NR == 1 { print; exit }')"

        if [[ -z "\$sample1_old" || -z "\$sample2_old" ]]; then
            echo "Unable to extract sample names from input VCFs" >&2
            exit 1
        fi

        printf '%s\\t%s\\n' "\$sample1_old" "\$merged_sample_name" > "\${workdir}/sample1.map"
        printf '%s\\t%s\\n' "\$sample2_old" "\$merged_sample_name" > "\${workdir}/sample2.map"

        bcftools reheader -s "\${workdir}/sample1.map" -o "\${workdir}/input1.reheader.vcf.gz" "\$vcf1"
        bcftools reheader -s "\${workdir}/sample2.map" -o "\${workdir}/input2.reheader.vcf.gz" "\$vcf2"
        index_vcf "\${workdir}/input1.reheader.vcf.gz"
        index_vcf "\${workdir}/input2.reheader.vcf.gz"

        bcftools view -h "\${workdir}/input1.reheader.vcf.gz" > "\${workdir}/header1.txt"
        bcftools view -h "\${workdir}/input2.reheader.vcf.gz" > "\${workdir}/header2.txt"

        awk \\
            -v sample_name="\$merged_sample_name" \\
            '
            function meta_key(line,    m) {
                if (match(line, /^##([^=]+)=<ID=([^,>]+)/, m)) {
                    return m[1] SUBSEP m[2]
                }
                return line
            }

            function remember(line, key) {
                key = meta_key(line)
                if (!(key in seen)) {
                    seen[key] = line
                    order[++n] = key
                } else if (seen[key] != line) {
                    printf("Header conflict for key %s; keeping first definition\\n", key) > "/dev/stderr"
                }
            }

            /^##fileformat=/ { next }
            /^#CHROM/ { next }
            { remember(\$0) }

            END {
                print "##fileformat=VCFv4.2"
                print "##INFO=<ID=SOURCES,Number=.,Type=String,Description=\\"Comma-separated list of input VCF sources contributing this record\\">"
                print "##INFO=<ID=SOURCE_COUNT,Number=1,Type=Integer,Description=\\"Number of input VCF sources contributing this record\\">"
                for (i = 1; i <= n; i++) {
                    print seen[order[i]]
                }
                print "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t" sample_name
            }
            ' "\${workdir}/header1.txt" "\${workdir}/header2.txt" > "\${workdir}/merged.header.vcf"

        bcftools view -H "\${workdir}/input1.reheader.vcf.gz" > "\${workdir}/body1.tsv"
        bcftools view -H "\${workdir}/input2.reheader.vcf.gz" > "\${workdir}/body2.tsv"

        awk \\
            -v source1="\$source1" \\
            -v source2="\$source2" \\
            '
            BEGIN {
                FS = OFS = "\\t"
            }

            function qual_num(q) {
                return (q == "." ? -1 : q + 0)
            }

            function filter_rank(f) {
                return (f == "PASS" ? 1 : 0)
            }

            function add_source(key, src, source_key) {
                source_key = key SUBSEP src
                if (!(source_key in source_seen)) {
                    source_seen[source_key] = 1
                    if (sources[key] == "") {
                        sources[key] = src
                    } else {
                        sources[key] = sources[key] "," src
                    }
                    source_count[key]++
                }
            }

            function should_replace(key, new_filter, new_qual) {
                old_filter = kept_filter[key]
                old_qual = kept_qual[key]

                if (filter_rank(new_filter) > filter_rank(old_filter)) {
                    return 1
                }
                if (filter_rank(new_filter) < filter_rank(old_filter)) {
                    return 0
                }
                if (qual_num(new_qual) > qual_num(old_qual)) {
                    return 1
                }
                return 0
            }

            function process_record(src, key) {
                key = \$1 SUBSEP \$2 SUBSEP \$4 SUBSEP \$5

                if (!(key in seen)) {
                    seen[key] = 1
                    order[++n] = key
                    record[key] = \$0
                    kept_filter[key] = \$7
                    kept_qual[key] = \$6
                } else if (should_replace(key, \$7, \$6)) {
                    record[key] = \$0
                    kept_filter[key] = \$7
                    kept_qual[key] = \$6
                }

                add_source(key, src)
            }

            FNR == NR {
                process_record(source1)
                next
            }

            {
                process_record(source2)
            }

            END {
                for (i = 1; i <= n; i++) {
                    key = order[i]
                    split(record[key], f, FS)

                    if (f[8] == "." || f[8] == "") {
                        f[8] = "SOURCES=" sources[key] ";SOURCE_COUNT=" source_count[key]
                    } else {
                        f[8] = f[8] ";SOURCES=" sources[key] ";SOURCE_COUNT=" source_count[key]
                    }

                    print f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10]
                }
            }
            ' "\${workdir}/body1.tsv" "\${workdir}/body2.tsv" > "\${workdir}/merged.body.vcf"

        cat "\${workdir}/merged.header.vcf" "\${workdir}/merged.body.vcf" > "\${workdir}/merged.unsorted.vcf"

        bcftools sort -Oz -o "\$output_vcf" "\${workdir}/merged.unsorted.vcf"
        index_vcf "\$output_vcf"
    }


    sample_name=${meta.sample}
    mutect2_final_vcf_file=${mutect2_final_vcf}
    rastair_final_vcf_file=${rastair_final_vcf}
    output_dir=./

    sample_prefix="\$(basename "\$sample_name")"

    for required_path in "\$mutect2_final_vcf_file" "\$rastair_final_vcf_file"; do
        if [[ ! -f "\$required_path" ]]; then
            echo "Missing input VCF: \$required_path" >&2
            exit 1
        fi
    done

    for cmd in bcftools awk bgzip tabix gzip mktemp; do
        require_cmd "\$cmd"
    done

    mkdir -p "\$output_dir"

    combined_final_vcf="\${output_dir}/\${sample_prefix}.combined.final.vcf.gz"

    tmp_root="\$(mktemp -d "\${TMPDIR:-/tmp}/combine-mutect2-rastair.XXXXXX")"
    cleanup() {
        rm -rf "\$tmp_root"
    }
    trap cleanup EXIT

    echo "Combining final Mutect2 and Rastair VCFs..."
    combine_union_vcfs \\
        "\$mutect2_final_vcf_file" \\
        "\$rastair_final_vcf_file" \\
        "\$combined_final_vcf" \\
        "\$sample_name" \\
        "\${tmp_root}/combine_union"

    echo "Done."
    echo "Combined final VCF: \${combined_final_vcf}"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | head -1 | sed 's/^bcftools //' )
END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | head -1 | sed 's/^bcftools //' )
END_VERSIONS
    """

}

process TUMOR_ONLY_FILTER {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/utils:0.0.2':
        'mskilab/utils:0.0.2' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    path(dbsnp)
    path(dbsnp_tbi)
    path(gnomAD_snv_db)
    path(gnomAD_snv_db_tbi)
    path(sage_germ_pon)
    path(sage_germ_pon_tbi)
    path(mills_gold_indel)
    path(mills_gold_indel_tbi)

    output:
    tuple val(meta), path("*.tumoronly.vcf.gz"), path("*.tumoronly.vcf.gz.tbi"), emit: tumor_only_filtered_vcf
    path "versions.yml", emit: versions, optional:true

    when:
    task.ext.when == null || task.ext.when
    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def output      = "${meta.id}.tumoronly.vcf.gz"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir -p filter_tumoronly/

    if [[ -e ${vcf} ]]; then
        echo "Filtering the tumor-only vcf by gnomAD..."
        bcftools isec -C -O z -o ${meta.id}.nognomAD.vcf.gz -p ./ ${vcf} ${gnomAD_snv_db} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.nognomAD.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.nognomAD.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by gnomAD. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.nognomAD.vcf.gz ]]; then
        echo "Now filtering by SAGE Germline PON ..."
        bcftools isec -C -O z -o ${meta.id}.nognomAD.noPON.vcf.gz -p ./ filter_tumoronly/${meta.id}.nognomAD.vcf.gz \\
        ${sage_germ_pon} && mv 0000.vcf.gz filter_tumoronly/${meta.id}.nognomAD.noPON.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.nognomAD.noPON.vcf.gz.tbi
    else
        echo "Cannot perform SNV filtering by SAGE Germline PON. Exiting..."
        exit 1;
    fi

    if [[ -e filter_tumoronly/${meta.id}.nognomAD.noPON.vcf.gz ]]; then
        echo "Now filtering by Mills and Gold Standard Indels ..."
        bcftools isec -C -O z -o ${meta.id}.nognomAD.noPON.noMGIndel.vcf.gz -p ./ \\
        filter_tumoronly/${meta.id}.nognomAD.noPON.vcf.gz ${mills_gold_indel} && \\
        mv 0000.vcf.gz filter_tumoronly/${meta.id}.nognomAD.noPON.noMGIndel.vcf.gz && \\
        mv 0000.vcf.gz.tbi filter_tumoronly/${meta.id}.nognomAD.noPON.noMGIndel.vcf.gz.tbi

        cp filter_tumoronly/${meta.id}.nognomAD.noPON.noMGIndel.vcf.gz ./${meta.id}.tumoronly.vcf.gz && \\
        cp filter_tumoronly/${meta.id}.nognomAD.noPON.noMGIndel.vcf.gz.tbi ./${meta.id}.tumoronly.vcf.gz.tbi

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
    touch ${meta.id}.tumoronly.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: ${VERSION}
END_VERSIONS
    """
}