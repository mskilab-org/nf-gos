process SV_CHIMERA_FILTER_ESVEE {
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
    def out_vcf = vcf.getName().replaceFirst(/\.vcf(\.gz|\.bgz)?$/, '.ffpe_filtered.vcf.gz')
    def tumor_id = meta.tumor_id ?: meta.sample ?: ''
    def normal_id = meta.normal_id ?: ''
    """
    # Resolve the tumor by exact VCF sample name, never by caller sample order.
    bcftools query -l ${vcf} > vcf.samples
    mapfile -t samples < vcf.samples
    n_samples=\${#samples[@]}
    tumor_id="${tumor_id}"
    tumor_idx=-1
    if [ -n "\${tumor_id}" ]; then
        for idx in "\${!samples[@]}"; do
            if [ "\${samples[\${idx}]}" = "\${tumor_id}" ]; then
                if [ "\${tumor_idx}" -ne -1 ]; then
                    echo "ERROR: ESVEE tumor sample '\${tumor_id}' matches multiple VCF samples." >&2
                    exit 1
                fi
                tumor_idx=\${idx}
            fi
        done
        if [ "\${tumor_idx}" -eq -1 ]; then
            echo "ERROR: ESVEE tumor sample '\${tumor_id}' does not match any VCF sample." >&2
            exit 1
        fi
    elif [ "\${n_samples}" -eq 1 ]; then
        tumor_idx=0
    else
        echo "ERROR: ESVEE filtering requires meta.tumor_id or meta.sample unless the VCF has exactly one sample; found \${n_samples} samples." >&2
        exit 1
    fi
    echo "samples=\${n_samples} tumor FORMAT index=\${tumor_idx}"

    # VF is ESVEE's total variant-fragment support. ESVEE site QUAL is not on
    # the GRIDSS FORMAT/QUAL scale: 30 is ESVEE's native WGS minQual threshold.
    # ASMID is not a GRIDSS INFO/AS assembly-support count, so no assembly term
    # is applied. Append advisory tags, preserving existing filters and all records.
    bcftools filter \\
        --soft-filter FFPE_SUPPORT \\
        --mode + \\
        -e "FORMAT/VF[\${tumor_idx}] <= 6 || QUAL < 30" \\
        -Oz -o support_tagged.vcf.gz \\
        ${vcf}

    bcftools index --tbi support_tagged.vcf.gz

    # FFPE_GEOM_CHIMERA is also advisory; use ESVEE VF for normal evidence.
    normal_arg=""
    if [ -n "${normal_id}" ] && [ "\${n_samples}" -gt 1 ]; then
        normal_arg="--normal-id ${normal_id}"
    fi

    python \${NEXTFLOW_BIN_DIR}/sv_ffpe_geometry_tag.py \\
        support_tagged.vcf.gz \\
        ${out_vcf} \\
        --caller esvee \\
        \${normal_arg} \\
        ${args}

    bcftools index --tbi ${out_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed -n 's/^bcftools //p' | head -1)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def out_vcf = vcf.getName().replaceFirst(/\.vcf(\.gz|\.bgz)?$/, '.ffpe_filtered.vcf.gz')
    """
    touch ${out_vcf}
    touch ${out_vcf}.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            bcftools: \$(echo \$(bcftools --version 2>&1) | sed -n 's/^bcftools //p' | head -1)
            samtools: \$(echo \$(samtools version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
