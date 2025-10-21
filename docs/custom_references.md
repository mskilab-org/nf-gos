# Using Custom Reference Genomes

This guide explains how to use custom reference genomes with nf-casereports, including non-canonical references, custom builds, or alternative assemblies.

## Table of Contents
- [Quick Start](#quick-start)
- [Understanding Reference Files](#understanding-reference-files)
- [Three Approaches](#three-approaches)
- [Common Use Cases](#common-use-cases)
- [Troubleshooting](#troubleshooting)
- [Performance & Resources](#performance--resources)

---

## Quick Start

**Most common scenario:** You have a custom FASTA file and want the pipeline to handle everything else:

```bash
nextflow run mskilab-org/nf-casereports \
   -profile singularity,nygc \
   --input samples.csv \
   --outdir results \
   --genome GATK.GRCh37 \
   --fasta /path/to/custom_genome.fasta \
   --regenerate_reference_indices true
```

That's it! The pipeline will:
1. Generate FASTA index (.fai)
2. Generate sequence dictionary (.dict)
3. Generate BWA index (for GRIDSS, SvABA)
4. Save all indices to `results/reference/`
5. Run the full analysis pipeline

---

## Understanding Reference Files

### Files Generated Automatically

When you use `--regenerate_reference_indices true`, these files are created from your FASTA:

| File | Purpose | Tool | Time (hg19) |
|------|---------|------|-------------|
| `.fai` | FASTA index | samtools faidx | ~5 min |
| `.dict` | Sequence dictionary | GATK | ~5 min |
| `bwa/` | BWA alignment index | bwa index | ~60-90 min |

**Total time:** ~90-120 minutes for human genome

### Files You Must Provide

These annotation files **cannot** be derived from FASTA and must be provided if needed:

- `--dbsnp` - dbSNP variants
- `--known_indels` - Known indels for BQSR
- `--known_snps` - Known SNPs for BQSR
- Tool-specific resources (see [Tool Resources](#tool-resources))

---

## Three Approaches

### Approach 1: Auto-Generate Everything (Recommended)

**When to use:** You have a FASTA file and want the simplest workflow.

```bash
nextflow run mskilab-org/nf-casereports \
   --input samples.csv \
   --outdir results \
   --fasta /refs/custom.fasta \
   --regenerate_reference_indices true \
   -profile singularity,nygc
```

**Pros:**
- ✅ Simplest workflow
- ✅ Guarantees matching indices
- ✅ Indices saved for reuse

**Cons:**
- ⏱️ Takes 90-120 minutes on first run

---

### Approach 2: Provide Pre-Built Indices

**When to use:** You already have indices or want to avoid regeneration time.

```bash
nextflow run mskilab-org/nf-casereports \
   --input samples.csv \
   --outdir results \
   --fasta /refs/genome.fasta \
   --fasta_fai /refs/genome.fasta.fai \
   --dict /refs/genome.dict \
   --bwa /refs/bwa_index/ \
   -profile singularity,nygc
```

**Pros:**
- ✅ Skips index generation
- ✅ Faster startup
- ✅ Use existing infrastructure

**Cons:**
- ⚠️ Must ensure indices match FASTA
- ⚠️ Must manage multiple files

---

### Approach 3: Build Reference Bundle Only

**When to use:** Pre-build indices once, reuse many times.

**Step 1: Build indices**
```bash
nextflow run mskilab-org/nf-casereports \
   --outdir reference_bundle \
   --fasta /refs/genome.fasta \
   --regenerate_reference_indices true \
   -profile singularity,nygc
```

Output:
```
reference_bundle/
└── reference/
    ├── dict/genome.dict
    ├── fai/genome.fasta.fai
    └── bwa/
        ├── genome.amb
        ├── genome.ann
        ├── genome.bwt
        ├── genome.pac
        └── genome.sa
```

**Step 2: Use in multiple runs**
```bash
nextflow run mskilab-org/nf-casereports \
   --input samples_cohort1.csv \
   --fasta /refs/genome.fasta \
   --fasta_fai reference_bundle/reference/fai/genome.fasta.fai \
   --dict reference_bundle/reference/dict/genome.dict \
   --bwa reference_bundle/reference/bwa/ \
   -profile singularity,nygc
```

---

## Common Use Cases

### Use Case 1: hg19 with Different Contigs

**Scenario:** You need hg19 but with only standard chromosomes (no alt contigs).

```bash
nextflow run mskilab-org/nf-casereports \
   --input samples.csv \
   --outdir results \
   --genome GATK.GRCh37 \
   --fasta /refs/hg19_no_alt.fasta \
   --regenerate_reference_indices true \
   -profile singularity,nygc
```

**Why `--regenerate_reference_indices`?**
- Without it: Pipeline uses old GATK.GRCh37 indices (with alt contigs) ❌
- With it: Pipeline generates new indices matching your FASTA ✅

---

### Use Case 2: CHM13 (T2T Assembly)

**Scenario:** Using the telomere-to-telomere CHM13 reference.

```bash
nextflow run mskilab-org/nf-casereports \
   --input samples.csv \
   --outdir results_chm13 \
   --fasta /refs/chm13v2.0.fasta \
   --regenerate_reference_indices true \
   -profile singularity,nygc
```

**Note:** CHM13 lacks many annotation resources. You'll need to provide or omit:
- `--dbsnp` (if available for CHM13)
- May need to skip certain annotation steps

---

### Use Case 3: Masked Reference

**Scenario:** You want to mask repetitive regions or problematic sequences.

```bash
# Step 1: Create masked FASTA
bedtools maskfasta -fi hg38.fasta -bed repeats.bed -fo hg38_masked.fasta

# Step 2: Run pipeline
nextflow run mskilab-org/nf-casereports \
   --input samples.csv \
   --fasta hg38_masked.fasta \
   --regenerate_reference_indices true \
   -profile singularity,nygc
```

---

### Use Case 4: Mouse or Non-Human Reference

**Scenario:** Running on mouse (mm10/GRCm38) or other species.

```bash
nextflow run mskilab-org/nf-casereports \
   --input mouse_samples.csv \
   --outdir results_mm10 \
   --fasta /refs/mm10.fasta \
   --regenerate_reference_indices true \
   -profile singularity,nygc
```

**Important:** Most tool-specific resources are human-only. You may need to:
- Provide species-specific annotation files
- Skip tools that lack non-human resources
- Use `--skip_tools` to disable incompatible tools

---

## Troubleshooting

### Error: "GRIDSS failed with reference mismatch"

**Cause:** Using custom FASTA without regenerating indices.

**Solution:**
```bash
# Add the regenerate flag
--regenerate_reference_indices true
```

---

### Error: "BWA index not found"

**Cause:** Pipeline couldn't find or generate BWA index.

**Check:**
1. Is `--regenerate_reference_indices true` set?
2. Is `--bwa /path/` provided if using pre-built?
3. Did index generation fail? Check `work/` logs.

---

### Warning: "Dict contigs don't match FASTA"

**Cause:** Provided dict doesn't match your FASTA.

**Solution:** Either:
1. Use `--regenerate_reference_indices true` to rebuild
2. Or provide matching dict: `--dict /path/to/matching.dict`

---

### Slow: "BWA indexing taking hours"

**Expected:** BWA indexing is single-threaded and takes ~60-90 minutes for human genome.

**Optimization:**
- Build indices once, reuse (Approach 3)
- Or use pre-built indices (Approach 2)

---

## Performance & Resources

### Index Generation Times (Human Genome)

| Process | Cores | Memory | Time (hg19) | Can parallelize? |
|---------|-------|--------|-------------|------------------|
| samtools faidx | 1 | 1 GB | ~5 min | No |
| GATK dict | 1 | 2 GB | ~5 min | No |
| BWA index | 1 | 8 GB | ~60-90 min | No |

**Total:** ~90-120 minutes, ~8 GB peak memory

### Disk Space Requirements

| Reference | FASTA | BWA Index | Total |
|-----------|-------|-----------|-------|
| hg19 | 3 GB | 5-6 GB | ~9 GB |
| hg38 | 3.2 GB | 5-6 GB | ~9 GB |
| CHM13 | 3.1 GB | 5-6 GB | ~9 GB |
| mm10 | 2.7 GB | 4-5 GB | ~7 GB |

---

## Tool Resources

Some tools require genome-specific resources beyond the basic indices:

### GRIDSS
- `--blacklist_gridss` - Blacklist regions (centromeres, telomeres)
- `--pon_gridss` - Panel of normals

### SvABA
- `--indel_mask` - Known problematic indel regions
- `--germ_sv_db` - Germline SV database

### PURPLE
- `--gc_profile` - GC content profile
- `--diploid_bed` - Expected diploid regions

### See Also
- [Tool-specific parameters](tool_parameters.md)
- [Creating custom annotation files](custom_annotations.md)

---

## Best Practices

1. **Always use `--regenerate_reference_indices true` when providing custom FASTA**
2. **Build reference bundle once, reuse many times** (Approach 3)
3. **Test with small dataset first** before full cohort
4. **Document your reference** - keep notes on what contigs/masking you used
5. **Archive generated indices** - save `${outdir}/reference/` for future runs

---

## Questions?

- Check [GitHub Issues](https://github.com/mskilab-org/nf-casereports/issues)
- Join [mskilab Slack](https://mskilab.org/slack)
- Read [nf-core documentation](https://nf-co.re/docs/)
