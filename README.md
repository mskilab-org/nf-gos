# NF-CaseReports (Nextflow - Case Reports Pipeline)

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**mskilab-org/nf-casereports** is a bioinformatics pipeline from [`mskilab-org`](https://www.mskilab.org/) for running [`JaBbA`](https://github.com/mskilab-org/JaBbA/), our algorithm for MIP based joint inference of copy number and rearrangement state in cancer whole genome sequence data. This pipeline runs all the pre-requisite tools (among others) and generates the necessary inputs for running JaBbA and loading into [gOS](https://github.com/mskilab-org/gOS), our clinical front-end. It is designed to take paired tumor-normal samples or tumor-only samples as input.

## Workflow Summary:
1. Align to Reference Genome (currently supports `BWA-MEM`, `BWA-MEM2`, and GPU accelerated `fq2bam`).
2. Quality Control (using [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`Picard CollectWGSMetrics`](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard), [`Picard CollectMultipleMetrics`](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard), and [`GATK4 EstimateLibraryComplexity`](https://gatk.broadinstitute.org/hc/en-us/articles/360037428891-EstimateLibraryComplexity-Picard))
4. Mark Duplicates (using [`GATK MarkDuplicates`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard))
5. Base recalibration (using [`GATK BaseRecalibrator`](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator))
6. Apply BQSR (using [`GATK ApplyBQSR`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR))
7. Perform structural variant calling (using [`GRIDSS`](https://github.com/PapenfussLab/gridss))
8. Perform pileups (using [`AMBER`](https://github.com/hartwigmedical/hmftools/blob/master/amber/README.md))
9. Generate raw coverages and correct for GC & Mappability bias (using [`fragCounter`](https://github.com/mskilab-org/fragCounter))
10. Remove biological and technical noise from coverage data. (using [`Dryclean`](https://github.com/mskilab-org/dryclean))
11. Perform segmentation using tumor/normal ratios of corrected read counts, (using the `CBS` (circular binary segmentation) algorithm)
12. Purity & ploidy estimation (using [`PURPLE`](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md))
13. Junction copy number estimation and event calling (using [`JaBbA`](https://github.com/mskilab-org/JaBbA/)
14. Call SNVs and indels (using [`SAGE`](https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md))
15. Annotate variants (using [`Snpeff`](https://pcingola.github.io/SnpEff/))
16. Assign mutational signatures (using [`SigProfiler`](https://github.com/AlexandrovLab/SigProfilerAssignment/))
17. Detect HRD (Homologous Recombination Deficiency) (using [`HRDetect`](https://github.com/Nik-Zainal-Group/signature.tools.lib))


## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

### Setting up the ***samplesheet.csv*** file for input:

You need to create a samplesheet with information regarding the samples you
want to run the pipeline on. You need to specify the path of your
**samplesheet** using the `--input` flag to specify the location. Make sure the
input `samplesheet.csv` file is a *comma-separated* file and contains the
headers discussed below. *It is highly recommended to provide the **absolute
path** for inputs inside the samplesheet rather than relative paths.*

For paired tumor-normal samples, use the same `patient` ID, but different
`sample` names. Indicate their respective tumor/normal `status`, where **1** in
the `status` field indicates a tumor sample, and **0** indicates a normal
sample. You may pass multiple `sample` IDs per patient, `nf-casereports` will
consider them as separate samples belonging to the same patient and output the
results accordingly.

Specify the desired output root directory using the `--outdir` flag.

The input samplesheet should look like this:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
TCXX49,XX,0,TCXX49_N,lane_1,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
```

Each row represents a pair of fastq files (paired end) for a single sample (in
this case a normal sample, status: 0). After the input file is ready, you can
run the pipeline using:

```bash
nextflow run mskilab-org/nf-jabba \
   -profile <docker|singularity|institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --genome <GATK.GRCh37/GATK.GRCh38>
```
> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow [`-params-file`](https://www.nextflow.io/blog/2020/cli-docs-release.html) option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

### Discussion of expected fields in input file and expected inputs for each `--step`

A typical sample sheet can populate with all or some of the column names as
shown below. The pipeline will use the information provided in the samplesheet
to parsimoniously run the steps of the pipeline to generate all remaining
outputs.

**N.B You do not need to supply all the columns in the table below. The table represents all the possible inputs that can be passed. If you are starting from BAMs just pass `bam` and `bai` columns. If you are starting from FASTQs, pass `fastq1` (and `fastq2` for paired reads). If you have already generated other outputs, you may pass them as well to prevent the pipeline from running tools for which you already have outputs.**

| Column Name         | Description                                                                                                                                  |
|---------------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| patient             | (required) Patient or Sample ID. This should differentiate each patient/sample. *Note*: Each patient can have multiple sample names.          |
| sample              | (required) Sample ID for each Patient. Should differentiate between tumor and normal (e.g `sample1_t` vs. `sample1_n`). Sample IDs should be unique to Patient IDs |
| lane                | If starting with FASTQ files, and if there are multiple lanes for each sample for each patient, mention lane name.                            |
| sex                 | If known, please provide the sex for the patient. For instance if **Male** type XY, else if **Female** type XX, otherwise put NA.             |
| status              | (required) This should indicate if your sample is **tumor** or **normal**. For **normal**, write 0, and for **tumor**, write 1.              |
| fastq_1             | Full path to FASTQ file read 1. The extension should be `.fastq.gz` or `.fq.gz`.                                                              |
| fastq_2             | Full path to FASTQ file read 2 (if paired reads). The extension should be `.fastq.gz` or `.fq.gz`.                                            |
| bam                 | Full path to BAM file. The extension should be `.bam`.                                                                                        |
| bai                 | Full path to BAM index file. The extension should be `.bam.bai`.                                                                              |
| hets                | Full path to sites.txt file.                                                                                                                  |
| amber_dir           | Full path to AMBER output directory.                                                                                                          |
| frag_cov            | Full path to the fragCounter coverage file.                                                                                                   |
| dryclean_cov        | Full path to the Dryclean corrected coverage file.                                                                                            |
| ploidy              | Ploidies for each sample.                                                                                                                     |
| seg                 | Full path to the CBS segmented file.                                                                                                          |
| nseg                | Full path to the CBS segmented file for normal samples.                                                                                       |
| vcf                 | Full path to the GRIDSS VCF file.                                                                                                             |
| vcf_tbi             | Full path to the GRIDSS VCF index file.                                                                                                       |
| jabba_rds           | Full path to the JaBbA RDS (`jabba.simple.rds`) file.                                                                                         |
| jabba_gg            | Full path to the JaBbA gGraph (`jabba.gg.rds`) file.                                                                                          |
| ni_balanced_gg      | Full path to the non-integer balanced gGraph (`non_integer.balanced.gg.rds`) file.                                                            |
| lp_phased_gg        | Full path to the LP phased gGraph (`lp_phased.balanced.gg.rds`) file.                                                                         |
| events              | Full path to the events file.                                                                                                                 |
| fusions             | Full path to the fusions file.                                                                                                                |
| snv_somatic_vcf     | Full path to the somatic SNV VCF file.                                                                                                        |
| snv_somatic_tbi     | Full path to the somatic SNV VCF index file.                                                                                                  |
| snv_germline_vcf    | Full path to the germline SNV VCF file.                                                                                                       |
| snv_germline_tbi    | Full path to the germline SNV VCF index file.                                                                                                 |
| variant_somatic_ann | Full path to the somatic SNV annotated VCF file.                                                                                              |
| variant_somatic_bcf | Full path to the somatic SNV BCF file.                                                                                                        |
| variant_germline_ann| Full path to the germline SNV annotated VCF file.                                                                                             |
| variant_germline_bcf| Full path to the germline SNV BCF file.                                                                                                       |
| snv_multiplicity    | Full path to the SNV multiplicity file.                                                                                                       |
| sbs_signatures      | Full path to the SBS signatures file.                                                                                                         |
| indel_signatures    | Full path to the indel signatures file.                                                                                                       |
| signatures_matrix   | Full path to the signatures matrix file.                                                                                                      |
| hrdetect            | Full path to the HRDetect file.                                                                                                               |

## Tumor-Only Samples

For tumor-only samples, simply add the flag `--tumor_only true` to the nextflow command. The pipeline will then run in tumor-only mode.

For more information regarding the pipeline usage and the inputs necesaary for each step, please follow the [Usage](docs//usage.md) documentation.

### Helpful Core Nextflow Commands:

#### `-resume`
If a process of the pipeline fails or is interrupted at some point, Nextflow can resume from that point without having to start over from the beginning. You must specify this in the `CLI` or on the `command-line` when restarting a pipeline. You can also supply a run name to resume a specific run using: `-resume` [run-name]. Use the `nextflow log` command to show previous run names.

#### `-profile`
Use this parameter for choosing a configuration profile. Profiles contain configuration presets for different computing environments.

Several generic profiles have been provided by default which instruct the pipeline to use software packaged using different methods. You can use this option to run the pipeline via containers (singularity/Docker) (**highly recommended**)

#### `-c`
You can mention custom configuration scripts to run the pipeline with using the `-c` flag and providing a path to the `.config` file. This is advised when you want to submit processes into an executor like `slurm` or `LSF`.

#### `-bg`
The Nextflow `-bg` flag launches the Nextflow pipeline as a background process. This allows you to detach or exit your terminal without interrupting the run. A log of the run will be saved inside a file upon completion. You can also use `screen` or `tmux` sessions to persist runs.

## Containers:
Every module in the pipeline has been containerized. Some modules are partially modified versions of [nf-core/modules](https://nf-co.re/modules), these modules use nf-core containers. Modules that use our lab packages and scripts were containerized into Docker images. These images can be found on our [DockerHub](https://hub.docker.com/repositories/mskilab).

> **Warning:**
> JaBbA depends on CPLEX MIP Optimizer to work. Because CPLEX is a proprietary software, it isn't included in the image and needs to be installed by the user.
> To add CPLEX:
>  1. Download CPLEX (Linux x86-64). (You may need to use the HTTP method.)
>  2. Pull image and run the container using:
> ```
> docker pull mskilab/jabba:latest
> docker run -it --rm --platform linux/amd64 mskilab/jabba:latest
> ```
>  3. Copy CPLEX binary into the container: docker cp /PATH/TO/DOWNLOADED_CPLEX.bin CONTAINER_ID:/opt/cplex_studio
>  4. Install CPLEX: /opt/cplex_studio/DOWNLOADED_CPLEX.bin (If you get a Permission denied error, run
>  chmod 777 /PATH/TO/DOWNLOADED_CPLEX.bin before copying it into the container.)
>  5. When prompted for an installation path, type /opt/cplex. This is what the CPLEX_DIR environmental variable is set to.
>  6. Save changes to a new image for future use:
>           Exit container (type exit or press Ctrl-D)
>           Run docker commit CONTAINER_ID NEW_IMAGE_ID


## Debugging any step/process:

To debug any step or process that failed, first check your current `execution_trace*.txt` file inside the `<outdir>/pipeline_info/` folder. There you'll find a `hash` number for that process. You can use that `hash` number to locate that process's working directory. This directory will contain multiple `.command.*` files that correspond to your run and contain valuable information that can help you debug your error. You can also run the `.command.sh` script to do a manual, isolated execution of the offending process for quick testing.

## Credits

`nf-casereports` was written by [`Shihab Dider`](https://github.com/shihabdider) and [`Tanubrata Dey`](https://github.com/tanubrata) and at the Perlmutter Cancer Center and the New York Genome Center.

We thank the following people for their extensive guidance in the development of this pipeline:
- [Marcin Imielinski](https://github.com/imielinski)
- [Joel Rosiene](https://github.com/jrosiene)

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **Most large structural variants in cancer genomes can be detected without long reads.**
> Choo, ZN., Behr, J.M., Deshpande, A. et al.
>
> _Nat Genet_ 2023 Nov 09. doi: [https://doi.org/10.1038/s41588-023-01540-6](https://doi.org/10.1038/s41588-023-01540-6)

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

