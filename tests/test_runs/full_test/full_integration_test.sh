# N.B: You should create a separate directory for each run or experiment and
# copy/run this script

# This is the path to the shared container images (NOT WRITABLE)
export NXF_SINGULARITY_LIBRARYDIR=/gpfs/commons/groups/imielinski_lab/data/pipeline/container_images_cache
# This is the path to your local container images (WRITABLE)
# Will be created if it doesn't exist
export NXF_SINGULARIY_CACHEDIR="${HOME}/local_singularity_cache"
mkdir -p $NXF_SINGULARIY_CACHEDIR

module unload java && module load java;
module load singularity/3.8.6;

mkdir hg19  # This is necessary for fragcounter to work

# Currently pointing to the 'mskilab' branch of the pipeline repo
# You should change this to your pipeline path if you want to run your own
# modified version of the pipeline
pipeline_dir=/gpfs/commons/home/sdider/Projects/nf-casereports

nygc_config=/gpfs/commons/home/sdider/Projects/nf-casereports/tests/test_runs/nygc.config

# You can change the flags/parameters below to suit your run.

# Use `nf-core launch` on the command line to launch a web based gui
# to modify any of the pipeline parameters from their defaults and then run

# Important parameters (see README.md for more details)
# --tools: tools listed are the minimal set necessary for JaBbA output
# --step: change this to 'alignment' if you are starting from fastqs
# -resume: always tries to resume if possible, otherwise will start from beginning
# -c nygc_config: include this flag/parameter if you want to run via slurm
# full run tools: --tools gridss,hetpileups,fragcounter,strelka,sage,snpeff,dryclean,ascat,cbs,jabba,sigprofilerassignment,allelic_cn,events,fusions \
#--tools strelka,sage,snpeff,sigprofilerassignment \
nextflow run "$pipeline_dir" \
    -c "$pipeline_dir/tests/nextflow_full_test.config" \
    -c $nygc_config \
    -params-file ./params_full_integration_test.json \
	-profile singularity	\
	-with-report "report_$(date +'%Y%m%d_%H%M%S').html"	\
	-with-trace	\
    -dump-channels \
	-resume
