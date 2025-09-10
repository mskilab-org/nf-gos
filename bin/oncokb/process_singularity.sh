#!/bin/bash -e

umask u=rwx,g=rwx,o=rx

tolower() {
    echo $1 | awk '{print tolower($0)}'
}


set -x
set -a
LIBDIR=.
VCF=$1
FUSIONS=$2
CNA=$3
MULTIPLICITY=$4
ASSEMBLY=$5
CORES=$6
DO_VEP=$(tolower $7)
VERBOSE=$(tolower $8)
TUMOR_TYPE=${9}
VEP_DIR=${10}
ONCOKB_GENES=${11}
GENCODE=${12}
OUTDIR=.
set +a
echo $TUMOR_TYPE
set +x


export TMPDIR=/tmp # as definied in the singularity container
export PATH=/root:scripts:${PATH}

## This should be the sample directory.
## It is a mount bindpoint
pushd .


files_only() {
    ls -p | grep -v /
}

parallel () {
    /opt/conda/envs/pact/bin/parallel --will-cite $@
}
export -f parallel

## OncoKB

# Test if ONCOKB_TOKEN is already in environment
if [ "${ONCOKB_TOKEN+set}" != set ]; then
    echo "ONCOKB_TOKEN not found, please pass it as an environment variable."
    exit 1
fi

if [ -e $VCF ] && [ ! $(wc -c <${VCF}) == 0 ]; then
    cleanup () {
        rm -f $tmpVcf
    }
    trap cleanup EXIT

    do_create_tmp=false
    case "$VCF" in
    *.gz | *.tgz )
        # it's gzipped
        echo "VCF is gzipped, creating tmp"
        tmpVcf="$(basename ${VCF%.*})"
        bgzip -d -c ${VCF} > ./${tmpVcf}
        VCF=${tmpVcf}
        do_create_tmp=true
        ;;
    esac

    is_vcf_ucsc=false
    bcftools view -h ${VCF} | grep "^##contig=" | grep "ID=chr" >/dev/null && is_vcf_ucsc=true
    if [ "${is_vcf_ucsc}" = "true" ]; then
        echo "VCF has UCSC style contigs (i.e. contains chr)"
        tmprename="$(mktemp -t tmp.XXXXXXXXXX.vcf)"
        ${HOME}/scripts/vcf_convert_ucsc_to_ncbi_contigs.py ${VCF} ${tmprename}
        mv ${tmprename} ./$(basename ${VCF})_rename
        VCF=./$(basename ${VCF})_rename
    fi

    vep_vcf_intermediate="${VCF%.*}.vep.vcf"
    if [ -e ${vep_vcf_intermediate} ]; then
        echo "VEP VCF intermediate found.. removing"
        rm -f ${vep_vcf_intermediate}
    fi

    if [ -e ./annotated.maf ]; then
        echo "Annotated MAF found.. removing"
        rm -f ./annotated.maf
    fi


    snames=$(bcftools query -l $VCF)
    tumor_id=$(parallel '[[ {} =~ ^TUMOR$|[Tt]umor|-T$ ]] && echo {} || echo ""' ::: $snames)
    normal_id=$(parallel '[[ {} =~ ^NORMAL$|[Nn]ormal|-T$ ]] && echo {} || echo ""' ::: $snames)
    if [ "$tumor_id" == "" ]; then
        tumor_id=$(bcftools query -l $VCF | awk 'NR==2') `# assume 2nd id is tumor if sample names don't conform`
    fi
    if [ "$normal_id" == "" ]; then
        normal_id=$(bcftools query -l $VCF | awk 'NR==1') `# assume 1st id is normal if the sample names don't conform`
    fi

    echo "Running vcf2maf"

    if [ ! -e ./annotated.maf ]
    then
        is_grch37=$( ( [ $(tolower $ASSEMBLY) == "hg19" ] || [ $(tolower $ASSEMBLY) == "grch37" ] ) && echo true || false )
        is_grch38=$( ( [ $(tolower $ASSEMBLY) == "hg38" ] || [ $(tolower $ASSEMBLY) == "grch38" ] ) && echo true || false )
        if $is_grch37; then
            ref_path=${VEP_DIR}/homo_sapiens/112_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
        elif $is_grch38; then
            ref_path=${VEP_DIR}/homo_sapiens/113_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        else
            echo "VEP Reference FASTA not found/supported!"
            exit 1
        fi

        perl ${HOME}/git/vcf2maf/vcf2maf.pl \
            --inhibit-vep \
            --vep-data=${VEP_DIR} \
            --ref-fasta ${ref_path} `# --ref-fasta ${VEP_DIR}/homo_sapiens/112_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz` \
            --vep-path /opt/conda/envs/pact/bin \
            --input-vcf ${VCF} \
            --vep-forks=1 \
            --output-maf ./annotated.maf \
            --vcf-tumor-id  ${tumor_id} \
            --vcf-normal-id ${normal_id} \
            || { echo "VCF2MAF failed!" && exit 1; }
    fi

    echo "Running OncoKB annotator"

    ## test line below tests if TUMOR_TYPE is not NA
    ## if not, then provide the argument, otherwise allow for default argument (Kevin: not clear what it is explicitly to me atm)
    python ${HOME}/git/oncokb-annotator/MafAnnotator.py \
        -i ./annotated.maf \
        -o ./oncokb_annotated.maf \
        -b ${ONCOKB_TOKEN} \
        -r GRCh37 \
        -d || \
        { echo "OncoKB Annotation failed!" && exit 1; }



    (
        Rscript "${HOME}/scripts/add_role_to_oncokb_maf.R" \
            --oncokb ${ONCOKB_GENES} \
            --oncokb_output ./oncokb_annotated.maf \
            --output_path ./ \
            --multiplicity ${MULTIPLICITY} \
            --cna ${CNA}

        # Output is "merged_oncokb.maf"

    )
fi

if [ -e ${FUSIONS} ] && [ ! $(wc -c <${FUSIONS}) == 0 ]; then

        # Below must be inside parentheses so the environment change happens within a subshell
        (
            Rscript ${HOME}/scripts/parse_fusions_for_oncokb.R --libdir ${LIBDIR} --fusions ${FUSIONS} --path ./input_fusions.tsv
        )
        set -x
        python ${HOME}/git/oncokb-annotator/FusionAnnotator.py \
            -i ./input_fusions.tsv \
            -o ./oncokb_fusions.tsv \
            $( test ! $(tolower ${TUMOR_TYPE}) = unknown && printf %s "-t ${TUMOR_TYPE}" )  \
            -b $ONCOKB_TOKEN \
            -d
        set +x

        Rscript ${HOME}/scripts/add_fusion_metadata_to_oncokb.R \
            --inter_file ./intermediate_fusions.rds \
            --oncokb ${ONCOKB_GENES} \
            --oncokb_output ./oncokb_fusions.tsv \
            --output_path ./

         # Output is merged_oncokb_fusions.tsv
    fi



if [ -e ${CNA} ] && [  ! $(wc -c <${CNA}) == 0 ]; then
    # Below must be inside parentheses so the environment change happens within a subshell
    (
        Rscript ${HOME}/scripts/parse_cna_for_oncokb.R \
            --gencode ${GENCODE} \
            --oncokb ${ONCOKB_GENES} \
            --jabba ${CNA} \
            --libdir ${LIBDIR}/OncoKB_Annotator \
            --output_path ./

    )
    python ${HOME}/git/oncokb-annotator/CnaAnnotator.py \
        -i ./oncokb_input_scna.txt \
        -o ./oncokb_cna.tsv \
        -b $ONCOKB_TOKEN \
        -d
    (
        Rscript ${HOME}/scripts/add_min_cn_to_cn_oncokb.R \
            --inter_file ./intermediate_scna.rds \
            --oncokb ${ONCOKB_GENES} \
            --oncokb_output ./oncokb_cna.tsv \
            --output_path ./
    )
    # Output is merged_oncokb_cna.tsv
fi


touch FINISHED
# chmod u=rwx,g=rwx,o=rx *

