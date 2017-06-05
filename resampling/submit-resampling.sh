#!/bin/env bash
#
# Co-expression network construction parameter optimization
#
# v6.0
# ----
#
# - GO enrichment now includes ancestor terms associated with each annotation
# - Updated TriTrypDB from 27 -> 31
# - Updated Bioconductor from 3.2 -> 3.4
#
# v5.0
# ----
# 
# - Fixed an issue with HTSeq causing many reads to be double-counted.
# - Updated human reference annotations from Ensembl 76 -> 83
# - Updated TriTrypDB from 25 -> 27
# - Updated Bioconductor from 3.1 -> 3.2
# - Added unsigned back in to list of params to optimize
#
# v4.0
# ----
#
# - Enabling filtering of non-DE genes (using a p-value cutoff of 0.1) 
#   for all datasets
# - For L. major and T. cruzi, multicopy genes are also excluded.
# - Including spearman correlation and biweight midcorrelation (using
#   suggested parameterization)
#
# v3.0
# ----
#
# - Added cor-dist and signed-dist similarity measures
#

# Qsub parameters
#export QSUB_PARAMS="-q workstation -l mem=46gb,ncpus=32,walltime=32:00:00"
export QSUB_PARAMS="-q throughput -l mem=36gb,ncpus=12,walltime=18:00:00"

# Base network construction settings to use
export WORKING_DIR=$(pwd)

# Maximum number of jobs to run when using qsub/parallel
export MAX_JOBS_QSUB=35

# Time to wait between submitting batches of jobs on qsub
export PAUSE_TIME_IN_SECS=600

# Counter start
start_i=1

# Number of resampled networks to construct
num_networks=500

# current date
datestr=$(date +'%Y%m%d%H%M%S')

# first eight chracters of queue name 
# used to filter results from qstat
export QUEUE=$(echo $QSUB_PARAMS| egrep -o "\-q ([a-z]+)" | cut -c 4- | cut -c -8)

#
# Settings choices:
#
#  lmmm - L. major infecting mouse
#  lmhs - L. major infecting human
#  tchs - T. cruzi infecting human
#  hslb - Human infected with L. braziliensis
#  hslm - Human infected with L. major
#  hstc - Human infected with T. cruzi
#  mmlm - Mouse infected with L. major
#  lmall - L. major all samples 
#  tcall - T. cruzi all samples 
#
if [[ "$1" == "lmmm-70" ]]; then
    # L. major / M. musculus
    echo "Optimizing network construction for L. major infecting M. musculus"
    export SETTINGS_FILE="${RESEARCH}/2015/02-coex-network-lmajor-infecting-mmusculus/settings/lmajor_intracellular-v5.0.Rmd"
    export JOB_NAME="lmajor_mouse_resampled_param_opt-$datestr"
    export SUBDIR="resampled/lmajor_infecting_mmusculus-net71-v5.0"
elif [[ "$1" == "lmhs-68" ]]; then
    # L. major / H. sapiens
    echo "Optimizing network construction for L. major infecting H. sapiens"
    export SETTINGS_FILE="${RESEARCH}/2015/01-coex-network-lmajor-infecting-hsapiens/settings/lmajor_intracellular-v5.0.Rmd"
    export JOB_NAME="lmajor_human_resampled_param_opt-$datestr"
    export SUBDIR="resampled/lmajor_infecting_hsapiens-net68-v5.0"
elif [[ "$1" == "lmall-762" ]]; then
    # L. major / H. sapiens
    echo "Optimizing network construction for L. major - All samples"
    export SETTINGS_FILE="${RESEARCH}/2015/14-coex-network-lmajor-all-samples/settings/lmajor_all_samples-v5.0.Rmd"
    export JOB_NAME="lmajor_all_resampled_param_opt-$datestr"
    export SUBDIR="resampled/lmajor_all_samples-net762-v5.0"
elif [[ "$1" == "mmlm" ]]; then
    # M. musculus / L. major
    echo "Optimizing network construction for M. musculus infected with L. major"
    export SETTINGS_FILE="${RESEARCH}/2015/12-coex-network-mmusculus-infected-with-lmajor/settings/mouse_intracellular-v5.0.Rmd"
    export JOB_NAME="mouse_lmajor_resampled_param_opt-$datestr"
    export SUBDIR="resampled/mmusculus_infected_with_lmajor-v5.0"
elif [[ "$1" == "hslm-101" ]]; then
    # H. sapiens / L. major
    echo "Optimizing network construction for H. sapiens infected with L. major"
    export SETTINGS_FILE="${RESEARCH}/2015/11-coex-network-hsapiens-infected-with-lmajor/settings/hsapiens_intracellular-v5.0.Rmd"
    export JOB_NAME="human_lmajor_resampled_param_opt-$datestr"
    export SUBDIR="resampled/hsapiens_infected_with_lmajor-net101-v5.0"
elif [[ "$1" == "hstc-446" ]]; then
    # H. sapiens / L. major
    echo "Optimizing network construction for H. sapiens infected with T. cruzi"
    export SETTINGS_FILE="${RESEARCH}/2015/13-coex-network-hsapiens-infected-with-tcruzi/settings/hsapiens_intracellular-v5.0.Rmd"
    export JOB_NAME="human_tcruzi_resampled_param_opt-$datestr"
    export SUBDIR="resampled/hsapiens_infected_with_tcruzi-net446-v5.0"
elif [[ "$1" == "hstc-777" ]]; then
    # H. sapiens / L. major
    echo "Optimizing network construction for H. sapiens infected with T. cruzi"
    export SETTINGS_FILE="${RESEARCH}/2015/13-coex-network-hsapiens-infected-with-tcruzi/settings/hsapiens_intracellular-v5.0b.Rmd"
    export JOB_NAME="human_tcruzi_resampled_param_opt-$datestr"
    export SUBDIR="resampled/hsapiens_infected_with_tcruzi-net777-v5.0"
elif [[ "$1" == "hstc-329" ]]; then
    # H. sapiens / L. major
    echo "Optimizing network construction for H. sapiens infected with T. cruzi"
    export SETTINGS_FILE="${RESEARCH}/2015/13-coex-network-hsapiens-infected-with-tcruzi/settings/hsapiens_intracellular-v5.0c.Rmd"
    export JOB_NAME="human_tcruzi_resampled_param_opt-$datestr"
    export SUBDIR="resampled/hsapiens_infected_with_tcruzi-net329-v5.0"
elif [[ "$1" == "tchs-294" ]]; then
    # T. cruzi / H. sapiens
    echo "Optimizing network construction for T. cruzi infecting H. sapiens"
    export SETTINGS_FILE="${RESEARCH}/2015/03-coex-network-tcruzi-infecting-hsapiens/settings/tcruzi_intracellular-v5.0.Rmd"
    export JOB_NAME="tchs_resampled_param_opt-$datestr"
    export SUBDIR="resampled/tcruzi_infecting_hsapiens-net294-v5.0"
elif [[ "$1" == "tchs-331" ]]; then
    # T. cruzi / H. sapiens
    echo "Optimizing network construction for T. cruzi infecting H. sapiens"
    export SETTINGS_FILE="${RESEARCH}/2015/03-coex-network-tcruzi-infecting-hsapiens/settings/tcruzi_intracellular-v5.0c.Rmd"
    export JOB_NAME="tchs_resampled_param_opt-$datestr"
    export SUBDIR="resampled/tcruzi_infecting_hsapiens-net331-v5.0"
elif [[ "$1" == "tchs-677" ]]; then
    # T. cruzi / H. sapiens
    echo "Optimizing network construction for T. cruzi infecting H. sapiens"
    export SETTINGS_FILE="${RESEARCH}/2015/03-coex-network-tcruzi-infecting-hsapiens/settings/tcruzi_intracellular-v5.0b.Rmd"
    export JOB_NAME="tchs_resampled_param_opt-$datestr"
    export SUBDIR="resampled/tcruzi_infecting_hsapiens-net677-v5.0"
elif [[ "$1" == "bodymap" ]]; then
    # T. cruzi / H. sapiens
    echo "Optimizing network construction for Illumina BodyMap"
    export SETTINGS_FILE="${RESEARCH}/2015/18-coex-network-illumina-bodymap/settings/illumina_bodymap2-v5.0.Rmd"
    export JOB_NAME="bodymap_resampled_param_opt-$datestr"
    export SUBDIR="resampled/illumina_bodymap-v5.0"
else
    echo 'Invalid option selected...'
    echo 'Usage: submit.sh <settings_choices>'
    exit
fi

# Log directory
LOGDIR="$(pwd)/log/${SUBDIR}"

# Output directory
OUTDIR="$(pwd)/output/${SUBDIR}"

# Set TMPDIR to avoid running out of space in /tmp
export TMPDIR=/cbcb/nelsayed-scratch/keith/tmp/

mkdir -p $LOGDIR
mkdir -p $OUTDIR
mkdir -p $TMPDIR

# Iteration counter
i=$start_i

# Minimum iteration (for continuing an earlier submission)
if [[ -z "$min_i" ]]; then
    min_i=$start_i
fi

# Maximum number of threads to use per job
export ALLOW_WGCNA_THREADS=12

# Loop over parameter ranges
for i in $(seq 1 $num_networks); do
    # Skip files outside of desired range
    if [[ $i -lt $min_i ]]; then
        let i=$i+1
        continue
    fi

    LOG="${LOGDIR}/${i}.log"

    # Output directory
    export OUTFILE="${OUTDIR}/${i}.out"

    # skip files that have already been processed
    if [[ -f "$OUTFILE" ]]; then
        let i=$i+1
        continue
    fi

    # build command
    export rscript=$(pwd)/coex_network_param_opt_resampling.R
    export cmd="Rscript $rscript $OUTFILE $SETTINGS_FILE"

    # Submit job (qsub)
    echo "Submitting job $i..."

    cat <<"EOF" | qsub $QSUB_PARAMS \
                    -N "param.opt.${i}" -V -m n \
                    -j eo -w $WORKING_DIR -e $LOG -
#------------------------------------------------------------------------------
# START QSUB SCRIPT
#------------------------------------------------------------------------------
echo $cmd
eval $cmd
EOF
#------------------------------------------------------------------------------
# END QSUB SCRIPT
#------------------------------------------------------------------------------
    let i=$i+1

    # Pause between queue submissions to avoid
    # overwhelming the queue
    if [[ $i -gt $min_i ]]; then
        # Get the number of currently queued/running jobs
        num_jobs=$(qstat -u $USER | grep param | grep $QUEUE | wc -l)

        while [[ $num_jobs -ge $MAX_JOBS_QSUB ]]; do
            printf "Pausing for %d seconds...\n" $PAUSE_TIME_IN_SECS
            sleep $PAUSE_TIME_IN_SECS
            num_jobs=$(qstat -u $USER | grep param | grep $QUEUE | wc -l)
        done
    fi
done
