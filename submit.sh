#!/bin/env bash
#
# Co-expression network construction parameter optimization
# Keith Hughitt (khughitt@umd.edu)
#
# v6.0
# ----
#
# - GO enrichment now includes ancestor terms associated with each annotation
# - Updated TriTrypDB from 27 -> 32
# - Updated Bioconductor from 3.2 -> 3.5
# - Updated human reference annotations from Ensembl 83 -> 88
# - Using Marbach TF regulon dataset in place of ENCODE regulon data for human
# - Added LeishCyc enrichment analysis for L. major
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

# Slurm parameters
#export SBATCH_PARAMS="-q workstation -l mem=46gb,ncpus=32,walltime=32:00:00"
#export SBATCH_PARAMS="-q workstation -l mem=46gb,ncpus=4,walltime=32:00:00"
export SBATCH_PARAMS="--qos=throughput --mem=3600 --ntasks-per-node=12 --time=0-18:00:00"

# Base network construction settings to use
export WORKING_DIR=$(pwd)

# Parameters to test
export NET_TYPE="signed"
export LOW_COUNT_THRESHOLD="1"
export CPM="false true"
export LOG2="false true"
export VOOM="false"
export QUANTILE_NORMALIZE="false true"
export BATCH_ADJUST="none limma combat"
export SIM_MEASURE="cor bicor spearman cor-dist"
export ADJ_POWER="1 2 3 4 5 6 7 8 9 10 11 12 13 14"
export TOPOLOGICAL_OVERLAP="false"
export MERGE_COR="1.0"

# current date
datestr=$(date +'%Y%m%d%H%M%S')

# first eight chracters of queue name 
# used to filter results from squeue
export QUEUE=$(echo $SBATCH_PARAMS| egrep -o "\-qos=([a-z]+)" | cut -c 6- | cut -c -8)

# Maximum number of jobs to run when using sbatch/parallel
export MAX_JOBS_PARALLEL=1

if [[ "$QUEUE" == "workstat" ]]; then
    export MAX_JOBS_SBATCH=5
else
    export MAX_JOBS_SBATCH=35
fi

# Time to wait between submitting batches of jobs on sbatch
export PAUSE_TIME_IN_SECS=600

# Counter start
start_i=1

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
#  bodymap - Illumina Bodymap
#  modencode-worm - ModEncode Worm
#  modencode-fly - ModEncode Fly
#
if [[ "$1" == "lmmm" ]]; then
    # L. major / M. musculus
    echo "Optimizing network construction for L. major infecting M. musculus"
    export SETTINGS_FILE="${RESEARCH}/2015/02-coex-network-lmajor-infecting-mmusculus/settings/lmajor_intracellular-v6.0.Rmd"
    export JOB_NAME="lmajor_mouse_param_opt-$datestr"
    export SUBDIR="normal/lmajor_infecting_mmusculus-v6.0"
elif [[ "$1" == "lmhs" ]]; then
    # L. major / H. sapiens
    echo "Optimizing network construction for L. major infecting H. sapiens"
    export SETTINGS_FILE="${RESEARCH}/2015/01-coex-network-lmajor-infecting-hsapiens/settings/lmajor_intracellular-v6.0.Rmd"
    export JOB_NAME="lmajor_human_opt-$datestr"
    export SUBDIR="normal/lmajor_infecting_hsapiens-v6.0"
elif [[ "$1" == "lmall" ]]; then
    # L. major / H. sapiens
    echo "Optimizing network construction for L. major - All samples"
    export SETTINGS_FILE="${RESEARCH}/2015/14-coex-network-lmajor-all-samples/settings/lmajor_all_samples-v6.0.Rmd"
    export JOB_NAME="lmajor_all_opt-$datestr"
    export SUBDIR="normal/lmajor_all_samples-v6.0"
elif [[ "$1" == "mmlm" ]]; then
    # M. musculus / L. major
    echo "Optimizing network construction for M. musculus infected with L. major"
    export SETTINGS_FILE="${RESEARCH}/2015/12-coex-network-mmusculus-infected-with-lmajor/settings/mouse_intracellular-v6.0.Rmd"
    export JOB_NAME="mm_lm_opt-$datestr"
    export SUBDIR="normal/mmusculus_infected_with_lmajor-v6.0"
elif [[ "$1" == "hslb" ]]; then
    # L. braziliensis / Human
    echo "Optimizing network construction for Human infected with L. braziliensis"
    export SETTINGS_FILE="${RESEARCH}/2015/17-coex-network-hsapiens-infected-with-lbraziliensis/settings/hsapiens_inf_with_lbraziliensis-v6.0.Rmd"
    export JOB_NAME="hslb_param_opt-$datestr"
    export SUBDIR="normal/hsapiens_infected_with_lbraziliensis-v6.0"
    export BATCH_ADJUST="none"
elif [[ "$1" == "hslm" ]]; then
    # H. sapiens / L. major
    echo "Optimizing network construction for H. sapiens infected with L. major"
    export SETTINGS_FILE="${RESEARCH}/2015/11-coex-network-hsapiens-infected-with-lmajor/settings/hsapiens_intracellular-v6.0.Rmd"
    export JOB_NAME="human_lmajor_param_opt-$datestr"
    export SUBDIR="normal/hsapiens_infected_with_lmajor-v6.0"
elif [[ "$1" == "hslm-full" ]]; then
    # H. sapiens / L. major (including TOM, signed/unsigned nets)
    echo "Optimizing network construction for H. sapiens infected with L. major"
    export SETTINGS_FILE="${RESEARCH}/2015/11-coex-network-hsapiens-infected-with-lmajor/settings/hsapiens_intracellular-v6.0.Rmd"
    export JOB_NAME="human_lmajor_param_opt-$datestr"
    export SUBDIR="normal/hsapiens_infected_with_lmajor-full-v6.0"
    export NET_TYPE="signed unsigned"
    export TOPOLOGICAL_OVERLAP="false true"
elif [[ "$1" == "hstc" ]]; then
    # H. sapiens / T. cruzi
    echo "Optimizing network construction for H. sapiens infected with T. cruzi"
    export SETTINGS_FILE="${RESEARCH}/2015/13-coex-network-hsapiens-infected-with-tcruzi/settings/hsapiens_intracellular-v6.0.Rmd"
    export JOB_NAME="human_tcruzi_param_opt-$datestr"
    export SUBDIR="normal/hsapiens_infected_with_tcruzi-v6.0"
elif [[ "$1" == "tchs" ]]; then
    # T. cruzi / H. sapiens
    echo "Optimizing network construction for T. cruzi infecting H. sapiens"
    export SETTINGS_FILE="${RESEARCH}/2015/03-coex-network-tcruzi-infecting-hsapiens/settings/tcruzi_intracellular-v6.0.Rmd"
    export JOB_NAME="tcruzi_hsapiens_param_opt-$datestr"
    export SUBDIR="normal/tcruzi_infecting_hsapiens-v6.0"
elif [[ "$1" == "bodymap" ]]; then
    # T. cruzi / H. sapiens
    echo "Optimizing network construction for Illumina BodyMap"
    export SETTINGS_FILE="${RESEARCH}/2015/18-coex-network-illumina-bodymap/settings/illumina_bodymap2-v6.0.Rmd"
    export JOB_NAME="bodymap_param_opt-$datestr"
    export SUBDIR="normal/illumina_bodymap-v6.0"
    export BATCH_ADJUST="none"
elif [[ "$1" == "modencode-fly" ]]; then
    # ModENCODE Fly
    echo "Optimizing network construction for ModENCODE Fly"
    export SETTINGS_FILE="${RESEARCH}/2016/02-coex-networks/01-recount-modencode-fly/settings/recount-modencode-fly-settings-v6.0.Rmd"
    export JOB_NAME="modencode_fly_param_opt-$datestr"
    export SUBDIR="normal/modencode_fly-v6.0"
    export BATCH_ADJUST="none"
elif [[ "$1" == "modencode-worm" ]]; then
    # ModENCODE Worm
    echo "Optimizing network construction for ModENCODE Worm"
    export SETTINGS_FILE="${RESEARCH}/2016/02-coex-networks/02-recount-modencode-worm/settings/recount-modencode-worm-settings-v6.0.Rmd"
    export JOB_NAME="modencode_worm_param_opt-$datestr"
    export SUBDIR="normal/modencode_worm-v6.0"
    export BATCH_ADJUST="none"
else
    echo 'Invalid option selected...'
    echo 'Usage: submit.sh <settings_choices>'
    exit
fi


# Log directory
#SUBDIR=$(date +"%Y%m%d%H%M%S")

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

# Check for cluster availability
# If not found, will fall back on GNU parallel
export USE_CLUSTER=$(type sbatch &>/dev/null && echo "true" || echo "false") 

# Maximum number of threads to use per job
if [[ $USE_CLUSTER == true ]]; then
    export ALLOW_WGCNA_THREADS=12
else
    export ALLOW_WGCNA_THREADS=4
fi

# Loop over parameter ranges
for p0 in $NET_TYPE; do
    for p1 in $LOW_COUNT_THRESHOLD; do
        for p2 in $CPM; do
            for p3 in $LOG2; do
                for p4 in $VOOM; do
                    # Voom automatically log2-CPM transforms data; skip any parameter
                    # combinations where voom is enabled and the data is not
                    # log- and CPM-transformed.

                    # If voom is enabled..
                    if [[ $p4 == true ]]; then
                        # Skip 
                        if [[ $p2 == false || $p3 == false ]]; then
                            continue
                        fi
                    fi
                    for p5 in $QUANTILE_NORMALIZE; do
                        for p6 in $BATCH_ADJUST; do
                            for p7 in $SIM_MEASURE; do
                                for p8 in $ADJ_POWER; do
                                    for p9 in $TOPOLOGICAL_OVERLAP; do
                                        for p10 in $MERGE_COR; do
                                            # Skip files outside of desired range
                                            if [[ $i -lt $min_i ]]; then
                                                let i=$i+1
                                                continue
                                            fi

                                            # Make parameters visible to sbatch
                                            export p0 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10

                                            LOG="${LOGDIR}/${i}.log"

                                            # Output directory
                                            export OUTFILE="${OUTDIR}/${i}.out"

                                            # skip files that have already been processed
                                            if [[ -f "$OUTFILE" ]]; then
                                                let i=$i+1
                                                continue
                                            fi

                                            # build command
                                            export rscript=$(pwd)/coex_network_param_opt.R
                                            export cmd="Rscript $rscript $OUTFILE $SETTINGS_FILE $p0 $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10"

                                            # Optional: skip certain parameter combinations
                                            SKIP_JOB=false 

                                            # 2015/06/30
                                            # To speed up processing of host networks, we will
                                            # temporarily skip some parameters
                                            #if [[ "$p1" != "1" || "$p7" != "cor" || "$p8" != "1" || "$p9" == true ]]; then
                                            #    SKIP_JOB=true
                                            #fi

                                            if [[ "$SKIP_JOB" == false ]]; then
                                                # Submit job (sbatch)
                                                echo "Submitting job $i..."

                                                if [[ "$USE_CLUSTER" == true ]]; then
                                                    #cat <<"EOF" | sbatch $SBATCH_PARAMS \
                                                    #                -N "param.opt.${i}" --mail-type=NONE \
                                                    #                -D $WORKING_DIR -o $LOG -
                                                    sbatch $SBATCH_PARAMS \
                                                                    -J "param.opt.${i}" --mail-type=NONE \
                                                                    -D $WORKING_DIR -o $LOG <<"EOF"
#!/usr/bin/env bash
#------------------------------------------------------------------------------
# START SBATCH SCRIPT
#------------------------------------------------------------------------------
echo $cmd
eval $cmd
EOF
#------------------------------------------------------------------------------
# END SBATCH SCRIPT
#------------------------------------------------------------------------------
                                                else
                                                    # Submit job (parallel)
                                                    # NOTE 2015/03/11: Attempts to capture STDOUT not working at the moment
                                                    #parallel --semaphore -j$MAX_JOBS_PARALLEL $cmd | tee $LOG 
                                                    #parallel --semaphore -j$MAX_JOBS_PARALLEL $cmd > $LOG
                                                    #parallel --semaphore -j$MAX_JOBS_PARALLEL $cmd
                                                    #sem --jobs $MAX_JOBS_PARALLEL --id $JOB_NAME -u '$cmd' | tee $LOG
                                                    sem --jobs $MAX_JOBS_PARALLEL --id $JOB_NAME -u '$cmd'

                                                    # Wait if running maximum number of jobs
                                                    num_running=$(($(ps -Af | grep -c parallel) - 1))
                                                    echo "Number of jobs running: " $num_running

                                                    #if [[ $num_running -ge $MAX_JOBS_PARALLEL ]]; then
                                                    #if [[ $num_running -ge 5 ]]; then
                                                    #    echo "Waiting for jobs to finish before continuing..."
                                                    #    parallel --wait
                                                    #fi
                                                fi
                                            fi
                                            let i=$i+1

                                            # Pause between queue submissions to avoid
                                            # overwhelming the queue
                                            if [[ "$USE_CLUSTER" == true ]]; then
                                                if [[ $i -gt $min_i ]]; then
                                                    # Get the number of currently queued/running jobs
                                                    num_jobs=$(squeue --user=$USER | grep param | grep $QUEUE | wc -l)

                                                    while [[ $num_jobs -ge $MAX_JOBS_SBATCH ]]; do
                                                        printf "Pausing for %d seconds...\n" $PAUSE_TIME_IN_SECS
                                                        sleep $PAUSE_TIME_IN_SECS
                                                        num_jobs=$(squeue --user=$USER | grep param | grep $QUEUE | wc -l)
                                                    done
                                                fi
                                            fi
                                        done
                                    done
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done

