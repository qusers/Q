#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH --mem-per-cpu=512  # more than enough for 25A sphere size FEP
#SBATCH -A ACCOUNT
#              d-hh:mm:ss
#SBATCH --time=TIME
#SBATCH -J JOBNAME
#SBATCH --array=1-TOTAL_JOBS
#SBATCH -o slurm.run%a.%N.%j.out

# Define your parameters
temperatures=(TEMP_VAR)
seeds=(RANDOM_SEEDS)
runs=${#seeds[@]}
restartfile=RESTARTFILE
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
fepfiles=(FEPS)

# Validate inputs
if [ -z "$runs" ] || [ "$runs" -eq 0 ]; then
    echo "Error: 'runs' variable is not set or is zero"
    exit 1
fi

if [ ${#temperatures[@]} -eq 0 ]; then
    echo "Error: No temperatures specified"
    exit 1
fi

# Calculate indices based on array task ID
TID=$((SLURM_ARRAY_TASK_ID - 1))  # Convert to 0-based indexing
temp_idx=$((TID / runs))
run_num=$((TID % runs + 1))

if [ $temp_idx -ge ${#temperatures[@]} ]; then
    echo "Error: Temperature index ($temp_idx) out of bounds"
    exit 1
fi

temperature=${temperatures[$temp_idx]}
seed=${seeds[$run_num-1]}
bridge_restart=$workdir/.fep1_${temperature}_${run_num}.re

## Load modules for qdynp
MODULES

## define qdynp location
QDYN

# Fail on the first error, so a failed mpirun can't exit 0 and be recorded COMPLETED. Set after
# `module`, which returns non-zero on some stacks just for unload warnings.
set -eo pipefail

starttime=$(date +%s)
starttime_readable=$(date)

# Emitted from an EXIT trap so the footer is written even when set -e aborts the job midway.
# Keep in sync with run.sh; parsed by QligFEP.IO.read_slurm_diagnostics.
write_footer() {
    rc=$?
    rm -f -- "$bridge_restart"
    endtime=$(date +%s)
    runtime=$((endtime - starttime))
    echo "#    EXPRESS LOG for jobid: $SLURM_JOB_ID"
    echo "#    Slurm tasks: $SLURM_NTASKS"
    echo "#    Starttime: $starttime_readable"
    echo "#    Endtime: $(date)"
    echo "#    Runtime: $((runtime / 3600))h:$((runtime % 3600 / 60))m:$((runtime % 60))s"
    echo "#    Random seed: ${seed:-unknown}"
    echo "#    Replicate Number: ${run_num:-unknown}"
    echo "#    Working Directory: ${workdir:-unknown}"
    echo "#    Exit status: $rc"
}
trap write_footer EXIT

# A dual-topology mutation runs as two consecutive FEP stages: stage 1 discharges the
# wild-type side chain while the mutant one is a soft-core ghost, stage 2 grows the mutant
# side chain in and removes the wild-type. Stage 2 starts from the final coordinates of
# stage 1, so the two cannot be reordered or run independently.
length=${#fepfiles[@]}
length=$((length-1))
for index in $(seq 0 $length); do
stage=$((index+1))
fepfile=${fepfiles[$index]}
fepdir=$workdir/FEP$stage
mkdir -p $fepdir
cd $fepdir || exit
tempdir=$fepdir/$temperature
mkdir -p $tempdir
cd $tempdir || exit

rundir=$tempdir/$run_num
mkdir -p $rundir
cd $rundir || exit

echo "Running FEP stage $stage in $rundir"
echo "Parameters T=$temperature, replicate=$run_num, seed=$seed"

cp $inputfiles/md_${stage}_*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
sed -i "s/SEED_VAR/$seed/" eq1.inp # change the random seed to custom
else
    # Stage 2 inherits stage 1's endpoint instead of equilibrating again.
    if [ ! -s "$bridge_restart" ]; then
        echo "ERROR: FEP1 bridge restart is missing or empty: $bridge_restart"
        exit 1
    fi
    mv -- "$bridge_restart" "$rundir/eq5.re"
fi
sed -i "s/T_VAR/$temperature/" *.inp
sed -i "s/FEP_VAR/$fepfile/" *.inp
if [ $index -lt 1 ]; then
#EQ_FILES
# Validate equilibration immediately. Its logs may then be archived or removed
# for quota without making a completed FEP1 appear to have failed hours later.
eq_incomplete=""
for log in eq*.log; do
    grep -q "terminated normally" "$log" 2>/dev/null \
        || eq_incomplete="$eq_incomplete $log"
done
if [ -n "$eq_incomplete" ]; then
    echo "ERROR: equilibration for replicate $run_num (T=$temperature, seed=$seed) is incomplete."
    echo "ERROR: logs without a normal qdyn termination:$eq_incomplete"
    exit 1
fi
#RUN_1_FILES
else
#RUN_2_FILES
fi

# qfep reads the energy files from lambda 1 down to 0; the list is built here because
# it differs per stage. Reverse lexicographic order is that order, since the file names
# carry lambda1 zero-padded (`sort -r` rather than `tac`, which BSD systems lack).
ls -1 md_${stage}_*.en | sort -r >> qfep.inp
timeout 3m QFEP < qfep.inp > qfep.out || [ $? -eq 124 ]

# qdyn's serial abort still exits 0, so a completed run is confirmed by "terminated normally"
# in its log, not by the exit code. Failing here also stops stage 2 from starting on a bad eq5.re.
incomplete=""
logs_to_check=(md_${stage}_*.log)
for log in "${logs_to_check[@]}"; do
    grep -q "terminated normally" "$log" 2>/dev/null || incomplete="$incomplete $log"
done
if [ -n "$incomplete" ]; then
    echo "ERROR: FEP$stage replicate $run_num (T=$temperature, seed=$seed) is incomplete."
    echo "ERROR: logs without a normal qdyn termination:$incomplete"
    exit 1
fi
# Preserve the sole restart needed by FEP2 before stage-local cleanup removes
# FEP1's restart chain. The bridge name is unique per temperature and replicate.
if [ $index -lt 1 ]; then
    cp -- "$restartfile" "$bridge_restart"
fi
# Cleanup must stay inside the stage loop. The working directory changes to the
# current FEP stage above, so cleaning after `done` would only clean FEP2.
#CLEANUP
done
