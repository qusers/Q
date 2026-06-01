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

# Number of forward/reverse switching pairs run per replicate.
neq_reps=NEQ_REPS

# Define your parameters
temperatures=(TEMP_VAR)
seeds=(RANDOM_SEEDS)
runs=${#seeds[@]}
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
fepfile=FEPS

# Debug prints
echo "Number of temperatures: ${#temperatures[@]}"
echo "Number of runs: $runs"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"

# Validate inputs
if [ -z "$runs" ] || [ "$runs" -eq 0 ]; then
    echo "Error: 'runs' variable is not set or is zero"
    exit 1
fi
if [ ${#temperatures[@]} -eq 0 ]; then
    echo "Error: No temperatures specified"
    exit 1
fi

# Calculate which temperature and replicate this array task corresponds to
TID=$((SLURM_ARRAY_TASK_ID - 1))  # 0-based
temp_idx=$((TID / runs))
run_num=$((TID % runs + 1))
temperature=${temperatures[$temp_idx]}
seed=${seeds[$run_num-1]}

## Load modules for qdyn_neq
MODULES

## define qdyn_neq location
QDYN_NEQ

starttime=$(date +%s)
starttime_readable=$(date)

rundir=$workdir/$temperature/$run_num
mkdir -p $rundir
cd $rundir || exit

echo "Running NEQ in $rundir"
echo "Parameters T=$temperature, replicate=$run_num, seed=$seed, neq_reps=$neq_reps"
echo

cp $inputfiles/eq*.inp .
cp $inputfiles/neq_*.inp .
cp $inputfiles/*.top .
cp $inputfiles/$fepfile .

# Custom random seed, temperature and FEP file for this realization
sed -i "s/SEED_VAR/$seed/" eq1.inp
sed -i "s/T_VAR/$temperature/" *.inp
sed -i "s/FEP_VAR/$fepfile/" *.inp

# Equilibration eq1 -> eq5 (no [lambda_scaling] section: plain equilibrium MD)
for i in 1 2 3 4 5; do
    time $qdyn_neq eq$i.inp > eq$i.log
done

# Seed both endpoint equilibrations from the equilibrated eq5 snapshot
cp eq5.re eq6_0_prev.re
cp eq5.re eq6_1_prev.re

# For each realization, advance a continuous endpoint equilibration (eq6) and
# fire a switch from each fresh snapshot. State 0 sweeps lambda 0->1 (reverse,
# logged as neq_0); state 1 sweeps lambda 1->0 (forward, logged as neq_1).
for rep in $(seq 0 $((neq_reps - 1))); do
    for s in 0 1; do
        # endpoint equilibration step (decorrelates successive switch starts)
        sed "s|RESTART_VAR|eq6_${s}_prev.re|; s|FINAL_VAR|eq6_${s}_${rep}.re|" \
            eq6_${s}.inp > eq6_${s}_run${rep}.inp
        time $qdyn_neq eq6_${s}_run${rep}.inp > eq6_${s}_${rep}.log
        cp eq6_${s}_${rep}.re eq6_${s}_prev.re
        # switching run: lambda is driven by the [lambda_scaling] section and the
        # accumulated work is written to the log (parsed by qligfep_neq_analyze)
        sed "s|RESTART_VAR|eq6_${s}_${rep}.re|; s|FINAL_VAR|neq_${s}_${rep}.re|" \
            neq_${s}.inp > neq_${s}_run${rep}.inp
        time $qdyn_neq neq_${s}_run${rep}.inp > neq_${s}_${rep}.log
    done
done

#CLEANUP

endtime=$(date +%s)
endtime_readable=$(date)
runtime=$((endtime - starttime))
hours=$(($runtime / 3600))
minutes=$(($runtime % 3600 / 60))
seconds=$(($runtime % 60))

echo "#    EXPRESS LOG for jobid: $SLURM_JOB_ID"
echo "#    Starttime: $starttime_readable"
echo "#    Endtime: $endtime_readable"
echo "#    Runtime: ${hours}h:${minutes}m:${seconds}s"
echo "#    Random seed: $seed"
echo "#    Replicate Number: $run_num"
echo "#    Working Directory: $workdir"
