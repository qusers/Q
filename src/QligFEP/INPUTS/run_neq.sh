#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH --mem-per-cpu=2000
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
ncores=$SLURM_NTASKS

# Debug prints
echo "Number of temperatures: ${#temperatures[@]}"
echo "Number of runs: $runs"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Cores available: $ncores"

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

## Load modules
MODULES

## define the MPI equilibration engine (qdynp) and serial switching engine (qdyn_neq)
QDYN
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

# 1) Equilibration eq1 -> eq5 with the MPI engine across all cores (fixed lambda)
for i in 1 2 3 4 5; do
    time mpirun -np $ncores --bind-to core $qdyn eq$i.inp > eq$i.log
done

# 2) Endpoint equilibration chain (MPI): one continuous trajectory per endpoint,
# saving a checkpoint per replicate to decorrelate the switch starting points.
# State 0 sits at lambda (0 1), state 1 at (1 0).
cp eq5.re eq6_0_prev.re
cp eq5.re eq6_1_prev.re
for rep in $(seq 0 $((neq_reps - 1))); do
    for s in 0 1; do
        sed "s|RESTARTFILE|eq6_${s}_prev.re|; s|FINALFILE|eq6_${s}_${rep}.re|" \
            eq6_${s}.inp > eq6_${s}_run${rep}.inp
        time mpirun -np $ncores --bind-to core $qdyn eq6_${s}_run${rep}.inp > eq6_${s}_${rep}.log
        cp eq6_${s}_${rep}.re eq6_${s}_prev.re
    done
done

# 3) Lambda switches with the serial engine, one per core. Each switch only needs
# its own eq6 checkpoint, so they are independent and run concurrently: mpirun
# launches one rank per switch, --bind-to core pins each to its own core, and the
# launcher routes each rank to its (input, log) pair. State 1 (1->0) is the forward
# work logged as neq_1; state 0 (0->1) is the reverse work logged as neq_0.
: > switch_list.txt
for rep in $(seq 0 $((neq_reps - 1))); do
    for s in 0 1; do
        sed "s|RESTARTFILE|eq6_${s}_${rep}.re|; s|FINALFILE|neq_${s}_${rep}.re|" \
            neq_${s}.inp > neq_${s}_run${rep}.inp
        echo "neq_${s}_run${rep}.inp neq_${s}_${rep}.log" >> switch_list.txt
    done
done

cat > neq_launch.sh <<EOF
#!/bin/bash
# one switch per MPI rank: read this rank's (input, log) from the list
idx=\$(( OMPI_COMM_WORLD_RANK + 1 + \${BATCH_OFFSET:-0} ))
line=\$(sed -n "\${idx}p" switch_list.txt)
[ -z "\$line" ] && exit 0
set -- \$line
$qdyn_neq "\$1" > "\$2"
EOF
chmod +x neq_launch.sh

nsw=$(wc -l < switch_list.txt)
for (( off=0; off<nsw; off+=ncores )); do
    remaining=$(( nsw - off ))
    np=$(( remaining < ncores ? remaining : ncores ))
    export BATCH_OFFSET=$off
    time mpirun -x BATCH_OFFSET --bind-to core -np $np ./neq_launch.sh
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
