#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH --mem-per-cpu=1700  # stay at/below the per-core memory so billing matches the core count
#SBATCH -A ACCOUNT
#              d-hh:mm:ss
#SBATCH --time=TIME
#SBATCH -J JOBNAME
#SBATCH --array=1-TOTAL_JOBS
#SBATCH -o slurm.run%a.%N.%j.out

starttime=$(date +%s)
starttime_readable=$(date)

# Emitted from an EXIT trap so the footer is written even when set -e aborts the job midway.
# Keep in sync with run.sh; parsed by QligFEP.IO.read_slurm_diagnostics.
write_footer() {
    rc=$?
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

## define the parallel equilibration engine (qdynp) and serial switching engine (qdyn)
QDYNP
QDYN_NEQ

# Fail on the first error, so a failed mpirun can't exit 0 and be recorded COMPLETED. Set after
# `module`, which returns non-zero on some stacks just for unload warnings.
set -eo pipefail

# Pin one rank per core: above two ranks mpirun maps by socket, which misplaces ranks when SLURM
# splits the allocated cores unevenly across sockets.
mpi_map="--map-by core --bind-to core"

rundir=$workdir/$temperature/$run_num
mkdir -p $rundir
cd $rundir || exit

echo "Running NEQ in $rundir"
echo "Parameters T=$temperature, replicate=$run_num, seed=$seed, neq_reps=$neq_reps"
echo

cp $inputfiles/eq*.inp .
cp $inputfiles/relax_*.inp .
cp $inputfiles/neq_*.inp .
cp $inputfiles/*.top .
cp $inputfiles/$fepfile .

# Custom random seed, temperature and FEP file for this realization
sed -i "s/SEED_VAR/$seed/" eq1.inp
sed -i "s/T_VAR/$temperature/" *.inp
sed -i "s/FEP_VAR/$fepfile/" *.inp

# 1) Equilibration eq1 -> eq5 with the MPI engine across all cores (fixed lambda)
for i in 1 2 3 4 5; do
    time mpirun -np $ncores $mpi_map $qdynp eq$i.inp > eq$i.log
done

# 2) Endpoint equilibration chain (MPI): one continuous trajectory per endpoint,
# saving a checkpoint per replicate to decorrelate the switch starting points.
# The first iteration is a longer one-time relaxation (relax_${s}.inp) that settles the
# nearly-decoupled ligand at each endpoint before any switching; later iterations use the
# shorter tEQ spacing (eq6_${s}.inp). State 0 sits at lambda (0 1), state 1 at (1 0).
cp eq5.re eq6_0_prev.re
cp eq5.re eq6_1_prev.re
for rep in $(seq 0 $((neq_reps - 1))); do
    for s in 0 1; do
        if [ "$rep" -eq 0 ]; then eqsrc=relax_${s}.inp; else eqsrc=eq6_${s}.inp; fi
        sed "s|RESTARTFILE|eq6_${s}_prev.re|; s|FINALFILE|eq6_${s}_${rep}.re|" \
            "$eqsrc" > eq6_${s}_run${rep}.inp
        time mpirun -np $ncores $mpi_map $qdynp eq6_${s}_run${rep}.inp > eq6_${s}_${rep}.log
        cp eq6_${s}_${rep}.re eq6_${s}_prev.re
    done
done

# 3) Lambda switches with the serial engine, one per core. Each switch only needs its own eq6
# checkpoint, so they are independent and run concurrently: mpirun launches one rank per switch,
# pins each to its own core, and the launcher routes each rank to its (input, log) pair.
# State 1 (1->0) is the forward work logged as neq_1; state 0 (0->1) is the reverse work
# logged as neq_0.
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
$qdyn "\$1" > "\$2"
EOF
chmod +x neq_launch.sh

nsw=$(wc -l < switch_list.txt)
for (( off=0; off<nsw; off+=ncores )); do
    remaining=$(( nsw - off ))
    np=$(( remaining < ncores ? remaining : ncores ))
    export BATCH_OFFSET=$off
    # A batch that loses ranks must not abort the others; the count below decides what to do.
    time mpirun -x BATCH_OFFSET $mpi_map -np $np ./neq_launch.sh || true
done

# qdyn's serial abort still exits 0, so a completed run is confirmed by "terminated normally" in
# its log, not by the exit code.

# Equilibration must finish in full -- with no eq6 checkpoint there is nothing to switch from.
eq_bad=""
for log in eq*.log; do
    grep -q "terminated normally" "$log" 2>/dev/null || eq_bad="$eq_bad $log"
done
if [ -n "$eq_bad" ]; then
    echo "ERROR: replicate $run_num (T=$temperature, seed=$seed): equilibration incomplete:$eq_bad"
    exit 1
fi

# A lost switch is tolerated -- analyze_neq discards it (analyze_neq._scan_leg). A partial loss is
# just sampling noise; only a total loss (e.g. mpirun never launched) is fatal.
nsw_ok=$(grep -l "terminated normally" neq_*.log 2>/dev/null | wc -l || true)
if [ "$nsw_ok" -eq 0 ]; then
    echo "ERROR: replicate $run_num (T=$temperature, seed=$seed): no switch produced work."
    exit 1
elif [ "$nsw_ok" -lt "$nsw" ]; then
    echo "WARNING: replicate $run_num: $((nsw - nsw_ok))/$nsw switches lost (tolerated)."
fi

#CLEANUP
