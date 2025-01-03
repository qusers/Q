#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH --mem-per-cpu=256  # more than enough for 25A sphere size FEP
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
restartfile=md_0000_1000.re
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
fepfiles=(FEPS)

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

if [ ${#seeds[@]} -eq 0 ]; then
    echo "Error: No seeds specified"
    exit 1
fi

# Calculate indices based on array task ID
NUMTEMPS=${#temperatures[@]}
TASKID=$((SLURM_ARRAY_TASK_ID * 2 - 2)) # substract 2 for 0-based indexing

# Calculate which temperature and run this job corresponds to
temp_idx=$((TASKID / runs))
run_num1=$((TASKID % runs + 1))
run_num2=$(((TASKID + 1) % runs + 1))

# Validate indices
if [ $temp_idx -ge ${#temperatures[@]} ]; then
    echo "Error: Temperature index ($temp_idx) out of bounds"
    exit 1
fi

if [ $((run_num - 1)) -ge ${#seeds[@]} ]; then
    echo "Error: Run number ($run_num) exceeds number of seeds"
    exit 1
fi

# Get the actual values for this job
temperature=${temperatures[$temp_idx]}

## Load modules for qdynp
MODULES

## define qdynp location
QDYN

starttime=$(date +%s)
starttime_readable=$(date)

length=${#fepfiles[@]}
length=$((length-1))
for index in $(seq 0 $length); do
fepfile=${fepfiles[$index]}
fepdir=$workdir/FEP$((index+1))
mkdir -p $fepdir
cd $fepdir || exit
tempdir=$fepdir/$temperature
mkdir -p $tempdir
cd $tempdir || exit

rundirs=("$tempdir/$run_num1" "$tempdir/$run_num2")
jobseeds=("${seeds[$run_num1-1]}" "${seeds[$run_num2-1]}")

echo "Running FEP$((index+1)) at $temperature K with seeds ${jobseeds[0]} and ${jobseeds[1]}"
echo "Run directory 1 is ${rundirs[0]}"
echo "Run directory 2 is ${rundirs[1]}"

# Create directories and copy files
for i in "${!rundirs[@]}"; do
dir="${rundirs[$i]}"
seed="${seeds[$i]}"

mkdir -p "$dir"
cd "$dir" || exit

cp "$inputfiles"/md*.inp .
cp "$inputfiles"/*.top .
cp "$inputfiles"/qfep.inp .
cp "$inputfiles"/$fepfile .

if [ $index -lt 1 ]; then
cp "$inputfiles"/eq*.inp .
sed -i "s/SEED_VAR/$seed/" eq1.inp # change the random seed to custom
else
    lastfep=FEP$index
    cp "$workdir/$lastfep/$temperature/$run_num/$finalMDrestart" ./eq5.re
fi

sed -i "s/T_VAR/$temperature/" *.inp
sed -i "s/FEP_VAR/$fepfile/" *.inp

cd - || exit
done

# Just to make sure twin jobs will run on different cores
(srun -n 8 --cpu-bind=rank bash -c 'echo "Job 1 cores:"; taskset -pc $$') &
(srun -n 8 --cpu-bind=rank bash -c 'echo "Job 2 cores:"; taskset -pc $$') &
wait

if [ $index -lt 1 ]; then
#EQ_FILES
fi
#RUN_FILES
timeout 3m QFEP < qfep.inp > qfep.out || [ $? -eq 124 ]
done
#CLEANUP

endtime=$(date +%s)
endtime_readable=$(date)
# Calculate runtime
runtime=$((endtime - starttime))

# Convert runtime to hours:minutes:seconds
hours=$(($runtime / 3600))
minutes=$(($runtime % 3600 / 60))
seconds=$(($runtime % 60))

echo "#    EXPRESS LOG for jobid: $SLURM_JOB_ID"
echo "#    Slurm tasks: $SLURM_NTASKS"
echo "#    Starttime: $starttime_readable"
echo "#    Endtime: $endtime_readable"
echo "#    Runtime: ${hours}h:${minutes}m:${seconds}s"
echo "#    Random seed: $seed"
echo "#    Replicate Number: $run"
echo "#    Working Directory: $workdir"