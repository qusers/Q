#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=512  # more than enough for 25A sphere size FEP
#              d-hh:mm:ss
#SBATCH --time=0-08:00:00
#SBATCH -p rome
#SBATCH -J p_jmc_28_jmc_30
#SBATCH --array=1-10
#SBATCH -o slurm.run%a.%N.%j.out

# Define your parameters
temperatures=(298)
seeds=(2924 25360 21448 14380 14188 28133 2816 22850 6601 3085)
runs=${#seeds[@]}
restartfile=md_0000_1000.re
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
fepfiles=(FEP1.fep)

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
TID=$((SLURM_ARRAY_TASK_ID - 1))  # Convert to 0-based indexing

# Calculate which temperature and run this job corresponds to
temp_idx=$((TID / runs))
run_num=$((TID % runs + 1))

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
seed=${seeds[$run_num-1]}

## Load modules for qdynp
module purge
module load 2023
module load AOCC/4.0.0-GCCcore-12.3.0
module load OpenMPI/4.1.5-GCC-12.3.0


## define qdynp location
qdyn=/zfsdata/data/david/Q/src/q6/bin/q6/qdynp

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

rundir=$tempdir/$run_num
mkdir -p $rundir
cd $rundir || exit

echo "Running job in $rundir"
echo "Parameters T=$temperature, replicate=$run_num, seed=$seed"
echo

echo "=== Process Binding Information ==="
echo "slurm tasks per node: $SLURM_TASKS_PER_NODE"
echo "slurm cpus per node: $SLURM_JOB_CPUS_PER_NODE"
echo

echo -e "\n=== CPU Model Information ==="
lscpu | grep -E "Model name|Architecture|CPU op|Thread|Core|Socket|NUMA|CPU(s)"
echo

echo -e "\n=== Available CPU List ==="
cpu_list=$(cat /proc/self/status | grep Cpus_allowed_list | awk '{print $2}')
echo "CPU list: $cpu_list"
echo

cp $inputfiles/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
sed -i "s/SEED_VAR/$seed/" eq1.inp # change the random seed to custom
else
    lastfep=FEP$index
    cp $workdir/$lastfep/$temperature/$run_num/$finalMDrestart $rundir/eq5.re
fi
sed -i "s/T_VAR/$temperature/" *.inp
sed -i "s/FEP_VAR/$fepfile/" *.inp
if [ $index -lt 1 ]; then
#EQ_FILES
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn eq1.inp > eq1.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn eq2.inp > eq2.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn eq3.inp > eq3.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn eq4.inp > eq4.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn eq5.inp > eq5.log
fi
#RUN_FILES
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0500_0500.inp > md_0500_0500.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0592_0408.inp > md_0592_0408.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0408_0592.inp > md_0408_0592.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0657_0343.inp > md_0657_0343.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0343_0657.inp > md_0343_0657.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0706_0294.inp > md_0706_0294.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0294_0706.inp > md_0294_0706.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0744_0256.inp > md_0744_0256.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0256_0744.inp > md_0256_0744.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0775_0225.inp > md_0775_0225.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0225_0775.inp > md_0225_0775.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0800_0200.inp > md_0800_0200.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0200_0800.inp > md_0200_0800.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0821_0179.inp > md_0821_0179.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0179_0821.inp > md_0179_0821.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0838_0162.inp > md_0838_0162.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0162_0838.inp > md_0162_0838.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0854_0146.inp > md_0854_0146.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0146_0854.inp > md_0146_0854.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0867_0133.inp > md_0867_0133.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0133_0867.inp > md_0133_0867.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0878_0122.inp > md_0878_0122.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0122_0878.inp > md_0122_0878.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0888_0112.inp > md_0888_0112.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0112_0888.inp > md_0112_0888.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0897_0103.inp > md_0897_0103.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0103_0897.inp > md_0103_0897.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0905_0095.inp > md_0905_0095.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0095_0905.inp > md_0095_0905.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0912_0088.inp > md_0912_0088.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0088_0912.inp > md_0088_0912.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0919_0081.inp > md_0919_0081.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0081_0919.inp > md_0081_0919.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0925_0075.inp > md_0925_0075.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0075_0925.inp > md_0075_0925.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0930_0070.inp > md_0930_0070.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0070_0930.inp > md_0070_0930.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0935_0065.inp > md_0935_0065.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0065_0935.inp > md_0065_0935.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0940_0060.inp > md_0940_0060.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0060_0940.inp > md_0060_0940.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0944_0056.inp > md_0944_0056.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0056_0944.inp > md_0056_0944.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0948_0052.inp > md_0948_0052.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0052_0948.inp > md_0052_0948.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0952_0048.inp > md_0952_0048.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0048_0952.inp > md_0048_0952.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0955_0045.inp > md_0955_0045.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0045_0955.inp > md_0045_0955.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0958_0042.inp > md_0958_0042.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0042_0958.inp > md_0042_0958.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0961_0039.inp > md_0961_0039.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0039_0961.inp > md_0039_0961.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0964_0036.inp > md_0964_0036.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0036_0964.inp > md_0036_0964.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0967_0033.inp > md_0967_0033.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0033_0967.inp > md_0033_0967.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0969_0031.inp > md_0969_0031.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0031_0969.inp > md_0031_0969.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0971_0029.inp > md_0971_0029.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0029_0971.inp > md_0029_0971.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0974_0026.inp > md_0974_0026.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0026_0974.inp > md_0026_0974.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0976_0024.inp > md_0976_0024.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0024_0976.inp > md_0024_0976.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0978_0022.inp > md_0978_0022.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0022_0978.inp > md_0022_0978.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0979_0021.inp > md_0979_0021.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0021_0979.inp > md_0021_0979.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0981_0019.inp > md_0981_0019.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0019_0981.inp > md_0019_0981.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0983_0017.inp > md_0983_0017.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0017_0983.inp > md_0017_0983.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0985_0015.inp > md_0985_0015.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0015_0985.inp > md_0015_0985.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0986_0014.inp > md_0986_0014.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0014_0986.inp > md_0014_0986.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0988_0013.inp > md_0988_0013.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0013_0988.inp > md_0013_0988.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0989_0011.inp > md_0989_0011.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0011_0989.inp > md_0011_0989.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0990_0010.inp > md_0990_0010.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0010_0990.inp > md_0010_0990.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0991_0009.inp > md_0991_0009.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0009_0991.inp > md_0009_0991.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0993_0007.inp > md_0993_0007.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0007_0993.inp > md_0007_0993.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0994_0006.inp > md_0994_0006.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0006_0994.inp > md_0006_0994.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0995_0005.inp > md_0995_0005.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0005_0995.inp > md_0005_0995.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0996_0004.inp > md_0996_0004.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0004_0996.inp > md_0004_0996.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0997_0003.inp > md_0997_0003.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0003_0997.inp > md_0003_0997.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0998_0002.inp > md_0998_0002.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0002_0998.inp > md_0002_0998.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0999_0001.inp > md_0999_0001.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0001_0999.inp > md_0001_0999.log

time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_1000_0000.inp > md_1000_0000.log
time mpirun -n $SLURM_NTASKS --bind-to core $qdyn md_0000_1000.inp > md_0000_1000.log

timeout 3m /zfsdata/data/david/Q/src/q6/bin/q6/qfep < qfep.inp > qfep.out || [ $? -eq 124 ]
done
#CLEANUP
rm -r *dcd *inp
#Cleaned dcd inp files


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
echo "#    Replicate Number: $run_num"
echo "#    Working Directory: $workdir"