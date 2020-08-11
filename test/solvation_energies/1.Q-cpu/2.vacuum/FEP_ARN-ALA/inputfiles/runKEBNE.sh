#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH -A SNIC2018-2-3 
#              d-hh:mm:ss
#SBATCH --time=0-04:00:00

export TMP=/tmp
export TEMP=/tmp
export TMPDIR=/tmp
## Load modules for qdynp
module load gompi/2017b


## define qdynp location
qdyn=/home/w/wije/pfs/software/q/bin/qdynp
fepfiles=(FEP1.fep)
temperature=298
run=10
finalMDrestart=md_0000_1000.re

workdir=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo
inputfiles=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo/inputfiles
length=${#fepfiles[@]}
length=$((length-1))
for index in $(seq 0 $length);do
fepfile=${fepfiles[$index]}
fepdir=$workdir/FEP$((index+1))
mkdir -p $fepdir
cd $fepdir
tempdir=$fepdir/$temperature
mkdir -p $tempdir
cd $tempdir

rundir=$tempdir/$run
mkdir -p $rundir
cd $rundir

cp $inputfiles/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .
cp $inputfiles/run_0500-1000.sh .
cp $inputfiles/run_0500-0000.sh .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp
else
lastfep=FEP$index
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re
fi

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
if [ $index -lt 1 ]; then
#time mpirun -np 16 $qdyn eq1.inp > eq1.log
#EQ_FILES
time mpirun -np 28 $qdyn eq1.inp > eq1.log
time mpirun -np 28 $qdyn eq2.inp > eq2.log
time mpirun -np 28 $qdyn eq3.inp > eq3.log
time mpirun -np 28 $qdyn eq4.inp > eq4.log
time mpirun -np 28 $qdyn eq5.inp > eq5.log
fi
#RUN_FILES
time mpirun -np 28 $qdyn md_1000_0000.inp > md_1000_0000.log
time mpirun -np 28 $qdyn md_0998_0002.inp > md_0998_0002.log
time mpirun -np 28 $qdyn md_0996_0004.inp > md_0996_0004.log
time mpirun -np 28 $qdyn md_0994_0006.inp > md_0994_0006.log
time mpirun -np 28 $qdyn md_0991_0009.inp > md_0991_0009.log
time mpirun -np 28 $qdyn md_0989_0011.inp > md_0989_0011.log
time mpirun -np 28 $qdyn md_0986_0014.inp > md_0986_0014.log
time mpirun -np 28 $qdyn md_0983_0017.inp > md_0983_0017.log
time mpirun -np 28 $qdyn md_0979_0021.inp > md_0979_0021.log
time mpirun -np 28 $qdyn md_0976_0024.inp > md_0976_0024.log
time mpirun -np 28 $qdyn md_0971_0029.inp > md_0971_0029.log
time mpirun -np 28 $qdyn md_0967_0033.inp > md_0967_0033.log
time mpirun -np 28 $qdyn md_0961_0039.inp > md_0961_0039.log
time mpirun -np 28 $qdyn md_0955_0045.inp > md_0955_0045.log
time mpirun -np 28 $qdyn md_0948_0052.inp > md_0948_0052.log
time mpirun -np 28 $qdyn md_0940_0060.inp > md_0940_0060.log
time mpirun -np 28 $qdyn md_0930_0070.inp > md_0930_0070.log
time mpirun -np 28 $qdyn md_0919_0081.inp > md_0919_0081.log
time mpirun -np 28 $qdyn md_0905_0095.inp > md_0905_0095.log
time mpirun -np 28 $qdyn md_0888_0112.inp > md_0888_0112.log
time mpirun -np 28 $qdyn md_0867_0133.inp > md_0867_0133.log
time mpirun -np 28 $qdyn md_0838_0162.inp > md_0838_0162.log
time mpirun -np 28 $qdyn md_0800_0200.inp > md_0800_0200.log
time mpirun -np 28 $qdyn md_0744_0256.inp > md_0744_0256.log
time mpirun -np 28 $qdyn md_0657_0343.inp > md_0657_0343.log
time mpirun -np 28 $qdyn md_0500_0500.inp > md_0500_0500.log
time mpirun -np 28 $qdyn md_0343_0657.inp > md_0343_0657.log
time mpirun -np 28 $qdyn md_0256_0744.inp > md_0256_0744.log
time mpirun -np 28 $qdyn md_0200_0800.inp > md_0200_0800.log
time mpirun -np 28 $qdyn md_0162_0838.inp > md_0162_0838.log
time mpirun -np 28 $qdyn md_0133_0867.inp > md_0133_0867.log
time mpirun -np 28 $qdyn md_0112_0888.inp > md_0112_0888.log
time mpirun -np 28 $qdyn md_0095_0905.inp > md_0095_0905.log
time mpirun -np 28 $qdyn md_0081_0919.inp > md_0081_0919.log
time mpirun -np 28 $qdyn md_0070_0930.inp > md_0070_0930.log
time mpirun -np 28 $qdyn md_0060_0940.inp > md_0060_0940.log
time mpirun -np 28 $qdyn md_0052_0948.inp > md_0052_0948.log
time mpirun -np 28 $qdyn md_0045_0955.inp > md_0045_0955.log
time mpirun -np 28 $qdyn md_0039_0961.inp > md_0039_0961.log
time mpirun -np 28 $qdyn md_0033_0967.inp > md_0033_0967.log
time mpirun -np 28 $qdyn md_0029_0971.inp > md_0029_0971.log
time mpirun -np 28 $qdyn md_0024_0976.inp > md_0024_0976.log
time mpirun -np 28 $qdyn md_0021_0979.inp > md_0021_0979.log
time mpirun -np 28 $qdyn md_0017_0983.inp > md_0017_0983.log
time mpirun -np 28 $qdyn md_0014_0986.inp > md_0014_0986.log
time mpirun -np 28 $qdyn md_0011_0989.inp > md_0011_0989.log
time mpirun -np 28 $qdyn md_0009_0991.inp > md_0009_0991.log
time mpirun -np 28 $qdyn md_0006_0994.inp > md_0006_0994.log
time mpirun -np 28 $qdyn md_0004_0996.inp > md_0004_0996.log
time mpirun -np 28 $qdyn md_0002_0998.inp > md_0002_0998.log
time mpirun -np 28 $qdyn md_0000_1000.inp > md_0000_1000.log
timeout 30s /home/w/wije/pfs/software/q/bin/qfep < qfep.inp > qfep.out
done
