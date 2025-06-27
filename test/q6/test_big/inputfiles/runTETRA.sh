#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -A naiss2023-3-5 
#              d-hh:mm:ss
#SBATCH --time=0-24:00:00
#SBATCH -J p_ejm_31_ejm_42
#SBATCH -o slurm.%N.%j.out # STDOUT

## Load modules for qdynp
#buildenv-gcc/2023a-eb


## define qdynp location
qdyn=/Q/src/q6/bin/q6/qdyn
fepfiles=(FEP1.fep)
temperature=298
run=1
finalMDrestart=md_0000_1000.re
seed="$[1 + $[RANDOM % 32767]]"
starttime=$(date +%s)
starttime_readable=$(date)
workdir=/proj/uucompbiochem/users/x_davfi/openFF-QligFEP/perturbations/Tyk2/2.protein/FEP_ejm_31_ejm_42
inputfiles=/proj/uucompbiochem/users/x_davfi/openFF-QligFEP/perturbations/Tyk2/2.protein/FEP_ejm_31_ejm_42/inputfiles
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

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
sed -i s/SEED_VAR/$seed/ eq1.inp # change the random seed to custom
else
lastfep=FEP$index
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re
fi

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
if [ $index -lt 1 ]; then
#EQ_FILES
time $qdyn eq1.inp > eq1.log
time $qdyn eq2.inp > eq2.log
time $qdyn eq3.inp > eq3.log
time $qdyn eq4.inp > eq4.log
time $qdyn eq5.inp > eq5.log
fi
#RUN_FILES
time $qdyn md_0500_0500.inp > md_0500_0500.log

time $qdyn md_0592_0408.inp > md_0592_0408.log
time $qdyn md_0408_0592.inp > md_0408_0592.log

time $qdyn md_0657_0343.inp > md_0657_0343.log
time $qdyn md_0343_0657.inp > md_0343_0657.log

time $qdyn md_0706_0294.inp > md_0706_0294.log
time $qdyn md_0294_0706.inp > md_0294_0706.log

time $qdyn md_0744_0256.inp > md_0744_0256.log
time $qdyn md_0256_0744.inp > md_0256_0744.log

time $qdyn md_0775_0225.inp > md_0775_0225.log
time $qdyn md_0225_0775.inp > md_0225_0775.log

time $qdyn md_0800_0200.inp > md_0800_0200.log
time $qdyn md_0200_0800.inp > md_0200_0800.log

time $qdyn md_0821_0179.inp > md_0821_0179.log
time $qdyn md_0179_0821.inp > md_0179_0821.log

time $qdyn md_0838_0162.inp > md_0838_0162.log
time $qdyn md_0162_0838.inp > md_0162_0838.log

time $qdyn md_0854_0146.inp > md_0854_0146.log
time $qdyn md_0146_0854.inp > md_0146_0854.log

time $qdyn md_0867_0133.inp > md_0867_0133.log
time $qdyn md_0133_0867.inp > md_0133_0867.log

time $qdyn md_0878_0122.inp > md_0878_0122.log
time $qdyn md_0122_0878.inp > md_0122_0878.log

time $qdyn md_0888_0112.inp > md_0888_0112.log
time $qdyn md_0112_0888.inp > md_0112_0888.log

time $qdyn md_0897_0103.inp > md_0897_0103.log
time $qdyn md_0103_0897.inp > md_0103_0897.log

time $qdyn md_0905_0095.inp > md_0905_0095.log
time $qdyn md_0095_0905.inp > md_0095_0905.log

time $qdyn md_0912_0088.inp > md_0912_0088.log
time $qdyn md_0088_0912.inp > md_0088_0912.log

time $qdyn md_0919_0081.inp > md_0919_0081.log
time $qdyn md_0081_0919.inp > md_0081_0919.log

time $qdyn md_0925_0075.inp > md_0925_0075.log
time $qdyn md_0075_0925.inp > md_0075_0925.log

time $qdyn md_0930_0070.inp > md_0930_0070.log
time $qdyn md_0070_0930.inp > md_0070_0930.log

time $qdyn md_0935_0065.inp > md_0935_0065.log
time $qdyn md_0065_0935.inp > md_0065_0935.log

time $qdyn md_0940_0060.inp > md_0940_0060.log
time $qdyn md_0060_0940.inp > md_0060_0940.log

time $qdyn md_0944_0056.inp > md_0944_0056.log
time $qdyn md_0056_0944.inp > md_0056_0944.log

time $qdyn md_0948_0052.inp > md_0948_0052.log
time $qdyn md_0052_0948.inp > md_0052_0948.log

time $qdyn md_0952_0048.inp > md_0952_0048.log
time $qdyn md_0048_0952.inp > md_0048_0952.log

time $qdyn md_0955_0045.inp > md_0955_0045.log
time $qdyn md_0045_0955.inp > md_0045_0955.log

time $qdyn md_0958_0042.inp > md_0958_0042.log
time $qdyn md_0042_0958.inp > md_0042_0958.log

time $qdyn md_0961_0039.inp > md_0961_0039.log
time $qdyn md_0039_0961.inp > md_0039_0961.log

time $qdyn md_0964_0036.inp > md_0964_0036.log
time $qdyn md_0036_0964.inp > md_0036_0964.log

time $qdyn md_0967_0033.inp > md_0967_0033.log
time $qdyn md_0033_0967.inp > md_0033_0967.log

time $qdyn md_0969_0031.inp > md_0969_0031.log
time $qdyn md_0031_0969.inp > md_0031_0969.log

time $qdyn md_0971_0029.inp > md_0971_0029.log
time $qdyn md_0029_0971.inp > md_0029_0971.log

time $qdyn md_0974_0026.inp > md_0974_0026.log
time $qdyn md_0026_0974.inp > md_0026_0974.log

time $qdyn md_0976_0024.inp > md_0976_0024.log
time $qdyn md_0024_0976.inp > md_0024_0976.log

time $qdyn md_0978_0022.inp > md_0978_0022.log
time $qdyn md_0022_0978.inp > md_0022_0978.log

time $qdyn md_0979_0021.inp > md_0979_0021.log
time $qdyn md_0021_0979.inp > md_0021_0979.log

time $qdyn md_0981_0019.inp > md_0981_0019.log
time $qdyn md_0019_0981.inp > md_0019_0981.log

time $qdyn md_0983_0017.inp > md_0983_0017.log
time $qdyn md_0017_0983.inp > md_0017_0983.log

time $qdyn md_0985_0015.inp > md_0985_0015.log
time $qdyn md_0015_0985.inp > md_0015_0985.log

time $qdyn md_0986_0014.inp > md_0986_0014.log
time $qdyn md_0014_0986.inp > md_0014_0986.log

time $qdyn md_0988_0013.inp > md_0988_0013.log
time $qdyn md_0013_0988.inp > md_0013_0988.log

time $qdyn md_0989_0011.inp > md_0989_0011.log
time $qdyn md_0011_0989.inp > md_0011_0989.log

time $qdyn md_0990_0010.inp > md_0990_0010.log
time $qdyn md_0010_0990.inp > md_0010_0990.log

time $qdyn md_0991_0009.inp > md_0991_0009.log
time $qdyn md_0009_0991.inp > md_0009_0991.log

time $qdyn md_0993_0007.inp > md_0993_0007.log
time $qdyn md_0007_0993.inp > md_0007_0993.log

time $qdyn md_0994_0006.inp > md_0994_0006.log
time $qdyn md_0006_0994.inp > md_0006_0994.log

time $qdyn md_0995_0005.inp > md_0995_0005.log
time $qdyn md_0005_0995.inp > md_0005_0995.log

time $qdyn md_0996_0004.inp > md_0996_0004.log
time $qdyn md_0004_0996.inp > md_0004_0996.log

time $qdyn md_0997_0003.inp > md_0997_0003.log
time $qdyn md_0003_0997.inp > md_0003_0997.log

time $qdyn md_0998_0002.inp > md_0998_0002.log
time $qdyn md_0002_0998.inp > md_0002_0998.log

time $qdyn md_0999_0001.inp > md_0999_0001.log
time $qdyn md_0001_0999.inp > md_0001_0999.log

time $qdyn md_1000_0000.inp > md_1000_0000.log
time $qdyn md_0000_1000.inp > md_0000_1000.log

timeout 30s /Q/src/q6/bin/q6/qfep < qfep.inp > qfep.out
done
#CLEANUP
rm -f *.dcd
#Cleaned .dcd files


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