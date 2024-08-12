#! /bin/bash

temperatures=(TEMP_VAR)
runs=RUN_VAR
restartfile=md_0000_1000.re
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
echo $wrkdir
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/RUNFILE
seeds=(SEED_VAR)

sed -i s/finalMDrestart=.*/finalMDrestart="$restartfile"/g $submitfile
sed -i s#workdir=.*#workdir="$workdir"#g $submitfile
sed -i s#inputfiles=.*#inputfiles="$inputfiles"#g $submitfile
for temp in ${temperatures[*]};do
sed -i s/temperature=.*/temperature="$temp"/g $submitfile
for i in $(seq 1 $runs);do
current_seed=${seeds[$i-1]}
sed -i s/run=.*/run="$i"/g $submitfile
sed -i s/seed=.*/seed=$current_seed/g $submitfile
sbatch $submitfile
done
done
