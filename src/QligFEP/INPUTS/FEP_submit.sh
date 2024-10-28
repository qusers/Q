#! /bin/bash

workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
echo $wrkdir
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/RUNFILE
sbatch $submitfile
