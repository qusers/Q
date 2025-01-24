#! /bin/bash

workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
echo $wrkdir
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/RUNFILE

if [ $# -eq 0 ]; then # if no arguments provided, submit the whole job array
    echo "No arguments provided - submitting full job array (all replicates)"
    sbatch $submitfile
else
    array_indices=$(echo "$@" | tr ' ' ',') # convert space-separated arguments to comma-separated
    echo "Submitting job array indexes: $array_indices"
    sbatch --array=$array_indices $submitfile
fi
