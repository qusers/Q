#! /bin/bash

workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
echo $workdir
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/RUNFILE

# No SLURM scheduler available: run the run script directly. The run script loops
# over every temperature/replicate combination internally, so there is no job
# array to dispatch and any command-line arguments are ignored.
echo "No scheduler - running the full set of replicates sequentially"
bash $submitfile
