#! /bin/bash

workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
echo $workdir
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/RUNFILE

# No SLURM scheduler available: run the run script directly. The run script runs
# exactly one replicate per invocation - the orchestrator starts a container per
# replicate to get several - so there is no job array to dispatch and any
# command-line arguments are ignored.
echo "No scheduler - running a single replicate"
bash $submitfile
