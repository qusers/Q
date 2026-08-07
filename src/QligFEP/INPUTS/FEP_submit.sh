#!/bin/sh
set -eu

workdir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
inputfiles="$workdir/inputfiles"
submitfile="$inputfiles/RUNFILE"
submission_record="$workdir/submission-jobid.txt"

if [ -s "$submission_record" ]; then
    echo "Already submitted as job $(cat "$submission_record")"
    exit 0
fi

if [ "$#" -eq 0 ]; then
    echo "No arguments provided - submitting full job array (all replicates)"
    job_id=$(sbatch --parsable "$submitfile")
else
    array_indices=$(printf '%s\n' "$*" | tr ' ' ',')
    echo "Submitting job array indexes: $array_indices"
    job_id=$(sbatch --parsable --array="$array_indices" "$submitfile")
fi

printf '%s\n' "$job_id" > "$submission_record.partial"
mv "$submission_record.partial" "$submission_record"
echo "Submitted Slurm job $job_id"
