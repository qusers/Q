#!/bin/sh
set -eu

workdir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
inputfiles="$workdir/inputfiles"
submitfile="$inputfiles/RUNFILE"
submission_record="$workdir/submission-jobid.txt"
submission_history="$workdir/submission-history.tsv"

# A repeated argument-free call is usually an accidental duplicate from a campaign
# loop. Explicit array indices intentionally bypass this guard to resume failed tasks.
if [ "$#" -eq 0 ] && [ -s "$submission_record" ]; then
    echo "Already submitted as job $(cat "$submission_record")"
    echo "Pass explicit array indices to resubmit selected replicates."
    exit 0
fi

if [ "$#" -eq 0 ]; then
    selection="all"
    echo "No arguments provided - submitting full job array (all replicates)"
    job_id=$(sbatch --parsable "$submitfile")
else
    array_indices=$(printf '%s\n' "$*" | tr ' ' ',')
    selection="$array_indices"
    echo "Submitting job array indexes: $array_indices"
    job_id=$(sbatch --parsable --array="$array_indices" "$submitfile")
fi

printf '%s\n' "$job_id" > "$submission_record.partial"
mv "$submission_record.partial" "$submission_record"
printf '%s\t%s\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$selection" "$job_id" >> "$submission_history"
echo "Submitted Slurm job $job_id"
