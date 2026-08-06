#! /bin/bash
#
# One entry point, two roles.
#
# Run from a shell, this submits the edge's replicates as a job array - what it
# has always done. Handed to sbatch by an orchestrator, it *is* the batch job and
# runs the single replicate it was allocated. That second role is what makes this
# file the entry point for the SLURM architectures in the same way the
# schedulerless submit scripts are for theirs: whoever drives the run always
# points at FEP_submit.sh and never has to know which script underneath it does
# the work.
#
# Without it, sbatch would submit this wrapper as the job and every array task
# would call sbatch again - one whole array per task - while the submitter was
# left holding the wrapper's job ID rather than the real one. The wrapper exits
# in seconds and writes no sentinel, so the run would report itself finished
# while the actual jobs were still queueing.
#
# Only directives the submitter cannot know belong in the header below. The
# partition, time limit, array bounds and job name are the caller's: by hand they
# come from the run script's own header, and under an orchestrator they come from
# its sbatch command line.
#SBATCH --nodes=NODES
#SBATCH --mem-per-cpu=512  # more than enough for 25A sphere size FEP
#SBATCH -o slurm.run%a.%N.%j.out

# Locating the edge takes a different answer in each role, because the two roles
# run this file from two different places.
#
# From a shell it sits in the edge directory, so its own path finds the tree.
#
# As the batch job it does not. SLURM does not execute the submitted file in
# place: it copies the script and runs that copy out of the node's spool
# directory, so BASH_SOURCE points at something like
# /var/spool/slurmd/job00123/slurm_script and $workdir/inputfiles is nowhere
# near the tree. The working directory is the thing that does survive, because
# --chdir sets it, and the run script derives its own working directory from it
# too — so taking it from anywhere else lets the two disagree about which edge
# they are running.
#
# Not SLURM_SUBMIT_DIR: it is documented as following --chdir but does not do so
# on every SLURM build, and where it does not it reports the directory sbatch
# was invoked from, which for a remote submitter is the login user's home.
if [ -n "$SLURM_JOB_ID" ]; then
    workdir="$PWD"
else
    workdir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
fi
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/RUNFILE

# Checked rather than left to `exec` so the failure names the tree it looked in.
# The run script owns the .failed sentinel, so nothing reaching this line has
# written one — an orchestrator sees only a job that vanished, and on a cluster
# without accounting storage that is indistinguishable from one still queued.
if [ ! -f "$submitfile" ]; then
    echo "FEP_submit.sh: no run script at $submitfile" >&2
    echo "FEP_submit.sh: expected an edge directory containing inputfiles/RUNFILE" >&2
    if [ -n "$SLURM_JOB_ID" ]; then
        echo "FEP_submit.sh: working directory is $PWD - was sbatch given --chdir=<edge>?" >&2
    fi
    exit 1
fi

if [ -n "$SLURM_JOB_ID" ]; then
    # Already running as the batch job. $SLURM_ARRAY_TASK_ID is inherited, and
    # the run script reads it to work out which replicate this task is.
    echo "Running as SLURM job $SLURM_JOB_ID, array task ${SLURM_ARRAY_TASK_ID:-1}"
    exec bash "$submitfile"
fi

# --chdir is passed explicitly because the run script derives its working
# directory from the current one, so submitting from anywhere else would run it
# against the wrong tree.
if [ $# -eq 0 ]; then # if no arguments provided, submit the whole job array
    echo "No arguments provided - submitting full job array (all replicates)"
    sbatch --chdir="$workdir" "$submitfile"
else
    array_indices=$(echo "$@" | tr ' ' ',') # convert space-separated arguments to comma-separated
    echo "Submitting job array indexes: $array_indices"
    sbatch --chdir="$workdir" --array=$array_indices "$submitfile"
fi
