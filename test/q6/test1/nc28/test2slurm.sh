#!/bin/bash
# Remember to change here the number of nodes to what is available in the queue
# you're calling. Perhaps you'll also need to issue the --exclusive option.
sbatch -A snic2018-2-3 -n 28 --exclusive -t 01:00:00 -J testqdynp run_test_mpi.sh

