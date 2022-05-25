#!/bin/bash
# This sends requests to start test1 with a changing number of cores.
#
# NOTE1: You may need to change the project ID to something relevant
# to your use.
#
# NOTE2: Don't forget to change run_test_mpi.sh to load the appropriate module.
#
for i in $( seq 4 4 32 )
do
    cd nc$i
    sbatch -A snic2018-2-3 -n $i --exclusive -t 01:00:00 -J testqdynp run_test_mpi.sh
    cd ..
done
