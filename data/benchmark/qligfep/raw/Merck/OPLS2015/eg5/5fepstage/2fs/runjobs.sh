#!/bin/bash
cwd=$(pwd)
for leg in */ ; do
    echo $(pwd)
    cd $leg
    for d in */ ; do
        cd $d
        echo "Submitting $leg$d ..."
	    ./FEP_submit.sh
        cd ..
        continue
    done
    cd $cwd
done
