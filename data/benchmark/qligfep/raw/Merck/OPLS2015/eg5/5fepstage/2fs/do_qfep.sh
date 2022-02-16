#!/bin/bash
cwd=$(pwd)
for leg in */ ; do
    echo $(pwd)
    cd $leg
    for d in */ ; do
        cd $d
        for rep in FEP1/298/*/; do
            cd $rep
            echo $(pwd)
            qfep < qfep.FEP1.inp > qfep1.out
            qfep < qfep.FEP2.inp > qfep2.out
            qfep < qfep.FEP3.inp > qfep3.out
            qfep < qfep.FEP4.inp > qfep4.out
            qfep < qfep.FEP5.inp > qfep5.out
            cd ../../..
            continue
        done
        cd ..
        continue
    done
    cd $cwd
done
