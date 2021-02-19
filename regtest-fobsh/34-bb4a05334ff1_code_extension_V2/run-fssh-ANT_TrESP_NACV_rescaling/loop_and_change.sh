#!/bin/bash


for i in `seq 0 9`;   #TO CHANGE!!!!
do
        echo $i
        cd run-fssh-$i
        sed -i -e 's/METHOD_RESCALING                      SIMPLE/METHOD_RESCALING                      NACV/g' run.inp
        cd ../
done

