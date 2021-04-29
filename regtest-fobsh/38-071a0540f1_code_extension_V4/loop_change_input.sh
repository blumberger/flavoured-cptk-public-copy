#!/bin/bash

for i in run-fssh-*;
do
        echo $i
        cd $i
            sed -i -e 's/NUMBER_ATOMS_PER_SITE/NUMBER_AOM_ATOMS_PER_SITE/g' run.inp
            #mv run.inp run_old.inp
            #cp ../run.inp .
        cd ../
done
