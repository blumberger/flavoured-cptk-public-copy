#!/bin/bash

for i in run-fssh-*;
do
        echo $i
        cd $i
        sed -i -e "/PPI_FILE_NAME/a \                        \PSIGMA_C_S_FILE_NAME             ..\/topologies\/S_psigma_carbon_sulfur.dat" run.inp
        sed -i -e "/PSIGMA_C_S_FILE_NAME/a \                        \PPI_C_S_FILE_NAME                ..\/topologies\/S_ppi_carbon_sulfur.dat" run.inp
        sed -i -e "/PPI_C_S_FILE_NAME/a \                        \PSIGMA_S_S_FILE_NAME             ..\/topologies\/S_psigma_sulfur_sulfur.dat" run.inp
        sed -i -e "/PSIGMA_S_S_FILE_NAME/a \                        \PPI_S_S_FILE_NAME                ..\/topologies\/S_ppi_sulfur_sulfur.dat" run.inp
        cd ../
done
