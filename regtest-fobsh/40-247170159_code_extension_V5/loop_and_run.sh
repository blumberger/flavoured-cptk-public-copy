#!/bin/bash

for i in run-fssh-*;
do
        echo $i
        cd $i
              #/scratch/sgiannini/flavoured-cptk/cp2k/exe/local/cp2k.sopt run.inp > run.log & 
              #/scratch/sgiannini/flavoured-cptk/cp2k/exe/local_hetero/cp2k.sopt run.inp > run.log & 
             /home/sam/flavoured-cptk/cp2k/exe/cp2k_local_dell/cp2k.sopt run.inp > run.log 
        cd ../
done

wait

### check the end ###
for i in run-fssh-*;
do
        echo $i
        cd $i
        if grep -q 'T I M I N G' run.log;
         then
            echo "Here are the Strings with the Pattern 'T I M I N G':"
            #echo -e "$(grep $PATTERN $FILE)\n"
            cd ../
         else
            echo "Error: The Pattern was NOT Found in $i"
            echo "Exiting..."
            #exit 0
            cd ../
            #mv run-fssh-$i NOT_FINISHED
        fi
done
