#!/bin/bash

## This is an old bash to check the outputs, however it's better to use the regtest_code.py attached 

file1=$1
file2=$2

echo $file1
echo $file2
echo " " 

#outfiles=(run-coeff-1.xyz run-pseudo-hamilt-1.xyz run-sh-1.log pos-init.xyz)
outfiles=(run-coeff-1.xyz run-pseudo-hamilt-1.xyz)
#directory=(run-fssh-0 run-fssh-BOMD run-fssh-ADIAB_INIT run-fssh-DLI_MTS1 run-fssh-ISOTROP run-fssh-LONG run-fssh-MTS1 run-fssh-MTS10 run-fssh-MTS1_EFH run-fssh-MTS1_FAST run-fssh-MTS1_FBD run-fssh-MTS1_NFF run-fssh-NO_REVERSAL run-fssh-NOSPT_NOREOR run-fssh-UNMOD_SH run-fssh-MTS1_NODECO run-fssh-MTS1_RESTART run-fssh-MTS1_SIMPLEST run-fssh-FROZEN_C_FIRST_OFF run-fssh-MTS1_FROZEN_H run-fssh-MTS1_NFF_FROZEN_C_FIRST_OFF)
directory=(run-fssh-MTS1 run-fssh-MTS1_EFH run-fssh-MTS1_FBD run-fssh-DLI_MTS1 run-fssh-MTS1_EFH run-fssh-MTS1_FBD run-fssh-MTS1_FAST)

for dir in ${directory[@]};
do
   for t in ${outfiles[@]};
   do 
      echo "file 1:" $file1/$dir/$t 
      echo "file 2:" $file2/$dir/$t
 
      outfile1=$file1/$dir/$t
      outfile2=$file2/$dir/$t  
      
      if cmp -s "$outfile1" "$outfile2"; then
         printf "PASSED!"
      else 
         printf "FAILED!"
      fi
      echo " " 
      echo " " 
   done
done


echo " " 

#for dir in $file1/run-fssh-*;
#do 
#   echo "DIRECTORY:" $dir 
#   for filename in $dir/run-*;
#   do 
#        echo "FILE:" $filename 
#        
#       # if cmp -s "$file1" "$file2"; then
#       #    printf 'The file "%s" is the same as "%s"\n' "$file1" "$file2"
#       # fi
#   done
#   
#   echo " "     
#done
