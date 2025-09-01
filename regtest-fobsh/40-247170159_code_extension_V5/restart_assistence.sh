#!/bin/bash

#version 1.1 
#update to read different starting points for different trajs

#This script reads the restart file produced by the usual FOB-SH algorithm (default restart file), and adds the necessary keywords for FOB-SH restart. 
#Please NOTE "TO CHANGE" comments and adapt to the case at hand.  

#for i in MTS1_RESTART;   #TO CHANGE!!!! starting trajectory index (e.g run-fssh-0), final traj index (e.g. run-fssh-300)
#do
#    echo $i
    coup_type=AOM
    run_take=run-fssh-MTS1_FAST_V3
    i=1
    #mv run-fssh-MTS1_RESTART run-fssh-$i
    mv $run_take run-fssh-$i
    cd run-fssh-$i
   
    #select the restart file name to use
    restart_file="run-1.restart.bak-1"
    
    updated_restart="run-1.restart.bak"
    
    #select step number  to do in total
    totstep=200 ###########TO CHANGE!! 
      
    step=$( grep "STEP_START_VAL" $restart_file | awk {'print $NF'} )
    states=$( grep "NUMBER_DIABATIC_STATES" $restart_file | awk {'print $NF'} )
    
    echo "n. states is: $states"
    echo "step considered is: $step "
    
    #one must be very careful with this line (choose the file to truncate, they must exist):
    cd ../
     #The following like MUST take the file names created by the initial FOB-SH run
    ./truncate2 $i $i $step adiab adiab_pop coeff couplings pop pseudo-hamilt site  ##### TO CHANGE IF NECESSARY!!! 
    cd run-fssh-$i
    
    #grep diab coefficient 
    grep -A $states "i = $step" run-coeff-1.xyz | tail -$states > DIAB_COEFF.include
    
    active_state=$( sed -n -e "/i = $step/,/Final/ p" run-sh-1.log | tail -1 | awk {'print $4'} )

    #check if &WAVEFUNCTION_RESTART is already present in the restart file from the previus step
   if grep -q WAVEFUNCTION_RESTART $restart_file; then 
      echo "&WAVEFUNCTION_RESTART IS ALREADY PRESENT"
      sed '/\&WAVEFUNCTION_RESTART/,/\&END WAVEFUNCTION_RESTART/d' $restart_file > temp.txt
      mv temp.txt $restart_file
   else 
       echo "&WAVEFUNCTION_RESTART NOT PRESENT"
   fi 

    #echo section in the new restart file
    wf_section="     &WAVEFUNCTION_RESTART
            RESTART_KEY  T
             @INCLUDE   DIAB_COEFF.include
            ACTIVE_STATE_RESTART  $active_state
          &END WAVEFUNCTION_RESTART"
    
    echo "$wf_section" > SECTION.txt 
    if [ $coup_type = TRESP ]; then 
       # get active state from run-sh.log with print more T of F
       sed '/&END TRESP\_COUPLINGS/r SECTION.txt' $restart_file > $updated_restart
    elif [ $coup_type = AOM ]; then 
       # get active state from run-sh.log with print more T of F
       sed '/&END AOM/r SECTION.txt' $restart_file > $updated_restart
    else 
       echo "NO METHOD FOUND FOR INSERTING WF"
    fi
    rm SECTION.txt
    
    #modify updated restart file with remaining number of steps
    line_to_change=$( grep " STEPS " $updated_restart ) 
    echo "line to change: $line_to_change"
    #remaining_steps=$(( $( grep " STEPS " $updated_restart | awk {'print $NF'} ) - $step ))
    remaining_steps=$(( $totstep - $step ))
    echo "step yet to do: $remaining_steps"
    
    sed -i -e "s/$line_to_change/STEPS    $remaining_steps/g" $updated_restart
    echo "####################################################################"
    cd ../
    mv run-fssh-$i $run_take 
#done


