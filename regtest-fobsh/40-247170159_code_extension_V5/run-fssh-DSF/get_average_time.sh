#!/bin/bash

grep "Alpha parameter" run.log
grep "Real Space Cutoff" run.log 
grep "G-space max. Miller index" run.log 

#grep "CPU TIME [s]" run.log | awk '{sum += $8} END {if (NR > 0) print sum / NR; else print "No data"}'
grep "CPU TIME" run.log | tail -1 
