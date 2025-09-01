#BSUB -n 1
#BSUB -J partestjob
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=32]"
#BSUB -W 1:00
#BSUB -q scafellpikeSKL

module load gcc9
/lustre/scafellpike/local/SCP011/rla62/lxs54-rla62/code_extension_V4_DSF/cp2k/exe/local/cp2k.sopt -i run.inp  -o run.out

