#!/bin/bash
#SBATCH -e slurm.err
#SBATCH --mem=200G

ndata=10
nsample=1306127
npd=$(((nsample + ndata-1)/ndata))


for ((i=0; i<ndata; i++)) 
do      
    awk -v i="$i" -v npd="$npd" '$2<=npd*(i+1)&&$2>npd*i {print (($2-1)%npd+1)"\t"$1"\t"$3}' 1M_matrix.mtx > matrix"$i".mtx
done 
