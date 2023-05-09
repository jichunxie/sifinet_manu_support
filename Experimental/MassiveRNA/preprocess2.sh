#!/bin/bash
#SBATCH -e slurm.err
#SBATCH --mem=200G

ndata=10
nsample=1306127
nfeature=27998
npd=$(((nsample + ndata-1)/ndata))


for ((i=0; i<9; i++)) 
do  
    nfeature=27998
    npd=$(((nsample + ndata-1)/ndata))	
    x=$(wc -l < matrix"$i".mtx)
    sed -e '1s/^/\%\%MatrixMarket matrix coordinate integer general\n'$npd'\t'$nfeature'\t'"$x"'\n/' matrix"$i".mtx > matrix"$i"_h.mtx
done 

x=$(($(wc -l < matrix9.mtx)-1))
a=$(((nsample-1)%npd+1))
sed -e '1s/.*/\%\%MatrixMarket matrix coordinate integer general\n'$a'\t'$nfeature'\t'"$x"'/' matrix9.mtx > matrix9_h.mtx 
