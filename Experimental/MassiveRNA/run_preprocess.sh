#!/bin/bash
#SBATCH -e slurm.err
#SBATCH --mem=200G
module load R/4.0.3-rhel8
python3 getdata.py 
sh download_data.sh
sh preprocess1.sh
sh preprocess2.sh
R CMD BATCH preprocess3.R
R CMD BATCH preprocess4.R  
R CMD BATCH split_data.R
R CMD BATCH get_genename.R
