#!/bin/bash
#SBATCH -e slurm.err
#SBATCH --mem=200G
module load R/4.0.3-rhel8
R CMD BATCH method.R
