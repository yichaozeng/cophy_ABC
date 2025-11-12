#!/bin/bash
#SBATCH --job-name=relative_error_sizes
#SBATCH --account=cromanpa
#SBATCH --partition=standard
#SBATCH --ntasks=94
#SBATCH --nodes=1
#SBATCH --time=40:00:00

module load R/4.4.0
Rscript -e "source('Desktop/cophy_ABC/R/relative_error_tree_size.R'); "