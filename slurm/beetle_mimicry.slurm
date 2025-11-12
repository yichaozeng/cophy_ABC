#!/bin/bash
#SBATCH --job-name=beetle_mimicry
#SBATCH --account=cromanpa
#SBATCH --partition=standard
#SBATCH --ntasks=94
#SBATCH --nodes=1
#SBATCH --time=4:00:00

module load R/4.4.0
Rscript -e "source('Desktop/cophy_ABC/R/beetle_data_analyses.R'); "