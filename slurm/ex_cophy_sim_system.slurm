#!/bin/bash
#SBATCH --job-name=ex_cophy_simulation_system
#SBATCH --account=cromanpa
#SBATCH --partition=standard
#SBATCH --ntasks=94
#SBATCH --nodes=1
#SBATCH --time=30:00:00

module load gnu8
module load R/4.4.0
Rscript -e "source('Desktop/cophy_ABC/R/cophy_sim_system.R')"