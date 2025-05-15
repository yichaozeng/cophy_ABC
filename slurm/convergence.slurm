#!/bin/bash
#SBATCH --job-name=convergence
#SBATCH --account=cromanpa
#SBATCH --partition=standard
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=5:00:00

module load gnu8
module load R/4.4.0
Rscript -e "source('Desktop/cophy_ABC/R/convergence_plot_compare_SS_choice.R'); "