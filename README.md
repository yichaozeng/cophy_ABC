# Macroevolutionary Rates of Species Interactions: Approximate Bayesian Inference from Cophylogenies

## Description
This project presents an Approximate Bayesian Computation (ABC) approach to estimating speciation rates from cophylogenetic system. The rates of four types of speciation can be inferred: host speciation ($\lambda_H$), symbiont speciation without host switching ($\lambda_S$), cospeciation ($\lambda_C$), and symbiont speciation with host switching ($\lambda_W$). This approach does not infer the rates of host extinction and symbiont extinction separately. Instead, it assumes that the relative extinctin rates ($\epsilon_H$ for the hosts and $\epsilon_S$ for the symbionts) are known. Under this assumption, the absolute extinction rates ($\mu_H$ for the hosts and $\mu_S$ for the symbionts) can be obtained by 

$$
\mu_H = \epsilon_H*(\lambda_H + \lambda_C)
$$

$$
\mu_S = \epsilon_S*(\lambda_S + \lambda_C + \lambda_W)
$$

For the user's information, some parameters are assigned a different name in the code: $\lambda_W$ - `exp_H`, $\epsilon_H$ - `mu_H_frac`, and $\epsilon_S$ - `mu_S_frac`.

The tutorial provided here is the full procedure to estimate speciation rates from a given cophylogeny, using both the BLenD curve and tree sizes are summary statistics. This procedure is illustrated with a re-analysis of the cophylogenetic dataset from Van Dam et al. (2024). Our convergence analyses can be reproduced by applying this procedure to simulated cophylogenies whose true parameter values are known and using the BLenD curve the tree sizes alone as summary statistics.

## Table of Contents
- [Setup](#setup)
- [Simulating cophylogenies](#simulating-cophylogenies)
- [ABC](#abc)
- [References](#references)
- [Contact](#contact)

## Setup
The provided code works on a cluster running on SLURM. The `batch_command` folder contain commands you can copy and paste into your SSH terminal. The `SLURM` folder contain scripts with which to submit your jobs. The `R` folder contain the R scripts you will need to run the analyses. Before starting, put all three folders in your work directory.

## Simulating cophylogenies
Cophylogenies are simulated with the `sim_cophyBD` function from the `treeducken` package (Dismukes & Heath 2022). We start by creating several folders named `ex_cophy_sims`, `ex_cophy_plots`, `ex_cophy_statistics`, `ex_cophy_checks`, `ex_cophy_convergence`, and `ex_cophy_results` in yout working directory. Then, run `R/cophy_sim_system.R`. The R script is written in a way that avoids new simulations to overwrite older ones, so you can run this R script simultanesouly on multiple jobs to generate large numbers of simulations in a relatively short period of time.

When your simulations are complete, run `R/archive_older_data.R` to put all simulations into a subfolder within `ex_cophy_sims`. The subfoler will be randomly named by a number; use this number as your `folder_id`, which will be need in subsequent steps.

If needed, the assumed relative extinction rates ($\epsilon_H$ for the hosts and $\epsilon_S$ for the symbionts) can be modifed in `R/cophy_sim_system.R`. The prior distributions of parameters can be modified in `R/cophy_sim.R`.

## ABC
The real data from Van Dam et al. (2024) is stored in `R/real_data/cophy_real.rds`. The host and symbiont trees meet these requirements:
1. Both contains a root edge (branch).
2. Both trees are of the same height (form the root to the tips, including the length of the root edge).
3. Both have been rescaled such that their height is equal to that of the cophylogenies simulated with `treeducken`.

To perform the ABC for this real cophylogenetic dataset with the simulations that we generated in the last step, specify your `folder_id` (you can have multiple of them corresponding to different assumptions or prior distributions) on Lines 22 & 23 in `R/beetle_data_analyses.R`. Then, run `R/beetle_data_analyses.R`. The parameter estimates will be stored in `cophy_ABC_results/real_para_est.rds`. Visualizations of parameter estimates and posterior predictive checks will be stored in `cophy_ABC_results/beetle_results.pdf`, with examples provided below.

Parameter estimates:
<p align="center">
    <img src=images/img1.png width="350">
</p>

Posterior predictive checks (where the line represents the real data, and the dots represent the accepted simulations):
<p align="center">
    <img src=images/img2.png width="700">
</p>


## References
Dismukes, W., & Heath, T. A. (2021). treeducken: An R package for simulating cophylogenetic systems. Methods in Ecology and Evolution, 12(8), 1358-1364.

Van Dam, M. H., Parisotto, A., Medina, M. N., Cabras, A. A., Gutiérrez-Trejo, N., Wilts, B. D., & Lam, A. W. (2024). Biogeography confounds the signal of cospeciation in Batesian mimicry. Current Biology, 34(23), 5554-5563.

## Citation
[To be updated when the preprint or the journal article gets published.]

## Contact
- Any questions or comments are welcome and should be sent to Yichao Zeng (yichaozeng44@gmail.com).
