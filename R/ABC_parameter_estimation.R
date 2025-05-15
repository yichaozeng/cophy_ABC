library(abc)
library(ape)
library(rlist)
library(dplyr)
library(pbapply)
library(RANN)
library(reshape2)
library(combinat)
library(pracma)
library(phytools)
library(treeducken)

# for parallelization
library(future)
library(doParallel)
library(parallel)
n.cores <- future::availableCores()
registerDoParallel(cores=n.cores - 3)

setwd('/home/u29/yichaozeng/Desktop')
source('cophy_ABC/R/functions.R')

# remove unnecessary SSs (if they are present)
SS$X <- NULL
SS$X.1 <- NULL
SS$X.2 <- NULL
SS$X.3 <- NULL
SS$NODF_z <- NULL
SS$mod_z <- NULL
SS$mantel_u_host <- NULL

# keep only the "good" simulations that are complete and whose SSs can be computed
para_ABC <- para_dat_separate

# remove rows where any variable is NA or (+ or -)Inf, keeping only the finites
SS_fin <- SS

# after parameters to generate observational data are selected, select only some of the SSs
SS_ABC <- SS_fin

# SS_real
if(exists('real_data_run') == 1){
  cophy_real <- list(list.load('cophy_ABC/R/real_data/cophy_real.rds'))
}else{
  cophy_real <- sim_cophyBD(hbr = lambda_H_real, hdr = mu_H_frac * (lambda_H_real + lambda_C_real), sbr = lambda_S_real, sdr = mu_S_frac * (lambda_S_real + lambda_C_real + exp_H_real), cosp_rate =lambda_C_real, host_exp_rate = exp_H_real, time_to_sim = 2, hs_mode = T, numbsim = 500)
}

# compute the summary statistics
tr_ht <- cophy_real[[1]][[1]]$root.edge + max(nodeHeights(cophy_real[[1]][[1]])) # tree height
breaks <- seq(from = -tr_ht, to = tr_ht, length.out = n_bin + 1) # length.out should be no. of bins + 1

SS_real <- c(
  SS_norm(cophy_real, breaks = breaks, tr_ht = tr_ht),
  SS_size(cophy_real)
)

# rescale the summary statistics
SS_real <- renormalize(SS_real)
SS_real <- rescale(SS_real, sd_vec)

# the real parameter values for the test run
if(exists('real_data_run') == 0){ # if this is a test run
  para_real <- para_ABC[1, ] # for formatting
  
  para_real$lambda_H <- lambda_H_real
  para_real$lambda_S <- lambda_S_real
  para_real$lambda_C <- lambda_C_real
  para_real$exp_H <- exp_H_real
  para_real$mu_H_frac
  para_real$mu_S_frac
  
  print(para_real)
}

para_sim <- para_ABC[,c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')]
SS_sim <- SS_ABC

para_sim_md <- para_sim
SS_sim_md <- SS_sim
