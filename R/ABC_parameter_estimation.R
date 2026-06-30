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
if(exists('cros_val') == 1){
  
  para_real <- list(
    para_ABC[cros_val_ids,1],
    para_ABC[cros_val_ids,2],
    para_ABC[cros_val_ids,3],
    para_ABC[cros_val_ids,4],
    para_ABC[cros_val_ids,5],
    para_ABC[cros_val_ids,6]
  )
  
  sizes_real <- list(
    SS_comb[cros_val_ids,301],
    SS_comb[cros_val_ids,302]
  )
  
  SS_real <- split(SS_ABC[cros_val_ids, ], seq_len(nrow(SS_ABC[cros_val_ids, ])))
  SS_real <- lapply(SS_real, function(x) as.numeric(x[1, ]))
  
  # similarly, SS_real abd sizes_real
  
}else if(exists('real_data_run') == 1){
  # read in real data
  # an important realization here is that scaling branch lengths of the cophylogeny changes neither the BLenD nor the size-based summary statistics. So there is no need to rescale the beetle cophylogeny here!
  cophy_real <- list(list.load('cophy_ABC/R/real_data/cophy_real.rds'))

  # here, have a list containing multiple vectors of SSs
  SS_real <- foreach(i = 1:length(cophy_real)) %dopar% {
    
    #print(i)
    unlist(c(
      SS_norm(list(cophy_real[[i]])),
      SS_size(list(cophy_real[[i]]))
    ))
    
  }
  
  # save the raw summary statistics of the real cophylogeny
  SS_real_raw <- SS_real
  
  # rescale the summary statistics
  SS_real <- lapply(SS_real, FUN = renormalize)
  SS_real <- lapply(SS_real, FUN = rescale, sd_vec = sd_vec)
  
  # the size measures of the "real" cophylogenies
  sizes_real <- SS_size(cophy_real)
  
}


para_sim <- para_ABC[,c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H', 'mu_H_frac', 'mu_S_frac')]
SS_sim <- SS_ABC

para_sim_md <- para_sim
SS_sim_md <- SS_sim
