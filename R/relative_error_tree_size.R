# here we plot the convergence for different parameterizations
library(treeducken)
library(reshape2)
library(ggpubr)
library(patchwork)
library(quantreg)

sim_time <- 2 # this is the time of simulations
n_bin <- 150 # this is the number of summary statistics for each of the three normalized distributions

panels_OG <- NULL
# panels_OG_percentage <- NULL # this records the percentage of convergence, defined as the percentage of unnormalized cophylogenies that results in parameter convergence
para_est <- NULL

folder_ids <-  list('27', '87', '78', '40', '90')[[5]] # 0/0, 0.3/0.3, 0.7/0, 0/0.7, 0.7/0.7
folder_ids_slash <-  list('27/', '87/', '78/', '40/', '90/')[[5]]
prefix <- "ex_"

# allow simulations from two or more folders to be loaded
SS_comb <- NULL
para_dat_separate_comb <- NULL

# this loop helps you read in the simulations if they are stored across multiple folders
for (pos in 1:length(folder_ids)) {
  
  folder_id <- folder_ids[pos]
  folder_id_slash <- folder_ids_slash[pos]
  
  setwd("/home/u29/yichaozeng/Desktop")
  source('cophy_ABC/R/parameter_distribution.R')
  
  setwd("/groups/cromanpa/yzeng")
  # re-read in the computed the SSs (optional)
  # SS <- read.csv(file = paste(prefix, "cophy_ABC_statistics/", folder_id_slash, "SS.csv", sep = ''))
  # para_dat_separate <- read.csv(file = paste(prefix, "cophy_ABC_statistics/", folder_id_slash, "para_dat_separate.csv", sep = ''))
  
  SS_comb <- rbind(SS_comb, SS)
  para_dat_separate_comb <- rbind(para_dat_separate_comb, para_dat_separate)
  
}

para_dat_separate <- para_dat_separate_comb

# (maybe optional) rescale the summary statistics
renormalize <- function(vec){
  return(c(
    vec[1:150] / sum(vec[1:150]),
    vec[151:300] / sum(vec[151:300]),
    log(vec[301:303]) # sizes are log-transformed
  ))
}

rescale <- function(vec, sd_vec){
  return(c(
    vec[1:150] / sd_vec[1], # BLenD
    vec[151:300] / sd_vec[2], # BLenD
    vec[301] / sd_vec[3], # host tree size
    vec[302] / sd_vec[4], # symbiont tree size
    vec[303] / sd_vec[5] # network size (not actually used in the inference)
  ))
}

SS <- as.data.frame(t(apply(SS_comb, 1, renormalize)))
sd_vec <- c(
  sum(apply(SS[, 1:150], 2, sd)), # BLenD
  sum(apply(SS[, 151:300], 2, sd)), # BLenD
  sd(SS[, 301]), # host tree size
  sd(SS[, 302]), # symbiont tree size
  sd(SS[, 303]) # network size (not actually used in the inference)
)
SS <- as.data.frame(t(apply(SS, 1, function(x) rescale(vec = x, sd_vec = sd_vec))))

mu_H_frac <- para[1,5]
mu_S_frac <- para[1,6]

rel_err <- NULL

# which summary statistics to use
# these include the three normalized distributions - the BLenD distribution, the distribution of difference in centrality between each host-symbiont pair, and the distribution of symbiont centrality (number of actual hosts devided by the number of potential hosts)
# #1:ncol(SS) if all SS are used

SS_sel <- list(
  c(301:302), # sizes only
  c(1:300), # BLenD only
  c(1:300, 301:302) # combined
)[[3]]

# # here, draw the true rates in a monte-carlo manner
# lambda_H_real <- runif(1, min = 0, max = 2)
# lambda_S_real <- runif(1, min = 0, max = 2)
# lambda_C_real <- runif(1, min = 0, max = 2)
# exp_H_real <- runif(1, min = 0, max = 2)

# here, draw the true rates in a grid-based manner
# lambda_H_real <- c(0.1, 1, 1.9)[1]
# lambda_S_real <- c(0.1, 1, 1.9)[1]
# lambda_C_real <- c(0.1, 1, 1.9)[1]
# exp_H_real <- c(0.1, 1, 1.9)[1]

# the "survival rates"
surv_H <- 1 - mu_H_frac
surv_S <- 1 - mu_S_frac

setwd("/home/u29/yichaozeng/Desktop")

# cross validation
n_test_set <- 1000
n_rep_cophy <- 100

cros_val <- 1

cros_val_ids <- 1:n_test_set * n_rep_cophy # the indices of the SS chosen as "true" SS

source('cophy_ABC/R/ABC_parameter_estimation.R')

source('cophy_ABC/R/convergence_checks.R')
setwd("/groups/cromanpa/yzeng")

# run this line if you want to use relative rather than absolute errors
res_dat_comb$lambda_H <- res_dat_comb$lambda_H / res_dat_comb$lambda_H_true
res_dat_comb$lambda_S <- res_dat_comb$lambda_S / res_dat_comb$lambda_S_true
res_dat_comb$lambda_C <- res_dat_comb$lambda_C / res_dat_comb$lambda_C_true
res_dat_comb$exp_H <- res_dat_comb$exp_H / res_dat_comb$exp_H_true

rel_err <- res_dat_comb

file_name <- paste(prefix, 'cophy_ABC_convergence/cross_validation_sizes',"_", as.character(mu_H_frac), "_", as.character(mu_S_frac), '.csv', sep = '')
write.csv(rel_err, file_name)
