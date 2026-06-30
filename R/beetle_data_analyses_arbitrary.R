# here we plot the inference results of the beetle mimicry data
setwd("/home/u29/yichaozeng/Desktop")
library(treeducken)
library(reshape2)
library(patchwork)
library(GGally)
library(ggridges)
library(ggpubr)
library(ggbreak)
library(dplyr)

sim_time <- 2 # this is the time of simulations
n_bin <- 150 # this is the Number of summary statistics for each of the three normalized distributions

# this variable differentiate the real run from the test runs
real_data_run <- 1

para_est <- NULL
SS_est <- NULL
panels_OG <- NULL

folder_ids <-  ''
folder_ids_slash <-  '/'
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

rel_err <- NULL

# which summary statistics to use
# these include the three normalized distributions - the BLenD distribution, the distribution of difference in centrality between each host-symbiont pair, and the distribution of symbiont centrality (Number of actual hosts devided by the Number of potential hosts)
# #1:ncol(SS) if all SS are used

SS_sel <- list(
  c(301:302), # sizes only
  c(1:300), # BLenD only
  c(1:300, 301:302) # combined
)[[3]]

setwd("/home/u29/yichaozeng/Desktop")
source('cophy_ABC/R/ABC_parameter_estimation.R')
source('cophy_ABC/R/actual_run.R')

para_est[[0 + 1]] <- list(para_sim_acc)
SS_est[[0 + 1]] <- list(SS_sim_acc)


test <- data.frame(
  H_S = log(para_sim_acc$lambda_H / para_sim_acc$lambda_S),
  H_C = log(para_sim_acc$lambda_H / para_sim_acc$lambda_C),
  H_W = log(para_sim_acc$lambda_H / para_sim_acc$exp_H),
  S_C = log(para_sim_acc$lambda_S / para_sim_acc$lambda_C),
  S_W = log(para_sim_acc$lambda_S / para_sim_acc$exp_H),
  C_W = log(para_sim_acc$lambda_C / para_sim_acc$exp_H)
)

boxplot(test)


test <- para_sim_acc
test$lambda_H <- log(test$lambda_H)
test$lambda_S <- log(test$lambda_S)
test$lambda_C <- log(test$lambda_C)
test$exp_H <- log(test$exp_H)

pairs(test)

library(tidyr)

test_long <- pivot_longer(
  test,
  cols = everything(),
  names_to = "comparison",
  values_to = "log_ratio"
)


library(ggplot2)

ggplot(test_long, aes(x = comparison, y = log_ratio)) +
  geom_violin(fill = "lightblue", color = "black") #+
  # geom_boxplot(width = 0.1, outlier.size = 0.5)
