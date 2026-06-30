# here we compile a table showing:
# (1) for each normalized cophylogeny (whose summary statistics are normalized across 500 replicates) the mean, standard deviation of each parameter estimate for each tested cophylogeny
# (2) for each non-normalized cophylogeny (each of the 500 replicates), the percentage of good convergence (defined as the percentage of cophylogenies whose residuals fall within the true value +- 1 event/unit time)
###################################################################

library(rlist)
setwd("/home/u29/yichaozeng/Desktop")

# read in the residuals and convergence percentages
folder_ids <-  c('86', '69', '17', '40', '76') # 0/0, 0.3/0.3, 0.7/0.3, 0.3/0.7, 0.7/0.7
epsilons <- c('0_0', '0.3_0.3', '0.7_0', '0_0.7', '0.7_0.7')

para_est <- NULL
conv_perc <- NULL

for (i in 1:length(folder_ids)) {
  para_est <- c(
    para_est,
    list(list.load(paste('ex_cophy_ABC_convergence/para_est_', as.character(epsilons[i]), '.rds', sep = ''))) # each of these corresponds to a different combination of epsilons
  )
  
  conv_perc <- c(
    conv_perc,
    list(list.load(paste('ex_cophy_ABC_convergence/panels_percentage_', as.character(epsilons[i]), '.rds', sep = ''))) # each of these corresponds to a different combination of epsilons
  )
}

# the indices differ depending on which summary statistics to look at
#ids <- 9:12; table_id <- 'sizes+BLenD' # both sizes and BLenD
#ids <- 1:4; table_id <- 'sizes' # sizes only
ids <- 5:8; table_id <- 'BLenD' # BLenD only

# specify the acceptance rate
# acc_rate <- '1/512'
acc_rate <- '1/4096'

# now we construct data frames to store the convergence percentage data
conv_tab <- matrix(NA, nrow = 20, ncol = 13)

colnames(conv_tab) <- c('epsilon_H',
                        'epsilon_S',
                        'lambda_H',
                        'lambda_S',
                        'lambda_C',
                        'lambda_W',
                        'lambda_H_hat',
                        'lambda_S_hat',
                        'lambda_C_hat',
                        'lambda_W_hat',
                        'conv_perc_1',
                        'conv_perc_2',
                        'lift'
                        )

conv_tab[, 'epsilon_H'] <- c(0, 0, 0, 0, 0.3, 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0, 0, 0, 0, 0.7, 0.7, 0.7, 0.7)
conv_tab[, 'epsilon_S'] <- c(0, 0, 0, 0, 0.3, 0.3, 0.3, 0.3, 0, 0, 0, 0, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)

conv_tab[, 'lambda_H'] <- rep(c(2.4, 0.3, 0.3, 0.3), 5)
conv_tab[, 'lambda_S'] <- rep(c(0.3, 2.4, 0.3, 0.3), 5)
conv_tab[, 'lambda_C'] <- rep(c(0.3, 0.3, 2.4, 0.3), 5)
conv_tab[, 'lambda_W'] <- rep(c(0.3, 0.3, 0.3, 2.4), 5)

# a function that computes the mean and standard deviation and format them as "mean (sd)"
mean_sd <- function(vec){
  return(
    paste(
      round(mean(vec), digits = 2),
      ' (',
      round(sd(vec), digits = 2),
      ')',
      sep = ''
    )
  )
}

# summarize the parameter estimates
para_est_mat <- NULL
for (epsilon_id in 1:5) {
  for (para_id in 1:4) {
    para_est_mat <- rbind(para_est_mat, apply(para_est[[epsilon_id]][[ids[para_id]]][[1]][para_est[[epsilon_id]][[ids[para_id]]][[1]]$run_name == acc_rate, 1:4], MARGIN = 2, FUN = mean_sd))
  }
}

# rearrange the convergence percentages
conv_perc_vec <- NULL
for (epsilon_id in 1:5) {
  conv_perc_vec <- c(conv_perc_vec, conv_perc[[epsilon_id]][ids])
}

# record these data in the matrixs
conv_tab[, c('lambda_H_hat', 'lambda_S_hat', 'lambda_C_hat', 'lambda_W_hat')] <- para_est_mat
conv_tab[, 'conv_perc_1'] <- unlist(lapply(conv_perc_vec, FUN = function(x) x[1]))
conv_tab[, 'conv_perc_2'] <- unlist(lapply(conv_perc_vec, FUN = function(x) x[2]))
conv_tab[, 'lift'] <- as.numeric(conv_tab[, 'conv_perc_2'])/0.25

# remove the old convergence percentage
conv_tab <- conv_tab[, -which(colnames(conv_tab) == 'conv_perc_1')]

# save the convergence plot
write.csv(conv_tab, file = paste('ex_cophy_ABC_convergence/convergence_', as.character(table_id), '.csv', sep = ''))
