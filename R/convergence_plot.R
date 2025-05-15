# here we plot the convergence for different parameterizations
setwd("/home/u29/yichaozeng/Desktop")

# specify the folder id
# these in turn correspond to simulations where the extinction fractions are 0.3/0.3, 0.3/0.7, 0.7/0.3, 0.7/0.7
ex_ind <- 1; SMC <- 0 # these are obsolete variables and do not mean anything

###############################################################
ext_id <- 1
folder_ids <- list('86', '36')[[ext_id]]
folder_ids_slash <- list('86/', '36/')[[ext_id]]
prefix <- "ex_"

###############################################################

# allow simulations from two or more folders to be loaded
SS_comb <- NULL
para_dat_separate_comb <- NULL

for (pos in 1:length(folder_ids)) {
  
  folder_id <- folder_ids[pos]
  folder_id_slash <- folder_ids_slash[pos]
  source('cophy_ABC/R/parameter_distribution.R')
  
  # read in the computed the SSs
  SS <- read.csv(file = paste(prefix, "cophy_ABC_statistics/", folder_id_slash, "SS.csv", sep = ''))
  para_dat_separate <- read.csv(file = paste(prefix, "cophy_ABC_statistics/", folder_id_slash, "para_dat_separate.csv", sep = ''))
  
  SS_comb <- rbind(SS_comb, SS)
  para_dat_separate_comb <- rbind(para_dat_separate_comb, para_dat_separate)
  
}

para_dat_separate <- para_dat_separate_comb

# (optional) rescale the summary statistics
sd_vec <- c(
  sum(apply(SS_comb[, 1:150], 2, sd)), # BLenD
  sum(apply(SS_comb[, 151:300], 2, sd)), # difference in centrality
  sum(apply(SS_comb[, 301:450], 2, sd)), # centrality of symbionts
  sd(SS_comb[, 451]), # host tree size
  sd(SS_comb[, 452]), # symbiont tree size
  sd(SS_comb[, 453]) # network size
)

rescale <- function(vec, sd_vec){
  return(c(
    vec[1:150] / sd_vec[1], # BLenD
    vec[151:300] / sd_vec[2], # difference in centrality
    vec[301:450] / sd_vec[3], # centrality of symbionts
    vec[451] / sd_vec[4], # host tree size
    vec[452] / sd_vec[5], # symbiont tree size
    vec[453] / sd_vec[6] # network size
  ))
}

renormalize <- function(vec){
  return(c(
    vec[1:150] / sum(vec[1:150]),
    vec[151:300] / sum(vec[151:300]),
    vec[301:450] / sum(vec[301:450]),
    vec[451:453]
  ))
}

SS <- as.data.frame(t(apply(SS_comb, 1, renormalize)))
SS <- as.data.frame(t(apply(SS, 1, function(x) rescale(vec = x, sd_vec = sd_vec))))
# SS <- SS_comb

mu_H_frac <- para[1,5]
mu_S_frac <- para[1,6]

panels_OG <- NULL
residuals_OG <- NULL
para_real_OG <- NULL

#rm(list = setdiff(ls(), c("panels_OG", "residuals_OG", "para_real_OG", "ext_id", "SMC")))

for (row_id in 2:5) {
  
  print(c(ext_id, row_id))
  
  # which summary statistics to use
  SS_sel <- 1:ncol(SS) # these include the three normalized distributions - the BLenD distribution, the distribution of difference in centrality between each host-symbiont pair, and the distribution of symbiont centrality (number of actual hosts devided by the number of potential hosts)
  
  # different true parameters
  lambda_H_real <- c(1, 2.5, 0.3, 0.3, 0.3)[row_id]
  lambda_S_real <- c(1, 0.3, 2.5, 0.3, 0.3)[row_id]
  lambda_C_real <- c(1, 0.3, 0.3, 2.5, 0.3)[row_id]
  exp_H_real <- c(1, 0.3, 0.3, 0.3, 2.5)[row_id]
  
  # different size (the number of association in the cophylogeny) ranges
  h_lower <- 59 #62 * 0.8
  h_upper <- Inf #62 * 1.2
  s_lower <- 29 #30 * 0.8
  s_upper <- Inf #30 * 1.2
  size_lower <- 99 #102 * 0.8
  size_upper <- Inf #102 * 1.2
  
  source('cophy_ABC/R/ABC_parameter_estimation.R')
  source('cophy_ABC/R/convergence_checks.R')
  
  residuals_OG[[length(residuals_OG) + 1]] <- list(res_dat_comb)
  para_real_OG <- rbind(para_real_OG, para_real)
  
}

# save the true parameters used
file_name <- paste(prefix, 'cophy_ABC_convergence/', folder_id_slash, "real_para.csv", sep = '')
write.csv(para_real_OG, file = file_name)

# save the convergence plot object
file_name <- paste(prefix, 'cophy_ABC_convergence/', folder_id_slash, "residuals.rds", sep = '')
list.save(residuals_OG, file_name)

# read in the dataset for plotting
residuals_OG <- list.load(file_name)


# arrange the panels
library(ggpubr)

for (row_id in 1:4) {
  
  res_dat_comb <- residuals_OG[[row_id]][[1]]
  

  summary <- res_dat_comb %>%
    group_by(run_name) %>%
    summarise(
      median_value = median(lambda_H),      # Mean (or use median(value) for median)
      lower = quantile(lambda_H, 0.95), # 5th percentile
      upper = quantile(lambda_H, 0.05)  # 95th percentile
    )
  lambda_H_plot <- ggplot(summary, aes(x = run_name, y = median_value)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #theme_minimal() +
    labs(y = "lambda_H", x = "Run")
  
  summary <- res_dat_comb %>%
    group_by(run_name) %>%
    summarise(
      median_value = median(lambda_S),      # Mean (or use median(value) for median)
      lower = quantile(lambda_S, 0.95), # 5th percentile
      upper = quantile(lambda_S, 0.05)  # 95th percentile
    )
  lambda_S_plot <- ggplot(summary, aes(x = run_name, y = median_value)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #theme_minimal() +
    labs(y = "lambda_S", x = "Run")
  
  summary <- res_dat_comb %>%
    group_by(run_name) %>%
    summarise(
      median_value = median(lambda_C),      # Mean (or use median(value) for median)
      lower = quantile(lambda_C, 0.95), # 5th percentile
      upper = quantile(lambda_C, 0.05)  # 95th percentile
    )
  lambda_C_plot <- ggplot(summary, aes(x = run_name, y = median_value)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #theme_minimal() +
    labs(y = "lambda_C", x = "Run")
  
  summary <- res_dat_comb %>%
    group_by(run_name) %>%
    summarise(
      median_value = median(exp_H),      # Mean (or use median(value) for median)
      lower = quantile(exp_H, 0.95), # 5th percentile
      upper = quantile(exp_H, 0.05)  # 95th percentile
    )
  exp_H_plot <- ggplot(summary, aes(x = run_name, y = median_value)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #theme_minimal() +
    labs(y = "exp_H", x = "Run")
  
  panels_OG[[length(panels_OG) + 1]] <- list(lambda_H_plot, lambda_S_plot, lambda_C_plot, exp_H_plot)
}

panels <- panels_OG
# for(i in 1:5){
#   for(j in 1:4){
#     panels[[i]][[j]] <- panels_OG[[i]][[j]] + xlim(-4,0) + scale_x_reverse()
#   }
# }

conv_plot <- ggarrange(panels[[1]][[1]], panels[[1]][[2]], panels[[1]][[3]], panels[[1]][[4]], 
                       panels[[2]][[1]], panels[[2]][[2]], panels[[2]][[3]], panels[[2]][[4]], 
                       panels[[3]][[1]], panels[[3]][[2]], panels[[3]][[3]], panels[[3]][[4]], 
                       panels[[4]][[1]], panels[[4]][[2]], panels[[4]][[3]], panels[[4]][[4]], 
                       # panels[[5]][[1]], panels[[5]][[2]], panels[[5]][[3]], panels[[5]][[4]],
                       # nrow = 5, ncol = 4
                       nrow = 4, ncol = 4
)

png(paste(prefix, 'cophy_ABC_convergence/', folder_id_slash, "convergence.png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
print(conv_plot)
dev.off()






library(ggplot2)
library(dplyr)

# Example data
set.seed(123)
df <- data.frame(
  group = rep(c("A", "B", "C"), each = 30),
  value = c(rnorm(30, mean = 10, sd = 2), 
            rnorm(30, mean = 15, sd = 3), 
            rnorm(30, mean = 20, sd = 4))
)

# Summarize data
summary_df <- df %>%
  group_by(group) %>%
  summarise(
    mean_value = mean(value),      # Mean (or use median(value) for median)
    lower = quantile(value, 0.05), # 5th percentile
    upper = quantile(value, 0.95)  # 95th percentile
  )

# Plot using geom_pointrange()
ggplot(summary_df, aes(x = group, y = mean_value)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), color = "blue") +
  theme_minimal() +
  labs(y = "Value", x = "Group", title = "Mean and 5th-95th Percentiles")

