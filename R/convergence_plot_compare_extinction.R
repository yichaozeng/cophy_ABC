# here we plot the convergence for different parameterizations
setwd("/home/u29/yichaozeng/Desktop")
library(treeducken)
library(reshape2)

sim_time <- 2 # this is the time of simulations
n_bin <- 150 # this is the number of summary statistics for each of the three normalized distributions

panels_OG <- NULL
panels_OG_percentage <- NULL # this recrods the percentage of convergence, defined as the percentage of unnormalized cophylogenies that results in parameter convergence

for(ext_id in 1:4){
  
  print(ext_id)
  
  # these, in turn, contain simulations for which extinction fractions are 0/0, 0.7/0, 0/0.7, 0.7/0.7
  folder_ids <-  list('41', '87', '28', '64')[[ext_id]] #list('86', '36', '0', '64')[[ext_id]]
  folder_ids_slash <-  list('41/', '87/', '28/', '64/')[[ext_id]] #list('86/', '36/', '0/', '64/')[[ext_id]]
  prefix <- "ex_"
  
  # allow simulations from two or more folders to be loaded
  SS_comb <- NULL
  para_dat_separate_comb <- NULL
  
  # this loop helps you read in the simulations if they are stored across multiple folders
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
  
  mu_H_frac <- para[1,5]
  mu_S_frac <- para[1,6]
  
  residuals_OG <- NULL
  para_real_OG <- NULL
  
  for (row_id in 1:4) {
    
    print(c(ext_id, row_id))
    
    # which summary statistics to use
    SS_sel <- c(1:150, 451:453) #1:ncol(SS) #these include the three normalized distributions - the BLenD distribution, the distribution of difference in centrality between each host-symbiont pair, and the distribution of symbiont centrality (number of actual hosts devided by the number of potential hosts)
    
    lambda_H_real <- c(2.4, 0.3, 0.3, 0.3)[row_id]
    lambda_S_real <- c(0.3, 2.4, 0.3, 0.3)[row_id]
    lambda_C_real <- c(0.3, 0.3, 2.4, 0.3)[row_id]
    exp_H_real <- c(0.3, 0.3, 0.3, 2.4)[row_id]
    
    source('cophy_ABC/R/ABC_parameter_estimation.R')
    source('cophy_ABC/R/convergence_checks.R')
    
    residuals_OG[[length(residuals_OG) + 1]] <- list(res_dat_comb)
    para_real_OG <- rbind(para_real_OG, para_real)
    
    panels_OG_percentage[[length(panels_OG_percentage) + 1]] <- conv_perc
    
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
    
    summary_comb <- NULL
    
    res_dat_comb <- residuals_OG[[row_id]][[1]]
    
    res_dat_comb_melt <- reshape2::melt(res_dat_comb, id.vars = "run_name", 
                      measure.vars = c("lambda_H", "lambda_S", "lambda_C", "exp_H"),
                      variable.name = "parameter", value.name = "residual"
                      )
    
    all_para_plot <- ggplot(res_dat_comb_melt, aes(x = run_name, y = residual, fill = parameter)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(y = "Residual", x = "Selectivity")
    
    panels_OG[[length(panels_OG) + 1]] <- all_para_plot
  }

}


panels <- panels_OG

conv_plot <- ggarrange(panels[[1]], panels[[2]], panels[[3]], panels[[4]], 
                       panels[[5]], panels[[6]], panels[[7]], panels[[8]], 
                       panels[[9]], panels[[10]], panels[[11]], panels[[12]], 
                       panels[[13]], panels[[14]], panels[[15]], panels[[16]], 
                       nrow = 4, ncol = 4
)
png(paste(prefix, 'cophy_ABC_convergence/', "convergence_comb.png", sep = ''), width = 15, height = 20, units = 'in', pointsize = 12, res = 300)
print(conv_plot)
dev.off()

# save the list of plots and the list of percetages
file_name <- paste(prefix, 'cophy_ABC_convergence/panels.rds', sep = '')
list.save(panels, file_name)

file_name <- paste(prefix, 'cophy_ABC_convergence/panels_percentage.rds', sep = '')
list.save(panels_OG_percentage, file_name)

temp <- unlist(panels_OG_percentage)
print(
  rbind(temp[1:4], temp[5:8], temp[9:12], temp[13:16])
)