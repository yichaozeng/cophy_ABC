# here we plot the convergence for different parameterizations
setwd("/home/u29/yichaozeng/Desktop")
library(treeducken)
library(reshape2)
library(ggpubr)
library(patchwork)

sim_time <- 2 # this is the time of simulations
n_bin <- 150 # this is the number of summary statistics for each of the three normalized distributions

panels_OG <- NULL
panels_OG_percentage <- NULL # this records the percentage of convergence, defined as the percentage of unnormalized cophylogenies that results in parameter convergence
para_est <- NULL

for(SS_choice_id in 1:3){
  
  print(SS_choice_id)
  
  folder_ids <-  list('86', '69', '17', '40', '76')[[5]] # 0/0, 0.3/0.3, 0.7/0, 0/0.7, 0.7/0.7
  folder_ids_slash <-  list('86/', '69/', '17/', '40/', '76/')[[5]]
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
  renormalize <- function(vec){
    return(c(
      vec[1:150] / sum(vec[1:150]),
      vec[151:300] / sum(vec[151:300]),
      vec[301:450] / sum(vec[301:450]),
      log(vec[451:453]) # sizes are log-transformed
    ))
  }
  
  rescale <- function(vec, sd_vec){
    return(c(
      vec[1:150] / sd_vec[1], # BLenD
      vec[151:300] / sd_vec[2], # difference in centrality
      vec[301:450] / sd_vec[3], # centrality of symbionts
      vec[451] / sd_vec[4], # host tree size
      vec[452] / sd_vec[5], # symbiont tree size
      vec[453] / sd_vec[6] # network size (not actually used in the inference)
    ))
  }
  
  SS <- as.data.frame(t(apply(SS_comb, 1, renormalize)))
  sd_vec <- c(
    sum(apply(SS[, 1:150], 2, sd)), # BLenD
    sum(apply(SS[, 151:300], 2, sd)), # difference in centrality
    sum(apply(SS[, 301:450], 2, sd)), # centrality of symbionts
    sd(SS[, 451]), # host tree size
    sd(SS[, 452]), # symbiont tree size
    sd(SS[, 453]) # network size (not actually used in the inference)
  )
  SS <- as.data.frame(t(apply(SS, 1, function(x) rescale(vec = x, sd_vec = sd_vec))))
  
  mu_H_frac <- para[1,5]
  mu_S_frac <- para[1,6]
  
  residuals_OG <- NULL
  para_real_OG <- NULL
  
  for (row_id in 1:4) {
    
    print(c(SS_choice_id, row_id))
    
    # which summary statistics to use
    # these include the three normalized distributions - the BLenD distribution, the distribution of difference in centrality between each host-symbiont pair, and the distribution of symbiont centrality (number of actual hosts devided by the number of potential hosts)
    # #1:ncol(SS) if all SS are used
    
    SS_sel <- list(
      c(451:452), # sizes only
      c(1:150), # BLenD only
      c(1:150, 451:452) # combined
    )[[SS_choice_id]]
    
    lambda_H_real <- c(2.4, 0.3, 0.3, 0.3)[row_id]
    lambda_S_real <- c(0.3, 2.4, 0.3, 0.3)[row_id]
    lambda_C_real <- c(0.3, 0.3, 2.4, 0.3)[row_id]
    exp_H_real <- c(0.3, 0.3, 0.3, 2.4)[row_id]
    
    source('cophy_ABC/R/ABC_parameter_estimation.R')
    source('cophy_ABC/R/convergence_checks.R')
    
    residuals_OG[[length(residuals_OG) + 1]] <- list(res_dat_comb)
    para_real_OG <- rbind(para_real_OG, para_real)
    
    panels_OG_percentage[[length(panels_OG_percentage) + 1]] <- conv_perc
    
    para_est[[length(para_est) + 1]] <- list(para_sim_acc)
    
  }
  
  # save the true parameters used
  file_name <- paste(prefix, 'cophy_ABC_convergence/', folder_id_slash, "real_para.csv", sep = '')
  write.csv(para_real_OG, file = file_name)
  
  # create the panels
  for (row_id in 1:4) {
    
    summary_comb <- NULL
    
    res_dat_comb <- residuals_OG[[row_id]][[1]]
    
    res_dat_comb_melt <- reshape2::melt(res_dat_comb, id.vars = "run_name", 
                      measure.vars = c("lambda_H", "lambda_S", "lambda_C", "exp_H"),
                      variable.name = "parameter", value.name = "residual"
                      )
    
    # The palette with grey:
    cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    all_para_plot <- ggplot(res_dat_comb_melt, aes(x = run_name, y = residual, fill = parameter)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
      geom_boxplot(outlier.size = 0.1) +
      scale_fill_manual(values = cbp1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      # theme_minimal() +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(legend.position = "none") +
      labs(y = "Residual", x = "Tolerance rate")
    
    panels_OG[[length(panels_OG) + 1]] <- all_para_plot
  }

}


panels <- panels_OG

# save the list of plots, list of parameter estimates, and the list of percetages
file_name <- paste(prefix, 'cophy_ABC_convergence/panels',"_", as.character(mu_H_frac), "_", as.character(mu_S_frac), '.rds', sep = '')
list.save(panels, file_name)

file_name <- paste(prefix, 'cophy_ABC_convergence/para_est',"_", as.character(mu_H_frac), "_", as.character(mu_S_frac), '.rds', sep = '')
list.save(para_est, file_name)

file_name <- paste(prefix, 'cophy_ABC_convergence/panels_percentage', "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), '.rds', sep = '')
list.save(panels_OG_percentage, file_name)

# quick preview of the percentages
# temp <- unlist(panels_OG_percentage)
# print(
#   t(rbind(temp[1:4], temp[5:8], temp[9:12]))
# )


############################################
# further modify the appearance of the the panels for public
for (p_id in 1:length(panels)) {
  panels[[p_id]] <- panels[[p_id]] +
    geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
    scale_fill_manual(
      values = cbp1,
      name = "Parameter estimates:",
      labels = c("lambda_H" = expression(hat(lambda)[H]), "lambda_S" = expression(hat(lambda)[S]), "lambda_C" = expression(hat(lambda)[C]), "exp_H" = expression(hat(lambda)[W])) 
      ) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14))
    
}


row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(lambda[H], " = ", 2.4), paste(lambda[S], " = ", lambda[C], " = ", lambda[W], " = ", 0.3))), angle = 90) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(lambda[S], " = ", 2.4), paste(lambda[H], " = ", lambda[C], " = ", lambda[W], " = ", 0.3))), angle = 90) + theme_void() 
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(lambda[C], " = ", 2.4), paste(lambda[H], " = ", lambda[S], " = ", lambda[W], " = ", 0.3))), angle = 90) + theme_void() 
row4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(lambda[W], " = ", 2.4), paste(lambda[H], " = ", lambda[S], " = ", lambda[C], " = ", 0.3))), angle = 90) + theme_void() 

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Tree sizes", size = 5) + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="BLenD", size = 5) + theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Tree sizes + BLenD", size = 5) + theme_void() 

x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Tolerance rate", size = 5) + theme_void()
y_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Residuals", size = 5, angle = 90) + theme_void() 


layoutplot <- "
#qqqqrrrrssss#
uaaaaeeeeiiiim
uaaaaeeeeiiiim
uaaaaeeeeiiiim
uaaaaeeeeiiiim
ubbbbffffjjjjn
ubbbbffffjjjjn
ubbbbffffjjjjn
ubbbbffffjjjjn
uccccggggkkkko
uccccggggkkkko
uccccggggkkkko
uccccggggkkkko
uddddhhhhllllp
uddddhhhhllllp
uddddhhhhllllp
uddddhhhhllllp
#tttttttttttt#
"
plotlist <- list(
  a = panels[[1]],
  b = panels[[2]],
  c = panels[[3]],
  d = panels[[4]],
  e = panels[[5]],
  f = panels[[6]],
  g = panels[[7]],
  h = panels[[8]],
  i = panels[[9]],
  j = panels[[10]],
  k = panels[[11]],
  l = panels[[12]],
  m = row1,
  n = row2,
  o = row3,
  p = row4,
  q = col1,
  r = col2,
  s = col3,
  t = x_lab,
  u = y_lab
  )

print('plotting')
#png(paste(prefix, 'cophy_ABC_convergence/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
pdf(paste(prefix, 'cophy_ABC_convergence/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".pdf", sep = ''), width = 10, height = 10, pointsize = 12)
print(wrap_plots(plotlist, guides = 'collect', design = layoutplot) &
  theme(legend.position = "bottom"),
  legend.text = element_text(size = 14))
dev.off()

