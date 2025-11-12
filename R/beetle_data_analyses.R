# here we plot the inference results of the beetle mimicry data
setwd("/home/u29/yichaozeng/Desktop")
library(treeducken)
library(reshape2)
library(patchwork)
library(GGally)
library(ggridges)
library(ggpubr)
library(dplyr)

sim_time <- 2 # this is the time of simulations
n_bin <- 150 # this is the Number of summary statistics for each of the three normalized distributions

# this variable differentiate the real run from the test runs
real_data_run <- 1

para_est <- NULL
SS_est <- NULL
panels_OG <- NULL

for(ext_id in 1:5){
  
  folder_ids <-  list('27', '87', '78', '40', '90')[[ext_id]] # 0/0, 0.3/0.3, 0.7/0, 0/0.7, 0.7/0.7
  folder_ids_slash <-  list('27/', '87/', '78/', '40/', '90/')[[ext_id]]
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
  
  para_est[[length(para_est) + 1]] <- list(para_sim_acc)
  SS_est[[length(SS_est) + 1]] <- list(SS_sim_acc)

  # create the panels
  rate_dat <- data.frame(
    rate = c(para_sim_acc$lambda_H, para_sim_acc$lambda_S, para_sim_acc$lambda_C, para_sim_acc$exp_H), #/ 13.75,
    event = c(rep('lambda_H', nrow(para_sim_acc)), rep('lambda_S', nrow(para_sim_acc)), rep('lambda_C', nrow(para_sim_acc)), rep('lambda_W', nrow(para_sim_acc)) )
  )
  rate_dat$event <- factor(rate_dat$event, levels = c('lambda_H', 'lambda_S', 'lambda_C', 'lambda_W'))
  
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # the panel for speciation rate estimates
  # rate_plot <-   ggplot(rate_dat, aes(rate)) +
  #   geom_density(aes(color = event), linewidth = 0.5, alpha = 0.2, kernel = 'gaussian') +
  #   xlim(0, 0.3) +
  #   labs(title = "Speciation rate estimates", x = "events/lineage/Myr", y = "density") +
  #   scale_color_manual(
  #     values = cbp1,
  #     name = "Events:",
  #     labels = c(
  #       "lambda_H" = expression(paste("Host speciation (", hat(lambda)[H], ")")),
  #       "lambda_S" = expression(paste("Symbiont speciation without host switching (", hat(lambda)[S], ")")),
  #       "lambda_C" = expression(paste("Cospeciation (", hat(lambda)[C], ")")),
  #       "lambda_W" = expression(paste("Symbiont speciation with host switching (", hat(lambda)[W], ")"))
  #       ) 
  #     ) +
  #   guides(color = guide_legend(ncol = 1)) +  
  #   theme_bw() +
  #   theme(panel.grid = element_blank()) +
  #   theme(legend.position = "bottom") +
  #   theme(plot.title = element_text(hjust = 0.5))
  
  scale_factor <- 27.5 / sim_time
  rate_dat$rate <- rate_dat$rate / scale_factor
  
  stats <- rate_dat %>%
    group_by(event) %>%
    summarise(
      mean = mean(rate),
      sd   = sd(rate)
    )
  
  stats <- stats %>%
    mutate(label = sprintf("mean = %.3f\nsd = %.3f", mean, sd))
  
  rate_plot <- ggplot(rate_dat, aes(x = rate, y = event, height = after_stat(density))) +
    geom_density_ridges(aes(fill = event), linewidth = 0.5, scale = 0.9, alpha = 0.2, stat = 'density') +
    scale_y_discrete(
      labels = c(
        "lambda_H" = expression(paste(hat(lambda)[H])),
        "lambda_S" = expression(paste(hat(lambda)[S])),
        "lambda_C" = expression(paste(hat(lambda)[C])),
        "lambda_W" = expression(paste(hat(lambda)[W]))
        ),
      expand = expansion(add = c(0.5, 1))
      ) +
    xlim(-1 / scale_factor, 8.3 / scale_factor) +
    labs(title = "Speciation rate estimates", x = "events/lineage/Myr", y = NULL) +
    scale_fill_manual(
      values = cbp1,
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(size = 14, face = "bold")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  rate_plot <- rate_plot +
    geom_text(
      data = stats,
      aes(x = mean + sd + 0.05, y = event, label = label),
      inherit.aes = FALSE,
      color = "black",
      size = 3,
      hjust = 0,
      vjust = -0.5
    )
  
  # the BLenD a panel of posterior predictive checks
  sim_x <- NULL
  sim_y <- NULL
  for (row_id in 1:nrow(SS_sim_acc)) {
    sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, 1:n_bin]) / sum(unlist(SS_sim_acc[row_id, 1:n_bin])) ) # revert the BLenD SSs to its original scale for plotting
    sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  }
  point_data <- data.frame(x = sim_x / tr_ht, y = sim_y / (2 / n_bin))
  line_data <- data.frame(x = (breaks[-1] + breaks[-length(breaks)]) / (2 * tr_ht), y = (SS_real[SS_sel[1:n_bin]] / sum(SS_real[SS_sel[1:n_bin]])) / (2 / n_bin)) # revert the BLenD SSs to its original scale for plotting
  blend_a_plot <- ggplot() +
    geom_point(data = point_data, aes(x = x, y = y), size = 0.5, color = "red", alpha = 0.3) +  # Scatter plot
    geom_line(data = line_data, aes(x = x, y = y), linewidth = 0.5, color = "black") +      # Line plot
    labs(title = "BLenD-a", x = expression(d[a]), y = "density") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # the BLenD b panel of posterior predictive checks
  sim_x <- NULL
  sim_y <- NULL
  for (row_id in 1:nrow(SS_sim_acc)) {
    sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, (n_bin+1):(2*n_bin)]) / sum(unlist(SS_sim_acc[row_id, (n_bin+1):(2*n_bin)])) ) # revert the BLenD SSs to its original scale for plotting
    sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  }
  point_data <- data.frame(x = sim_x / tr_ht, y = sim_y / (2 / n_bin))
  line_data <- data.frame(x = (breaks[-1] + breaks[-length(breaks)]) / (2 * tr_ht), y = (SS_real[SS_sel[(n_bin+1):(2*n_bin)]] / sum(SS_real[SS_sel[(n_bin+1):(2*n_bin)]])) / (2 / n_bin)) # revert the BLenD SSs to its original scale for plotting
  blend_b_plot <- ggplot() +
    geom_point(data = point_data, aes(x = x, y = y), size = 0.5, color = "red", alpha = 0.3) +  # Scatter plot
    geom_line(data = line_data, aes(x = x, y = y), linewidth = 0.5, color = "black") +      # Line plot
    labs(title = "BLenD-b", x = expression(delta[b]), y = "density") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))
    
  # the host tree size panels of posterior predictive checks
  sim_y <- NULL
  for (row_id in 1:length(ids_sim_acc)) {
    sim_y <- c(sim_y, unlist(SS_comb[ids_sim_acc[row_id], (n_bin * 2 + 1)]))
  }
  point_data <- data.frame(y = sim_y)
  line_data <- SS_real_raw[[1]][SS_sel[(n_bin * 2 + 1)]]
  tree_size_h_plot <- ggplot(point_data, aes(x = "", y = y)) +  # Empty string for x-axis
    geom_jitter(width = 0.2, color = "red", size = 1, alpha = 0.3) +
    geom_hline(yintercept = line_data, color = "black", linewidth = 0.5) +
    ylim(10, 100) +
    labs(title = "Hosts", x = '', y = "Number") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )

  # the symbiont tree size panels of posterior predictive checks
  sim_y <- NULL
  for (row_id in 1:length(ids_sim_acc)) {
    sim_y <- c(sim_y, unlist(SS_comb[ids_sim_acc[row_id], (n_bin * 2 + 2)]))
  }
  point_data <- data.frame(y = sim_y)
  line_data <- SS_real_raw[[1]][SS_sel[(n_bin * 2 + 2)]]
  tree_size_s_plot <- ggplot(point_data, aes(x = "", y = y)) +  # Empty string for x-axis
    geom_jitter(width = 0.2, color = "red", size = 1, alpha = 0.3) +
    geom_hline(yintercept = line_data, color = "black", linewidth = 0.5) +
    ylim(10, 100) +
    labs(title = "Symbionts", x = '', y = "Number") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  panels_OG[[length(panels_OG) + 1]] <- list(rate_plot, blend_a_plot, blend_b_plot, tree_size_h_plot, tree_size_s_plot)

}

panels <- panels_OG

# save the list of plots, list of parameter estimates, and the list of percentages
file_name <- paste(prefix, 'cophy_ABC_results/real_panels.rds', sep = '')
list.save(panels, file_name)

file_name <- paste(prefix, 'cophy_ABC_results/real_para_est.rds', sep = '')
list.save(para_est, file_name)

row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(paste(epsilon[H], " = ", 0, ",  ", epsilon[S], " = ", 0)), size = 5, angle = 90) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(paste(epsilon[H], " = ", 0.3, ",  ", epsilon[S], " = ", 0.3)), size = 5, angle = 90) + theme_void() 
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(paste(epsilon[H], " = ", 0.7, ",  ", epsilon[S], " = ", 0)), size = 5, angle = 90) + theme_void() 
row4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(paste(epsilon[H], " = ", 0, ",  ", epsilon[S], " = ", 0.7)), size = 5, angle = 90) + theme_void() 
row5 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(paste(epsilon[H], " = ", 0.7, ",  ", epsilon[S], " = ", 0.7)), size = 5, angle = 90) + theme_void() 


layoutplot <- "
uaaaaaabbbbbcccccdd
uaaaaaabbbbbcccccdd
uaaaaaabbbbbcccccDD
uaaaaaabbbbbcccccDD
veeeeeefffffggggghh
veeeeeefffffggggghh
veeeeeefffffgggggHH
veeeeeefffffgggggHH
wiiiiiijjjjjkkkkkll
wiiiiiijjjjjkkkkkll
wiiiiiijjjjjkkkkkLL
wiiiiiijjjjjkkkkkLL
xmmmmmmnnnnnooooopp
xmmmmmmnnnnnooooopp
xmmmmmmnnnnnoooooPP
xmmmmmmnnnnnoooooPP
yqqqqqqrrrrrssssstt
yqqqqqqrrrrrssssstt
yqqqqqqrrrrrsssssTT
yqqqqqqrrrrrsssssTT
"
plotlist <- list(
  a = panels[[1]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  b = panels[[1]][[2]],
  c = panels[[1]][[3]],
  d = panels[[1]][[4]],
  D = panels[[1]][[5]],
  e = panels[[2]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  f = panels[[2]][[2]],
  g = panels[[2]][[3]],
  h = panels[[2]][[4]],
  H = panels[[2]][[5]],
  i = panels[[3]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  j = panels[[3]][[2]],
  k = panels[[3]][[3]],
  l = panels[[3]][[4]],
  L = panels[[3]][[5]],
  m = panels[[4]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  n = panels[[4]][[2]],
  o = panels[[4]][[3]],
  p = panels[[4]][[4]],
  P = panels[[4]][[5]],
  q = panels[[5]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  r = panels[[5]][[2]],
  s = panels[[5]][[3]],
  t = panels[[5]][[4]],
  'T' = panels[[5]][[5]],
  u = row1,
  v = row2,
  w = row3,
  x = row4,
  y = row5
)

print('plotting')
setwd("/groups/cromanpa/yzeng")
#png(paste(prefix, 'cophy_ABC_results/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
pdf(paste(prefix, 'cophy_ABC_results/', "beetle_results.pdf", sep = ''), width = 15, height = 15, pointsize = 12)
print(wrap_plots(plotlist, guides = 'collect', design = layoutplot) &
        theme(legend.position = 'none'))
dev.off()




# here is a modified version with cosmetic changes for Fig 4
layoutplot <- "
#aaaaaaaaaaaaa#bbbbbbbb#cccccccc#XXXX
meeeeeeeeeeeeenffffffffoggggggggwzzzz
meeeeeeeeeeeeenffffffffoggggggggwzzzz
meeeeeeeeeeeeenffffffffoggggggggwzzzz
meeeeeeeeeeeeenffffffffoggggggggwzzzz
meeeeeeeeeeeeenffffffffoggggggggwzzzz
meeeeeeeeeeeeenffffffffogggggggg#YYYY
meeeeeeeeeeeeenffffffffoggggggggWZZZZ
meeeeeeeeeeeeenffffffffoggggggggWZZZZ
meeeeeeeeeeeeenffffffffoggggggggWZZZZ
meeeeeeeeeeeeenffffffffoggggggggWZZZZ
meeeeeeeeeeeeenffffffffoggggggggWZZZZ
#iiiiiiiiiiiii#jjjjjjjj#kkkkkkkk#####
"

# layoutplot <- "
# #aaaaaaaaaaaaa#bbbbbbbb#cccccccc#XX
# meeeeeeeeeeeeenffffffffoggggggggWzz
# meeeeeeeeeeeeenffffffffoggggggggWzz
# meeeeeeeeeeeeenffffffffoggggggggWzz
# meeeeeeeeeeeeenffffffffoggggggggWzz
# meeeeeeeeeeeeenffffffffoggggggggWzz
# meeeeeeeeeeeeenffffffffogggggggg#YY
# meeeeeeeeeeeeenffffffffoggggggggWZZ
# meeeeeeeeeeeeenffffffffoggggggggWZZ
# meeeeeeeeeeeeenffffffffoggggggggWZZ
# meeeeeeeeeeeeenffffffffoggggggggWZZ
# #iiiiiiiiiiiii#jjjjjjjj#kkkkkkkkWZZ
# "

a <- ggplot() + annotate(geom = 'text', x=1, y=1, label="b) Speciation rate estimates", size = 5.2) + theme_void()
b <- ggplot() + annotate(geom = 'text', x=1, y=1, label="c) BLenD-a", size = 5.2) + theme_void()
c <- ggplot() + annotate(geom = 'text', x=1, y=1, label="d) BLenD-b", size = 5.2) + theme_void()

i <- ggplot() + annotate(geom = 'text', x=1, y=1, label="events/lineage/Myr", size = 5) + theme_void()
j <- ggplot() + annotate(geom = 'text', x=1, y=1, label=expression(paste(d[a])), size = 5) + theme_void()

k <- ggplot() + annotate(geom = 'text', x=1, y=1, label=expression(paste(delta[b])), size = 5) + theme_void()

# m <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void()
n <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void() 
o <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void() 

e <- panels[[2]][[1]] + labs(title = NULL, x = NULL, y = NULL) + xlim(-1 / scale_factor, 7.5 / scale_factor) #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian') +
  # scale_color_manual(
  #   values = cbp1,
  #   name = NULL,
  #   labels = c(
  #     "lambda_H" = expression(paste(hat(lambda)[H])),
  #     "lambda_S" = expression(paste(hat(lambda)[S])),
  #     "lambda_C" = expression(paste(hat(lambda)[C])),
  #     "lambda_W" = expression(paste(hat(lambda)[W]))
  #   ) 
  # ) +
  # theme(legend.position.inside = c(1,1))
f <- panels[[2]][[2]] + labs(title = NULL, x = NULL, y = NULL)
g <- panels[[2]][[3]] + labs(title = NULL, x = NULL, y = NULL)
# h <- panels[[2]][[4]] + labs(title = NULL, x = NULL, y = NULL)

X <- ggplot() + annotate(geom = 'text', x=1, y=1, label="e) Hosts", size = 5.2) + theme_void()
Y <- ggplot() + annotate(geom = 'text', x=1, y=1, label="f) Symbionts", size = 5.2) + theme_void()
z <- panels[[2]][[4]] + labs(title = NULL, x = NULL, y = NULL)
Z <- panels[[2]][[5]] + labs(title = NULL, x = NULL, y = NULL)
w <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Number", size = 5, angle = 90) + theme_void() 
W <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Number", size = 5, angle = 90) + theme_void() 

plotlist <- list(
  a = a,
  b = b,
  c = c,
  e = e + theme(
    legend.position = c(0.85,0.6),
    legend.spacing.y = unit(0.3, "cm"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  ),
  f = f + theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  ), #+
    #ylim(0,0.23), #+
    #annotate("text", x = -0, y = 0.23, label = "BLenD", size = 4),
  g = g + theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  ), #+
    #annotate("text", x = 1.5, y = 3.8, label = "Tree sizes", size = 4),
  # h = h + theme(
  #   legend.position = "bottom",
  #   panel.grid = element_blank(),
  #   axis.title = element_text(size = 14),
  #   plot.title = element_text(size = 16)
  # ),
  i = i,
  j = j,
  k = k,
  # l = l,
  # m = m,
  n = n,
  o = o,
  X = X,
  Y = Y,
  z = z,
  Z = Z,
  w = w,
  W = W
)

setwd("/groups/cromanpa/yzeng")
pdf(paste(prefix, 'cophy_ABC_results/', "beetle_results_Fig4.pdf", sep = ''), width = 12, height = 4, pointsize = 12)
print(wrap_plots(plotlist, guides = 'keep', design = layoutplot) &
        theme(
          legend.position = 'none' #c(0.85,0.7),
          # legend.text = element_text(size = 12),
          # legend.spacing.y = unit(0.3, "cm")
          # panel.grid = element_blank(),
          # axis.title = element_text(size = 12 * 1.5),
          # plot.title = element_text(size = 14 * 1.5)
          ) &
        guides(color = guide_legend(override.aes = list(linewidth = 0.7))))
dev.off()
