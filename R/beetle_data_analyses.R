# here we plot the inference results of the beetle mimicry data
setwd("/home/u29/yichaozeng/Desktop")
library(treeducken)
library(reshape2)
library(patchwork)
library(GGally)
library(ggridges)
library(ggpubr)

sim_time <- 2 # this is the time of simulations
n_bin <- 150 # this is the number of summary statistics for each of the three normalized distributions

# this variable differentiate the real run from the test runs
real_data_run <- 1

para_est <- NULL
SS_est <- NULL
panels_OG <- NULL

for(ext_id in 1:5){
  
  folder_ids <-  list('86', '69', '17', '40', '76')[[ext_id]] # 0/0, 0.3/0.3, 0.7/0, 0/0.7, 0.7/0.7
  folder_ids_slash <-  list('86/', '69/', '17/', '40/', '76/')[[ext_id]]
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
    
  SS_sel <- list(
    c(451:452), # Tree sizes only
    c(1:150), # BLenD only
    c(1:150, 451:452) # combined
  )[[3]]
  
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
  #   labs(title = "Speciation rate estimates", x = "rates", y = "density") +
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
    xlim(0, 4) +
    labs(title = "Speciation rate estimates", x = "rates", y = NULL) +
    scale_fill_manual(
      values = cbp1,
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(size = 14, face = "bold")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # the BLenD panel of posterior predictive checks
  sim_x <- NULL
  sim_y <- NULL
  for (row_id in 1:nrow(SS_sim_acc)) {
    sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, 1:n_bin]) / sum(unlist(SS_sim_acc[row_id, 1:n_bin])) ) # revert the BLenD SSs to its original scale for plotting
    sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  }
  point_data <- data.frame(x = sim_x / tr_ht, y = sim_y / (2 / n_bin))
  line_data <- data.frame(x = (breaks[-1] + breaks[-length(breaks)]) / (2 * tr_ht), y = (SS_real[SS_sel[1:n_bin]] / sum(SS_real[SS_sel[1:n_bin]])) / (2 / n_bin)) # revert the BLenD SSs to its original scale for plotting
  blend_plot <- ggplot() +
    geom_point(data = point_data, aes(x = x, y = y), size = 0.1, color = "black", alpha = 0.3) +  # Scatter plot
    geom_line(data = line_data, aes(x = x, y = y), linewidth = 1, color = "red") +      # Line plot
    labs(title = "Model fit (BLenD)", x = expression(delta), y = "density") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # # the host tree size panels of posterior predictive checks
  # sim_y <- NULL
  # for (row_id in 1:nrow(SS_sim_acc)) {
  #   sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, (n_bin + 1)]))
  # }
  # point_data <- data.frame(y = sim_y)
  # line_data <- SS_real[SS_sel[(n_bin + 1)]]
  # tree_size_h_plot <- ggplot(point_data, aes(x = "", y = y)) +  # Empty string for x-axis
  #   geom_jitter(width = 0.2, size = 2, alpha = 0.3) +
  #   geom_hline(yintercept = line_data, color = "red", linewidth = 1) +
  #   ylim(2, 4) +
  #   labs(title = "Tree size", x = 'host', y = "log(N. of tips)") +
  #   theme_bw() +
  #   theme(panel.grid = element_blank()) +
  #   theme(plot.title = element_text(hjust = 0.5))
  # 
  # # the symbiont tree size panels of posterior predictive checks
  # sim_y <- NULL
  # for (row_id in 1:nrow(SS_sim_acc)) {
  #   sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, (n_bin + 2)]))
  # }
  # point_data <- data.frame(y = sim_y)
  # line_data <- SS_real[SS_sel[(n_bin + 2)]]
  # tree_size_s_plot <- ggplot(point_data, aes(x = "", y = y)) +  # Empty string for x-axis
  #   geom_jitter(width = 0.2, size = 2, alpha = 0.3) +
  #   geom_hline(yintercept = line_data, color = "red", linewidth = 1) +
  #   ylim(2, 4) +
  #   labs(title = "Tree size", x = 'symbiont', y = "log(N. of tips)") +
  #   theme_bw() +
  #   theme(panel.grid = element_blank()) +
  #   theme(plot.title = element_text(hjust = 0.5))
  
  # the combined panel showing the host and symbiont tree sizes
  sim_y_h <- NULL
  for (row_id in 1:nrow(SS_sim_acc)) {
    sim_y_h <- c(sim_y_h, unlist(SS_sim_acc[row_id, (n_bin + 1)]))
  }
  
  sim_y_s <- NULL
  for (row_id in 1:nrow(SS_sim_acc)) {
    sim_y_s <- c(sim_y_s, unlist(SS_sim_acc[row_id, (n_bin + 2)]))
  }
  
  size_dat <- data.frame(
    party = c(rep('host', length(sim_y_h)), rep('symbiont', length(sim_y_s))),
    size = c(sim_y_h, sim_y_s)
  )
  
  size_dat_summary <- aggregate(. ~ party, mean, data=size_dat)
  size_dat_summary$size <- c(SS_real[SS_sel[(n_bin + 1)]], SS_real[SS_sel[(n_bin + 2)]])
  
  tree_size_plot <- ggplot(size_dat, aes(x=party, y=size)) +  
    geom_jitter(width = 0.2, size = 2, alpha = 0.3) +
    geom_crossbar(data=size_dat_summary, aes(ymin = size, ymax = size),
                  linewidth=0.6,col="red", width = .8) +
    ylim(2, 4) +
    labs(title = "Model fit (tree sizes)", x = 'tree', y = "log(N. of tips)") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))
  
    
    panels_OG[[length(panels_OG) + 1]] <- list(rate_plot, blend_plot, tree_size_plot)

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
uaaaaaabbbbbbccc
uaaaaaabbbbbbccc
uaaaaaabbbbbbccc
uaaaaaabbbbbbccc
veeeeeeffffffggg
veeeeeeffffffggg
veeeeeeffffffggg
veeeeeeffffffggg
wiiiiiijjjjjjkkk
wiiiiiijjjjjjkkk
wiiiiiijjjjjjkkk
wiiiiiijjjjjjkkk
xmmmmmmnnnnnnooo
xmmmmmmnnnnnnooo
xmmmmmmnnnnnnooo
xmmmmmmnnnnnnooo
yqqqqqqrrrrrrsss
yqqqqqqrrrrrrsss
yqqqqqqrrrrrrsss
yqqqqqqrrrrrrsss
"
plotlist <- list(
  a = panels[[1]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  b = panels[[1]][[2]],
  c = panels[[1]][[3]],
  # d = panels[[1]][[4]],
  e = panels[[2]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  f = panels[[2]][[2]],
  g = panels[[2]][[3]],
  # h = panels[[2]][[4]],
  i = panels[[3]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  j = panels[[3]][[2]],
  k = panels[[3]][[3]],
  # l = panels[[3]][[4]],
  m = panels[[4]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  n = panels[[4]][[2]],
  o = panels[[4]][[3]],
  # p = panels[[4]][[4]],
  q = panels[[5]][[1]], #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian'),
  r = panels[[5]][[2]],
  s = panels[[5]][[3]],
  # t = panels[[5]][[4]],
  u = row1,
  v = row2,
  w = row3,
  x = row4,
  y = row5
)

print('plotting')
#png(paste(prefix, 'cophy_ABC_results/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
pdf(paste(prefix, 'cophy_ABC_results/', "beetle_results.pdf", sep = ''), width = 11, height = 15, pointsize = 12)
print(wrap_plots(plotlist, guides = 'collect', design = layoutplot) &
        theme(legend.position = 'none'))
dev.off()




# here is a modified version with cosmetic changes for Fig 4
layoutplot <- "
#aaaaaaaaaaaaa#bbbbbbbbbb#ccccc
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
meeeeeeeeeeeeenffffffffffoggggg
#iiiiiiiiiiiii#jjjjjjjjjj#kkkkk
"

a <- ggplot() + annotate(geom = 'text', x=1, y=1, label="b) Speciation rate estimates", size = 5.2) + theme_void()
b <- ggplot() + annotate(geom = 'text', x=1, y=1, label="c) Model fit", size = 5.2) + theme_void()
c <- ggplot() + annotate(geom = 'text', x=1, y=1, label="d) Model fit", size = 5.2) + theme_void()
# b <- ggplot() + annotate(geom = 'text', x=1, y=1, label="c) Model fit (BLenD & Tree sizes)", size = 5.2) + theme_void()

i <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Rates", size = 5) + theme_void()
j <- ggplot() + annotate(geom = 'text', x=1, y=1, label=expression(paste(delta, " (BLenD)")), size = 5) + theme_void()
# k <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Host", size = 5) + theme_void()
# l <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Symbiont", size = 5) + theme_void()
k <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Tree size", size = 5) + theme_void()

# m <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void()
n <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void() 
o <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "log(N. of tips)", size = 5, angle = 90) + theme_void() 

e <- panels[[2]][[1]] + labs(title = NULL, x = NULL, y = NULL) #+ geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian') +
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
  o = o
)

pdf(paste(prefix, 'cophy_ABC_results/', "beetle_results_Fig4.pdf", sep = ''), width = 10, height = 4, pointsize = 12)
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








# here is a modified version with cosmetic changes for Fig 4, optimized for compact spacing
layoutplot <- "
#aaaaaaaaaaa#bbbbbbbbbbb
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
meeeeeeeeeeenfffffffffff
#iiiiiiiiiii#jjjjjjjjjjj
"

a <- ggplot() + annotate(geom = 'text', x=1, y=1, label="b) Speciation rate estimates", size = 5.2) + theme_void()
b <- ggplot() + annotate(geom = 'text', x=1, y=1, label="c) Model fit", size = 5.2) + theme_void()

i <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Rates", size = 5) + theme_void()
j <- ggplot() + annotate(geom = 'text', x=1, y=1, label=expression(delta), size = 5) + theme_void()

m <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void()
n <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void() 

e <- panels[[2]][[1]] + labs(title = 'b) Speciation rate estimates', x = 'Rates', y = 'Density') + geom_density(aes(color = event), linewidth = 1, alpha = 0.2, kernel = 'gaussian') +
  scale_color_manual(
    values = cbp1,
    name = NULL,
    labels = c(
      "lambda_H" = expression(paste(hat(lambda)[H])),
      "lambda_S" = expression(paste(hat(lambda)[S])),
      "lambda_C" = expression(paste(hat(lambda)[C])),
      "lambda_W" = expression(paste(hat(lambda)[W]))
    ) 
  ) +
  theme(legend.position.inside = c(1,1))

f <- panels[[2]][[2]] +
  labs(title = 'c) Model fit', x = expression(delta), y = 'Density') +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  ) +
  annotate("text", x = -0.2, y = 0.2, label = "BLenD", size = 5)

g <- panels[[2]][[3]] + labs(title = NULL, x = 'Tree', y = 'log(N. of tips)') +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12)
  ) +
  annotate("text", x = 1.5, y = 3.9, label = "Tree sizes", size = 5)

plotlist <- list(
  a = a,
  b = b,
  c = c,
  e = e + theme(
    legend.position = c(0.85,0.6),
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  ),
  f = f +
    annotation_custom(ggplotGrob(g), xmin = 0.2, xmax = 1, 
                                         ymin = 0.05, ymax = 0.21)
  ,
  i = i,
  j = j,
  m = m,
  n = n
)

pdf(paste(prefix, 'cophy_ABC_results/', "beetle_results_Fig4_compact.pdf", sep = ''), width = 12, height = 5, pointsize = 12)
# print(wrap_plots(plotlist, guides = 'keep', design = layoutplot) &
#         theme(
#           legend.position = "none",
#           legend.text = element_text(size = 12)
#           # panel.grid = element_blank(),
#           # axis.title = element_text(size = 12 * 1.5),
#           # plot.title = element_text(size = 14 * 1.5)
#         ))
print(ggarrange(plotlist$e, plotlist$f))
dev.off()

