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

para_est[[length(para_est) + 1]] <- list(para_sim_acc)
SS_est[[length(SS_est) + 1]] <- list(SS_sim_acc)

# create the panels
rate_spec <- data.frame(
  rate = c(para_sim_acc$lambda_H, para_sim_acc$lambda_S, para_sim_acc$lambda_C, para_sim_acc$exp_H),
  event = c(rep('lambda_H', nrow(para_sim_acc)), rep('lambda_S', nrow(para_sim_acc)), rep('lambda_C', nrow(para_sim_acc)), rep('lambda_W', nrow(para_sim_acc)) )
)
rate_spec$event <- factor(rate_spec$event, levels = c('lambda_H', 'lambda_S', 'lambda_C', 'lambda_W'))

rate_ext <- data.frame(
  rate = c(para_sim_acc$mu_H_frac, para_sim_acc$mu_S_frac), #/ 13.75,
  event = c(rep('mu_H_frac', nrow(para_sim_acc)), rep('mu_S_frac', nrow(para_sim_acc)) )
)
rate_ext$event <- factor(rate_ext$event, levels = c('mu_H_frac', 'mu_S_frac'))

spec_cl <- c(
  "#4E79A7",
  "#F28E2B",
  "#E15759",
  "#76B7B2"
)

ext_cl <- c(
  "#B07AA1",
  "#59A14F"
)

# here we convert the speciation rates back to events/Myr
rate_spec$rate <- rate_spec$rate * sim_time / 27.5 # essentially this gives you the number of each type of events per lineage

# rate_spec$rate <- rate_spec$rate * sim_time # essentially this gives you the number of each type of events per lineage

spec_stats <- rate_spec %>%
  group_by(event) %>%
  summarise(
    mean = mean(rate),
    sd   = sd(rate)
  )
spec_stats <- spec_stats %>%
  mutate(label = sprintf("%.2f ± %.2f", mean, sd))

ext_stats <- rate_ext %>%
  group_by(event) %>%
  summarise(
    mean = mean(rate),
    sd   = sd(rate)
  )
ext_stats <- ext_stats %>%
  mutate(label = sprintf("%.2f ± %.2f", mean, sd))


spec_plot <- ggplot(rate_spec, aes(x = rate, y = event, height = after_stat(density))) +
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
  # xlim(, ) +
  labs(title = "Speciation rate estimates", x = "events/lineage/time", y = NULL) +
  scale_fill_manual(
    values = spec_cl,
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5))

spec_plot <- spec_plot +
  labs(title = NULL, x =  "events/lineage/Myr", y = NULL) +
  xlim(-1* sim_time / 27.5, 6* sim_time / 27.5) +
  geom_text(
    data = spec_stats,
    aes(x = mean + sd, y = event, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 3,
    hjust = 0,
    vjust = -0.5
  )

ext_plot <- ggplot(rate_ext, aes(x = rate, y = event, height = after_stat(density))) +
  geom_density_ridges(aes(fill = event), linewidth = 0.5, scale = 0.45, alpha = 0.2, stat = 'density') +
  scale_y_discrete(
    labels = c(
      "mu_H_frac" = expression(paste(hat(epsilon)[H])),
      "mu_S_frac" = expression(paste(hat(epsilon)[S]))
    ),
    expand = expansion(add = c(0.5, 1))
  ) +
  xlim(-0.2, 1.2) +
  labs(title = "extinction rate estimates", x = "", y = NULL) +
  scale_fill_manual(
    values = ext_cl,
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5))

ext_plot <- ext_plot +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_x_continuous(
    limits = c(-0.2, 1.2),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
  ) +
  geom_text(
    data = ext_stats,
    aes(x = mean + sd * -0.5 + 0, y = event, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 3,
    hjust = 0,
    vjust = -0.5
  )

# the BLenD a panel of posterior predictive checks
sim_x <- NULL
sim_y <- NULL
group <- NULL
breaks <- seq(from = -1, to = 1, length.out = n_bin + 1) 
for (row_id in 1:nrow(SS_sim_acc)) {
  sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, 1:n_bin]) / sum(unlist(SS_sim_acc[row_id, 1:n_bin])) ) # revert the BLenD SSs to its original scale for plotting
  sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  group <- c(group, rep(as.character(row_id), length(breaks)-1))
}
point_data <- data.frame(x = sim_x, y = sim_y / (2 / n_bin), group)
line_data <- data.frame(x = (breaks[-1] + breaks[-length(breaks)]) / 2, y = (SS_real[SS_sel[1:n_bin]] / sum(SS_real[SS_sel[1:n_bin]])) / (2 / n_bin)) # revert the BLenD SSs to its original scale for plotting
blend_a_plot <- ggplot() +
  geom_line(data = point_data, aes(x = x, y = y, group = group), size = 0.5, color = "black", alpha = 0.3) +  # Scatter plot
  geom_line(data = line_data, aes(x = x, y = y), linewidth = 1, color = "red") +      # Line plot
  xlim(-1, 1) +
  ylim(0, 7) +
  labs(title = NULL, x = NULL, y = NULL) +
  # labs(title = "BLenD-a", x = expression(d[a]), y = "density") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

# the BLenD b panel of posterior predictive checks
sim_x <- NULL
sim_y <- NULL
group <- NULL
for (row_id in 1:nrow(SS_sim_acc)) {
  sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, (n_bin+1):(2*n_bin)]) / sum(unlist(SS_sim_acc[row_id, (n_bin+1):(2*n_bin)])) ) # revert the BLenD SSs to its original scale for plotting
  sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  group <- c(group, rep(as.character(row_id), length(breaks)-1))
}
point_data <- data.frame(x = sim_x, y = sim_y / (2 / n_bin))
line_data <- data.frame(x = (breaks[-1] + breaks[-length(breaks)]) / 2, y = (SS_real[SS_sel[(n_bin+1):(2*n_bin)]] / sum(SS_real[SS_sel[(n_bin+1):(2*n_bin)]])) / (2 / n_bin)) # revert the BLenD SSs to its original scale for plotting
blend_b_plot <- ggplot() +
  geom_line(data = point_data, aes(x = x, y = y, group = group), size = 0.5, color = "black", alpha = 0.3) +  # Scatter plot
  geom_line(data = line_data, aes(x = x, y = y), linewidth = 1, color = "red") +      # Line plot
  xlim(-1, 1) +
  ylim(0, 7) +
  labs(title = NULL, x = NULL, y = NULL) +
  # labs(title = "BLenD-b", x = expression(delta[b]), y = "density") +
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
  geom_jitter(width = 0.2, color = "black", size = 2, alpha = 0.3) +
  geom_hline(yintercept = line_data, color = "red", linewidth = 1) +
  ylim(10, 100) +
  labs(title = NULL, x = NULL, y = NULL) +
  # labs(title = "Hosts", x = '', y = "Number") +
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
  geom_jitter(width = 0.2, color = "black", size = 2, alpha = 0.3) +
  geom_hline(yintercept = line_data, color = "red", linewidth = 1) +
  ylim(10, 100) +
  labs(title = NULL, x = NULL, y = NULL) +
  # labs(title = "Symbionts", x = '', y = "Number") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

# panels_OG[[length(panels_OG) + 1]] <- list(spec_plot, ext_plot, blend_a_plot, blend_b_plot, tree_size_h_plot, tree_size_s_plot)
panels_OG[[0 + 1]] <- list(spec_plot, ext_plot, blend_a_plot, blend_b_plot, tree_size_h_plot, tree_size_s_plot)

panels <- panels_OG

# save the list of plots, list of parameter estimates, and the list of percentages
file_name <- paste(prefix, 'cophy_ABC_results/real_panels.rds', sep = '')
list.save(panels, file_name)

file_name <- paste(prefix, 'cophy_ABC_results/real_para_est.rds', sep = '')
list.save(para_est, file_name)

# For Fig 4
# layoutplot <- "
# aaaaaaaaaa#ccccccccccc#dddddddddddXXXYYY
# eeeeeeeeeengggggggggggohhhhhhhhhhhwzzWZZ
# eeeeeeeeeengggggggggggohhhhhhhhhhhwzzWZZ
# eeeeeeeeeengggggggggggohhhhhhhhhhhwzzWZZ
# eeeeeeeeeengggggggggggohhhhhhhhhhhwzzWZZ
# eeeeeeeeeengggggggggggohhhhhhhhhhhwzzWZZ
# eeeeeeeeeengggggggggggohhhhhhhhhhhwzzWZZ
# eeeeeeeeeengggggggggggohhhhhhhhhhhwzzWZZ
# bbbbbbbbbbngggggggggggohhhhhhhhhhhwzzWZZ
# ffffffffffngggggggggggohhhhhhhhhhhwzzWZZ
# ffffffffffngggggggggggohhhhhhhhhhhwzzWZZ
# ffffffffffngggggggggggohhhhhhhhhhhwzzWZZ
# ###########jjjjjjjjjjj#kkkkkkkkkkk######
# "

layoutplot <- "
#####################ccccc
####################nggggg
####################nggggg
####################nggggg
####################nggggg
####################nggggg
#####################jjjjj
#####################ddddd
####################ohhhhh
####################ohhhhh
####################ohhhhh
####################ohhhhh
####################ohhhhh
#####################kkkkk
aaaaaaaaaabbbbbbbbbbXXXYYY
eeeeeeeeeeffffffffffwzzWZZ
eeeeeeeeeeffffffffffwzzWZZ
eeeeeeeeeeffffffffffwzzWZZ
eeeeeeeeeeffffffffffwzzWZZ
eeeeeeeeeeffffffffffwzzWZZ
"

a <- ggplot() + annotate(geom = 'text', x=1, y=1, label="b) Speciation rate", size = 5.2) + theme_void()
b <- ggplot() + annotate(geom = 'text', x=1, y=1, label="c) Relative extinction rate", size = 5.2) + theme_void()
c <- ggplot() + annotate(geom = 'text', x=1, y=1, label="d) BLenD-a", size = 5.2) + theme_void()
d <- ggplot() + annotate(geom = 'text', x=1, y=1, label="e) BLenD-b", size = 5.2) + theme_void()

e <- panels[[1]][[1]]
f <- panels[[1]][[2]]
g <- panels[[1]][[3]]
h <- panels[[1]][[4]]

# i <- ggplot() + annotate(geom = 'text', x=1, y=1, label="events/lineage/Myr", size = 5) + theme_void()
j <- ggplot() + annotate(geom = 'text', x=1, y=1, label=expression(paste(d[a])), size = 5) + theme_void()
k <- ggplot() + annotate(geom = 'text', x=1, y=1, label=expression(paste(delta[b])), size = 5) + theme_void()

n <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void() 
o <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Density", size = 5, angle = 90) + theme_void() 

X <- ggplot() + annotate(geom = 'text', x=1, y=1, label="f) Hosts", size = 5.2) + theme_void()
Y <- ggplot() + annotate(geom = 'text', x=1, y=1, label="g) Symbionts", size = 5.2) + theme_void()
z <- panels[[1]][[5]]
Z <- panels[[1]][[6]]
w <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Number", size = 5, angle = 90) + theme_void() 
W <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Number", size = 5, angle = 90) + theme_void() 

plotlist <- list(
  a = a,
  b = b,
  c = c,
  d = d,
  e = e + theme(
    legend.position = c(0.85,0.6),
    legend.spacing.y = unit(0.3, "cm"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  ),
  f = f + theme(
    legend.position = c(0.85,0.6),
    legend.spacing.y = unit(0.3, "cm"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  ),
  g = g + theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  ), #+
    #ylim(0,0.23), #+
    #annotate("text", x = -0, y = 0.23, label = "BLenD", size = 4),
  h = h + theme(
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
  # i = i,
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
pdf(paste(prefix, 'cophy_ABC_results/', "beetle_results_Fig4_20260415.pdf", sep = ''), width = 10, height = 10, pointsize = 12)
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

