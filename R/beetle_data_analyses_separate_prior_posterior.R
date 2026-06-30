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
library(viridis)

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

# read in the data frame coontaining the prior
file_name <- paste("ex_", 'cophy_ABC_convergence/cross_validation_sizes', '.csv', sep = '')
rel_err <- read.csv(file_name)

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


# create the panels
# averaging function
aver <- function(vec, every){ # averaging function for every X rows
  temp <- tapply(vec, (seq_along(vec)-1) %/% every, mean)
  return(temp)
}

# here we process the rates (for both the prior and posterior distributions)
# for the posterior
para_sim_acc$mu_H <- para_sim_acc$mu_H_frac * (para_sim_acc$lambda_H + para_sim_acc$lambda_C)
para_sim_acc$mu_S <- para_sim_acc$mu_S_frac * (para_sim_acc$lambda_S + para_sim_acc$lambda_C + para_sim_acc$exp_H)

r_host <- para_sim_acc$lambda_H + para_sim_acc$lambda_C - para_sim_acc$mu_H
r_symb <- para_sim_acc$lambda_S + para_sim_acc$lambda_C + para_sim_acc$exp_H - para_sim_acc$mu_S
para_sim_acc$r_cophy <- 2 * r_host * r_symb / (r_host + r_symb)

dat_joint_posterior_relative <- data.frame(
  lambda_H_frac = para_sim_acc$lambda_H / para_sim_acc$r_cophy,
  lambda_S_frac = para_sim_acc$lambda_S / para_sim_acc$r_cophy,
  lambda_C_frac = para_sim_acc$lambda_C / para_sim_acc$r_cophy,
  lambda_W_frac = para_sim_acc$exp_H / para_sim_acc$r_cophy,
  mu_H_frac = para_sim_acc$mu_H / para_sim_acc$r_cophy,
  mu_S_frac = para_sim_acc$mu_S / para_sim_acc$r_cophy,
  group = rep('posterior', nrow(para_sim_acc))
)

dat_joint_posterior_absolute <- data.frame(
  lambda_H = para_sim_acc$lambda_H,
  lambda_S = para_sim_acc$lambda_S,
  lambda_C = para_sim_acc$lambda_C,
  lambda_W = para_sim_acc$exp_H,
  mu_H = para_sim_acc$mu_H,
  mu_S = para_sim_acc$mu_S,
  group = rep('posterior', nrow(para_sim_acc))
)

# for the prior
para_sim_md$mu_H <- para_sim_md$mu_H_frac * (para_sim_md$lambda_H + para_sim_md$lambda_C)
para_sim_md$mu_S <- para_sim_md$mu_S_frac * (para_sim_md$lambda_S + para_sim_md$lambda_C + para_sim_md$exp_H)

para_sim_md$r_cophy <- para_sim_md$lambda_H + para_sim_md$lambda_C - para_sim_md$mu_H + para_sim_md$lambda_S + para_sim_md$lambda_C + para_sim_md$exp_H - para_sim_md$mu_S

dat_joint_prior_relative <- data.frame(
  lambda_H_frac = para_sim_md$lambda_H / para_sim_md$r_cophy,
  lambda_S_frac = para_sim_md$lambda_S / para_sim_md$r_cophy,
  lambda_C_frac = para_sim_md$lambda_C / para_sim_md$r_cophy,
  lambda_W_frac = para_sim_md$exp_H / para_sim_md$r_cophy,
  mu_H_frac = para_sim_md$mu_H / para_sim_md$r_cophy,
  mu_S_frac = para_sim_md$mu_S / para_sim_md$r_cophy,
  group = rep('prior', nrow(para_sim_acc))
)[1:1000 * 1200,] # subsetting to reduce plotting pressure

dat_joint_prior_absolute <- data.frame(
  lambda_H = para_sim_md$lambda_H,
  lambda_S = para_sim_md$lambda_S,
  lambda_C = para_sim_md$lambda_C,
  lambda_W = para_sim_md$exp_H,
  mu_H = para_sim_md$mu_H,
  mu_S = para_sim_md$mu_S,
  group = rep('prior', nrow(para_sim_acc))
)[1:1000 * 1200,] # subsetting to reduce plotting pressure

# now we plot the joint distribution
dat_relative <- rbind(dat_joint_prior_relative, dat_joint_posterior_relative)
dat_absolute <- rbind(dat_joint_prior_absolute, dat_joint_posterior_absolute)

dat_relative$group <- factor(dat_relative$group, levels = c('prior', 'posterior'))
dat_absolute$group <- factor(dat_absolute$group, levels = c('prior', 'posterior'))

pair_relative <- GGally::ggpairs(
  log10(dat_relative[,1:6]),
  columnLabels = c(
    'hat(lambda)[H] / hat(r)',
    'hat(lambda)[S] / hat(r)',
    'hat(lambda)[C] / hat(r)',
    'hat(lambda)[W] / hat(r)',
    'hat(mu)[H] / hat(r)',
    'hat(mu)[S] / hat(r)'
  ),
  labeller = "label_parsed",
  lower = list(continuous = "points"),
  upper = list(continuous = "blank"),
  diag = list(continuous = wrap("densityDiag", alpha = 0.6)),
  mapping = aes(color = dat_relative$group),
)+
  scale_color_manual(
    values = c(
      "black", "#B1CAE1"
    )
  )+
  scale_fill_manual(
    values = c(
      "black", "#B1CAE1"
    )
  )
# pair_relative

pair_absolute <- GGally::ggpairs(
  log10(dat_absolute[,1:6]),
  columnLabels = c(
    'hat(lambda)[H]',
    'hat(lambda)[S]',
    'hat(lambda)[C]',
    'hat(lambda)[W]',
    'hat(mu)[H]',
    'hat(mu)[S]'
  ),
  labeller = "label_parsed",
  lower = list(continuous = "points"),
  upper = list(continuous = "blank"),
  diag = list(continuous = wrap("densityDiag", alpha = 0.6)),
  mapping = aes(color = dat_absolute$group),
)+
  scale_color_manual(
    values = c(
      "black", "#F28E2B"
    )
  )+
  scale_fill_manual(
    values = c(
      "black", "#F28E2B"
    )
  )
# pair_absolute


setwd("/groups/cromanpa/yzeng")

pdf(paste(prefix, 'cophy_ABC_results/', "pair_relative.pdf", sep = ''), width = 10, height = 10, pointsize = 12)
print(pair_relative)
dev.off()

pdf(paste(prefix, 'cophy_ABC_results/', "pair_absolute.pdf", sep = ''), width = 10, height = 10, pointsize = 12)
print(pair_absolute)
dev.off()


# now we plot the marginal distribution

rate_relative <- data.frame(
  rate = c(dat_joint_posterior_relative$lambda_H_frac, dat_joint_posterior_relative$lambda_S_frac, dat_joint_posterior_relative$lambda_C_frac, dat_joint_posterior_relative$lambda_W_frac, dat_joint_posterior_relative$mu_H_frac, dat_joint_posterior_relative$mu_S_frac),
  event = c(rep('lambda_H_frac', nrow(dat_joint_posterior_relative)), rep('lambda_S_frac', nrow(dat_joint_posterior_relative)), rep('lambda_C_frac', nrow(dat_joint_posterior_relative)), rep('lambda_W_frac', nrow(dat_joint_posterior_relative)), rep('mu_H_frac', nrow(dat_joint_posterior_relative)), rep('mu_S_frac', nrow(dat_joint_posterior_relative)))
)
rate_relative$event <- factor(rate_relative$event, levels = c('lambda_H_frac', 'lambda_S_frac', 'lambda_C_frac', 'lambda_W_frac', 'mu_H_frac', 'mu_S_frac'))

rate_absolute <- data.frame(
  rate = c(dat_joint_posterior_absolute$lambda_H, dat_joint_posterior_absolute$lambda_S, dat_joint_posterior_absolute$lambda_C, dat_joint_posterior_absolute$lambda_W, dat_joint_posterior_absolute$mu_H, dat_joint_posterior_absolute$mu_S),
  event = c(rep('lambda_H', nrow(dat_joint_posterior_absolute)), rep('lambda_S', nrow(dat_joint_posterior_absolute)), rep('lambda_C', nrow(dat_joint_posterior_absolute)), rep('lambda_W', nrow(dat_joint_posterior_absolute)), rep('mu_H', nrow(dat_joint_posterior_absolute)), rep('mu_S', nrow(dat_joint_posterior_absolute)))
)
rate_absolute$event <- factor(rate_absolute$event, levels = c('lambda_H', 'lambda_S', 'lambda_C', 'lambda_W', 'mu_H', 'mu_S'))



# rate_spec <- data.frame(
#   rate = c(para_sim_acc$lambda_H, para_sim_acc$lambda_S, para_sim_acc$lambda_C, para_sim_acc$exp_H),
#   event = c(rep('lambda_H', nrow(para_sim_acc)), rep('lambda_S', nrow(para_sim_acc)), rep('lambda_C', nrow(para_sim_acc)), rep('lambda_W', nrow(para_sim_acc)) )
# )
# rate_spec$event <- factor(rate_spec$event, levels = c('lambda_H', 'lambda_S', 'lambda_C', 'lambda_W'))
# 
# rate_spec_ratio <- data.frame(
#   ratio = c(para_sim_acc$lambda_H/para_sim_acc$sum_lambda_host, para_sim_acc$lambda_S/para_sim_acc$sum_lambda_symb, para_sim_acc$lambda_C/para_sim_acc$sum_lambda_host, para_sim_acc$lambda_C/para_sim_acc$sum_lambda_symb, para_sim_acc$exp_H/para_sim_acc$sum_lambda_symb),
#   event = c(rep('lambda_H', nrow(para_sim_acc)), rep('lambda_S', nrow(para_sim_acc)), rep('lambda_C_host', nrow(para_sim_acc)), rep('lambda_C_symb', nrow(para_sim_acc)), rep('lambda_W', nrow(para_sim_acc)) )
# )
# rate_spec_ratio$event <- factor(rate_spec_ratio$event, levels = c('lambda_H', 'lambda_S', 'lambda_C_host', 'lambda_C_symb', 'lambda_W'))
# 
# rate_ext <- data.frame(
#   rate = c(para_sim_acc$mu_H_frac, para_sim_acc$mu_S_frac), #/ 13.75,
#   event = c(rep('mu_H_frac', nrow(para_sim_acc)), rep('mu_S_frac', nrow(para_sim_acc)) )
# )
# rate_ext$event <- factor(rate_ext$event, levels = c('mu_H_frac', 'mu_S_frac'))

# relative_cl <- c(
#   "#4E79A7",
#   "#F28E2B",
#   "#E15759",
#   "#76B7B2"
# )
# 
# relative_cl <- c(
#   "#B07AA1",
#   "#59A14F"
# )

# here we convert the speciation rates back to events/Myr
rate_absolute_arbitrary <- rate_absolute # rates in events/lineage/Myr
rate_absolute$rate <- rate_absolute$rate * sim_time / 27.5 # essentially this gives you the number of each type of events per lineage

# rate_spec$rate <- rate_spec$rate * sim_time # essentially this gives you the number of each type of events per lineage

relative_stats <- rate_relative %>%
  group_by(event) %>%
  summarise(
    mean = mean(rate),
    sd   = sd(rate)
  )
relative_stats <- relative_stats %>%
  mutate(label = sprintf("%.2f ± %.2f", mean, sd))

absolute_stats <- rate_absolute %>%
  group_by(event) %>%
  summarise(
    mean = mean(rate),
    sd   = sd(rate)
  )
absolute_stats <- absolute_stats %>%
  mutate(label = sprintf("%.2f ± %.2f", mean, sd))


relative_plot <- ggplot(rate_relative, aes(x = event, y = rate)) +
  
  geom_violin(
    aes(fill = event),
    alpha = 0.2,
    linewidth = 0.5,
    trim = FALSE
  ) +
  
  geom_boxplot(
    width = 0.12,
    outlier.size = 0.8,
    linewidth = 0.4,
    alpha = 0.8
  ) +
  
  scale_y_log10() +
  
  # coord_flip() +
  
  scale_x_discrete(
    labels = c(
      "lambda_H_frac" = expression(frac(hat(lambda)[H], hat(r))),
      "lambda_S_frac" = expression(frac(hat(lambda)[S], hat(r))),
      "lambda_C_frac" = expression(frac(hat(lambda)[C], hat(r))),
      "lambda_W_frac" = expression(frac(hat(lambda)[W], hat(r))),
      "mu_H_frac" = expression(frac(hat(mu)[H], hat(r))),
      "mu_S_frac" = expression(frac(hat(mu)[S], hat(r)))
    )
  ) +
  
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  
  # scale_fill_manual(values = spec_cl) +
  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  # theme(axis.text.x = element_text(
  #   angle = 45,
  #   hjust = 1
  # )) +
  theme(axis.text.x = element_text(size = 12, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5))





absolute_plot <- ggplot(rate_absolute, aes(x = event, y = rate)) +
  
  geom_violin(
    aes(fill = event),
    alpha = 0.2,
    linewidth = 0.5,
    trim = FALSE
  ) +
  
  geom_boxplot(
    width = 0.12,
    outlier.size = 0.8,
    linewidth = 0.4,
    alpha = 0.8
  ) +
  
  scale_y_log10() +
  
  # coord_flip() +
  
  scale_x_discrete(
    labels = c(
      "lambda_H" = expression(hat(lambda)[H]),
      "lambda_S" = expression(hat(lambda)[S]),
      "lambda_C" = expression(hat(lambda)[C]),
      "lambda_W" = expression(hat(lambda)[W]),
      "mu_H" = expression(hat(mu)[H]),
      "mu_S" = expression(hat(mu)[S])
    )
  ) +
  
  labs(
    title = NULL,
    x = NULL,
    y = "events / lineage / Myr"
  ) +
  
  # scale_fill_manual(values = spec_cl) +
  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  # theme(axis.text.x = element_text(
  #   angle = 45,
  #   hjust = 1
  # )) +
  theme(axis.text.x = element_text(size = 12, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5))






# the BLenD a panel of posterior predictive checks
breaks <- seq(from = -1, to = 1, length.out = n_bin + 1)
sim_x <- NULL
sim_y <- NULL
group <- NULL
for (row_id in 1:12000 * 100) { # use only 1000 samples from the prior
  sim_y <- c(sim_y, unlist(SS_sim_md[row_id, 1:n_bin]) / sum(unlist(SS_sim_md[row_id, 1:n_bin])) ) # revert the BLenD SSs to its original scale for plotting
  sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  group <- c(group, rep(as.character(row_id), length(breaks)-1))
}
prior_data <- data.frame(x = sim_x, y = sim_y / (2 / n_bin), group)
sim_x <- NULL
sim_y <- NULL
group <- NULL
for (row_id in 1:nrow(SS_sim_acc)) {
  sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, 1:n_bin]) / sum(unlist(SS_sim_acc[row_id, 1:n_bin])) ) # revert the BLenD SSs to its original scale for plotting
  sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  group <- c(group, rep(as.character(row_id), length(breaks)-1))
}
posterior_data <- data.frame(x = sim_x, y = sim_y / (2 / n_bin), group)
line_data <- data.frame(x = (breaks[-1] + breaks[-length(breaks)]) / 2, y = (SS_real[SS_sel[1:n_bin]] / sum(SS_real[SS_sel[1:n_bin]])) / (2 / n_bin)) # revert the BLenD SSs to its original scale for plotting
blend_a_plot <- ggplot() +
  geom_line(data = prior_data, aes(x = x, y = y, group = group), size = 1, color = "gray", alpha = 0.05) +
  geom_line(data = posterior_data, aes(x = x, y = y, group = group), size = 0.5, color = "red", alpha = 1) +
  geom_line(data = line_data, aes(x = x, y = y), linewidth = 0.8, color = "blue") +      # Line plot
  xlim(-1, 1) +
  ylim(0, 40) +
  labs(title = NULL, x = NULL, y = NULL) +
  # labs(title = "BLenD-a", x = expression(d[a]), y = "density") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

# the BLenD b panel of posterior predictive checks
sim_x <- NULL
sim_y <- NULL
group <- NULL
for (row_id in 1:12000 * 100) { # use only 1000 samples from the prior
  sim_y <- c(sim_y, unlist(SS_sim_md[row_id, (n_bin+1):(2*n_bin)]) / sum(unlist(SS_sim_md[row_id, (n_bin+1):(2*n_bin)])) ) # revert the BLenD SSs to its original scale for plotting
  sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  group <- c(group, rep(as.character(row_id), length(breaks)-1))
}
prior_data <- data.frame(x = sim_x, y = sim_y / (2 / n_bin), group)
sim_x <- NULL
sim_y <- NULL
group <- NULL
for (row_id in 1:nrow(SS_sim_acc)) {
  sim_y <- c(sim_y, unlist(SS_sim_acc[row_id, (n_bin+1):(2*n_bin)]) / sum(unlist(SS_sim_acc[row_id, (n_bin+1):(2*n_bin)])) ) # revert the BLenD SSs to its original scale for plotting
  sim_x <- c(sim_x, (breaks[-1] + breaks[-length(breaks)]) / 2)
  group <- c(group, rep(as.character(row_id), length(breaks)-1))
}
posterior_data <- data.frame(x = sim_x, y = sim_y / (2 / n_bin), group)
line_data <- data.frame(x = (breaks[-1] + breaks[-length(breaks)]) / 2, y = (SS_real[SS_sel[(n_bin+1):(2*n_bin)]] / sum(SS_real[SS_sel[(n_bin+1):(2*n_bin)]])) / (2 / n_bin)) # revert the BLenD SSs to its original scale for plotting
blend_b_plot <- ggplot() +
  geom_line(data = prior_data, aes(x = x, y = y, group = group), size = 1, color = "gray", alpha = 0.05) +
  geom_line(data = posterior_data, aes(x = x, y = y, group = group), size = 0.5, color = "red", alpha = 1) +
  geom_line(data = line_data, aes(x = x, y = y), linewidth = 0.8, color = "blue") +      # Line plot
  xlim(-1, 1) +
  ylim(0, 40) +
  labs(title = NULL, x = NULL, y = NULL) +
  # labs(title = "BLenD-a", x = expression(d[a]), y = "density") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))
  
# the host tree size panels of posterior predictive checks
sim_y <- NULL
for (row_id in ids_sim_acc) {
  sim_y <- c(sim_y, unlist(SS_comb[row_id, (n_bin * 2 + 1)]))
}
posterior_data <- data.frame(y = sim_y)
sim_y <- NULL
for (row_id in 1:12000 * 100) {
  sim_y <- c(sim_y, unlist(SS_comb[row_id, (n_bin * 2 + 1)]))
}
prior_data <- data.frame(y = sim_y)
line_data <- SS_real_raw[[1]][SS_sel[(n_bin * 2 + 1)]]
tree_size_h_plot <- ggplot() +  # Empty string for x-axis
  geom_jitter(data = prior_data, mapping = aes(x = "", y = y), width = 0.2, color = "gray", size = 0.5, alpha = 0.2) +
  geom_jitter(data = posterior_data, mapping = aes(x = "", y = y), width = 0.2, color = "red", size = 2, alpha = 1) +
  geom_hline(yintercept = line_data, color = "blue", linewidth = 0.8) +
  # ylim(10, 100) +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_y_log10() +
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
for (row_id in ids_sim_acc) {
  sim_y <- c(sim_y, unlist(SS_comb[row_id, (n_bin * 2 + 2)]))
}
posterior_data <- data.frame(y = sim_y)
sim_y <- NULL
for (row_id in 1:12000 * 100) {
  sim_y <- c(sim_y, unlist(SS_comb[row_id, (n_bin * 2 + 2)]))
}
prior_data <- data.frame(y = sim_y)
line_data <- SS_real_raw[[1]][SS_sel[(n_bin * 2 + 2)]]
tree_size_s_plot <- ggplot() +  # Empty string for x-axis
  geom_jitter(data = prior_data, mapping = aes(x = "", y = y), width = 0.2, color = "gray", size = 0.5, alpha = 0.2) +
  geom_jitter(data = posterior_data, mapping = aes(x = "", y = y), width = 0.2, color = "red", size = 2, alpha = 1) +
  geom_hline(yintercept = line_data, color = "blue", linewidth = 0.8) +
  # ylim(10, 100) +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_y_log10() +
  # labs(title = "Hosts", x = '', y = "Number") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

# panels_OG[[length(panels_OG) + 1]] <- list(spec_plot, ext_plot, blend_a_plot, blend_b_plot, tree_size_h_plot, tree_size_s_plot)
panels_OG[[0 + 1]] <- list(relative_plot, absolute_plot, blend_a_plot, blend_b_plot, tree_size_h_plot, tree_size_s_plot)

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

a <- ggplot() + annotate(geom = 'text', x=1, y=1, label="b) Relative rates", size = 5.2) + theme_void()
b <- ggplot() + annotate(geom = 'text', x=1, y=1, label="c) Absolute rates", size = 5.2) + theme_void()
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
pdf(paste(prefix, 'cophy_ABC_results/', "beetle_results_Fig4_20260607.pdf", sep = ''), width = 10, height = 10, pointsize = 12)
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

