library(ggpubr)
library(patchwork)
library(viridis)
library(dplyr)
library(quantreg)
library(splines)


setwd("/groups/cromanpa/yzeng")

err_plots <- NULL
ratio_plots <- NULL

file_name <- paste("ex_", 'cophy_ABC_convergence/cross_validation_sizes', '.csv', sep = '')
rel_err <- read.csv(file_name)

# here, we reverse the relative errors back to absolute errors
rel_err$lambda_H_abs <- rel_err$lambda_H * rel_err$lambda_H_true
rel_err$lambda_S_abs <- rel_err$lambda_S * rel_err$lambda_S_true
rel_err$lambda_C_abs <- rel_err$lambda_C * rel_err$lambda_C_true
rel_err$exp_H_abs <- rel_err$exp_H * rel_err$exp_H_true
rel_err$mu_H_frac_abs <- rel_err$mu_H_frac * rel_err$mu_H_frac_true
rel_err$mu_S_frac_abs <- rel_err$mu_S_frac * rel_err$mu_H_frac_true

# averaging function for every X rows
aver <- function(vec, every){
  temp <- tapply(vec, (seq_along(vec)-1) %/% every, mean)
  return(temp)
}

# "standard deviation"ing function for every X rows
stan_dev <- function(vec, every){
  temp <- tapply(vec, (seq_along(vec)-1) %/% every, sd)
  return(temp)
}

row1 <- as.vector(rel_err[1, 9:14])
every <- sum( apply(rel_err[, 9:14], 1, function(row) identical(as.numeric(row), as.numeric(row1))) )
# every <- 1 # this nullify the averaging

# for posterior
rel_err_aver_posterior <- data.frame(
  lambda_H_abs_aver = aver(rel_err$lambda_H_abs, every),
  lambda_S_abs_aver = aver(rel_err$lambda_S_abs, every),
  lambda_C_abs_aver = aver(rel_err$lambda_C_abs, every),
  exp_H_abs_aver = aver(rel_err$exp_H_abs, every),
  mu_H_frac_abs_aver = aver(rel_err$mu_H_frac_abs, every),
  mu_S_frac_abs_aver = aver(rel_err$mu_S_frac_abs, every),

  lambda_H_abs_stan_dev = stan_dev(rel_err$lambda_H_abs, every),
  lambda_S_abs_stan_dev = stan_dev(rel_err$lambda_S_abs, every),
  lambda_C_abs_stan_dev = stan_dev(rel_err$lambda_C_abs, every),
  exp_H_abs_stan_dev = stan_dev(rel_err$exp_H_abs, every),
  mu_H_frac_abs_stan_dev = stan_dev(rel_err$mu_H_frac_abs, every),
  mu_S_frac_abs_stan_dev = stan_dev(rel_err$mu_S_frac_abs, every),

  h_tree_size = aver(rel_err$h_tree_size, every),
  s_tree_size = aver(rel_err$s_tree_size, every)#,
  # net_size = aver(rel_err$net_size, every),
  # blend_1_size = aver(rel_err$blend_1_size, every),
  # blend_2_size = aver(rel_err$blend_2_size, every)
)

# for prior
lambda_H_prior <- aver(rel_err$lambda_H_true, every)
lambda_S_prior <- aver(rel_err$lambda_S_true, every)
lambda_C_prior <- aver(rel_err$lambda_C_true, every)
exp_H_prior <- aver(rel_err$exp_H_true, every)
mu_H_frac_prior <- aver(rel_err$mu_H_frac_true, every)
mu_S_frac_prior <- aver(rel_err$mu_S_frac_true, every)

rel_err_aver_prior <- data.frame(
  lambda_H_abs_aver = mean(lambda_H_prior) - lambda_H_prior,
  lambda_S_abs_aver = mean(lambda_S_prior) - lambda_S_prior,
  lambda_C_abs_aver = mean(lambda_C_prior) - lambda_C_prior,
  exp_H_abs_aver = mean(exp_H_prior) - exp_H_prior,
  mu_H_frac_abs_aver = mean(mu_H_frac_prior) - mu_H_frac_prior,
  mu_S_frac_abs_aver = mean(mu_S_frac_prior) - mu_S_frac_prior,

  lambda_H_abs_stan_dev = sd(lambda_H_prior),
  lambda_S_abs_stan_dev = sd(lambda_S_prior),
  lambda_C_abs_stan_dev = sd(lambda_C_prior),
  exp_H_abs_stan_dev = sd(exp_H_prior),
  mu_H_frac_abs_stan_dev = sd(mu_H_frac_prior),
  mu_S_frac_abs_stan_dev = sd(mu_S_frac_prior),

  h_tree_size = aver(rel_err$h_tree_size, every),
  s_tree_size = aver(rel_err$s_tree_size, every)#,
  # net_size = aver(rel_err$net_size, every),
  # blend_1_size = aver(rel_err$blend_1_size, every),
  # blend_2_size = aver(rel_err$blend_2_size, every)
)

for (parameter in 1:6) {
  
  box_data_prior_post <- NULL
  
  for (temp_id in 1:2) {
    
    rel_err_aver <- list(rel_err_aver_prior, rel_err_aver_posterior)[[temp_id ]]
    
    # the harmonic means of the two tree sizes
    x <- log10( 2* rel_err_aver$h_tree_size * rel_err_aver$s_tree_size / (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) )
    
    y1 <- abs(rel_err_aver$lambda_H_abs_aver) / rel_err_aver$lambda_H_abs_stan_dev
    y2 <- abs(rel_err_aver$lambda_S_abs_aver) / rel_err_aver$lambda_S_abs_stan_dev
    y3 <- abs(rel_err_aver$lambda_C_abs_aver) / rel_err_aver$lambda_C_abs_stan_dev
    y4 <- abs(rel_err_aver$exp_H_abs_aver) / rel_err_aver$exp_H_abs_stan_dev
    y5 <- abs(rel_err_aver$mu_H_frac_abs_aver) / rel_err_aver$mu_H_frac_abs_stan_dev
    y6 <- abs(rel_err_aver$mu_S_frac_abs_aver) / rel_err_aver$mu_S_frac_abs_stan_dev
    
    y <- rbind(y1, y2, y3, y4, y5, y6)[parameter,]
    df <- data.frame(x=x, y=y)
    
    # discretize cophylogeny size and add variable indicating prior or posterior
    df$x[df$x <= 1] <- '0-10'
    df$x[df$x > 1 & df$x <= 2] <- '11-100'
    df$x[df$x > 2] <- '>101'
    df$x <- factor(df$x, levels = c('0-10', '11-100', '>101'))
    
    df$prior_posterior <- rep(c('prior', 'posterior')[temp_id], length(df$x))
    df$prior_posterior <- factor(df$prior_posterior, levels = c('prior', 'posterior'))
    
    # stack the data frames for the prior and posterior
    box_data_prior_post <- rbind(box_data_prior_post, df)
    
  }
  
  title_char <- list(
    expression(hat(lambda)[H]),
    expression(hat(lambda)[S]),
    expression(hat(lambda)[C]),
    expression(hat(lambda)[W]),
    expression(hat(epsilon)[H]),
    expression(hat(epsilon)[S])
  )[[parameter]]
  
  err_plots[[length(err_plots) + 1]] <- ggplot(box_data_prior_post, aes(x=x, y=log(y), fill=prior_posterior)) + 
    geom_split_violin()
  
}

# arrange the panels
x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Cophylogeny size", size = 5) + theme_void()
y_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "z-score (unsigned)", size = 5, angle = 90) + theme_void()

plot_list_ratio <- list(
  a = err_plots[[1]],
  b = err_plots[[2]],
  c = err_plots[[3]],
  d = err_plots[[4]],
  e = err_plots[[5]],
  f = err_plots[[6]],
  x = x_lab,
  y = y_lab
)











library(ggpubr)
library(patchwork)
library(viridis)
library(dplyr)
library(quantreg)
library(splines)

setwd("/groups/cromanpa/yzeng")

err_plots <- NULL

file_name <- paste("ex_", 'cophy_ABC_convergence/cross_validation_sizes', '.csv', sep = '')
rel_err <- read.csv(file_name)

# averaging function for every X rows
aver <- function(vec, every){
  temp <- tapply(vec, (seq_along(vec)-1) %/% every, mean)
  return(temp)
}

row1 <- as.vector(rel_err[1,9:14])
every <- sum( apply(rel_err[,9:14], 1, function(row) identical(as.numeric(row), as.numeric(row1))) )
# every <- 1 # this nullify the averaging

# for posterior
rel_err_aver_posterior <- data.frame(
    lambda_H = aver(rel_err$lambda_H, every),
    lambda_S = aver(rel_err$lambda_S, every),
    lambda_C = aver(rel_err$lambda_C, every),
    exp_H = aver(rel_err$exp_H, every),
    mu_H_frac = aver(rel_err$mu_H_frac, every),
    mu_S_frac = aver(rel_err$mu_S_frac, every),
    h_tree_size = aver(rel_err$h_tree_size, every),
    s_tree_size = aver(rel_err$s_tree_size, every)#,
    # net_size = aver(rel_err$net_size, every),
    # blend_1_size = aver(rel_err$blend_1_size, every),
    # blend_2_size = aver(rel_err$blend_2_size, every)
)

# for prior
lambda_H_prior <- aver(rel_err$lambda_H_true, every)
lambda_S_prior <- aver(rel_err$lambda_S_true, every)
lambda_C_prior <- aver(rel_err$lambda_C_true, every)
exp_H_prior <- aver(rel_err$exp_H_true, every)
mu_H_frac_prior <- aver(rel_err$mu_H_frac_true, every)
mu_S_frac_prior <- aver(rel_err$mu_S_frac_true, every)

rel_err_aver_prior <- data.frame(
  lambda_H = (mean(lambda_H_prior) - lambda_H_prior) / lambda_H_prior,
  lambda_S = (mean(lambda_S_prior) - lambda_S_prior) / lambda_S_prior,
  lambda_C = (mean(lambda_C_prior) - lambda_C_prior) / lambda_C_prior,
  exp_H = (mean(exp_H_prior) - exp_H_prior) / exp_H_prior,
  mu_H_frac = (mean(mu_H_frac_prior) - mu_H_frac_prior) / mu_H_frac_prior,
  mu_S_frac = (mean(mu_S_frac_prior) - mu_S_frac_prior) / mu_S_frac_prior,
  
  h_tree_size = aver(rel_err$h_tree_size, every),
  s_tree_size = aver(rel_err$s_tree_size, every)#,
  # net_size = aver(rel_err$net_size, every),
  # blend_1_size = aver(rel_err$blend_1_size, every),
  # blend_2_size = aver(rel_err$blend_2_size, every)
)

for (parameter in 1:6) {
  
  preds_prior_post <- NULL
  labels_prior_post <- NULL
  
  for (temp_id in 1:2) {
    
    rel_err_aver <- list(rel_err_aver_prior, rel_err_aver_posterior)[[temp_id ]]
    
    # the harmonic means of the two tree sizes
    x <- log10( 2* rel_err_aver$h_tree_size * rel_err_aver$s_tree_size / (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) )
    
    # relative errors
    y1 <- log(abs(rel_err_aver$lambda_H))
    y2 <- log(abs(rel_err_aver$lambda_S))
    y3 <- log(abs(rel_err_aver$lambda_C))
    y4 <- log(abs(rel_err_aver$exp_H))
    y5 <- log(abs(rel_err_aver$mu_H_frac))
    y6 <- log(abs(rel_err_aver$mu_S_frac))
    
    # # absolute errors
    # y1 <- log(abs(rel_err_aver$lambda_H) * lambda_H_prior)
    # y2 <- log(abs(rel_err_aver$lambda_S) * lambda_S_prior)
    # y3 <- log(abs(rel_err_aver$lambda_C) * lambda_C_prior)
    # y4 <- log(abs(rel_err_aver$exp_H) * exp_H_prior)
    # y5 <- log(abs(rel_err_aver$mu_H_frac) * mu_H_frac_prior)
    # y6 <- log(abs(rel_err_aver$mu_S_frac) * mu_S_frac_prior)
    
    y <- rbind(y1, y2, y3, y4, y5, y6)[parameter,]
    
    df <- data.frame(x=x, y=y)
    # discretize cophylogeny size and add variable indicating prior or posterior
    df$x[df$x <= 1] <- '0-10'
    df$x[df$x > 1 & df$x <= 2] <- '11-100'
    df$x[df$x > 2] <- '>101'
    df$x <- factor(df$x, levels = c('0-10', '11-100', '>101'))
    
    df$prior_posterior <- rep(c('prior', 'posterior')[temp_id], length(df$x))
    df$prior_posterior <- factor(df$prior_posterior, levels = c('prior', 'posterior'))
    
    # stack the data frames for the prior and posterior
    box_data_prior_post <- rbind(box_data_prior_post, df)

  }
  
  title_char <- list(
    expression(hat(lambda)[H]),
    expression(hat(lambda)[S]),
    expression(hat(lambda)[C]),
    expression(hat(lambda)[W]),
    expression(hat(epsilon)[H]),
    expression(hat(epsilon)[S])
  )[[parameter]]
  
  err_plots[[length(err_plots) + 1]] <- ggplot(box_data_prior_post, aes(x=x, y=y, fill=prior_posterior)) + 
    geom_split_violin()
  
}

# arrange the panels
x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Cophylogeny size", size = 5) + theme_void()
y_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Relative error (log-transformed)", size = 5, angle = 90) + theme_void()

plot_list_error <- list(
  a = err_plots[[1]],
  b = err_plots[[2]],
  c = err_plots[[3]],
  d = err_plots[[4]],
  e = err_plots[[5]],
  f = err_plots[[6]],
  x = x_lab,
  y = y_lab
)












# run this after "plot_error_precision_ratio_cophylo_sizes.R" and "plot_accuracy_relative_error_cophylo_sizes.R"
plot_list_fig_4 <- list(
  A = plot_list_error$a,
  B = plot_list_error$b,
  C = plot_list_error$c,
  D = plot_list_error$d,
  E = plot_list_error$e,
  'F' = plot_list_error$f,
  X = plot_list_error$x,
  Y = plot_list_error$y,
  
  a = plot_list_ratio$a,
  b = plot_list_ratio$b,
  c = plot_list_ratio$c,
  d = plot_list_ratio$d,
  e = plot_list_ratio$e,
  f = plot_list_ratio$f,
  x = plot_list_ratio$x,
  y = plot_list_ratio$y
)

layoutplot_fig_4 <- "
YAABBCCDDEEFF
YAABBCCDDEEFF
YAABBCCDDEEFF
YAABBCCDDEEFF
YAABBCCDDEEFF
YAABBCCDDEEFF
#XXXXXXXXXXXX
yaabbccddeeff
yaabbccddeeff
yaabbccddeeff
yaabbccddeeff
yaabbccddeeff
yaabbccddeeff
#xxxxxxxxxxxx
"

print('plotting')
#png(paste(prefix, 'cophy_ABC_convergence/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
pdf(paste('ex_cophy_ABC_convergence/', "fig_4_errors_ratios_box_20260331", ".pdf", sep = ''), width = 11, height = 9, pointsize = 12)
print(wrap_plots(plot_list_fig_4, guides = 'collect', design = layoutplot_fig_4) &
        theme(legend.position = "left",
              legend.text = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
              legend.title = element_text(size = 9),
              legend.key.size = unit(0.5, "lines"),
              legend.spacing.y = unit(0.1, "cm")
        ))
dev.off()
