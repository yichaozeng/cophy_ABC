library(ggpubr)
library(patchwork)
library(viridis)
library(dplyr)
library(quantreg)
library(splines)

setwd("/groups/cromanpa/yzeng")

err_plots <- NULL
ratio_plots <- NULL

err_plots <- NULL
ratio_plots_ext_scen <- NULL

file_name <- paste("ex_", 'cophy_ABC_convergence/cross_validation_sizes', '.csv', sep = '')
rel_err <- read.csv(file_name)

# note that "lambda_H" is the relative error of untransformed lambda_H, not the rate estimate. Same for other parameters (all lambdas, exp_H, mu_H_frac, and mu_S_frac)
# we rename these variable to avoid confusion
colnames(rel_err)[2:7] <- c('lambda_H_raw_rel', 'lambda_S_raw_rel', 'lambda_C_raw_rel', 'exp_H_raw_rel', 'mu_H_frac_raw_rel', 'mu_S_frac_raw_rel')

# here, we reverse the relative errors back to absolute errors
rel_err$lambda_H_raw_abs <- rel_err$lambda_H_raw_rel * rel_err$lambda_H_true
rel_err$lambda_S_raw_abs <- rel_err$lambda_S_raw_rel * rel_err$lambda_S_true
rel_err$lambda_C_raw_abs <- rel_err$lambda_C_raw_rel * rel_err$lambda_C_true
rel_err$exp_H_raw_abs <- rel_err$exp_H_raw_rel * rel_err$exp_H_true
rel_err$mu_H_frac_raw_abs <- rel_err$mu_H_frac_raw_rel * rel_err$mu_H_frac_true
rel_err$mu_S_frac_raw_abs <- rel_err$mu_S_frac_raw_rel * rel_err$mu_S_frac_true

# then, we reverse the absolute errors back to estimates
rel_err$lambda_H_est <- rel_err$lambda_H_raw_abs + rel_err$lambda_H_true
rel_err$lambda_S_est <- rel_err$lambda_S_raw_abs + rel_err$lambda_S_true
rel_err$lambda_C_est <- rel_err$lambda_C_raw_abs + rel_err$lambda_C_true
rel_err$exp_H_est <- rel_err$exp_H_raw_abs + rel_err$exp_H_true
rel_err$mu_H_frac_est <- rel_err$mu_H_frac_raw_abs + rel_err$mu_H_frac_true
rel_err$mu_S_frac_est <- rel_err$mu_S_frac_raw_abs + rel_err$mu_S_frac_true

# now we remove the relative and absolute errors of untrnasformed parameters to avoid confusion
rel_err <- rel_err[, -c((2:7), c(17:22))]

# reparameterization
# a function for computing the harmonic mean
div_sum <- function(x, y){
  # x + y
  2 * x * y / (x + y)
}
# # a function for computing the harmonic mean
# div_sum <- function(x, y){
#   x + y
# }
# estimates
rel_err$mu_H_est <- rel_err$mu_H_frac_est * (rel_err$lambda_H_est + rel_err$lambda_C_est)
rel_err$mu_S_est <- rel_err$mu_S_frac_est * (rel_err$lambda_S_est + rel_err$lambda_C_est + rel_err$exp_H_est)

rel_err$lambda_H_frac_est <- rel_err$lambda_H_est / div_sum(rel_err$lambda_H_est + rel_err$lambda_C_est - rel_err$mu_H_est, rel_err$lambda_S_est + rel_err$lambda_C_est + rel_err$exp_H_est - rel_err$mu_S_est)
rel_err$lambda_S_frac_est <- rel_err$lambda_S_est / div_sum(rel_err$lambda_H_est + rel_err$lambda_C_est - rel_err$mu_H_est, rel_err$lambda_S_est + rel_err$lambda_C_est + rel_err$exp_H_est - rel_err$mu_S_est)
rel_err$lambda_C_frac_est <- rel_err$lambda_C_est / div_sum(rel_err$lambda_H_est + rel_err$lambda_C_est - rel_err$mu_H_est, rel_err$lambda_S_est + rel_err$lambda_C_est + rel_err$exp_H_est - rel_err$mu_S_est)
rel_err$exp_H_frac_est <- rel_err$exp_H_est / div_sum(rel_err$lambda_H_est + rel_err$lambda_C_est - rel_err$mu_H_est, rel_err$lambda_S_est + rel_err$lambda_C_est + rel_err$exp_H_est - rel_err$mu_S_est)
rel_err$mu_H_frac_est <- rel_err$mu_H_est / div_sum(rel_err$lambda_H_est + rel_err$lambda_C_est - rel_err$mu_H_est, rel_err$lambda_S_est + rel_err$lambda_C_est + rel_err$exp_H_est - rel_err$mu_S_est)
rel_err$mu_S_frac_est <- rel_err$mu_S_est / div_sum(rel_err$lambda_H_est + rel_err$lambda_C_est - rel_err$mu_H_est, rel_err$lambda_S_est + rel_err$lambda_C_est + rel_err$exp_H_est - rel_err$mu_S_est)

# true values
rel_err$mu_H_true <- rel_err$mu_H_frac_true * (rel_err$lambda_H_true + rel_err$lambda_C_true)
rel_err$mu_S_true <- rel_err$mu_S_frac_true * (rel_err$lambda_S_true + rel_err$lambda_C_true + rel_err$exp_H_true)

rel_err$lambda_H_frac_true <- rel_err$lambda_H_true / div_sum(rel_err$lambda_H_true + rel_err$lambda_C_true - rel_err$mu_H_true, rel_err$lambda_S_true + rel_err$lambda_C_true + rel_err$exp_H_true - rel_err$mu_S_true)
rel_err$lambda_S_frac_true <- rel_err$lambda_S_true / div_sum(rel_err$lambda_H_true + rel_err$lambda_C_true - rel_err$mu_H_true, rel_err$lambda_S_true + rel_err$lambda_C_true + rel_err$exp_H_true - rel_err$mu_S_true)
rel_err$lambda_C_frac_true <- rel_err$lambda_C_true / div_sum(rel_err$lambda_H_true + rel_err$lambda_C_true - rel_err$mu_H_true, rel_err$lambda_S_true + rel_err$lambda_C_true + rel_err$exp_H_true - rel_err$mu_S_true)
rel_err$exp_H_frac_true <- rel_err$exp_H_true / div_sum(rel_err$lambda_H_true + rel_err$lambda_C_true - rel_err$mu_H_true, rel_err$lambda_S_true + rel_err$lambda_C_true + rel_err$exp_H_true - rel_err$mu_S_true)
rel_err$mu_H_frac_true <- rel_err$mu_H_true / div_sum(rel_err$lambda_H_true + rel_err$lambda_C_true - rel_err$mu_H_true, rel_err$lambda_S_true + rel_err$lambda_C_true + rel_err$exp_H_true - rel_err$mu_S_true)
rel_err$mu_S_frac_true <- rel_err$mu_S_true / div_sum(rel_err$lambda_H_true + rel_err$lambda_C_true - rel_err$mu_H_true, rel_err$lambda_S_true + rel_err$lambda_C_true + rel_err$exp_H_true - rel_err$mu_S_true)

# absolute errors
rel_err$lambda_H_frac_abs <- rel_err$lambda_H_frac_est - rel_err$lambda_H_frac_true
rel_err$lambda_S_frac_abs <- rel_err$lambda_S_frac_est - rel_err$lambda_S_frac_true
rel_err$lambda_C_frac_abs <- rel_err$lambda_C_frac_est - rel_err$lambda_C_frac_true
rel_err$exp_H_frac_abs <- rel_err$exp_H_frac_est - rel_err$exp_H_frac_true
rel_err$mu_H_frac_abs <- rel_err$mu_H_frac_est - rel_err$mu_H_frac_true
rel_err$mu_S_frac_abs <- rel_err$mu_S_frac_est - rel_err$mu_S_frac_true
rel_err$lambda_H_abs <- rel_err$lambda_H_est - rel_err$lambda_H_true
rel_err$lambda_S_abs <- rel_err$lambda_S_est - rel_err$lambda_S_true
rel_err$lambda_C_abs <- rel_err$lambda_C_est - rel_err$lambda_C_true
rel_err$exp_H_abs <- rel_err$exp_H_est - rel_err$exp_H_true
rel_err$mu_H_abs <- rel_err$mu_H_est - rel_err$mu_H_true
rel_err$mu_S_abs <- rel_err$mu_S_est - rel_err$mu_S_true

# relative errors
rel_err$lambda_H_frac_rel <- rel_err$lambda_H_frac_abs / rel_err$lambda_H_frac_true
rel_err$lambda_S_frac_rel <- rel_err$lambda_S_frac_abs / rel_err$lambda_S_frac_true
rel_err$lambda_C_frac_rel <- rel_err$lambda_C_frac_abs / rel_err$lambda_C_frac_true
rel_err$exp_H_frac_rel <- rel_err$exp_H_frac_abs / rel_err$exp_H_frac_true
rel_err$mu_H_frac_rel <- rel_err$mu_H_frac_abs / rel_err$mu_H_frac_true
rel_err$mu_S_frac_rel <- rel_err$mu_S_frac_abs / rel_err$mu_S_frac_true
rel_err$lambda_H_rel <- rel_err$lambda_H_abs / rel_err$lambda_H_true
rel_err$lambda_S_rel <- rel_err$lambda_S_abs / rel_err$lambda_S_true
rel_err$lambda_C_rel <- rel_err$lambda_C_abs / rel_err$lambda_C_true
rel_err$exp_H_rel <- rel_err$exp_H_abs / rel_err$exp_H_true
rel_err$mu_H_rel <- rel_err$mu_H_abs / rel_err$mu_H_true
rel_err$mu_S_rel <- rel_err$mu_S_abs / rel_err$mu_S_true

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

row1 <- as.vector(rel_err[1, 5:10])
every <- sum( apply(rel_err[, 5:10], 1, function(row) identical(as.numeric(row), as.numeric(row1))) )
# every <- 1 # this nullify the averaging

# for posterior
rel_err_aver <- data.frame(
  lambda_H_frac_abs_aver = aver(rel_err$lambda_H_frac_abs, every),
  lambda_S_frac_abs_aver = aver(rel_err$lambda_S_frac_abs, every),
  lambda_C_frac_abs_aver = aver(rel_err$lambda_C_frac_abs, every),
  exp_H_frac_abs_aver = aver(rel_err$exp_H_frac_abs, every),
  mu_H_frac_abs_aver = aver(rel_err$mu_H_frac_abs, every),
  mu_S_frac_abs_aver = aver(rel_err$mu_S_frac_abs, every),

  lambda_H_abs_aver = aver(rel_err$lambda_H_abs, every),
  lambda_S_abs_aver = aver(rel_err$lambda_S_abs, every),
  lambda_C_abs_aver = aver(rel_err$lambda_C_abs, every),
  exp_H_abs_aver = aver(rel_err$exp_H_abs, every),
  mu_H_abs_aver = aver(rel_err$mu_H_abs, every),
  mu_S_abs_aver = aver(rel_err$mu_S_abs, every),

  lambda_H_frac_abs_stan_dev = stan_dev(rel_err$lambda_H_frac_abs, every),
  lambda_S_frac_abs_stan_dev = stan_dev(rel_err$lambda_S_frac_abs, every),
  lambda_C_frac_abs_stan_dev = stan_dev(rel_err$lambda_C_frac_abs, every),
  exp_H_frac_abs_stan_dev = stan_dev(rel_err$exp_H_frac_abs, every),
  mu_H_frac_abs_stan_dev = stan_dev(rel_err$mu_H_frac_abs, every),
  mu_S_frac_abs_stan_dev = stan_dev(rel_err$mu_S_frac_abs, every),

  lambda_H_abs_stan_dev = stan_dev(rel_err$lambda_H_abs, every),
  lambda_S_abs_stan_dev = stan_dev(rel_err$lambda_S_abs, every),
  lambda_C_abs_stan_dev = stan_dev(rel_err$lambda_C_abs, every),
  exp_H_abs_stan_dev = stan_dev(rel_err$exp_H_abs, every),
  mu_H_abs_stan_dev = stan_dev(rel_err$mu_H_abs, every),
  mu_S_abs_stan_dev = stan_dev(rel_err$mu_S_abs, every),

  h_tree_size = aver(rel_err$h_tree_size, every),
  s_tree_size = aver(rel_err$s_tree_size, every)#,
  # net_size = aver(rel_err$net_size, every),
  # blend_1_size = aver(rel_err$blend_1_size, every),
  # blend_2_size = aver(rel_err$blend_2_size, every)
)

# # for prior
# lambda_H_frac_prior <- aver(rel_err$lambda_H_frac_true, every)
# lambda_S_frac_prior <- aver(rel_err$lambda_S_frac_true, every)
# lambda_C_frac_prior <- aver(rel_err$lambda_C_frac_true, every)
# exp_H_frac_prior <- aver(rel_err$exp_H_frac_true, every)
# mu_H_frac_prior <- aver(rel_err$mu_H_frac_true, every)
# mu_S_frac_prior <- aver(rel_err$mu_S_frac_true, every)
# 
# lambda_H_prior <- aver(rel_err$lambda_H_true, every)
# lambda_S_prior <- aver(rel_err$lambda_S_true, every)
# lambda_C_prior <- aver(rel_err$lambda_C_true, every)
# exp_H_prior <- aver(rel_err$exp_H_true, every)
# mu_H_prior <- aver(rel_err$mu_H_true, every)
# mu_S_prior <- aver(rel_err$mu_S_true, every)
# 
# rel_err_aver <- data.frame(
#   lambda_H_frac_abs_aver = mean(lambda_H_frac_prior) - lambda_H_frac_prior,
#   lambda_S_frac_abs_aver = mean(lambda_S_frac_prior) - lambda_S_frac_prior,
#   lambda_C_frac_abs_aver = mean(lambda_C_frac_prior) - lambda_C_frac_prior,
#   exp_H_frac_abs_aver = mean(exp_H_frac_prior) - exp_H_frac_prior,
#   mu_H_frac_abs_aver = mean(mu_H_frac_prior) - mu_H_frac_prior,
#   mu_S_frac_abs_aver = mean(mu_S_frac_prior) - mu_S_frac_prior,
# 
#   lambda_H_abs_aver = mean(lambda_H_prior) - lambda_H_prior,
#   lambda_S_abs_aver = mean(lambda_S_prior) - lambda_S_prior,
#   lambda_C_abs_aver = mean(lambda_C_prior) - lambda_C_prior,
#   exp_H_abs_aver = mean(exp_H_prior) - exp_H_prior,
#   mu_H_abs_aver = mean(mu_H_prior) - mu_H_prior,
#   mu_S_abs_aver = mean(mu_S_prior) - mu_S_prior,
# 
#   lambda_H_frac_abs_stan_dev = sd(lambda_H_frac_prior),
#   lambda_S_frac_abs_stan_dev = sd(lambda_S_frac_prior),
#   lambda_C_frac_abs_stan_dev = sd(lambda_C_frac_prior),
#   exp_H_frac_abs_stan_dev = sd(exp_H_frac_prior),
#   mu_H_frac_abs_stan_dev = sd(mu_H_frac_prior),
#   mu_S_frac_abs_stan_dev = sd(mu_S_frac_prior),
# 
#   lambda_H_abs_stan_dev = sd(lambda_H_prior),
#   lambda_S_abs_stan_dev = sd(lambda_S_prior),
#   lambda_C_abs_stan_dev = sd(lambda_C_prior),
#   exp_H_abs_stan_dev = sd(exp_H_prior),
#   mu_H_abs_stan_dev = sd(mu_H_prior),
#   mu_S_abs_stan_dev = sd(mu_S_prior),
# 
#   h_tree_size = aver(rel_err$h_tree_size, every),
#   s_tree_size = aver(rel_err$s_tree_size, every)#,
#   # net_size = aver(rel_err$net_size, every),
#   # blend_1_size = aver(rel_err$blend_1_size, every),
#   # blend_2_size = aver(rel_err$blend_2_size, every)
# )

# the smaller of two tree sizes
# x <- log10(rel_err_aver$h_tree_size * (rel_err_aver$h_tree_size <= rel_err_aver$s_tree_size) + rel_err_aver$s_tree_size * (rel_err_aver$h_tree_size > rel_err_aver$s_tree_size))
# or the harmonic means of the two tree sizes
x <- log10( 2* rel_err_aver$h_tree_size * rel_err_aver$s_tree_size / (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) )

# compute the z-scores
y1 <- abs(rel_err_aver$lambda_H_frac_abs_aver) / rel_err_aver$lambda_H_frac_abs_stan_dev
y2 <- abs(rel_err_aver$lambda_S_frac_abs_aver) / rel_err_aver$lambda_S_frac_abs_stan_dev
y3 <- abs(rel_err_aver$lambda_C_frac_abs_aver) / rel_err_aver$lambda_C_frac_abs_stan_dev
y4 <- abs(rel_err_aver$exp_H_frac_abs_aver) / rel_err_aver$exp_H_frac_abs_stan_dev
y5 <- abs(rel_err_aver$mu_H_frac_abs_aver) / rel_err_aver$mu_H_frac_abs_stan_dev
y6 <- abs(rel_err_aver$mu_S_frac_abs_aver) / rel_err_aver$mu_S_frac_abs_stan_dev
y7 <- abs(rel_err_aver$lambda_H_abs_aver) / rel_err_aver$lambda_H_abs_stan_dev
y8 <- abs(rel_err_aver$lambda_S_abs_aver) / rel_err_aver$lambda_S_abs_stan_dev
y9 <- abs(rel_err_aver$lambda_C_abs_aver) / rel_err_aver$lambda_C_abs_stan_dev
y10 <- abs(rel_err_aver$exp_H_abs_aver) / rel_err_aver$exp_H_abs_stan_dev
y11 <- abs(rel_err_aver$mu_H_abs_aver) / rel_err_aver$mu_H_abs_stan_dev
y12 <- abs(rel_err_aver$mu_S_abs_aver) / rel_err_aver$mu_S_abs_stan_dev

for (parameter in 1:12) {
  
  y <- rbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12)[parameter,]
  df <- data.frame(x=x, y=y)
  # remove NAs and infinite values if present
  df <- df[!is.na(df$y),]
  df <- df[is.finite(df$y),]
  
  # mod <- lm(y ~ x, data = df)
  # df_pred <- cbind(df, predict(mod, interval = "prediction", level = 0.90))
  
  qs <- c(0.25, 0.5, 0.75)
  
  # Create a sequence of x values for smooth prediction
  x_seq <- seq(min(df$x), max(df$x), length.out = 200)
  
  # Function to estimate local sample density
  local_density <- sapply(x_seq, function(x) {
    sum(abs(df$x - x) < 0.3)  # window size = 0.3
  })
  
  # Fit smooth quantile regressions using splines
  preds <- lapply(qs, function(q) {
    fit <- rq(y ~ bs(x, df = 3), tau = q, data = df)  # spline fit
    data.frame(
      x = x_seq,
      y = predict(fit, newdata = data.frame(x = x_seq)),
      quantile = q,
      n_local = local_density
    )
  }) %>%
    bind_rows()
  
  # keep only prediction based on high densities
  preds <- preds[preds$n_local >= 90, ]
  
  # Extract end points for labeling
  labels <- preds %>%
    group_by(quantile) %>%
    filter(x == max(x)) %>%
    mutate(label = paste0(quantile * 100, "%"))
  
  title_char <- list(
    
    expression(italic(paste(frac(lambda[H], r[cophy])))),
    expression(italic(paste(frac(lambda[S], r[cophy])))),
    expression(italic(paste(frac(lambda[C], r[cophy])))),
    expression(italic(paste(frac(lambda[W], r[cophy])))),
    expression(italic(paste(frac(mu[H], r[cophy])))),
    expression(italic(paste(frac(mu[S], r[cophy])))),
    
    expression(italic(paste(lambda[H]))),
    expression(italic(paste(lambda[S]))),
    expression(italic(paste(lambda[C]))),
    expression(italic(paste(lambda[W]))),
    expression(italic(paste(mu[H]))),
    expression(italic(paste(mu[S])))
  )[[parameter]]
  
  err_plots[[length(err_plots) + 1]] <- ggplot(df, aes(x = x, y = y)) +
    # geom_point(size = 1, alpha = 1, color = point_color) +
    geom_hex(bins = 20) +
    # geom_point(size = 1.8, alpha = 0.5)
    scale_fill_viridis(option="viridis") +
    # scale_fill_viridis(option="viridis", limits = c(0, 22)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_line(
      data = preds,
      aes(y = y, group = quantile),
      linewidth = 0.7,
      color = 'white'
    ) +
    geom_text(
      data = labels,
      aes(y = y, label = label),
      hjust = -0.1, size = 3, color = 'white'
    ) +
    xlim(0, 3.3) +
    # ylim(-4, 2) +
    ylim(0, 3) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(panel.background = element_rect(fill = "grey75")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = title_char, x = NULL, y = NULL)

}


# arrange the panels
x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label=expression(log[10](Cophylogeny~size)), size = 5) + theme_void()
y_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Posterior z-score (absolute value)", size = 5, angle = 90) + theme_void()

plot_list_ratio <- list(
  x = x_lab,
  y = y_lab,
  a = err_plots[[1]],
  b = err_plots[[2]],
  c = err_plots[[3]],
  d = err_plots[[4]],
  e = err_plots[[5]],
  f = err_plots[[6]],
  
  g = err_plots[[7]],
  h = err_plots[[8]],
  i = err_plots[[9]],
  j = err_plots[[10]],
  k = err_plots[[11]],
  l = err_plots[[12]]
)

layoutplot <- "
yaabbcc
yaabbcc
yaabbcc
yaabbcc
yddeeff
yddeeff
yddeeff
yddeeff
ygghhii
ygghhii
ygghhii
ygghhii
yjjkkll
yjjkkll
yjjkkll
yjjkkll
#xxxxxx
"

print('plotting')
#png(paste(prefix, 'cophy_ABC_convergence/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
pdf(paste('ex_cophy_ABC_convergence/', "ratios", ".pdf", sep = ''), width = 12, height = 15, pointsize = 12)
print(wrap_plots(plot_list_ratio, guides = 'keep', design = layoutplot) &
        theme(legend.position = "right",
          legend.text = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.5, "lines"),
          legend.spacing.y = unit(0.1, "cm")
        ))
dev.off()
