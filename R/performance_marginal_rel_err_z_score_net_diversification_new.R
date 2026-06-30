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
# # a function for computing the algorithmic mean
# div_sum <- function(x, y){
#   # x + y
#   (x + y) / 2
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
aver <- function(vec, every = length(vec)){ # for absolute rates
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
rel_err_aver_posterior <- data.frame(

    lambda_H_frac = aver(rel_err$lambda_H_frac_rel, every),
    lambda_S_frac = aver(rel_err$lambda_S_frac_rel, every),
    lambda_C_frac = aver(rel_err$lambda_C_frac_rel, every),
    exp_H_frac = aver(rel_err$exp_H_frac_rel, every),
    mu_H_frac = aver(rel_err$mu_H_frac_rel, every),
    mu_S_frac = aver(rel_err$mu_S_frac_rel, every),
    
    lambda_H = aver(rel_err$lambda_H_rel, every),
    lambda_S = aver(rel_err$lambda_S_rel, every),
    lambda_C = aver(rel_err$lambda_C_rel, every),
    exp_H = aver(rel_err$exp_H_rel, every),
    mu_H = aver(rel_err$mu_H_rel, every),
    mu_S = aver(rel_err$mu_S_rel, every),
    
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
    s_tree_size = aver(rel_err$s_tree_size, every)
)

# for prior
lambda_H_frac_prior <- aver(rel_err$lambda_H_frac_true, every)
lambda_S_frac_prior <- aver(rel_err$lambda_S_frac_true, every)
lambda_C_frac_prior <- aver(rel_err$lambda_C_frac_true, every)
exp_H_frac_prior <- aver(rel_err$exp_H_frac_true, every)
mu_H_frac_prior <- aver(rel_err$mu_H_frac_true, every)
mu_S_frac_prior <- aver(rel_err$mu_S_frac_true, every)

lambda_H_prior <- aver(rel_err$lambda_H_true, every)
lambda_S_prior <- aver(rel_err$lambda_S_true, every)
lambda_C_prior <- aver(rel_err$lambda_C_true, every)
exp_H_prior <- aver(rel_err$exp_H_true, every)
mu_H_prior <- aver(rel_err$mu_H_true, every)
mu_S_prior <- aver(rel_err$mu_S_true, every)

# transform the data frame for the 1st_prior-level inference
rel_err_aver_prior <- data.frame(
  lambda_H_frac = (mean(lambda_H_frac_prior) - lambda_H_frac_prior) / lambda_H_frac_prior,
  lambda_S_frac = (mean(lambda_S_frac_prior) - lambda_S_frac_prior) / lambda_S_frac_prior,
  lambda_C_frac = (mean(lambda_C_frac_prior) - lambda_C_frac_prior) / lambda_C_frac_prior,
  exp_H_frac = (mean(exp_H_frac_prior) - exp_H_frac_prior) / exp_H_frac_prior,
  mu_H_frac = (mean(mu_H_frac_prior) - mu_H_frac_prior) / mu_H_frac_prior,
  mu_S_frac = (mean(mu_S_frac_prior) - mu_S_frac_prior) / mu_S_frac_prior,
  
  lambda_H = (mean(lambda_H_prior) - lambda_H_prior) / lambda_H_prior,
  lambda_S = (mean(lambda_S_prior) - lambda_S_prior) / lambda_S_prior,
  lambda_C = (mean(lambda_C_prior) - lambda_C_prior) / lambda_C_prior,
  exp_H = (mean(exp_H_prior) - exp_H_prior) / exp_H_prior,
  mu_H = (mean(mu_H_prior) - mu_H_prior) / mu_H_prior,
  mu_S = (mean(mu_S_prior) - mu_S_prior) / mu_S_prior,
  
  lambda_H_frac_abs_aver = mean(lambda_H_frac_prior) - lambda_H_frac_prior,
  lambda_S_frac_abs_aver = mean(lambda_S_frac_prior) - lambda_S_frac_prior,
  lambda_C_frac_abs_aver = mean(lambda_C_frac_prior) - lambda_C_frac_prior,
  exp_H_frac_abs_aver = mean(exp_H_frac_prior) - exp_H_frac_prior,
  mu_H_frac_abs_aver = mean(mu_H_frac_prior) - mu_H_frac_prior,
  mu_S_frac_abs_aver = mean(mu_S_frac_prior) - mu_S_frac_prior,

  lambda_H_abs_aver = mean(lambda_H_prior) - lambda_H_prior,
  lambda_S_abs_aver = mean(lambda_S_prior) - lambda_S_prior,
  lambda_C_abs_aver = mean(lambda_C_prior) - lambda_C_prior,
  exp_H_abs_aver = mean(exp_H_prior) - exp_H_prior,
  mu_H_abs_aver = mean(mu_H_prior) - mu_H_prior,
  mu_S_abs_aver = mean(mu_S_prior) - mu_S_prior,

  lambda_H_frac_abs_stan_dev = sd(lambda_H_frac_prior),
  lambda_S_frac_abs_stan_dev = sd(lambda_S_frac_prior),
  lambda_C_frac_abs_stan_dev = sd(lambda_C_frac_prior),
  exp_H_frac_abs_stan_dev = sd(exp_H_frac_prior),
  mu_H_frac_abs_stan_dev = sd(mu_H_frac_prior),
  mu_S_frac_abs_stan_dev = sd(mu_S_frac_prior),

  lambda_H_abs_stan_dev = sd(lambda_H_prior),
  lambda_S_abs_stan_dev = sd(lambda_S_prior),
  lambda_C_abs_stan_dev = sd(lambda_C_prior),
  exp_H_abs_stan_dev = sd(exp_H_prior),
  mu_H_abs_stan_dev = sd(mu_H_prior),
  mu_S_abs_stan_dev = sd(mu_S_prior),
  
  h_tree_size = aver(rel_err$h_tree_size, every),
  s_tree_size = aver(rel_err$s_tree_size, every)
)

for (parameter in 1:24) { # 1:12 - relative errors; 13:24 - z-scores
  
  preds_prior_post <- NULL
  labels_prior_post <- NULL
  
  ribbon_prior_post <- NULL
  
  for (temp_id in 1:2) {
    
    rel_err_aver <- list(rel_err_aver_prior, rel_err_aver_posterior)[[temp_id ]]
    
    # the harmonic means of the two tree sizes
    x <- log10( 2* rel_err_aver$h_tree_size * rel_err_aver$s_tree_size / (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) )
    
    # # the algorithmic means of the two tree sizes
    # x <- log10( (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) / 2 )
    
    # relative errors
    y1 <- log(abs(rel_err_aver$lambda_H_frac))
    y2 <- log(abs(rel_err_aver$lambda_S_frac))
    y3 <- log(abs(rel_err_aver$lambda_C_frac))
    y4 <- log(abs(rel_err_aver$exp_H_frac))
    
    y5 <- log(abs(rel_err_aver$mu_H_frac))
    y6 <- log(abs(rel_err_aver$mu_S_frac))
    
    y7 <- log(abs(rel_err_aver$lambda_H))
    y8 <- log(abs(rel_err_aver$lambda_S))
    y9 <- log(abs(rel_err_aver$lambda_C))
    y10 <- log(abs(rel_err_aver$exp_H))

    y11 <- log(abs(rel_err_aver$mu_H))
    y12 <- log(abs(rel_err_aver$mu_S))
    
    # z-scores
    y13 <- abs(rel_err_aver$lambda_H_frac_abs_aver) / rel_err_aver$lambda_H_frac_abs_stan_dev
    y14 <- abs(rel_err_aver$lambda_S_frac_abs_aver) / rel_err_aver$lambda_S_frac_abs_stan_dev
    y15 <- abs(rel_err_aver$lambda_C_frac_abs_aver) / rel_err_aver$lambda_C_frac_abs_stan_dev
    y16 <- abs(rel_err_aver$exp_H_frac_abs_aver) / rel_err_aver$exp_H_frac_abs_stan_dev
    y17 <- abs(rel_err_aver$mu_H_frac_abs_aver) / rel_err_aver$mu_H_frac_abs_stan_dev
    y18 <- abs(rel_err_aver$mu_S_frac_abs_aver) / rel_err_aver$mu_S_frac_abs_stan_dev
    y19 <- abs(rel_err_aver$lambda_H_abs_aver) / rel_err_aver$lambda_H_abs_stan_dev
    y20 <- abs(rel_err_aver$lambda_S_abs_aver) / rel_err_aver$lambda_S_abs_stan_dev
    y21 <- abs(rel_err_aver$lambda_C_abs_aver) / rel_err_aver$lambda_C_abs_stan_dev
    y22 <- abs(rel_err_aver$exp_H_abs_aver) / rel_err_aver$exp_H_abs_stan_dev
    y23 <- abs(rel_err_aver$mu_H_abs_aver) / rel_err_aver$mu_H_abs_stan_dev
    y24 <- abs(rel_err_aver$mu_S_abs_aver) / rel_err_aver$mu_S_abs_stan_dev
    
    y <- rbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, y21, y22, y23, y24)[parameter,]
    
    df <- data.frame(x=x, y=y)
    # mod <- lm(y ~ x, data = df)
    # df_pred <- cbind(df, predict(mod, interval = "prediction", level = 0.90))
    
    # qs <- c(0.25, 0.5, 0.75)
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
    preds_prior_post[[length(preds_prior_post) + 1]] <- preds
    
    ribbon <- data.frame(
      x = preds$x[preds$quantile == qs[1]],
      ymin = preds$y[preds$quantile == qs[1]],
      ymax = preds$y[preds$quantile == qs[length(qs)]]
    )
    ribbon_prior_post[[length(ribbon_prior_post) + 1]] <- ribbon
    
    # Extract end points for labeling
    labels <- preds %>%
      group_by(quantile) %>%
      filter(x == range(x)[temp_id]) %>%
      mutate(label = paste0(quantile * 100, "%"))
    labels_prior_post[[length(labels_prior_post) + 1]] <- labels
    
  }
  
  title_char <- list(
    
    expression(italic(paste(frac(lambda[H], r)))),
    expression(italic(paste(frac(lambda[S], r)))),
    expression(italic(paste(frac(lambda[C], r)))),
    expression(italic(paste(frac(lambda[W], r)))),
    expression(italic(paste(frac(mu[H], r)))),
    expression(italic(paste(frac(mu[S], r)))),
    
    expression(italic(paste(lambda[H]))),
    expression(italic(paste(lambda[S]))),
    expression(italic(paste(lambda[C]))),
    expression(italic(paste(lambda[W]))),
    expression(italic(paste(mu[H]))),
    expression(italic(paste(mu[S]))),
    
    expression(italic(paste(frac(lambda[H], r)))),
    expression(italic(paste(frac(lambda[S], r)))),
    expression(italic(paste(frac(lambda[C], r)))),
    expression(italic(paste(frac(lambda[W], r)))),
    expression(italic(paste(frac(mu[H], r)))),
    expression(italic(paste(frac(mu[S], r)))),
    
    expression(italic(paste(lambda[H]))),
    expression(italic(paste(lambda[S]))),
    expression(italic(paste(lambda[C]))),
    expression(italic(paste(lambda[W]))),
    expression(italic(paste(mu[H]))),
    expression(italic(paste(mu[S])))
  )[[parameter]]
  
  err_plots[[length(err_plots) + 1]] <- ggplot(df, aes(x = x, y = y)) +
    # geom_hex(bins = 20) +
    # geom_point(size = 1.8, alpha = 0.9)
    # geom_smooth(data = subset(df, x >= quantile(df$x, 0) & x <= quantile(df$x, 1)), method = "loess", color = "white", linewidth = 0.9, se = FALSE) +
    # scale_fill_viridis(option="magma") +
    # scale_fill_viridis(option="magma", limits = c(0, 32)) +
    # geom_ribbon(aes(ymin = lwr, ymax = upr),
    #             fill = "chocolate", alpha = 0.75) +
    # geom_point(aes(y = y)) +
    # geom_line(aes(y = fit), colour = "chocolate", size = 1) +
    #geom_smooth(aes(y = y), method = "loess", se = TRUE, level = 0.95, color = "chocolate", fill = "lightchocolate") +
    geom_ribbon(
      data = ribbon_prior_post[[1]],
      aes(x = x, ymin = ymin, ymax = ymax),
      fill = "black",
      alpha = 0.75,
      inherit.aes = FALSE
    ) +
    geom_ribbon(
      data = ribbon_prior_post[[2]],
      aes(x = x, ymin = ymin, ymax = ymax),
      fill = c("#F28E2B", "#B1CAE1", "#F28E2B", "#B1CAE1")[ceiling(parameter / 6)],
      alpha = 0.6, #0.4,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = preds_prior_post[[2]],
      aes(y = y, group = quantile),
      linewidth = 0.7,
      color = c("#F28E2B", "#4E79A7", "#F28E2B", "#4E79A7")[ceiling(parameter / 6)]
    ) +
    geom_line(
      data = preds_prior_post[[1]],
      aes(y = y, group = quantile),
      linewidth = 0.5,
      color = 'black',
      linetype = 'dotted'
    ) +
    geom_text(
      data = labels_prior_post[[1]],
      aes(y = y, label = label),
      hjust = 1.2, size = 2, color = 'black',
      fontface = "bold"
    ) +
    geom_text(
      data = labels_prior_post[[2]],
      aes(y = y, label = label),
      hjust = -0.1, size = 2, color = c("#F28E2B", "#4E79A7", "#F28E2B", "#4E79A7")[ceiling(parameter / 6)],
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = 1.4, y = c(2.1, 2.1, 1.9, 1.9)[ceiling(parameter / 6)],
      label = title_char,
      hjust = 0.5, vjust = c(0.8, 0.5, 0.8, 0.5)[ceiling(parameter / 6)],
      size = c(4, 5, 4, 5)[ceiling(parameter / 6)]
    ) +
    scale_x_continuous(
      limits = c(-0.7, 3.3),
      breaks = c(0, 1, 2, 3),
      # labels = expression(10^0, 10^0.5, 10^1, 10^1.5, 10^2, 10^2.5, 10^3)
    ) +
    # ylim(-3, 2.5) +
    ylim(c(-2.5, -2.5, 0, 0)[ceiling(parameter / 6)], c(2.2, 2.2, 2, 2)[ceiling(parameter / 6)]) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill = "white")) +
    labs(title = NULL, x = NULL, y = NULL) #+ 
  # annotate(
  #   "text", 
  #   x = 0, y = 5.5, 
  #   label = paste0("Slope=", round(slope, 2), ", ",
  #                  "r=", round(r, 2)),
  #   hjust = 0,
  #   size = 3.5,
  #   color = 'white'
  # )
  
}

# arrange the panels
x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label= expression(log[10](Cophylogeny~size)), size = 5) + theme_void()
Y_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "z-score (absolute value)", size = 5, angle = 90) + theme_void()
y_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(ln(Relative~error)), size = 5, angle = 90) + theme_void()

plot_list_error <- list(

  A = err_plots[[1]],
  B = err_plots[[2]],
  C = err_plots[[3]],
  D = err_plots[[4]],
  E = err_plots[[5]],
  'F' = err_plots[[6]],
  
  a = err_plots[[7]],
  b = err_plots[[8]],
  c = err_plots[[9]],
  d = err_plots[[10]],
  e = err_plots[[11]],
  f = err_plots[[12]],
  
  G = err_plots[[13]],
  H = err_plots[[14]],
  I = err_plots[[15]],
  J = err_plots[[16]],
  K = err_plots[[17]],
  L = err_plots[[18]],
  
  g = err_plots[[19]],
  h = err_plots[[20]],
  i = err_plots[[21]],
  j = err_plots[[22]],
  k = err_plots[[23]],
  l = err_plots[[24]],
  
  x = x_lab,
  Y = Y_lab,
  y = y_lab
)















# run this after "plot_error_precision_ratio_cophylo_sizes.R" and "plot_accuracy_relative_error_cophylo_sizes.R"
plot_list_fig_4 <- list(
  
  'T' = ggplot() +
    annotate("text", x = 0, y = 0, label = "(a) Error", size = 5, fontface = "bold") +
    theme_void(),
  't' = ggplot() +
    annotate("text", x = 0, y = 0, label = "(b) Overconfidence", size = 5, fontface = "bold") +
    theme_void(),
  
  A = plot_list_error$A,
  B = plot_list_error$B,
  C = plot_list_error$C,
  D = plot_list_error$D,
  E = plot_list_error$E,
  'F' = plot_list_error$F,
  
  a = plot_list_error$a,
  b = plot_list_error$b,
  c = plot_list_error$c,
  d = plot_list_error$d,
  e = plot_list_error$e,
  f = plot_list_error$f,
  
  G = plot_list_error$G,
  H = plot_list_error$H,
  I = plot_list_error$I,
  J = plot_list_error$J,
  K = plot_list_error$K,
  L = plot_list_error$L,
  
  g = plot_list_error$g,
  h = plot_list_error$h,
  i = plot_list_error$i,
  j = plot_list_error$j,
  k = plot_list_error$k,
  l = plot_list_error$l,
  
  X = plot_list_error$x,
  x = plot_list_error$x,
  Y = plot_list_error$Y,
  y = plot_list_error$y
)

# layoutplot_fig_4 <- "
# yAaBbCc
# yDdEeFf
# YGgHhIi
# YJjKkLl
# xxxxxxx
# "

# layoutplot_fig_4 <- "
# yAa#DdYGg#Jj
# yBb#EeYHh#Kk
# yCc#FfYIi#Ll
# xxxxxxxxxxxx
# "


layoutplot_fig_4 <- "
#TTTTTTTTTTTT
yAABBCCDDEEFF
yAABBCCDDEEFF
yAABBCCDDEEFF
yAABBCCDDEEFF
yAABBCCDDEEFF
yAABBCCDDEEFF
yaabbccddeeff
yaabbccddeeff
yaabbccddeeff
yaabbccddeeff
yaabbccddeeff
yaabbccddeeff
#XXXXXXXXXXXX
#tttttttttttt
YGGHHIIJJKKLL
YGGHHIIJJKKLL
YGGHHIIJJKKLL
YGGHHIIJJKKLL
YGGHHIIJJKKLL
YGGHHIIJJKKLL
Ygghhiijjkkll
Ygghhiijjkkll
Ygghhiijjkkll
Ygghhiijjkkll
Ygghhiijjkkll
Ygghhiijjkkll
#xxxxxxxxxxxx
"

print('plotting')
#png(paste(prefix, 'cophy_ABC_convergence/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
pdf(paste('ex_cophy_ABC_convergence/', "fig_4_errors_ratios_regression_20260606", ".pdf", sep = ''), width = 11*0.9, height = 13*0.45*2, pointsize = 12)
print(wrap_plots(plot_list_fig_4, guides = 'collect', design = layoutplot_fig_4) &
        theme(legend.position = "left",
              legend.text = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
              legend.title = element_text(size = 9),
              legend.key.size = unit(0.5, "lines"),
              legend.spacing.y = unit(0.1, "cm")
        ))
dev.off()
