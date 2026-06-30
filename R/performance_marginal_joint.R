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

# note that "lambda_H" is the relative error of lambda_H, not the rate estimate

# here, we reverse the relative errors back to absolute errors
rel_err$lambda_H_abs <- rel_err$lambda_H * rel_err$lambda_H_true
rel_err$lambda_S_abs <- rel_err$lambda_S * rel_err$lambda_S_true
rel_err$lambda_C_abs <- rel_err$lambda_C * rel_err$lambda_C_true
rel_err$exp_H_abs <- rel_err$exp_H * rel_err$exp_H_true
rel_err$mu_H_frac_abs <- rel_err$mu_H_frac * rel_err$mu_H_frac_true
rel_err$mu_S_frac_abs <- rel_err$mu_S_frac * rel_err$mu_S_frac_true

# then, we reverse the absolute errors back to estimates
rel_err$lambda_H_est <- rel_err$lambda_H_abs + rel_err$lambda_H_true
rel_err$lambda_S_est <- rel_err$lambda_S_abs + rel_err$lambda_S_true
rel_err$lambda_C_est <- rel_err$lambda_C_abs + rel_err$lambda_C_true
rel_err$exp_H_est <- rel_err$exp_H_abs + rel_err$exp_H_true
rel_err$mu_H_frac_est <- rel_err$mu_H_frac_abs + rel_err$mu_H_frac_true
rel_err$mu_S_frac_est <- rel_err$mu_S_frac_abs + rel_err$mu_S_frac_true

# now we reparameterize
# estimates
rel_err$ratio_S_H_est <- rel_err$lambda_S_est / rel_err$lambda_H_est
rel_err$ratio_C_H_est <- rel_err$lambda_C_est / rel_err$lambda_H_est
rel_err$ratio_W_H_est <- rel_err$exp_H_est / rel_err$lambda_H_est # exp_H == lambda_W
rel_err$mu_H_est <- rel_err$mu_H_frac_est * (rel_err$lambda_H_est + rel_err$lambda_C_est)
rel_err$mu_S_est <- rel_err$mu_S_frac_est * (rel_err$lambda_S_est + rel_err$lambda_C_est + rel_err$exp_H_est)
# true values
rel_err$ratio_S_H_true <- rel_err$lambda_S_true / rel_err$lambda_H_true
rel_err$ratio_C_H_true <- rel_err$lambda_C_true / rel_err$lambda_H_true
rel_err$ratio_W_H_true <- rel_err$exp_H_true / rel_err$lambda_H_true # exp_H == lambda_W
rel_err$mu_H_true <- rel_err$mu_H_frac_true * (rel_err$lambda_H_true + rel_err$lambda_C_true)
rel_err$mu_S_true <- rel_err$mu_S_frac_true * (rel_err$lambda_S_true + rel_err$lambda_C_true + rel_err$exp_H_true)
# absolute errors
rel_err$ratio_S_H_abs <- rel_err$ratio_S_H_est - rel_err$ratio_S_H_true
rel_err$ratio_C_H_abs <- rel_err$ratio_C_H_est - rel_err$ratio_C_H_true
rel_err$ratio_W_H_abs <- rel_err$ratio_W_H_est - rel_err$ratio_W_H_true
rel_err$mu_H_abs <- rel_err$mu_H_est - rel_err$mu_H_true
rel_err$mu_S_abs <- rel_err$mu_S_est - rel_err$mu_S_true
# relative errors
rel_err$ratio_S_H <- rel_err$ratio_S_H_abs / rel_err$ratio_S_H_true
rel_err$ratio_C_H <- rel_err$ratio_C_H_abs / rel_err$ratio_C_H_true
rel_err$ratio_W_H <- rel_err$ratio_W_H_abs / rel_err$ratio_W_H_true
rel_err$mu_H <- rel_err$mu_H_abs / rel_err$mu_H_true
rel_err$mu_S <- rel_err$mu_S_abs / rel_err$mu_S_true

# averaging function for every X rows
aver <- function(vec, every){
  temp <- tapply(vec, (seq_along(vec)-1) %/% every, mean)
  return(temp)
}

# compute the effective dimensionality for every n rows
eff_dim_by_chunk <- function(df, n, method = "pearson") {
  
  # number of rows
  nr <- nrow(df)
  
  # split indices into chunks of size n
  idx <- split(seq_len(nr), ceiling(seq_len(nr) / n))
  
  # compute correlation matrix for each chunk
  cor_list <- lapply(idx, function(i) {
    cor(df[i, , drop = FALSE], method = method)
  })
  
  #  # based on Shannon entropy
  p_list <- lapply(cor_list, function(cor){
    eigen(cor)$values / sum(eigen(cor)$values)
  })

  n_dim_list <- lapply(p_list, function(p){
    prod( (p ^ (-p) ) )
  })
  
  # or
  # # based on quadratic entropy
  # eigen_list <- lapply(cor_list, function(cor){ # the set of eigenvalues (spectrum)
  #   eigen(cor)$values
  # })
  # n_dim_list <- lapply(eigen_list, function(eig_val){
  #   sum(eig_val)^2 / sum(eig_val^2)
  # })
  
  return(unlist(n_dim_list))
}


# "standard deviation"ing function for every X rows
stan_dev <- function(vec, every){
  temp <- tapply(vec, (seq_along(vec)-1) %/% every, sd)
  return(temp)
}

row1 <- as.vector(rel_err[1, 9:16])
every <- sum( apply(rel_err[, 9:16], 1, function(row) identical(as.numeric(row), as.numeric(row1))) )
# every <- 1 # this nullify the averaging

# transform the data frame for the 1st-level inference
logit <- function(p) { # a function for logit trnsformation
  log(p / (1 - p))
}
df_1st_post <- data.frame(
  ratio_S_H_est = log(rel_err$ratio_S_H_est),
  ratio_C_H_est = log(rel_err$ratio_C_H_est),
  ratio_W_H_est = log(rel_err$ratio_W_H_est),
  mu_H_frac_est = logit(rel_err$mu_H_frac_est),
  mu_S_frac_est = logit(rel_err$mu_S_frac_est)
)

df_2nd_post <- data.frame(
  lambda_H_est = log(rel_err$lambda_H_est),
  lambda_S_est = log(rel_err$lambda_S_est),
  lambda_C_est = log(rel_err$lambda_C_est),
  exp_H_est = log(rel_err$exp_H_est),
  mu_H_est = log(rel_err$mu_H_est),
  mu_S_est = log(rel_err$mu_S_est)
)


# for posterior
rel_err_aver_posterior <- data.frame(

    ratio_S_H = aver(rel_err$ratio_S_H, every),
    ratio_C_H = aver(rel_err$ratio_C_H, every),
    ratio_W_H = aver(rel_err$ratio_W_H, every),
    
    mu_H_frac = aver(rel_err$mu_H_frac, every),
    mu_S_frac = aver(rel_err$mu_S_frac, every),
    
    lambda_H = aver(rel_err$lambda_H, every),
    lambda_S = aver(rel_err$lambda_S, every),
    lambda_C = aver(rel_err$lambda_C, every),
    exp_H = aver(rel_err$exp_H, every),
    
    mu_H = aver(rel_err$mu_H, every),
    mu_S = aver(rel_err$mu_S, every),
    
    h_tree_size = aver(rel_err$h_tree_size, every),
    s_tree_size = aver(rel_err$s_tree_size, every),

    n_dim_1st = eff_dim_by_chunk(df_1st_post, n = every),
    n_dim_2nd = eff_dim_by_chunk(df_2nd_post, n = every)
)

# for prior
ratio_S_H_prior <- aver(rel_err$ratio_S_H_true, every)
ratio_C_H_prior <- aver(rel_err$ratio_C_H_true, every)
ratio_W_H_prior <- aver(rel_err$ratio_W_H_true, every)
mu_H_frac_prior <- aver(rel_err$mu_H_frac_true, every)
mu_S_frac_prior <- aver(rel_err$mu_S_frac_true, every)

lambda_H_prior <- aver(rel_err$lambda_H_true, every)
lambda_S_prior <- aver(rel_err$lambda_S_true, every)
lambda_C_prior <- aver(rel_err$lambda_C_true, every)
exp_H_prior <- aver(rel_err$exp_H_true, every)
mu_H_prior <- aver(rel_err$mu_H_true, every)
mu_S_prior <- aver(rel_err$mu_S_true, every)

# transform the data frame for the 1st_prior-level inference
df_1st_prior <- data.frame(
  ratio_S_H_true = aver(log(rel_err$ratio_S_H_true), every),
  ratio_C_H_true = aver(log(rel_err$ratio_C_H_true), every),
  ratio_W_H_true = aver(log(rel_err$ratio_W_H_true), every),
  mu_H_frac_true = aver(logit(rel_err$mu_H_frac_true), every),
  mu_S_frac_true = aver(logit(rel_err$mu_S_frac_true), every)
)

df_2nd_prior <- data.frame(
  lambda_H_true = aver(log(rel_err$lambda_H_true), every),
  lambda_S_true = aver(log(rel_err$lambda_S_true), every),
  lambda_C_true = aver(log(rel_err$lambda_C_true), every),
  exp_H_true = aver(log(rel_err$exp_H_true), every),
  mu_H_true = aver(log(rel_err$mu_H_true), every),
  mu_S_true = aver(log(rel_err$mu_S_true), every)
)

rel_err_aver_prior <- data.frame(
  ratio_S_H = (mean(ratio_S_H_prior) - ratio_S_H_prior) / ratio_S_H_prior,
  ratio_C_H = (mean(ratio_C_H_prior) - ratio_C_H_prior) / ratio_C_H_prior,
  ratio_W_H = (mean(ratio_W_H_prior) - ratio_W_H_prior) / ratio_W_H_prior,

  mu_H_frac = (mean(mu_H_frac_prior) - mu_H_frac_prior) / mu_H_frac_prior,
  mu_S_frac = (mean(mu_S_frac_prior) - mu_S_frac_prior) / mu_S_frac_prior,
  
  lambda_H = (mean(lambda_H_prior) - lambda_H_prior) / lambda_H_prior,
  lambda_S = (mean(lambda_S_prior) - lambda_S_prior) / lambda_S_prior,
  lambda_C = (mean(lambda_C_prior) - lambda_C_prior) / lambda_C_prior,
  exp_H = (mean(exp_H_prior) - exp_H_prior) / exp_H_prior,
  
  mu_H = (mean(mu_H_prior) - mu_H_prior) / mu_H_prior,
  mu_S = (mean(mu_S_prior) - mu_S_prior) / mu_S_prior,
  
  h_tree_size = aver(rel_err$h_tree_size, every),
  s_tree_size = aver(rel_err$s_tree_size, every),
  
  n_dim_1st = rep(eff_dim_by_chunk(df_1st_prior, n = nrow(df_1st_prior)), nrow(df_1st_prior)),
  n_dim_2nd = rep(eff_dim_by_chunk(df_2nd_prior, n = nrow(df_2nd_prior)), nrow(df_2nd_prior))
)

for (parameter in 1:11) {
  
  preds_prior_post <- NULL
  labels_prior_post <- NULL
  
  ribbon_prior_post <- NULL
  
  for (temp_id in 1:2) {
    
    rel_err_aver <- list(rel_err_aver_prior, rel_err_aver_posterior)[[temp_id ]]
    
    # the harmonic means of the two tree sizes
    x <- log10( 2* rel_err_aver$h_tree_size * rel_err_aver$s_tree_size / (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) )
    
    # relative errors
    y1 <- log(abs(rel_err_aver$ratio_S_H))
    y2 <- log(abs(rel_err_aver$ratio_C_H))
    y3 <- log(abs(rel_err_aver$ratio_W_H))
    
    y4 <- log(abs(rel_err_aver$mu_H_frac))
    y5 <- log(abs(rel_err_aver$mu_S_frac))
    
    y6 <- log(abs(rel_err_aver$lambda_H))
    y7 <- log(abs(rel_err_aver$lambda_S))
    y8 <- log(abs(rel_err_aver$lambda_C))
    y9 <- log(abs(rel_err_aver$exp_H))

    y10 <- log(abs(rel_err_aver$mu_H))
    y11 <- log(abs(rel_err_aver$mu_S))
    
    y <- rbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11)[parameter,]
    
    df <- data.frame(x=x, y=y)
    # mod <- lm(y ~ x, data = df)
    # df_pred <- cbind(df, predict(mod, interval = "prediction", level = 0.90))
    
    # qs <- c(0.2, 0.5, 0.80)
    qs <- c(0.2, 0.5, 0.80)
    
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
    
    expression(italic(frac(lambda[S], lambda[H]))),
    expression(italic(frac(lambda[C], lambda[H]))),
    expression(italic(frac(lambda[W], lambda[H]))),
    
    expression(italic(paste(frac(mu[H], lambda[H] + lambda[C])))),
    expression(italic(paste(frac(mu[S], lambda[S] + lambda[C] + lambda[W])))),
    
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
      fill = c("#B1CAE1", "#F28E2B")[(parameter > 5) + 1],
      alpha = 0.6, #0.4,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = preds_prior_post[[2]],
      aes(y = y, group = quantile),
      linewidth = 0.7,
      color = c("#4E79A7", "#F28E2B")[(parameter > 5) + 1]
    ) +
    geom_line(
      data = preds_prior_post[[1]],
      aes(y = y, group = quantile),
      linewidth = 0.5,
      color = 'black',
      linetype = 'dotted'
    ) +
    # geom_text(
    #   data = labels_prior_post[[1]],
    #   aes(y = y, label = label),
    #   hjust = 1.2, size = 3, color = 'black',
    #   fontface = "bold"
    # ) +
    # geom_text(
    #   data = labels_prior_post[[2]],
    #   aes(y = y, label = label),
    #   hjust = -0.1, size = 3, color = c("#4E79A7", "#F28E2B")[(parameter > 5) + 1],
    #   fontface = "bold"
    # ) +
    annotate(
      "text",
      x = 1.4, y = c(4, 2.8)[(parameter>5) + 1],
      label = title_char,
      hjust = 0.5, vjust = 0.5,
      size = c(4, 5)[(parameter>5) + 1]
    ) +
    annotate(
      "text",
      x = 1.4, y = 3.1,
      label = list(
        NULL,
        NULL,
        NULL,
        expression(paste(' ( = ', italic(epsilon[H]), ' )')),
        expression(paste(' ( = ', italic(epsilon[S]), ' )')),
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
      )[[parameter]],
      hjust = 0.5, vjust = 0.5,
      size = c(4, 5)[(parameter>5) + 1]
    )+
    scale_x_continuous(
      # limits = c(-0, 3),
      breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
      # labels = expression(10^0, 10^0.5, 10^1, 10^1.5, 10^2, 10^2.5, 10^3)
    ) +
    # ylim(-3, 2.5) +
    ylim(-2.7, 4.3 - (parameter > 5)) +
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
x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Cophylogeny size (log10 scale)", size = 5) + theme_void()
y_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Relative error (log scale)", size = 5, angle = 90) + theme_void()

plot_list_error <- list(

  B = err_plots[[1]],
  C = err_plots[[2]],
  D = err_plots[[3]],
  E = err_plots[[4]],
  'F' = err_plots[[5]],
  
  a = err_plots[[6]],
  b = err_plots[[7]],
  c = err_plots[[8]],
  d = err_plots[[9]],
  e = err_plots[[10]],
  f = err_plots[[11]],
  x = x_lab,
  y = y_lab
)



# run this after "plot_error_precision_ratio_cophylo_sizes.R" and "plot_accuracy_relative_error_cophylo_sizes.R"
plot_list_fig_4 <- list(
  
  B = plot_list_error$B,
  C = plot_list_error$C,
  D = plot_list_error$D,
  E = plot_list_error$E,
  'F' = plot_list_error$F,
  
  X = plot_list_error$x,
  Y = plot_list_error$y,
  
  a = plot_list_error$a,
  b = plot_list_error$b,
  c = plot_list_error$c,
  d = plot_list_error$d,
  e = plot_list_error$e,
  f = plot_list_error$f,
  
  x = plot_list_error$x,
  y = plot_list_error$y
)

layoutplot_fig_4 <- "
Y##BBCCDDEEFF
Y##BBCCDDEEFF
Y##BBCCDDEEFF
Y##BBCCDDEEFF
Y##BBCCDDEEFF
Y##BBCCDDEEFF
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
pdf(paste('ex_cophy_ABC_convergence/', "fig_4_errors_ratios_regression_20260426", ".pdf", sep = ''), width = 11*0.9, height = 9*0.9, pointsize = 12)
print(wrap_plots(plot_list_fig_4, guides = 'collect', design = layoutplot_fig_4) &
        theme(legend.position = "left",
              legend.text = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
              legend.title = element_text(size = 9),
              legend.key.size = unit(0.5, "lines"),
              legend.spacing.y = unit(0.1, "cm")
        ))
dev.off()






# here we plot the effective dimensionallity as a function of cophylogeny size
x <- log10( 2* rel_err_aver_posterior$h_tree_size * rel_err_aver_posterior$s_tree_size / (rel_err_aver_posterior$h_tree_size + rel_err_aver_posterior$s_tree_size) )
plot(x, rel_err_aver_posterior$n_dim_1st / rel_err_aver_prior$n_dim_1st, ylim = c(0.45, 1.2))
plot(x, rel_err_aver_posterior$n_dim_2nd / rel_err_aver_prior$n_dim_2nd, ylim = c(0.45, 1.2))


plot(x, rel_err_aver_posterior$n_dim_1st / 5, ylim = c(0.45, 1.2))
plot(x, rel_err_aver_posterior$n_dim_2nd / 6, ylim = c(0.45, 1.2))

plot(x, (rel_err_aver_posterior$n_dim_1st - rel_err_aver_prior$n_dim_1st) / 5, ylim = c(-0.5, 0.5))
plot(x, (rel_err_aver_posterior$n_dim_2nd - rel_err_aver_prior$n_dim_2nd) / 6, ylim = c(-0.5, 0.5))



# now we plot the pairwise correlations in the joint distribution of parameters
# a function for z-standardize in chunks
z_score <- function(v){
  (v - mean(v)) / sd(v)
}
z_by_n <- function(x, n) {
  split_x <- split(x, ceiling(seq_along(x) / n))
  unlist(lapply(split_x, z_score ))
}

# a function for logit trnsformation
logit <- function(p) {
  log(p / (1 - p))
}

# for 1st-level inference
joint_z_1st_post <- data.frame(
  ratio_S_H = z_by_n(log(rel_err$ratio_S_H_est), every),
  ratio_C_H = z_by_n(log(rel_err$ratio_C_H_est), every),
  ratio_W_H = z_by_n(log(rel_err$ratio_W_H_est), every),
  mu_H_frac = z_by_n(logit(rel_err$mu_H_frac_est), every),
  mu_S_frac = z_by_n(logit(rel_err$mu_S_frac_est), every)
)
joint_z_1st_post$group <- rep("Posterior", nrow(joint_z_1st_post))

joint_z_1st_prior <- data.frame(
  ratio_S_H = z_score(aver(log(rel_err$ratio_S_H_true), every)),
  ratio_C_H = z_score(aver(log(rel_err$ratio_C_H_true), every)),
  ratio_W_H = z_score(aver(log(rel_err$ratio_W_H_true), every)),
  mu_H_frac = z_score(aver(logit(rel_err$mu_H_frac_true), every)),
  mu_S_frac = z_score(aver(logit(rel_err$mu_S_frac_true), every))
)
joint_z_1st_prior$group <- rep("Prior", nrow(joint_z_1st_prior))


GGally::ggpairs(
  rbind(
    joint_z_1st_post[1:1000,],
    joint_z_1st_prior
    ),
  columns = 1:5,
  aes(color = group),
  # lower = list(continuous = "points"),
  lower = list(continuous = "density"),
  diag = list(continuous = "blank")
  )





# for 2nd-level inference
joint_z_2nd_post <- data.frame(
  lambda_H = z_by_n(log(rel_err$lambda_H_est), every),
  lambda_S = z_by_n(log(rel_err$lambda_S_est), every),
  lambda_C = z_by_n(log(rel_err$lambda_C_est), every),
  lambda_W = z_by_n(log(rel_err$exp_H_est), every),
  mu_H = z_by_n(log(rel_err$mu_H_est), every),
  mu_S = z_by_n(log(rel_err$mu_S_est), every)
)
joint_z_2nd_post$group <- rep("Posterior", nrow(joint_z_2nd_post))

joint_z_2nd_prior <- data.frame(
  lambda_H = z_score(aver(log(rel_err$lambda_H_true), every)),
  lambda_S = z_score(aver(log(rel_err$lambda_S_true), every)),
  lambda_C = z_score(aver(log(rel_err$lambda_C_true), every)),
  lambda_W = z_score(aver(log(rel_err$exp_H_true), every)),
  mu_H = z_score(aver(log(rel_err$mu_H_true), every)),
  mu_S = z_score(aver(log(rel_err$mu_S_true), every))
)
joint_z_2nd_prior$group <- rep("Prior", nrow(joint_z_2nd_prior))


GGally::ggpairs(
  rbind(
    joint_z_2nd_post[1:1000,],
    joint_z_2nd_prior
  ),
  columns = 1:6,
  aes(color = group),
  # lower = list(continuous = "points"),
  lower = list(continuous = "density"),
  diag = list(continuous = "blank")
)
