library(ggpubr)
library(patchwork)
library(viridis)
library(dplyr)
library(quantreg)
library(splines)


setwd("/groups/cromanpa/yzeng")

err_plots <- NULL
ratio_plots <- NULL

rel_err_list <- NULL

for (ext_scen in 1:5) {
  
  err_plots_ext_scen <- NULL
  ratio_plots_ext_scen <- NULL
  
  mu_H_frac <- c(0, 0.3, 0.7, 0,   0.7)[ext_scen]
  mu_S_frac <- c(0, 0.3, 0,   0.7, 0.7)[ext_scen]
  
  file_name <- paste("ex_", 'cophy_ABC_convergence/cross_validation_sizes',"_", as.character(mu_H_frac), "_", as.character(mu_S_frac), '.csv', sep = '')
  rel_err <- read.csv(file_name)
  
  # here, we reverse the relative errors back to absolute errors
  rel_err$lambda_H_abs <- rel_err$lambda_H * rel_err$lambda_H_true
  rel_err$lambda_S_abs <- rel_err$lambda_S * rel_err$lambda_S_true
  rel_err$lambda_C_abs <- rel_err$lambda_C * rel_err$lambda_C_true
  rel_err$exp_H_abs <- rel_err$exp_H * rel_err$exp_H_true
  
  rel_err_list[[length(rel_err_list) + 1]] <- rel_err #save the data frame for relative errors
  
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
  
  row1 <- as.vector(rel_err[1,7:11])
  every <- sum( apply(rel_err[,7:11], 1, function(row) identical(as.numeric(row), as.numeric(row1))) )
  # every <- 1 # this nullify the averaging
  
  # for posterior
  rel_err_aver_posterior <- data.frame(
    lambda_H_abs_aver = aver(rel_err$lambda_H_abs, every),
    lambda_S_abs_aver = aver(rel_err$lambda_S_abs, every),
    lambda_C_abs_aver = aver(rel_err$lambda_C_abs, every),
    exp_H_abs_aver = aver(rel_err$exp_H_abs, every),
    
    lambda_H_abs_stan_dev = stan_dev(rel_err$lambda_H_abs, every),
    lambda_S_abs_stan_dev = stan_dev(rel_err$lambda_S_abs, every),
    lambda_C_abs_stan_dev = stan_dev(rel_err$lambda_C_abs, every),
    exp_H_abs_stan_dev = stan_dev(rel_err$exp_H_abs, every),
    
    h_tree_size = aver(rel_err$h_tree_size, every),
    s_tree_size = aver(rel_err$s_tree_size, every)
  )
  
  # for prior
  lambda_H_prior <- aver(rel_err$lambda_H_true, every)
  lambda_S_prior <- aver(rel_err$lambda_S_true, every)
  lambda_C_prior <- aver(rel_err$lambda_C_true, every)
  exp_H_prior <- aver(rel_err$exp_H_true, every)

  rel_err_aver_prior <- data.frame(
    lambda_H_abs_aver = mean(lambda_H_prior) - lambda_H_prior,
    lambda_S_abs_aver = mean(lambda_S_prior) - lambda_S_prior,
    lambda_C_abs_aver = mean(lambda_C_prior) - lambda_C_prior,
    exp_H_abs_aver = mean(exp_H_prior) - exp_H_prior,

    lambda_H_abs_stan_dev = sd(lambda_H_prior),
    lambda_S_abs_stan_dev = sd(lambda_S_prior),
    lambda_C_abs_stan_dev = sd(lambda_C_prior),
    exp_H_abs_stan_dev = sd(exp_H_prior),

    h_tree_size = aver(rel_err$h_tree_size, every),
    s_tree_size = aver(rel_err$s_tree_size, every)
  )
  
  for (lambda in 1:4) {
    
    preds_prior_post <- NULL
    labels_prior_post <- NULL
    
    for (temp_id in 1:2) {
      
      rel_err_aver <- list(rel_err_aver_prior, rel_err_aver_posterior)[[temp_id ]]
      
      # the harmonic means of the two tree sizes
      x <- log10( 2* rel_err_aver$h_tree_size * rel_err_aver$s_tree_size / (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) )
      
      y1 <- log( abs(rel_err_aver$lambda_H_abs_aver) / rel_err_aver$lambda_H_abs_stan_dev )
      y2 <- log( abs(rel_err_aver$lambda_S_abs_aver) / rel_err_aver$lambda_S_abs_stan_dev )
      y3 <- log( abs(rel_err_aver$lambda_C_abs_aver) / rel_err_aver$lambda_C_abs_stan_dev )
      y4 <- log( abs(rel_err_aver$exp_H_abs_aver) / rel_err_aver$exp_H_abs_stan_dev )
      
      y <- rbind(y1, y2, y3, y4)[lambda,]
      df <- data.frame(x=x, y=y)
      # mod <- lm(y ~ x, data = df)
      # df_pred <- cbind(df, predict(mod, interval = "prediction", level = 0.90))
      
      qs <- c(0.1, 0.5, 0.9)
      
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
      preds <- preds[preds$n_local >= 150, ]
      preds_prior_post[[length(preds_prior_post) + 1]] <- preds
      
      # Extract end points for labeling
      labels <- preds %>%
        group_by(quantile) %>%
        filter(x == range(x)[temp_id]) %>%
        mutate(label = paste0(quantile * 100, "%"))
      labels_prior_post[[length(labels_prior_post) + 1]] <- labels
      
    }
    
    err_plots_ext_scen[[length(err_plots_ext_scen) + 1]] <- ggplot(df, aes(x = x, y = y)) +
      # geom_point(size = 1, alpha = 1, color = point_color) +
      # geom_hex(bins = 20) +
      # geom_point(size = 1.8, alpha = 0.5)
      # scale_fill_viridis(option="viridis") +
      # scale_fill_viridis(option="viridis", limits = c(0, 22)) +
      # geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = 0.5) +
      geom_line(
        data = preds_prior_post[[1]],
        aes(y = y, group = quantile),
        linewidth = 0.5,
        linetype = "dashed",
        color = 'black'
      ) +
      geom_line(
        data = preds_prior_post[[2]],
        aes(y = y, group = quantile),
        linewidth = 0.9,
        color = 'darkorange'
      ) +
      geom_text(
        data = labels_prior_post[[1]],
        aes(y = y, label = label),
        hjust = 1.2, size = 3, color = 'black',
        fontface = "bold"
      ) +
      geom_text(
        data = labels_prior_post[[2]],
        aes(y = y, label = label),
        hjust = -0.1, size = 3, color = 'darkorange1',
        fontface = "bold"
      ) +
      xlim(-0.5, 3.3) +
      ylim(-3, 1.2) +
      # ylim(0, 3) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(panel.background = element_rect(fill = "white")) +
      labs(x = NULL, y = NULL)
    
  }
  
  err_plots[[length(err_plots) + 1]] <- err_plots_ext_scen
  ratio_plots[[length(ratio_plots) + 1]] <- ratio_plots_ext_scen
  
}

# arrange the panels
x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Cophylogeny size", size = 5) + theme_void()
y_lab_err <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "z-score (unsigned)", size = 5, angle = 90) + theme_void()
y_lab_acc <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Accuracy (%)", size = 5, angle = 90) + theme_void()

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0), paste(epsilon[S], " = ", 0))), angle = 0) + theme_void()
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0.3), paste(epsilon[S], " = ", 0.3))), angle = 0) + theme_void()
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0), paste(epsilon[S], " = ", 0.7))), angle = 0) + theme_void()
col4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0.7), paste(epsilon[S], " = ", 0))), angle = 0) + theme_void()
col5 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0.7), paste(epsilon[S], " = ", 0.7))), angle = 0) + theme_void()

row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(hat(lambda)[H]), size = 5, angle = 0) + theme_void()
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(hat(lambda)[S]), size = 5, angle = 0) + theme_void()
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(hat(lambda)[C]), size = 5, angle = 0) + theme_void()
row4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(hat(lambda)[W]), size = 5, angle = 0) + theme_void()

# without row and column labels

# portrait arrangement
# layoutplot <- "
# #qqrrssttuu#
# wAAEEIIMMQQm
# wAAEEIIMMQQm
# wAAEEIIMMQQm
# wBBFFJJNNRRn
# wBBFFJJNNRRn
# wBBFFJJNNRRn
# wCCGGKKOOSSo
# wCCGGKKOOSSo
# wCCGGKKOOSSo
# wDDHHLLPPTTp
# wDDHHLLPPTTp
# wDDHHLLPPTTp
# #vvvvvvvvvv#
# "

# landscape arrangement
layoutplot <- "
#mmmnnnoooppp##
wAAABBBCCCDDDqq
wAAABBBCCCDDDqq
wEEEFFFGGGHHHrr
wEEEFFFGGGHHHrr
wIIIJJJKKKLLLss
wIIIJJJKKKLLLss
wMMMNNNOOOPPPtt
wMMMNNNOOOPPPtt
wQQQRRRSSSTTTuu
wQQQRRRSSSTTTuu
#vvvvvvvvvvvv##
"

panels <- err_plots
y_lab <- y_lab_err

plot_list_ratio <- list(
  A = panels[[1]][[1]],
  B = panels[[1]][[2]],
  C = panels[[1]][[3]],
  D = panels[[1]][[4]],
  E = panels[[2]][[1]],
  'F' = panels[[2]][[2]],
  G = panels[[2]][[3]],
  H = panels[[2]][[4]],
  I = panels[[4]][[1]],
  J = panels[[4]][[2]],
  K = panels[[4]][[3]],
  L = panels[[4]][[4]],
  M = panels[[3]][[1]],
  N = panels[[3]][[2]],
  O = panels[[3]][[3]],
  P = panels[[3]][[4]],
  Q = panels[[5]][[1]],
  R = panels[[5]][[2]],
  S = panels[[5]][[3]],
  'T' = panels[[5]][[4]],
  m = row1,
  n = row2,
  o = row3,
  p = row4,
  q = col1,
  r = col2,
  s = col3,
  t = col4,
  u = col5,
  v = x_lab,
  w = y_lab
)












library(ggpubr)
library(patchwork)
library(viridis)
library(dplyr)
library(quantreg)
library(splines)

setwd("/groups/cromanpa/yzeng")

err_plots <- NULL

rel_err_list <- NULL

for (ext_scen in 1:5) {
  
  err_plots_ext_scen <- NULL
  
  mu_H_frac <- c(0, 0.3, 0.7, 0,   0.7)[ext_scen]
  mu_S_frac <- c(0, 0.3, 0,   0.7, 0.7)[ext_scen]
  
  file_name <- paste("ex_", 'cophy_ABC_convergence/cross_validation_sizes',"_", as.character(mu_H_frac), "_", as.character(mu_S_frac), '.csv', sep = '')
  rel_err <- read.csv(file_name)
  
  rel_err_list[[length(rel_err_list) + 1]] <- rel_err #save the data frame for relative errors
  
  # averaging function for every X rows
  aver <- function(vec, every){
    temp <- tapply(vec, (seq_along(vec)-1) %/% every, mean)
    return(temp)
  }
  
  row1 <- as.vector(rel_err[1,7:11])
  every <- sum( apply(rel_err[,7:11], 1, function(row) identical(as.numeric(row), as.numeric(row1))) )
  # every <- 1 # this nullify the averaging
  
  # for posterior
  rel_err_aver_posterior <- data.frame(
    lambda_H = aver(rel_err$lambda_H, every),
    lambda_S = aver(rel_err$lambda_S, every),
    lambda_C = aver(rel_err$lambda_C, every),
    exp_H = aver(rel_err$exp_H, every),
    h_tree_size = aver(rel_err$h_tree_size, every),
    s_tree_size = aver(rel_err$s_tree_size, every)
  )
  
  # for prior
  lambda_H_prior <- aver(rel_err$lambda_H_true, every)
  lambda_S_prior <- aver(rel_err$lambda_S_true, every)
  lambda_C_prior <- aver(rel_err$lambda_C_true, every)
  exp_H_prior <- aver(rel_err$exp_H_true, every)

  rel_err_aver_prior <- data.frame(
    lambda_H = (mean(lambda_H_prior) - lambda_H_prior) / lambda_H_prior,
    lambda_S = (mean(lambda_S_prior) - lambda_S_prior) / lambda_S_prior,
    lambda_C = (mean(lambda_C_prior) - lambda_C_prior) / lambda_C_prior,
    exp_H = (mean(exp_H_prior) - exp_H_prior) / exp_H_prior,

    h_tree_size = aver(rel_err$h_tree_size, every),
    s_tree_size = aver(rel_err$s_tree_size, every)
  )
  
  for (lambda in 1:4) {
    
    preds_prior_post <- NULL
    labels_prior_post <- NULL
    
    for (temp_id in 1:2) {
      
      rel_err_aver <- list(rel_err_aver_prior, rel_err_aver_posterior)[[temp_id ]]
      
      # the harmonic means of the two tree sizes
      x <- log10( 2* rel_err_aver$h_tree_size * rel_err_aver$s_tree_size / (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) )
      
      y1 <- log(abs(rel_err_aver$lambda_H))
      y2 <- log(abs(rel_err_aver$lambda_S))
      y3 <- log(abs(rel_err_aver$lambda_C))
      y4 <- log(abs(rel_err_aver$exp_H))
      
      y <- rbind(y1, y2, y3, y4)[lambda,]
      
      df <- data.frame(x=x, y=y)
      # mod <- lm(y ~ x, data = df)
      # df_pred <- cbind(df, predict(mod, interval = "prediction", level = 0.90))
      
      qs <- c(0.1, 0.5, 0.9)
      
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
      preds <- preds[preds$n_local >= 150, ]
      preds_prior_post[[length(preds_prior_post) + 1]] <- preds
      
      # Extract end points for labeling
      labels <- preds %>%
        group_by(quantile) %>%
        filter(x == range(x)[temp_id]) %>%
        mutate(label = paste0(quantile * 100, "%"))
      labels_prior_post[[length(labels_prior_post) + 1]] <- labels
      
    }
    
    err_plots_ext_scen[[length(err_plots_ext_scen) + 1]] <- ggplot(df, aes(x = x, y = y)) +
      # geom_hex(bins = 20) +
      # geom_point(size = 1.8, alpha = 0.5)
      # geom_smooth(data = subset(df, x >= quantile(df$x, 0) & x <= quantile(df$x, 1)), method = "loess", color = "white", linewidth = 0.9, se = FALSE) +
      # scale_fill_viridis(option="magma") +
      # scale_fill_viridis(option="magma", limits = c(0, 32)) +
      # geom_ribbon(aes(ymin = lwr, ymax = upr),
      #             fill = "chocolate", alpha = 0.2) +
      # geom_point(aes(y = y)) +
      # geom_line(aes(y = fit), colour = "chocolate", size = 1) +
      #geom_smooth(aes(y = y), method = "loess", se = TRUE, level = 0.95, color = "chocolate", fill = "lightchocolate") +
      geom_line(
        data = preds_prior_post[[1]],
        aes(y = y, group = quantile),
        linewidth = 0.5,
        color = 'black',
        linetype = 'dashed'
      ) +
      geom_line(
        data = preds_prior_post[[2]],
        aes(y = y, group = quantile),
        linewidth = 0.9,
        color = 'blue'
      ) +
      geom_text(
        data = labels_prior_post[[1]],
        aes(y = y, label = label),
        hjust = 1.2, size = 3, color = 'black',
        fontface = "bold"
      ) +
      geom_text(
        data = labels_prior_post[[2]],
        aes(y = y, label = label),
        hjust = -0.1, size = 3, color = 'blue',
        fontface = "bold"
      ) +
      xlim(-0.5, 3.3) +
      ylim(-4, 3.3) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(panel.background = element_rect(fill = "white")) +
      labs(x = NULL, y = NULL) #+ 
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
  
  err_plots[[length(err_plots) + 1]] <- err_plots_ext_scen
}

# arrange the panels
x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Cophylogeny size", size = 5) + theme_void()
y_lab_err <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Relative error", size = 5, angle = 90) + theme_void()
y_lab_acc <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Accuracy (%)", size = 5, angle = 90) + theme_void()

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0), paste(epsilon[S], " = ", 0))), angle = 0) + theme_void()
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0.3), paste(epsilon[S], " = ", 0.3))), angle = 0) + theme_void()
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0), paste(epsilon[S], " = ", 0.7))), angle = 0) + theme_void()
col4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0.7), paste(epsilon[S], " = ", 0))), angle = 0) + theme_void()
col5 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(atop(paste(epsilon[H], " = ", 0.7), paste(epsilon[S], " = ", 0.7))), angle = 0) + theme_void()

row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(hat(lambda)[H]), size = 5, angle = 0) + theme_void()
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(hat(lambda)[S]), size = 5, angle = 0) + theme_void()
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(hat(lambda)[C]), size = 5, angle = 0) + theme_void()
row4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label = expression(hat(lambda)[W]), size = 5, angle = 0) + theme_void()

# without row and column labels

# portrait arrangement
# layoutplot <- "
# #qqrrssttuu#
# wAAEEIIMMQQm
# wAAEEIIMMQQm
# wAAEEIIMMQQm
# wBBFFJJNNRRn
# wBBFFJJNNRRn
# wBBFFJJNNRRn
# wCCGGKKOOSSo
# wCCGGKKOOSSo
# wCCGGKKOOSSo
# wDDHHLLPPTTp
# wDDHHLLPPTTp
# wDDHHLLPPTTp
# #vvvvvvvvvv#
# "

# landscape arrangement
layoutplot <- "
#mmmnnnoooppp##
wAAABBBCCCDDDqq
wAAABBBCCCDDDqq
wEEEFFFGGGHHHrr
wEEEFFFGGGHHHrr
wIIIJJJKKKLLLss
wIIIJJJKKKLLLss
wMMMNNNOOOPPPtt
wMMMNNNOOOPPPtt
wQQQRRRSSSTTTuu
wQQQRRRSSSTTTuu
#vvvvvvvvvvvv##
"

panels <- err_plots
y_lab <- y_lab_err

plot_list_error <- list(
  A = panels[[1]][[1]],
  B = panels[[1]][[2]],
  C = panels[[1]][[3]],
  D = panels[[1]][[4]],
  E = panels[[2]][[1]],
  'F' = panels[[2]][[2]],
  G = panels[[2]][[3]],
  H = panels[[2]][[4]],
  I = panels[[4]][[1]],
  J = panels[[4]][[2]],
  K = panels[[4]][[3]],
  L = panels[[4]][[4]],
  M = panels[[3]][[1]],
  N = panels[[3]][[2]],
  O = panels[[3]][[3]],
  P = panels[[3]][[4]],
  Q = panels[[5]][[1]],
  R = panels[[5]][[2]],
  S = panels[[5]][[3]],
  'T' = panels[[5]][[4]],
  m = row1,
  n = row2,
  o = row3,
  p = row4,
  q = col1,
  r = col2,
  s = col3,
  t = col4,
  u = col5,
  v = x_lab,
  w = y_lab
)












# run this after "plot_error_precision_ratio_cophylo_sizes.R" and "plot_accuracy_relative_error_cophylo_sizes.R"
plot_list_fig_4 <- list(
  V = plot_list_error$m,
  X = plot_list_error$n,
  Y = plot_list_error$o,
  Z= plot_list_error$p,
  
  w = plot_list_error$w,
  
  e = plot_list_error$E,
  f = plot_list_error$'F',
  g = plot_list_error$G,
  h = plot_list_error$H,
  r = plot_list_error$r,
  
  i = plot_list_error$I,
  j = plot_list_error$J,
  k = plot_list_error$K,
  l = plot_list_error$L,
  s = plot_list_error$s,
  
  m = plot_list_error$M,
  n = plot_list_error$N,
  o = plot_list_error$O,
  p = plot_list_error$P,
  t = plot_list_error$t,
  
  W = plot_list_ratio$w,
  
  E = plot_list_ratio$E,
  'F' = plot_list_ratio$'F',
  G = plot_list_ratio$G,
  H = plot_list_ratio$H,
  R = plot_list_ratio$r,
  
  I = plot_list_ratio$I,
  J = plot_list_ratio$J,
  K = plot_list_ratio$K,
  L = plot_list_ratio$L,
  S = plot_list_ratio$s,
  
  M = plot_list_ratio$M,
  N = plot_list_ratio$N,
  O = plot_list_ratio$O,
  P = plot_list_ratio$P,
  'T' = plot_list_ratio$t,
  
  v = plot_list_error$v
)

layoutplot_fig_4 <- "
#VVVXXXYYYZZZ##
weeefffggghhhrr
weeefffggghhhrr
weeefffggghhhrr
wiiijjjkkklllss
wiiijjjkkklllss
wiiijjjkkklllss
wmmmnnnoooppptt
wmmmnnnoooppptt
wmmmnnnoooppptt
###############
WEEEFFFGGGHHHRR
WEEEFFFGGGHHHRR
WEEEFFFGGGHHHRR
WIIIJJJKKKLLLSS
WIIIJJJKKKLLLSS
WIIIJJJKKKLLLSS
WMMMNNNOOOPPPTT
WMMMNNNOOOPPPTT
WMMMNNNOOOPPPTT
#vvvvvvvvvvvv##
"

print('plotting')
#png(paste(prefix, 'cophy_ABC_convergence/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
pdf(paste('ex_cophy_ABC_convergence/', "fig_4_errors_ratios", ".pdf", sep = ''), width = 10, height = 12, pointsize = 12)
print(wrap_plots(plot_list_fig_4, guides = 'collect', design = layoutplot_fig_4) &
        theme(legend.position = "left",
              legend.text = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
              legend.title = element_text(size = 9),
              legend.key.size = unit(0.5, "lines"),
              legend.spacing.y = unit(0.1, "cm")
        ))
dev.off()