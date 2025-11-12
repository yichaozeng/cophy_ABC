library(ggpubr)
library(patchwork)
library(viridis)
library(dplyr)
library(quantreg)
library(splines)

setwd("/groups/cromanpa/yzeng")

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

err_plots <- NULL
acc_plots <- NULL

rel_err_list <- NULL

for (ext_scen in 1:5) {
  
  err_plots_ext_scen <- NULL
  acc_plots_ext_scen <- NULL
  
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
  rel_err_aver <- data.frame(
    lambda_H = aver(rel_err$lambda_H, every),
    lambda_S = aver(rel_err$lambda_S, every),
    lambda_C = aver(rel_err$lambda_C, every),
    exp_H = aver(rel_err$exp_H, every),
    h_tree_size = aver(rel_err$h_tree_size, every),
    s_tree_size = aver(rel_err$s_tree_size, every)#,
    # net_size = aver(rel_err$net_size, every),
    # blend_1_size = aver(rel_err$blend_1_size, every),
    # blend_2_size = aver(rel_err$blend_2_size, every)
  )
  
  # # for prior
  # lambda_H_prior <- aver(rel_err$lambda_H_true, every)
  # lambda_S_prior <- aver(rel_err$lambda_S_true, every)
  # lambda_C_prior <- aver(rel_err$lambda_C_true, every)
  # exp_H_prior <- aver(rel_err$exp_H_true, every)
  # 
  # rel_err_aver <- data.frame(
  #   lambda_H = (mean(lambda_H_prior) - lambda_H_prior) / lambda_H_prior,
  #   lambda_S = (mean(lambda_S_prior) - lambda_S_prior) / lambda_S_prior,
  #   lambda_C = (mean(lambda_C_prior) - lambda_C_prior) / lambda_C_prior,
  #   exp_H = (mean(exp_H_prior) - exp_H_prior) / exp_H_prior,
  # 
  #   h_tree_size = aver(rel_err$h_tree_size, every),
  #   s_tree_size = aver(rel_err$s_tree_size, every)#,
  #   # net_size = aver(rel_err$net_size, every),
  #   # blend_1_size = aver(rel_err$blend_1_size, every),
  #   # blend_2_size = aver(rel_err$blend_2_size, every)
  # )
  
  # x <- log10(rel_err_aver$h_tree_size * (rel_err_aver$h_tree_size <= rel_err_aver$s_tree_size) + rel_err_aver$s_tree_size * (rel_err_aver$h_tree_size > rel_err_aver$s_tree_size)) # the smaller of two tree sizes
  
  # or the harmonic means of the two tree sizes
  x <- log10( 2* rel_err_aver$h_tree_size * rel_err_aver$s_tree_size / (rel_err_aver$h_tree_size + rel_err_aver$s_tree_size) ) # the smaller of two tree sizes
  
  # x <- log10(rel_err_aver$h_tree_size)
  # x <- log10(rel_err_aver$s_tree_size)
  # x <- log10(rel_err_aver$net_size)
  # x <- log10(rel_err_aver$blend_1_size)
  # x <- log10(rel_err_aver$blend_2_size)
  
  y1 <- log(abs(rel_err_aver$lambda_H))
  y2 <- log(abs(rel_err_aver$lambda_S))
  y3 <- log(abs(rel_err_aver$lambda_C))
  y4 <- log(abs(rel_err_aver$exp_H))
  
  # y1 <- abs(rel_err_aver$lambda_H)
  # y2 <- abs(rel_err_aver$lambda_S)
  # y3 <- abs(rel_err_aver$lambda_C)
  # y4 <- abs(rel_err_aver$exp_H)
  
  # y1 <- rel_err_aver$lambda_H
  # y2 <- rel_err_aver$lambda_S
  # y3 <- rel_err_aver$lambda_C
  # y4 <- rel_err_aver$exp_H
  
  for (lambda in 1:4) {
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
    
    # Extract end points for labeling
    labels <- preds %>%
      group_by(quantile) %>%
      filter(x == max(x)) %>%
      mutate(label = paste0(quantile * 100, "%"))
    
    err_plots_ext_scen[[length(err_plots_ext_scen) + 1]] <- ggplot(df, aes(x = x, y = y)) +
      geom_hex(bins = 20) +
      # geom_point(size = 1.8, alpha = 0.5)
      # geom_smooth(data = subset(df, x >= quantile(df$x, 0) & x <= quantile(df$x, 1)), method = "loess", color = "white", linewidth = 0.7, se = FALSE) +
      # scale_fill_viridis(option="magma") +
      scale_fill_viridis(option="magma", limits = c(0, 32)) +
      # geom_ribbon(aes(ymin = lwr, ymax = upr),
      #             fill = "blue", alpha = 0.2) +
      # geom_point(aes(y = y)) +
      # geom_line(aes(y = fit), colour = "blue", size = 1) +
      #geom_smooth(aes(y = y), method = "loess", se = TRUE, level = 0.95, color = "blue", fill = "lightblue") +
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
      ylim(-8, 6) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(panel.background = element_rect(fill = "grey75")) +
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
    
    acc_thres <- 0.3 # the relative error threshold
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
    hist1 <- hist(x, breaks = breaks, plot = F)
    hist2 <- hist(x[abs(y) <= acc_thres], breaks = breaks, plot = F)
    
    size = NULL
    for (size_class in 1:(length(breaks) - 1)) {
      size <- c(size, paste(breaks[size_class], '-', breaks[size_class+1]))
    }
    df <- data.frame(
      size = size,
      accuracy = hist2$counts / hist1$counts * 100
    )
    df$size <- factor(df$size, levels = size)
    
    # Bar plot with y-axis range (0, 100)
    acc_plots_ext_scen[[length(acc_plots_ext_scen) + 1]] <- ggplot(df, aes(x = size, y = accuracy)) +
      geom_bar(stat = "identity", fill = cbp1[lambda]) +
      scale_y_continuous(limits = c(0, 100)) +
      # theme_minimal()
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(x = NULL, y = NULL)
    
  }
  
  err_plots[[length(err_plots) + 1]] <- err_plots_ext_scen
  acc_plots[[length(acc_plots) + 1]] <- acc_plots_ext_scen

}

# arrange the panels
x_lab <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Cophylogeny size", size = 5) + theme_void()
y_lab_err <- ggplot() + annotate(geom = 'text', x=1, y=1, label = "Posterior relative error", size = 5, angle = 90) + theme_void()
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

for (fig_id in 2:1) {
  
  panels <- list(err_plots, acc_plots)[[fig_id]]
  y_lab <- list(y_lab_err, y_lab_acc)[[fig_id]]

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
  
  
  print('plotting')
  #png(paste(prefix, 'cophy_ABC_convergence/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
  pdf(paste('ex_cophy_ABC_convergence/', c("errors_dont", 'accuracy')[fig_id], ".pdf", sep = ''), width = 12, height = 12, pointsize = 12)
  print(wrap_plots(plot_list_error, guides = 'keep', design = layoutplot) &
          theme(legend.position = "right",
            legend.text = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
            legend.title = element_text(size = 9),
            legend.key.size = unit(0.5, "lines"),
            legend.spacing.y = unit(0.1, "cm")
          ))
  dev.off()
  
}





# # here we assign the signed relative errors to each true parameter value
# for (ext_scen in 1:5) {
# 
#   print(ext_scen)
# 
#   mu_H_frac <- c(0, 0.3, 0.7, 0,   0.7)[ext_scen]
#   mu_S_frac <- c(0, 0.3, 0,   0.7, 0.7)[ext_scen]
# 
#   surv_H <- 1 - mu_H_frac
#   surv_S <- 1 - mu_S_frac
# 
#   df_true_para <- NULL
# 
#   # re-create the simulated parameters
#   for(lambda_H_real in c(1.9, 1.45, 1, 0.55, 0.1)){
#     for (lambda_S_real in c(1.9, 1.45, 1, 0.55, 0.1)) {
#       for (lambda_C_real in c(1.9, 1.45, 1, 0.55, 0.1)){
#         for (exp_H_real in c(1.9, 1.45, 1, 0.55, 0.1)){
# 
#           # here we set an upper limit for the sum of all four lambdas
#           if(lambda_H_real*surv_H + lambda_S_real*surv_S + lambda_C_real*(surv_H+surv_S)/2 + exp_H_real*surv_S <= 5){
# 
#             temp <- rel_err_list[[ext_scen]]
#             row1 <- as.vector(temp[1,7:11])
#             every <- sum( apply(temp[,7:11], 1, function(row) identical(as.numeric(row), as.numeric(row1))) )
# 
#             df_true_para <- rbind(
#               df_true_para,
#               cbind(
#                 rep(lambda_H_real, every),
#                 rep(lambda_S_real, every),
#                 rep(lambda_C_real, every),
#                 rep(exp_H_real, every)
#               )
#             )
# 
#           }
# 
#         }
#       }
#     }
#   }
# 
#   rel_err_list[[ext_scen]]$lambda_H_true <- df_true_para[,1]
#   rel_err_list[[ext_scen]]$lambda_S_true <- df_true_para[,2]
#   rel_err_list[[ext_scen]]$lambda_C_true <- df_true_para[,3]
#   rel_err_list[[ext_scen]]$exp_H_true <- df_true_para[,4]
# 
# }
# 
# # save them
# for (ext_scen in 1:5) {
# 
#   print(ext_scen)
# 
#   mu_H_frac <- c(0, 0.3, 0.7, 0,   0.7)[ext_scen]
#   mu_S_frac <- c(0, 0.3, 0,   0.7, 0.7)[ext_scen]
# 
#   file_name <- paste("ex_", 'cophy_ABC_convergence/relative_errors_sizes_truth_included',"_", as.character(mu_H_frac), "_", as.character(mu_S_frac), '.csv', sep = '')
#   write.csv(rel_err_list[[ext_scen]], file_name)
# 
# }

