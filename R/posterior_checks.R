# this is where we plot the posterior predictive checks
library(ggplot2)
library(ggpubr)

posterior_dat <- cbind(para_sim, SS_sim)
posterior_dat$acc_rej <- abc_cophy$region

posterior_dat <- posterior_dat[order(posterior_dat$acc_rej, decreasing = F), ]

SS_real_check <- SS_ABC[1,]

p_1_1 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=ntip1, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$ntip1, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_2 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=ntip2, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$ntip2, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_3 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=sackin_1_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$sackin_1_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_4 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=sackin_2_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$sackin_1_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_5 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=gamma_H, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$gamma_H, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_6 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=gamma_S, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$gamma_S, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_7 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=conn, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$conn, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_8 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=cent_mean, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$cent_mean, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_9 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=S_to_H_mean, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$S_to_H_mean, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_10 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=phylo_div_host, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$phylo_div_host, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_11 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=phylo_div_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$phylo_div_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_12 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=NODF_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$NODF_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_13 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=NODF_host, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$NODF_host, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_14 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=mod_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$mod_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_15 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=mantel_u_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$mantel_u_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_16 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif1, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif1, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_17 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif2, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif2, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_18 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif3, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif3, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_19 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif4, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif4, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_20 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif5, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif5, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_21 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif6, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif6, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_22 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif7, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif7, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_23 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif8, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif8, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_24 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif9, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif9, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_25 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif10, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif10, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_26 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif11, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif11, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_27 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif12, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif12, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_28 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif13, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif13, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_29 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif14, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif14, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_30 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif15, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif15, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_31 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif16, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif16, linetype = "dashed", linewidth = 0.5, color = "black")

p_1_32 <- ggplot(data = posterior_dat, aes(x=lambda_H, y=motif17, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif17, linetype = "dashed", linewidth = 0.5, color = "black")

post_plot_1 <- ggarrange(
  p_1_1, p_1_2, p_1_3, p_1_4, p_1_5, p_1_6, p_1_7, p_1_8, p_1_9, p_1_10, p_1_11, p_1_12, p_1_13, p_1_14, p_1_15,
  p_1_16, p_1_17, p_1_18, p_1_19, p_1_20, p_1_21, p_1_22, p_1_23, p_1_24, p_1_25, p_1_26, p_1_27, p_1_28, p_1_29, p_1_30, p_1_31, p_1_32,
  ncol = 8, nrow = 4,
  common.legend = TRUE, legend = "bottom" # Optional: shared legend
)





p_2_1 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=ntip1, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$ntip1, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_2 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=ntip2, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$ntip2, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_3 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=sackin_1_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$sackin_1_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_4 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=sackin_2_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$sackin_1_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_5 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=gamma_H, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$gamma_H, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_6 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=gamma_S, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$gamma_S, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_7 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=conn, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$conn, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_8 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=cent_mean, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$cent_mean, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_9 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=S_to_H_mean, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$S_to_H_mean, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_10 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=phylo_div_host, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$phylo_div_host, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_11 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=phylo_div_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$phylo_div_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_12 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=NODF_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$NODF_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_13 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=NODF_host, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$NODF_host, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_14 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=mod_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$mod_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_15 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=mantel_u_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$mantel_u_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_16 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif1, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif1, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_17 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif2, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif2, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_18 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif3, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif3, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_19 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif4, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif4, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_20 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif5, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif5, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_21 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif6, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif6, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_22 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif7, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif7, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_23 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif8, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif8, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_24 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif9, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif9, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_25 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif10, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif10, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_26 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif11, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif11, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_27 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif12, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif12, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_28 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif13, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif13, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_29 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif14, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif14, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_30 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif15, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif15, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_31 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif16, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif16, linetype = "dashed", linewidth = 0.5, color = "black")

p_2_32 <- ggplot(data = posterior_dat, aes(x=lambda_S, y=motif17, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif17, linetype = "dashed", linewidth = 0.5, color = "black")

post_plot_2 <- ggarrange(
  p_2_1, p_2_2, p_2_3, p_2_4, p_2_5, p_2_6, p_2_7, p_2_8, p_2_9, p_2_10, p_2_11, p_2_12, p_2_13, p_2_14, p_2_15,
  p_2_16, p_2_17, p_2_18, p_2_19, p_2_20, p_2_21, p_2_22, p_2_23, p_2_24, p_2_25, p_2_26, p_2_27, p_2_28, p_2_29, p_2_30, p_2_31, p_2_32,
  ncol = 8, nrow = 4,
  common.legend = TRUE, legend = "bottom" # Optional: shared legend
)





p_3_1 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=ntip1, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$ntip1, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_2 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=ntip2, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$ntip2, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_3 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=sackin_1_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$sackin_1_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_4 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=sackin_2_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$sackin_1_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_5 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=gamma_H, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$gamma_H, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_6 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=gamma_S, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$gamma_S, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_7 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=conn, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$conn, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_8 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=cent_mean, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$cent_mean, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_9 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=S_to_H_mean, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$S_to_H_mean, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_10 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=phylo_div_host, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$phylo_div_host, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_11 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=phylo_div_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$phylo_div_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_12 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=NODF_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$NODF_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_13 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=NODF_host, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$NODF_host, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_14 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=mod_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$mod_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_15 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=mantel_u_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$mantel_u_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_16 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif1, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif1, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_17 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif2, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif2, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_18 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif3, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif3, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_19 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif4, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif4, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_20 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif5, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif5, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_21 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif6, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif6, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_22 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif7, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif7, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_23 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif8, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif8, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_24 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif9, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif9, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_25 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif10, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif10, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_26 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif11, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif11, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_27 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif12, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif12, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_28 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif13, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif13, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_29 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif14, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif14, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_30 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif15, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif15, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_31 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif16, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif16, linetype = "dashed", linewidth = 0.5, color = "black")

p_3_32 <- ggplot(data = posterior_dat, aes(x=lambda_C, y=motif17, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif17, linetype = "dashed", linewidth = 0.5, color = "black")

post_plot_3 <- ggarrange(
  p_3_1, p_3_2, p_3_3, p_3_4, p_3_5, p_3_6, p_3_7, p_3_8, p_3_9, p_3_10, p_3_11, p_3_12, p_3_13, p_3_14, p_3_15,
  p_3_16, p_3_17, p_3_18, p_3_19, p_3_20, p_3_21, p_3_22, p_3_23, p_3_24, p_3_25, p_3_26, p_3_27, p_3_28, p_3_29, p_3_30, p_3_31, p_3_32,
  ncol = 8, nrow = 4,
  common.legend = TRUE, legend = "bottom" # Optional: shared legend
)






p_4_1 <- ggplot(data = posterior_dat, aes(x=exp_H, y=ntip1, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$ntip1, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_2 <- ggplot(data = posterior_dat, aes(x=exp_H, y=ntip2, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$ntip2, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_3 <- ggplot(data = posterior_dat, aes(x=exp_H, y=sackin_1_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$sackin_1_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_4 <- ggplot(data = posterior_dat, aes(x=exp_H, y=sackin_2_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$sackin_1_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_5 <- ggplot(data = posterior_dat, aes(x=exp_H, y=gamma_H, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$gamma_H, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_6 <- ggplot(data = posterior_dat, aes(x=exp_H, y=gamma_S, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$gamma_S, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_7 <- ggplot(data = posterior_dat, aes(x=exp_H, y=conn, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$conn, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_8 <- ggplot(data = posterior_dat, aes(x=exp_H, y=cent_mean, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$cent_mean, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_9 <- ggplot(data = posterior_dat, aes(x=exp_H, y=S_to_H_mean, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$S_to_H_mean, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_10 <- ggplot(data = posterior_dat, aes(x=exp_H, y=phylo_div_host, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$phylo_div_host, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_11 <- ggplot(data = posterior_dat, aes(x=exp_H, y=phylo_div_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$phylo_div_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_12 <- ggplot(data = posterior_dat, aes(x=exp_H, y=NODF_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$NODF_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_13 <- ggplot(data = posterior_dat, aes(x=exp_H, y=NODF_host, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$NODF_host, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_14 <- ggplot(data = posterior_dat, aes(x=exp_H, y=mod_raw, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$mod_raw, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_15 <- ggplot(data = posterior_dat, aes(x=exp_H, y=mantel_u_symb, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$mantel_u_symb, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_16 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif1, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif1, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_17 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif2, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif2, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_18 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif3, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif3, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_19 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif4, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif4, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_20 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif5, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif5, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_21 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif6, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif6, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_22 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif7, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif7, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_23 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif8, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif8, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_24 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif9, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif9, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_25 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif10, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif10, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_26 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif11, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif11, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_27 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif12, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif12, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_28 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif13, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif13, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_29 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif14, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif14, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_30 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif15, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif15, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_31 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif16, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif16, linetype = "dashed", linewidth = 0.5, color = "black")

p_4_32 <- ggplot(data = posterior_dat, aes(x=exp_H, y=motif17, color = acc_rej)) +
  geom_point() +
  geom_hline(yintercept = as.data.frame(SS_real_check)$motif17, linetype = "dashed", linewidth = 0.5, color = "black")

post_plot_4 <- ggarrange(
  p_4_1, p_4_2, p_4_3, p_4_4, p_4_5, p_4_6, p_4_7, p_4_8, p_4_9, p_4_10, p_4_11, p_4_12, p_4_13, p_4_14, p_4_15,
  p_4_16, p_4_17, p_4_18, p_4_19, p_4_20, p_4_21, p_4_22, p_4_23, p_4_24, p_4_25, p_4_26, p_4_27, p_4_28, p_4_29, p_4_30, p_4_31, p_4_32,
  ncol = 8, nrow = 4,
  common.legend = TRUE, legend = "bottom" # Optional: shared legend
)





png(paste(plot_folder, 'post_check_lambda_H', sep = ''), width = 14, height = 5, units = 'in', pointsize = 12, res = 300)
print(post_plot_1)
dev.off()

png(paste(plot_folder, 'post_check_lambda_S', sep = ''), width = 14, height = 5, units = 'in', pointsize = 12, res = 300)
print(post_plot_2)
dev.off()

png(paste(plot_folder, 'post_check_lambda_C', sep = ''), width = 14, height = 5, units = 'in', pointsize = 12, res = 300)
print(post_plot_3)
dev.off()

png(paste(plot_folder, 'post_check_exp_H', sep = ''), width = 14, height = 5, units = 'in', pointsize = 12, res = 300)
print(post_plot_4)
dev.off()
