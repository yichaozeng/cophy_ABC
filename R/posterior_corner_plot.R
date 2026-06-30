ggplot(rate_spec_ratio, aes(x = event, y = ratio)) +
  geom_violin(aes(fill = event), linewidth = 0.5, alpha = 0.4, trim = F) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.5) +
  scale_x_discrete(
    labels = c(
      "lambda_H" = expression(paste(hat(lambda)[H])),
      "lambda_S" = expression(paste(hat(lambda)[S])),
      "lambda_C" = expression(paste(hat(lambda)[C])),
      "lambda_W" = expression(paste(hat(lambda)[W]))
    )#,
    # expand = expansion(add = c(0.5, 1))
  ) +
  # xlim(, ) +
  labs(title = "Speciation rate estimates", y = "events/lineage/Myr", x = NULL) +
  scale_fill_manual(
    values = spec_cl,
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5))




car::scatterplotMatrix(para_sim_acc[,1:6], ellipse = TRUE)


pca <- prcomp(para_sim_acc[,1:6], scale = TRUE)
biplot(pca)


# cornor plot
# do this to plot relative rates
# do this to plot absolute rates
para_sim_acc$mu_H_frac <- para_sim_acc$mu_H_frac * (para_sim_acc$lambda_H + para_sim_acc$lambda_C)
para_sim_acc$mu_S_frac <- para_sim_acc$mu_S_frac * (para_sim_acc$lambda_S + para_sim_acc$lambda_C + para_sim_acc$exp_H)

logit <- function(p) {
  log(p / (1 - p))
}

para_sim_acc$lambda_H <- log(para_sim_acc$lambda_H)
para_sim_acc$lambda_S <- log(para_sim_acc$lambda_S)
para_sim_acc$lambda_C <- log(para_sim_acc$lambda_C)
para_sim_acc$exp_H <- log(para_sim_acc$exp_H)
para_sim_acc$mu_H_frac <- logit(para_sim_acc$mu_H_frac)
para_sim_acc$mu_S_frac <- logit(para_sim_acc$mu_S_frac)

GGally::ggpairs(
  para_sim_acc[, 1:6],
  lower = list(continuous = "points"),
  upper = list(continuous = "density"),
  diag = list(continuous = "densityDiag")
)

GGally::ggpairs(
  dat_joint, aes(color = group),
  lower = list(continuous = "points"),
  upper = list(continuous = "density"),
  # diag = list(continuous = "densityDiag")
  diag = list(continuous = "barDiag")
)
