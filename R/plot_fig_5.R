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
WEEEFFFGGGHHHRR
WEEEFFFGGGHHHRR
WEEEFFFGGGHHHRR
WIIIJJJKKKLLLSS
WIIIJJJKKKLLLSS
WIIIJJJKKKLLLSS
WMMMNNNOOOPPPTT
WMMMNNNOOOPPPTT
WMMMNNNOOOPPPTT
weeefffggghhhrr
weeefffggghhhrr
weeefffggghhhrr
wiiijjjkkklllss
wiiijjjkkklllss
wiiijjjkkklllss
wmmmnnnoooppptt
wmmmnnnoooppptt
wmmmnnnoooppptt
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