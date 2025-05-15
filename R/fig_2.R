library(treeducken)
library(phytools)
library(reshape)
library(rlist)

source('cophy_ABC/R/functions.R')

n_bin <- 300

# here we generate the speciation and extinction rates randomly
epsilon_H <- runif(1, min = 0.3, max = 0.3)
epsilon_S <- runif(1, min = 0.3, max = 0.3) 

cophys <- NULL
paras <- NULL

for (i in 1:2) {
  lambda_H <- c(0.9, 1.1)[i]
  lambda_S <- c(0.3, 1.2)[i]
  lambda_C <- c(0.4, 0.7)[i]
  lambda_W <- c(1.7, 0.3)[i]
  
  # p1 <- runif(1,0,100)
  # p2 <- runif(1,0,100)
  # p3 <- runif(1,0,100)
  # p4 <- runif(1,0,100)
  # 
  # lambda_H <- round(3.3 * p1 / (p1 + p2 + p3 + p4), digits = 1)
  # lambda_S <- round(3.3 * p2 / (p1 + p2 + p3 + p4), digits = 1)
  # lambda_C <- round(3.3 * p3 / (p1 + p2 + p3 + p4), digits = 1)
  # lambda_W <- round(3.3 * p4 / (p1 + p2 + p3 + p4), digits = 1)
  
  mu_H <- epsilon_H * (lambda_H + lambda_C)
  mu_S <- epsilon_S * (lambda_S + lambda_C + lambda_W)
  
  cophys <- c(
    cophys,
    list(sim_cophyBD(hbr = lambda_H, hdr = mu_H, sbr = lambda_S, sdr = mu_S, host_exp_rate = lambda_W, cosp_rate = lambda_C, time_to_sim = 2, numbsim = 500, hs_mode = T))
  )
  
  paras <- c(
    paras,
    list(c(lambda_H, lambda_S, lambda_C, lambda_W))
  )
}


# compute the SSs
tr_ht <- cophys[[1]][[1]][[1]]$root.edge + max(nodeHeights(cophys[[1]][[1]][[1]])) # tree height
breaks <- seq(from = -tr_ht, to = tr_ht, length.out = n_bin + 1) # length.out should be no. of bins + 1

samples_1 <- SS_norm(cophys[[1]], breaks = breaks, tr_ht = tr_ht)[1:n_bin]
samples_2 <- SS_norm(cophys[[2]], breaks = breaks, tr_ht = tr_ht)[1:n_bin]

# or
# cophy_1 <- sim_cophyBD(0.3, 0, 0.3, 0, 0.3, 0.3, 2, 500, hs_mode = T)
# cophy_2 <- sim_cophyBD(0.6, 0, 0.1, 0, 0.1, 0.6, 2, 500, hs_mode = T)
# 
# tr_ht <- cophy_1[[1]][[1]]$root.edge + max(nodeHeights(cophy_1[[1]][[1]])) # tree height
# breaks <- seq(from = -tr_ht, to = tr_ht, length.out = n_bin + 1) # length.out should be no. of bins + 1
# 
# samples_1 <- SS_norm(cophy_1, breaks = breaks, tr_ht = tr_ht)[1:n_bin] / tr_ht
# samples_2 <- SS_norm(cophy_2, breaks = breaks, tr_ht = tr_ht)[1:n_bin] / tr_ht

# plotting
x <- (breaks[-1] + breaks[-length(breaks)]) / (2 * tr_ht)
y1 <- samples_1 / (2 / n_bin)
y2 <- samples_2 / (2 / n_bin)

upper <- y1 * (y1 > y2) + y2 * (y1<=y2)
lower <- y2 * (y1 > y2) + y1 * (y1<=y2)

# Plot the curves
plot(x, y1, type = "l", col = "blue", lwd = 3, ylim = c(min(lower), max(upper)), xlim = c(-1,1), xlab = 'BL1 - BL2', ylab = "density", main = 'BLenD distribution')
lines(x, y2, col = "red", lwd = 3, lty = 2)  # Add lower bound

# Fill the area between the curves
polygon(c(x, rev(x)), c(lower, rev(upper)), col = rgb(0, 0, 0, 0.3, alpha = 0.15), border = NA)



# ggplot2 version of the same plot
library(ggplot2)

df <- data.frame(
  x = x,
  y1 = y1,
  y2 = y2,
  lower = lower,
  upper = upper
)

ribbon_df <- data.frame(
  x = c(x, rev(x)),
  y = c(lower, rev(upper))
)

plot_curves <- ggplot(df, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'black', alpha = 0.1, inherit.aes = TRUE) +
  geom_line(aes(y = y1), color = "red", size = 0.7) +
  geom_line(aes(y = y2), color = "black", linetype = I(2), size = 0.7) +
  labs(
    x = expression(paste(delta, ' (BLenD)', sep = '')),
    y = "Density",
    title = "b) BLenD curve"
  ) +
  xlim(-1, 1) +
  ylim(min(lower), max(upper)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0))

print(paras)

# save the plot the the parameters
file_name <- paste('ex_cophy_ABC_fig2/fig2_curves.rds', sep = '')
list.save(plot_curves, file_name)

file_name <- paste('ex_cophy_ABC_fig2/fig2_curves.csv', sep = '')
write.csv(rbind(paras[[1]], paras[[2]]), file = file_name)

print('plotting')
#png(paste(prefix, 'cophy_ABC_results/', "convergence_comb", "_", as.character(mu_H_frac), "_", as.character(mu_S_frac), ".png", sep = ''), width = 10, height = 10, units = 'in', pointsize = 12, res = 300)
pdf(paste('ex_cophy_ABC_fig2/', "fig2_curves.pdf", sep = ''), width = 4.5/1.2, height = 3.5/1.2, pointsize = 12)
print(plot_curves)
dev.off()




