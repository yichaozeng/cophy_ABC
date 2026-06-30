library(treeducken)
library(phytools)
library(reshape)

setwd('/home/u29/yichaozeng/Desktop')
source('cophy_ABC/R/functions.R')

n_bin <- 150

cophy_real <- sim_cophyBD(hbr = 2.4, hdr = 0.3 * (2.3 + 0.3), sbr = 0.3, sdr = 0.3 * (0.3 + 0.3 + 0.3), cosp_rate = 0.3, host_exp_rate = 0.3, time_to_sim = 2, hs_mode = T, numbsim = 100)

# compute the summary statistics
tr_ht <- cophy_real[[1]][[1]]$root.edge + max(nodeHeights(cophy_real[[1]][[1]])) # tree height
breaks <- seq(from = -tr_ht, to = tr_ht, length.out = n_bin + 1) # length.out should be no. of bins + 1

SS_real <- c(
  SS_norm(cophy_real, breaks = breaks, tr_ht = tr_ht),
  SS_size(cophy_real)
)


plot(1:150, SS_real[1:150])
