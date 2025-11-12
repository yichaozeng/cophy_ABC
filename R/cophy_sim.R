# create a extra variable when extinctions rates are allowed to vary independently: ex_ind <- 1
ex_ind <- 1

n_bin <- 150

# run the simulations, one at a time, with the Slurm scheduler
library(treeducken)
library(rlist)
library(MASS)
library(reshape2)
library(pracma)
library(phytools)

setwd("/groups/cromanpa/yzeng")

setwd("/home/u29/yichaozeng/Desktop")
source('cophy_ABC/R/functions.R')
setwd("/groups/cromanpa/yzeng")

# process the folder id
if(exists('ex_ind') == F){
  folder_id_slash <- paste(folder_id, '/', sep = '')
  prefix <- NULL
}else if(exists('ex_ind') == T){
  folder_id <- NULL
  folder_id_slash <- NULL
  prefix <- "ex_"
}

# initialize the status variable (timeout, network too small, or success)
stat <- "timeout"
draw <- NULL
cophy <- NULL
SS <- NULL

# time parameter
time <- 2

# acceptable size ranges
h_lower <- -Inf
h_upper <- Inf
s_lower <- -Inf
s_upper <- Inf
size_lower <- -Inf
size_upper <- Inf

# rates drawn from the prior distributions in a Monte-Carlo fashion
lambda_H <- rexp(1, 0.5) #rgamma(n=1, 2, 1) #rexp(n=1)
lambda_S <- rexp(1, 0.5) #rgamma(n=1, 2, 1) #rexp(n=1)
lambda_C <- rexp(1, 1.2) #rgamma(n=1, 2, 1) #rexp(n=1)
exp_H <- rexp(1, 0.5) #rgamma(n=1, 2, 1) #rexp(n=1)

ext_frac <- read.csv(file = 'ex_cophy_ABC_statistics/ext_frac.csv')

mu_H_frac <- ext_frac$mu_H_frac    
mu_S_frac <- ext_frac$mu_S_frac     

# compute the the extinction rates
mu_H <- mu_H_frac * (lambda_H + lambda_C)
mu_S <- mu_S_frac * (lambda_S + lambda_C + exp_H)

# run the simulation
draw <- c(lambda_H, lambda_S, lambda_C, exp_H, mu_H_frac, mu_S_frac)

sim <- list(stat, draw, SS)

sim_id <- as.integer(runif(1, min=0, max=99999))
file_name <- paste(prefix, "cophy_ABC_sims/", as.character(folder_id_slash), "sim_", as.character(sim_id), ".rds", sep = '')
folder_path <- paste(prefix, "cophy_ABC_sims/", as.character(folder_id), sep = '')

# this is the time-consuming step
cophy <- sim_cophyBD(hbr = lambda_H, hdr = mu_H, sbr = lambda_S, sdr = mu_S, cosp_rate =lambda_C, host_exp_rate = exp_H, time_to_sim = time, hs_mode = T, numbsim = 100)

# compute the SSs
tr_ht <- cophy[[1]][[1]]$root.edge + max(nodeHeights(cophy[[1]][[1]])) # tree height

# if using Option 1 for SS_norm
# breaks <- seq(from = -tr_ht, to = tr_ht, length.out = n_bin + 1) # length.out should be no. of bins + 1

# if using Option 2 for SS_norm
breaks <- seq(from = -tr_ht, to = tr_ht, length.out = n_bin + 1) # length.out should be no. of bins + 1

# here it is possible to consider separately the size of the cophylogeny when computing the SSs
ids_size <- rep(0, length(cophy))
for (cophy_id in 1:length(cophy)) {
  
  cophy_sing <- cophy[[cophy_id]]
  
  tree1 <- cophy_sing[[1]] # the host tree
  tree2 <- cophy_sing[[2]] # the symbiont tree
  net <- cophy_sing[[3]] # the network
  
  # first, prune host tree to keep only extant tips
  tree1 <- keep.tip(tree1, rownames(net))
  tree2 <- keep.tip(tree2, colnames(net))
  
  ids_size[cophy_id] <- (length(tree1$tip.label) > h_lower) *
                        (length(tree1$tip.label) < h_upper) *
                        (length(tree2$tip.label) > s_lower) *
                        (length(tree2$tip.label) < s_upper) *
                        (sum(net) > size_lower) *
                        (sum(net) < size_upper)
}

# obtain all file names in the folder
all_file_names <- list.files(path = folder_path)
# remove names belonging to folders
#all_file_names <- all_file_names[unlist(lapply(all_file_names, function(x) nchar(x)>4))]

while (sum(all_file_names == file_name) > 0) { # if a file with the same name exists
  sim_id <- as.integer(runif(1, min=0, max=99999))

  file_name <- paste(prefix, "cophy_ABC_sims/", as.character(folder_id_slash), "sim_", as.character(sim_id), ".rds", sep = '')
  folder_path <- paste(prefix, "cophy_ABC_sims/", as.character(folder_id), sep = '')
}

if(sum(ids_size) > 0){ # if there are at least some cophylogenies within the desired size range
  
  SS <- list(
    SS_norm(cophy[which(ids_size == 1)], breaks = breaks, tr_ht = tr_ht),
    # also record the sizes (tree1, tree1, net)
    SS_size(cophy[which(ids_size == 1)])
  )
  
  # record that the simulation is complete
  stat <- "simulation_complete"
  
  # store the status, drawn ext_frac, and the resulting cophy objective
  sim <- list(stat, draw, SS)
  
  list.save(sim, file_name)
  
}
