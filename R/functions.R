# a function for getting the branch length of a tip
get_BL <- function(tree, tip_name){
  return(
    tree$edge.length[tree$edge[,2] == which(tree$tip.label == tip_name)]
  )
}

# SS_norm is a function that computes the average difference in divergence between the host and symbiont for all pairs of species-to-species pairs
SS_norm <- function(cophys, breaks, tr_ht){
  
  SS_1 <- NULL
  SS_2 <- NULL
  SS_3 <- NULL

  for (cophy_id in 1:length(cophys)) {
    
    cophy <- cophys[[cophy_id]]
    
    # all summary statistics can be computed from these three objects
    tree1 <- cophy[[1]] # the host tree
    tree2 <- cophy[[2]] # the symbiont tree
    net <- cophy[[3]] # the network
    
    # first, prune host tree to keep only extant tips
    tree1 <- keep.tip(tree1, rownames(net))
    tree2 <- keep.tip(tree2, colnames(net))
    
    # find all species-to-species pairs
    net_long <- melt(net)
    S_to_S_pairs <- net_long[net_long$value == 1, 1:2]
    
    # compute the branch lengths
    S_to_S_pairs$BL1 <- apply(S_to_S_pairs, MARGIN = 1, function(x) get_BL(tree1, x[1]))
    S_to_S_pairs$BL2 <- apply(S_to_S_pairs, MARGIN = 1, function(x) get_BL(tree2, x[2]))
    
    # the difference
    S_to_S_pairs$div_time_diff <- S_to_S_pairs$BL1 - S_to_S_pairs$BL2
    
    # compute the number of connections
    S_to_S_pairs$conn1 <- apply(S_to_S_pairs, MARGIN = 1, function(x) sum(net[x[1], ]) / length(net[x[1], ]) )
    S_to_S_pairs$conn2 <- apply(S_to_S_pairs, MARGIN = 1, function(x) sum(net[, x[2]]) / length(net[, x[2]]) )
    
    # the difference
    S_to_S_pairs$conn_diff <- S_to_S_pairs$conn1 - S_to_S_pairs$conn2
    
    # the centrality distribution of the symbionts species
    symb_cent <- apply(net, 2, mean)
    
    # obtain the counts (binned-based histogram approach)
    # counts <- hist(S_to_S_pairs$div_time_diff, breaks = breaks, plot = F)$counts
    
    # or
    # obtain the counts (density-based approach with smoothing)
    dens <- density(S_to_S_pairs$div_time_diff, from = min(breaks), to = max(breaks))
    counts <- approx(dens$x, dens$y, xout = (breaks[-1] + breaks[-length(breaks)])/2, method = 'linear')$y
    
    dens_conn <- density(S_to_S_pairs$conn_diff, from = -1, to = 1)
    counts_conn <- approx(dens_conn$x, dens_conn$y, xout = (breaks[-1] + breaks[-length(breaks)])/2 / tr_ht, method = 'linear')$y
    
    dens_symb_cent <- density(symb_cent, from = -1, to = 1)
    counts_symb_cent <- approx(dens_symb_cent$x, dens_symb_cent$y, xout = (breaks[-1] + breaks[-length(breaks)])/2 / tr_ht, method = 'linear')$y
    
    SS_1 <- rbind(SS_1, counts)
    SS_2 <- rbind(SS_2, counts_conn)
    SS_3 <- rbind(SS_3, counts_symb_cent)
    
  }
  
  SS_1 <- apply(SS_1, 2, median)
  SS_2 <- apply(SS_2, 2, median)
  SS_3 <- apply(SS_3, 2, median)

  return(c(
    SS_1 / sum(SS_1), # normalized counts such that the area under the curve is equal to 1
    SS_2 / sum(SS_2),
    SS_3 / sum(SS_3)
  ))
  
}


# SS_size is a function that computes the average host tree size, symbiont tree size, and number of associations
SS_size <- function(cophys){
  
  size_tree1 <- NULL
  size_tree2 <- NULL
  size_net <- NULL
  
  for (cophy_id in 1:length(cophys)) {
    
    cophy <- cophys[[cophy_id]]
    
    # all summary statistics can be computed from these three objects
    tree1 <- cophy[[1]] # the host tree
    tree2 <- cophy[[2]] # the symbiont tree
    net <- cophy[[3]] # the network
    
    size_tree1 <- c(size_tree1, length(tree1$tip.label))
    size_tree2 <- c(size_tree2, length(tree2$tip.label))
    size_net <- c(size_net, sum(net))
    
  }
  
  return(c(
    median(size_tree1),
    median(size_tree2),
    median(size_net)
  ))
  
}


# a function for ABC based on taxicab distance
ABC_taxi <- function(target, param, sumstat, tol){
  distances_raw <- t(apply(sumstat, 1, function(row) row - target))
  dist_taxi <- t(apply(distances_raw, 1, function(row) sum(abs(row))))
  
  rej_thres <- quantile(dist_taxi , tol)
  indices <- dist_taxi <= rej_thres
  
  return(
    list(
      unadj.values = param[dist_taxi <= rej_thres, ],
      adj.values = param[dist_taxi <= rej_thres, ],
      indices = indices,
      ss = sumstat[indices, ],
      dist = dist_taxi[indices]
    )
  )
}