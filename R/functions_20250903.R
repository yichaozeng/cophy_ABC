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
    # Option 1: use the values of delta themselves
    blend_1 <- S_to_S_pairs$BL1 - S_to_S_pairs$BL2
    
    # Option 2: use the differences in the values of delta
    temp_vec <- S_to_S_pairs$BL1 - S_to_S_pairs$BL2
    
    # between any two values of delta
    # blend_2 <- combn(temp_vec, 2, function(y) diff(y))
    # or, between two closest values of delta
    temp_vec <- sort(temp_vec)
    blend_2 <- temp_vec[-1] - temp_vec[-length(temp_vec)]
    
    
    blend_2 <- c(blend_2, -blend_2) # this gives you all the permutations instead of combinations
    
    breaks_1 <- breaks
    breaks_2 <- 2 * breaks
    
    # obtain the counts (density-based approach with smoothing)
    dens_1 <- density(blend_1, from = min(breaks_1), to = max(breaks_1))
    counts_1 <- approx(dens_1$x, dens_1$y, xout = (breaks_1[-1] + breaks_1[-length(breaks_1)])/2, method = 'linear')$y
    
    dens_2 <- density(blend_2, from = min(breaks_2), to = max(breaks_2))
    counts_2 <- approx(dens_2$x, dens_2$y, xout = (breaks_2[-1] + breaks_2[-length(breaks_2)])/2, method = 'linear')$y
    
    SS_1 <- rbind(SS_1, counts_1)
    SS_2 <- rbind(SS_2, counts_2)
    
  }
  
  return(
    #list(
      #SS_1,
      SS_2
   #)
  )
  
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
  
  # again, this is a controversial step where only the medians are used
  # return(c(
  #   median(size_tree1),
  #   median(size_tree2),
  #   median(size_net)
  # ))
  
  return(list(
    size_tree1,
    size_tree2,
    size_net
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