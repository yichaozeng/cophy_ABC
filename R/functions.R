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
    
    # Option 1: use the values of delta, where delta is the difference between associated branches
    # here we only use the pairwise differences * 1/2 (such that the range is still between -tr_ht and tr_ht)
    blend_1 <- (S_to_S_pairs$BL1 - S_to_S_pairs$BL2) / tr_ht # the original BLenD
    blend_1 <- combn(blend_1, 2, function(y) abs(diff(y))) # the differences (all positives)
    blend_1 <- c(blend_1, -blend_1)  # the differences (positives and negatives)
    blend_1 <- blend_1/2
    
    # Option 2: use the values of delta, where delta is the difference between the lengths of host branches of sister symbiont species
    # go over every tip in the symbiont tree to find the root node id of each
    daug_moth_mat <- NULL # a matrix used to record the daughter-mother relationships
    for (tip_id in 1:length(tree2$tip.label)) {
      # the id and root node id of this tip
      daug_moth_mat <- rbind(
        daug_moth_mat,
        c(tip_id, tree2$edge[,1][which(tree2$edge[,2] == tip_id)])
      )
    }
    sist_sp_mat <- NULL # a matrix used to record the sister species relationships
    # record the species pairs that share the same root node
    for (root_id in unique(daug_moth_mat[,2])) {
      if(length(which(daug_moth_mat[,2] == root_id)) > 1){ # if there are at least two species sharing the root node
        sist_sp_mat <- rbind(
          sist_sp_mat,
          which(daug_moth_mat[,2] == root_id)
        )
      }
    }
    # now, go over the matrix of sister species to find the branch lengths of their host tips and their differences (between all possible pairs)
    blend_2 <- NULL
    for (pair_id in 1:nrow(sist_sp_mat)) {
      symb_sp1 <- sist_sp_mat[pair_id, 1]
      symb_sp1_host <- which(net[, symb_sp1] == 1)
      symb_sp1_host_BL <- sapply(symb_sp1_host, function(x) get_BL(tree1, tree1$tip.label[x]))
      
      symb_sp2 <- sist_sp_mat[pair_id, 2]
      symb_sp2_host <- which(net[, symb_sp2] == 1)
      symb_sp2_host_BL <- sapply(symb_sp2_host, function(x) get_BL(tree1, tree1$tip.label[x]))
      
      diff <- as.numeric(outer(symb_sp1_host_BL, symb_sp2_host_BL, `-`))
      
      # since the assignment of sp1 and sp2 are arbitrary, we take the opposite of diff as well
      blend_2 <- c(
        blend_2,
        c(diff, -diff)
      )
    }
    blend_2 <- blend_2 / tr_ht
    

    # obtain the counts (density-based approach with smoothing)
    dens_1 <- density(blend_1, from = -1, to = 1)
    counts_1 <- approx(dens_1$x, dens_1$y, xout = (breaks[-1] + breaks[-length(breaks)])/(2*tr_ht), method = 'linear')$y
    
    dens_2 <- density(blend_2, from = -1, to = 1)
    counts_2 <- approx(dens_2$x, dens_2$y, xout = (breaks[-1] + breaks[-length(breaks)])/(2*tr_ht), method = 'linear')$y
    
    SS_1 <- rbind(SS_1, counts_1)
    SS_2 <- rbind(SS_2, counts_2)
    
  }
  
  return(
    list(
      SS_1, #blend 1
      SS_2  #blend 2
   )
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
    
    # first, prune host tree to keep only extant tips
    tree1 <- keep.tip(tree1, rownames(net))
    tree2 <- keep.tip(tree2, colnames(net))
    
    # fill the tree size vectors
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

# a function for rejection-based ABC based on taxicab distance
ABC_taxi <- function(target, param, sumstat, tol){
  distances_raw <- t(apply(sumstat, 1, function(row) row - target))
  dist_taxi <- t(apply(distances_raw, 1, function(row) sum(abs(row))))
  
  # rej_thres <- quantile(dist_taxi , tol)
  # indices <- dist_taxi <= rej_thres
  
  # or an alternative version (rank-based) that can handle simulations of equal distance to the data
  n_acc <- ceiling(tol * length(dist_taxi))
  indices <- order(dist_taxi)[1:n_acc]

  return(
    list(
      unadj.values = param[indices, ], #unadj.values = param[dist_taxi <= rej_thres, ],
      adj.values = param[indices, ], #adj.values = param[dist_taxi <= rej_thres, ],
      indices = indices,
      ss = sumstat[indices, ],
      dist = dist_taxi[indices]
    )
  )
}