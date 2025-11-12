# for parallelization
library(future)
library(doParallel)
library(parallel)
n.cores <- future::availableCores()
registerDoParallel(cores=n.cores - 3)

# this is where we check whether the posterior distribution goes toward convergence with decreasing tolerance (epsilon)
perc <- 1 #the initial tolerance

para_sim_md <- para_sim
SS_sim_md <- SS_sim

para_sim_acc <- NULL # accepted values

# for(run_id in 1:14){
for(run_id in c(1,2,5,7,9,11)){
  
  run_name <- c('1', '1/2', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256', '1/512', '1/1024', '1/2048', '1/4096', '1/8192')[run_id]
  print(run_name)
  
  perc <- c(1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096, 1/8192)[run_id]
  
  # rejection based on Taxicab distance
  par_est_list <- foreach(cophy_id = 1: length(cophy_real)) %dopar% {
    
    cophy_real_sing <- cophy_real[cophy_id]
    # compute the summary statistics
    tr_ht <- cophy_real_sing[[1]][[1]]$root.edge + max(nodeHeights(cophy_real_sing[[1]][[1]])) # tree height
    breaks <- seq(from = -tr_ht, to = tr_ht, length.out = n_bin + 1) # length.out should be no. of bins + 1
    
    SS_real <- c(
      SS_norm(cophy_real_sing, breaks = breaks, tr_ht = tr_ht),
      SS_size(cophy_real_sing)
    )
    
    # rescale the summary statistics
    SS_real <- renormalize(SS_real)
    SS_real <- rescale(SS_real, sd_vec)
    
    para_real <- para_ABC[1, ] # for formatting
    
    para_real$lambda_H <- lambda_H_real
    para_real$lambda_S <- lambda_S_real
    para_real$lambda_C <- lambda_C_real
    para_real$exp_H <- exp_H_real
    para_real$mu_H_frac
    para_real$mu_S_frac
    
    # remove the true SSs and para from the two data frames
    para_sim <- para_ABC[,c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')] #para_sim <- para_ABC[-ids, ,c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')]
    SS_sim <- SS_ABC #SS_sim <- SS_ABC[-ids, ]
    
    para_sim_md <- para_sim
    SS_sim_md <- SS_sim
    
    # the ABC
    abc_cophy <- ABC_taxi(target = SS_real[SS_sel], param = para_sim_md[, c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')], sumstat = as.matrix(SS_sim_md)[, SS_sel], tol = perc) # 1:30 are column indices of the normalized BLenD
    
    as.data.frame(abc_cophy$unadj.values)
    
  }
  
  new_dat <- do.call(rbind, par_est_list)
  new_dat$run_name <- rep(run_name, nrow(new_dat))
  
  para_sim_acc <- rbind(para_sim_acc, new_dat)
    
}


# compute the residuals for plotting
res_dat_comb <- para_sim_acc
res_dat_comb$lambda_H <- res_dat_comb$lambda_H - para_real$lambda_H
res_dat_comb$lambda_S <- res_dat_comb$lambda_S - para_real$lambda_S
res_dat_comb$lambda_C <- res_dat_comb$lambda_C - para_real$lambda_C
res_dat_comb$exp_H <- res_dat_comb$exp_H - para_real$exp_H
res_dat_comb$run_name <- factor(res_dat_comb$run_name, levels = c('1', '1/2', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256', '1/512', '1/1024', '1/2048', '1/4096', '1/8192', 'NN'))


##############################################
# calculate the percentage of "good convergences"
conv_res <- foreach(cophy_id = 1: length(cophy_real)) %dopar% {
  
  cophy_real_sing <- cophy_real[cophy_id]
  # compute the summary statistics
  tr_ht <- cophy_real_sing[[1]][[1]]$root.edge + max(nodeHeights(cophy_real_sing[[1]][[1]])) # tree height
  breaks <- seq(from = -tr_ht, to = tr_ht, length.out = n_bin + 1) # length.out should be no. of bins + 1
  
  SS_real <- c(
    SS_norm(cophy_real_sing, breaks = breaks, tr_ht = tr_ht),
    SS_size(cophy_real_sing)
  )
  
  # rescale the summary statistics
  SS_real <- renormalize(SS_real)
  SS_real <- rescale(SS_real, sd_vec)
  
  para_real <- para_ABC[1, ] # for formatting
  
  para_real$lambda_H <- lambda_H_real
  para_real$lambda_S <- lambda_S_real
  para_real$lambda_C <- lambda_C_real
  para_real$exp_H <- exp_H_real
  para_real$mu_H_frac
  para_real$mu_S_frac
  
  # remove the true SSs and para from the two data frames
  para_sim <- para_ABC[,c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')] #para_sim <- para_ABC[-ids, ,c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')]
  SS_sim <- SS_ABC #SS_sim <- SS_ABC[-ids, ]
  
  para_sim_md <- para_sim
  SS_sim_md <- SS_sim
  
  # the ABC
  abc_cophy <- ABC_taxi(target = SS_real[SS_sel], param = para_sim_md[, c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')], sumstat = as.matrix(SS_sim_md)[, SS_sel], tol = perc) # 1:30 are column indices of the normalized BLenD
  
  # compute the residuls
  # compute the residuls
  res_dat_comb <- as.data.frame(abc_cophy$unadj.values)
  res_dat_comb$lambda_H <- res_dat_comb$lambda_H - para_real$lambda_H
  res_dat_comb$lambda_S <- res_dat_comb$lambda_S - para_real$lambda_S
  res_dat_comb$lambda_C <- res_dat_comb$lambda_C - para_real$lambda_C
  res_dat_comb$exp_H <- res_dat_comb$exp_H - para_real$exp_H
  
  yes_no <- (abs(median(res_dat_comb$lambda_H)) < 1) * 
    (abs(median(res_dat_comb$lambda_S)) < 1) * 
    (abs(median(res_dat_comb$lambda_C)) < 1) * 
    (abs(median(res_dat_comb$exp_H)) < 1)
  
  yes_no
  
}

conv_perc <- mean(unlist(conv_res))