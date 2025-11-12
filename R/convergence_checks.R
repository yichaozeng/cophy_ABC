# for parallelization
library(future)
library(doParallel)
library(parallel)
n.cores <- future::availableCores()
# registerDoParallel(cores=n.cores - 3)
registerDoParallel(1)

# this is where we check whether the posterior distribution goes toward convergence with decreasing tolerance (epsilon)
perc <- 1 #the initial tolerance

para_sim_md <- para_sim
SS_sim_md <- SS_sim

para_sim_acc <- NULL # accepted values

# for(run_id in 1:14){
for(run_id in c(10)){
  
  run_name <- '1/40000' # the tolerance rate
  print(run_name)

  perc <- 1/40000
  
  # rejection based on Taxicab distance of normalized BLED density
  
  ### a loop over all cophylogenies used as true cophylogenies
  ##################################################
  new_dat <- foreach(i = 1:length(SS_real)) %dopar% {
    
    if(exists('cros_val') == 1){
      print(i)
    }
    
    abc_cophy <- ABC_taxi(target = SS_real[[i]][SS_sel], param = para_sim_md[-cros_val_ids[i], c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')], sumstat = as.matrix(SS_sim_md)[-cros_val_ids[i], SS_sel], tol = perc)
    
    new_dat_sing <- as.data.frame(abc_cophy$unadj.values)
    new_dat_sing$run_name <- rep(run_name, nrow(new_dat_sing))
    
    # here, add the size measures to the date frame
    new_dat_sing$h_tree_size <- rep(sizes_real[[1]][i], nrow(new_dat_sing))
    new_dat_sing$s_tree_size <- rep(sizes_real[[2]][i], nrow(new_dat_sing))
    
    # here add the true parameters to the data frame
    new_dat_sing$lambda_H_true <- rep(para_real[[1]][i], nrow(new_dat_sing))
    new_dat_sing$lambda_S_true <- rep(para_real[[2]][i], nrow(new_dat_sing))
    new_dat_sing$lambda_C_true <- rep(para_real[[3]][i], nrow(new_dat_sing))
    new_dat_sing$exp_H_true <- rep(para_real[[4]][i], nrow(new_dat_sing))

    new_dat_sing
  }

  ##################################################
  new_dat <- do.call(rbind, new_dat)
  
  para_sim_acc <- rbind(para_sim_acc, new_dat)
  
    
}


# compute the residuals for plotting
if(exists('cros_val') == 1){
  res_dat_comb <- para_sim_acc
  res_dat_comb$lambda_H <- res_dat_comb$lambda_H - res_dat_comb$lambda_H_true
  res_dat_comb$lambda_S <- res_dat_comb$lambda_S - res_dat_comb$lambda_S_true
  res_dat_comb$lambda_C <- res_dat_comb$lambda_C - res_dat_comb$lambda_C_true
  res_dat_comb$exp_H <- res_dat_comb$exp_H - res_dat_comb$exp_H_true
}else if(exists('real_data_run') == 1){
  res_dat_comb <- para_sim_acc
  res_dat_comb$lambda_H <- res_dat_comb$lambda_H - para_real$lambda_H
  res_dat_comb$lambda_S <- res_dat_comb$lambda_S - para_real$lambda_S
  res_dat_comb$lambda_C <- res_dat_comb$lambda_C - para_real$lambda_C
  res_dat_comb$exp_H <- res_dat_comb$exp_H - para_real$exp_H
}
