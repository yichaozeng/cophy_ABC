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
for(run_id in c(1,4,7,10,13)){
  
  run_name <- c('1', '1/2', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256', '1/512', '1/1024', '1/2048', '1/4096', '1/8192')[run_id]
  print(run_name)
  
  perc <- c(1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096, 1/8192)[run_id]
  
  # rejection based on Taxicab distance of normalized BLED density
  abc_cophy <- ABC_taxi(target = SS_real[SS_sel], param = para_sim_md[, c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')], sumstat = as.matrix(SS_sim_md)[, SS_sel], tol = perc) # 1:30 are column indices of the normalized BLenD

  
  new_dat <- as.data.frame(abc_cophy$unadj.values)
  new_dat$run_name <- rep(run_name, nrow(new_dat))
  
  para_sim_acc <- rbind(para_sim_acc, new_dat)
    
}

# the parameter estimate
para_sim_acc$run_name <- factor(para_sim_acc$run_name, levels = c('1', '1/2', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256', '1/512', '1/1024', '1/2048', '1/4096', '1/8192', 'NN'))


# compute the residuals for plotting
res_dat_comb <- para_sim_acc
res_dat_comb$lambda_H <- res_dat_comb$lambda_H - para_real$lambda_H
res_dat_comb$lambda_S <- res_dat_comb$lambda_S - para_real$lambda_S
res_dat_comb$lambda_C <- res_dat_comb$lambda_C - para_real$lambda_C
res_dat_comb$exp_H <- res_dat_comb$exp_H - para_real$exp_H
res_dat_comb$run_name <- factor(res_dat_comb$run_name, levels = c('1', '1/2', '1/4', '1/8', '1/16', '1/32', '1/64', '1/128', '1/256', '1/512', '1/1024', '1/2048', '1/4096', '1/8192', 'NN'))


##############################################
# compute the percentage of convergence
library(future)
library(doParallel)
library(parallel)
n.cores <- future::availableCores()
registerDoParallel(cores=n.cores - 3)


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
  
  # compute the residuals
  res_dat_comb <- as.data.frame(abc_cophy$unadj.values)
  res_dat_comb$lambda_H_res <- res_dat_comb$lambda_H - para_real$lambda_H
  res_dat_comb$lambda_S_res <- res_dat_comb$lambda_S - para_real$lambda_S
  res_dat_comb$lambda_C_res <- res_dat_comb$lambda_C - para_real$lambda_C
  res_dat_comb$exp_H_res <- res_dat_comb$exp_H - para_real$exp_H
  
  # yes if the estimates fall within the true value + (-1, +1)
  yes_no_1 <- (abs(median(res_dat_comb$lambda_H_res)) < 1) * 
    (abs(median(res_dat_comb$lambda_S_res)) < 1) * 
    (abs(median(res_dat_comb$lambda_C_res)) < 1) * 
    (abs(median(res_dat_comb$exp_H_res)) < 1)
  
  # yes if the maximum rate can be correctly identified
  yes_no_2 <-  which.max(c(
    median(res_dat_comb$lambda_H),
    median(res_dat_comb$lambda_S),
    median(res_dat_comb$lambda_C),
    median(res_dat_comb$exp_H)
  )) == which.max(c(
    para_real$lambda_H,
    para_real$lambda_S,
    para_real$lambda_C,
    para_real$exp_H
  ))
  
  c(yes_no_1, yes_no_2)
  
}

# conv_perc <- mean(unlist(conv_res))

temp1 <- unlist(lapply(conv_res, function(x) x[1]))
temp2 <- unlist(lapply(conv_res, function(x) x[2]))
conv_perc <- list(
  mean(temp1),
  mean(temp2)
)

