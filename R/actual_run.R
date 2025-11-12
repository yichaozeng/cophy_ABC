para_sim_md <- para_sim
SS_sim_md <- SS_sim

perc <- 1/40000

# rejection based on Taxicab distance of normalized BLenD density
SS_real <- SS_real[[1]]
abc_cophy <- ABC_taxi(target = SS_real[SS_sel], param = para_sim_md[, c('lambda_H', 'lambda_S', 'lambda_C', 'exp_H')], sumstat = as.matrix(SS_sim_md)[, SS_sel], tol = perc) # 1:30 are column indices of the normalized BLenD

para_sim_acc <- as.data.frame(abc_cophy$unadj.values) # parameter values of the accepted simulations
SS_sim_acc <- as.data.frame(abc_cophy$ss) # summary statistics of the accepted simulations

ids_sim_acc <- abc_cophy$indices # indices of the accepted simulations
