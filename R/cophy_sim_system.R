# create a extra variable when extinctions rates are allowed to vary independently: ex_ind <- 1

# specify the extinctin fractions
setwd("/home/u29/yichaozeng/Desktop")

ext_frac <- data.frame(
  #lambda_H <- rexp(n=1), #rgamma(n=1, 2, 1) #rexp(n=1)
  #lambda_S <- rexp(n=1), #rgamma(n=1, 2, 1) #rexp(n=1)
  #lambda_C <- rexp(n=1), #rgamma(n=1, 2, 1) #rexp(n=1)
  #exp_H <- rexp(n=1), #rgamma(n=1, 2, 1) #rexp(n=1)
  mu_H_frac = runif(1, min = 0, max = 0), # here you can specify the upper limit of the extinction fraction - a high extinction fraction is perhaps more likely to lead to poor inference performance.
  mu_S_frac = runif(1, min = 0.7, max = 0.7) # same as above.
)

write.csv(ext_frac, file = 'ex_cophy_ABC_statistics/ext_frac.csv')


# run the simulations, one at a time, with the Slurm scheduler
library(future)
library(doParallel)
library(parallel)
n.cores <- future::availableCores()
registerDoParallel(cores=n.cores - 3)

foreach(i = 1:5000) %dopar% {
  system("Rscript /home/u29/yichaozeng/Desktop/cophy_ABC/R/cophy_sim.R", timeout = 600)
  print(i)
}
