# create a extra variable when extinctions rates are allowed to vary independently: ex_ind <- 1

# specify the extinctin fractions
setwd("/groups/cromanpa/yzeng")

# run the simulations, one at a time, with the Slurm scheduler
library(future)
library(doParallel)
library(parallel)
n.cores <- future::availableCores()
registerDoParallel(cores=n.cores - 3)

foreach(i = 1:5000) %dopar% {
  
  setwd("/home/u29/yichaozeng/Desktop")
  system("Rscript /home/u29/yichaozeng/Desktop/cophy_ABC/R/cophy_sim.R", timeout = 600)
  setwd("/groups/cromanpa/yzeng")
  
  print(i)
}
