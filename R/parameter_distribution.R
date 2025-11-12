library(future)
library(doParallel)
library(parallel)
n.cores <- future::availableCores()
registerDoParallel(cores=n.cores - 3)

library(rlist)

setwd("/home/u29/yichaozeng/Desktop")
source('cophy_ABC/R/functions.R')

setwd("/groups/cromanpa/yzeng")
SS <- NULL
para_dat_separate <- NULL

# specify the name of the folder where plots will be stored
plot_folder <- paste(prefix, 'cophy_ABC_plots/', folder_id_slash, sep = '')

# specify the name of the folder where simulations will be read
folder_path <- paste(prefix, "cophy_ABC_sims/", folder_id, sep = '')

all_file_names <- as.list(list.files(path = folder_path))
# remove names belonging to folders
all_file_names <- all_file_names[unlist(lapply(all_file_names, function(x) nchar(x)>8))][1:12000]


SS_para_list <- foreach(i = 1: length(all_file_names)) %dopar% {
  #SS_para_list <- foreach(i = 1: 10) %dopar% {
  #print(i)
  
  SS_row <- NA
  para_row <- NA
  group_row <- NA
  
  sim <- tryCatch({
    list.load(file = paste(folder_path, '/', all_file_names[[i]], sep = ''))
  }, error = function(e) {
    NA
  })
  
  if(sum(!is.na(sim)) > 0){
    if(sim[[1]] != 'timeout'){
      
      SS_row <- sim[[3]]
      para_row <- sim[[2]]
      group_row <- sim[[1]]
      
    }
  }
  
  list(
    SS_row,
    para_row,
    group_row
  )
}

SS <- NULL
para <- NULL
group <- NULL

# here the lengths of para and group need to match that of SS
for(row_id in 1:length(SS_para_list)){
  #print(row_id)
  if(sum(!is.na(SS_para_list[[row_id]]))>0){
    if(SS_para_list[[row_id]][[3]] == 'simulation_complete' & length(SS_para_list[[row_id]][[1]]) > 0){
      
      # below, nrow is the number of replicate runs (see numbsim in cophy_sim.R)
      ele <- SS_para_list[[row_id]][[2]]
      para <- rbind(para, 
                    matrix(ele, nrow = 100, ncol = length(ele), byrow = TRUE)
                    )
      
      ele <- SS_para_list[[row_id]][[3]]
      group <- rbind(group, 
                     matrix(ele, nrow = 100, ncol = length(ele), byrow = TRUE)
                     )
    }
  }
}

SS <- foreach(row_id = 1:length(SS_para_list)) %dopar% {
  #print(row_id)
  if(sum(!is.na(SS_para_list[[row_id]]))>0){
    if(SS_para_list[[row_id]][[3]] == 'simulation_complete' & length(SS_para_list[[row_id]][[1]]) > 0){
      
      all_SS <- SS_para_list[[row_id]][[1]]
      
      cbind(all_SS[[1]][[1]], all_SS[[1]][[2]], all_SS[[2]][[1]], all_SS[[2]][[2]], all_SS[[2]][[3]])
      
    }
  }
}

SS <- do.call(rbind, SS)

# transform into a data frame
para_dat_separate <- data.frame(
  lambda_H = para[,1],
  lambda_S = para[,2],
  lambda_C = para[,3],
  exp_H = para[,4],
  mu_H_frac = para[,5],
  mu_S_frac = para[,6],
  group = group
)
para_dat_separate$group[para_dat_separate$group == 'simulation_complete'] <- "finished"

row.names(SS) <- NULL
write.csv(SS, file = paste(prefix, "cophy_ABC_statistics/", folder_id, "/SS.csv", sep = ''), row.names = F)
write.csv(para_dat_separate, file = paste(prefix, "cophy_ABC_statistics/", folder_id, "/para_dat_separate.csv", sep = ''), row.names = F)
