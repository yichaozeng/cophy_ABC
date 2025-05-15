setwd("/home/u29/yichaozeng/Desktop")

# process the folder id
folder_id <- NULL
folder_id_slash <- NULL
prefix <- "ex_"

# for each folder
#for (lab2 in c('_0/', '/', '_SMC/', '_SMC_2/', '_SMC_3/', '_SMC_4/')) {
for (lab2 in c('/')) {
  for (lab1 in c('_plots', '_sims', '_statistics', '_checks', '_convergence', '_results')){
    
    folder_path <- paste(prefix, "cophy_ABC", lab1, lab2, folder_id, sep = '')
    print('Archiving:')
    print(folder_path)
    
    all_file_names <- as.list(list.files(path = folder_path))
    # remove names belonging to folders
    all_file_names_folder_only <- all_file_names[unlist(lapply(all_file_names, function(x) file.info(paste(folder_path, x, sep = ''))$isdir))]
    all_file_names_file_only <- all_file_names[unlist(lapply(all_file_names, function(x) file.info(paste(folder_path, x, sep = ''))$isdir)) == F]
    
    # create a folder
    # randomly pick a folder name
    if(exists('new_folder_name') == F){
      new_folder_name <- as.character(0)
      while (sum(all_file_names_folder_only == new_folder_name) == 1) { # if a folder with the same name already exists
        new_folder_name <-as.character(sample(0:99, size = 1))
      }
    }
    new_folder_path <- paste(folder_path, new_folder_name, '/', sep = '')
    dir.create(new_folder_path)
    
    # move the files into the folder
    if(lab1 != '_convergence'){
      lapply(all_file_names_file_only, function(x) file.rename(paste(folder_path, x, sep = ''), paste(new_folder_path, x, sep = '')))
    }
    
  }
}
