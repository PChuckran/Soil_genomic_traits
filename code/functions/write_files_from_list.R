write_files_from_list <- function(xx, path = 'data/input/'){
## also write data as CSVs for Pete
lapply(xx, function(x){
  #browser()
  ## only keep data files, these all start with three-letters and an underscore
  files_to_save <- x[grepl('^[a-z]{3}_', names(x))]
  file_names <- names(files_to_save)
  
  ### a for loop approach
  # for (file in 1:length(files_to_save)){
  # 
  #   write.csv(x = files_to_save[file], file = paste0(file_names[file], '.csv'), row.names=FALSE)
  #   
  # }
  
  lapply(file_names, function(file){
    #browser()
    ## the name of the particular table in the list
    file_name <- names(files_to_save[file])
    write.csv(x = files_to_save[[file]], file = paste0(path, file_name, '.csv'), row.names=FALSE)
    
  })
  
})
}
