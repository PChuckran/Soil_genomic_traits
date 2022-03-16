packages <- c('tidyverse', 'jsonlite','neonUtilities', 'httr')

# trigger for whether to retrieve data now, or to load data from existing file directory
get_data = FALSE

## method for loading all packages, and installing ones that are missing if this is launched on a different rig
lapply(packages, function(x){
  #browser()
  ## if the required package is present (TRUE), load it up
  if (require(x, character.only = TRUE)==TRUE){
    library(x, character.only=TRUE)
  } else {
    install.packages(x)
  }
})

#### load custom functions, if any
file_sources <- list.files(file.path('code', 'functions'), pattern="*.R$",
                            full.names=TRUE, ignore.case=TRUE)

file_sources <- base::tryCatch(invisible(sapply(file_sources, source, .GlobalEnv)),
                                error = function(e) {
                                  base::stop()
                                })

NEON_TOKEN = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJwZmMyNUBuYXUuZWR1Iiwic2NvcGUiOiJyYXRlOnB1YmxpYyIsImlzcyI6Imh0dHBzOi8vZGF0YS5uZW9uc2NpZW5jZS5vcmcvIiwiZXhwIjoxNzcxODgzNTE4LCJpYXQiOjE2MTQyMDM1MTgsImVtYWlsIjoicGZjMjVAbmF1LmVkdSJ9.p9fxkL6DbB6UOjc3YR2xw4nnsFilURlh3_-bme5gLu3IUbMVVUdSzrs_oXhgLNsgq9-p2G_lU15-dHwyofS2bQ"

## split workflow, if you want to retrieve all data again (because time has passed and there may be more available) run get_data==TRUE
### else, don't waste time on getting data and just load from the local folder
if (get_data == TRUE){
  
  ## list of soil related data products
  soil_products <- c(
    #dpid1 <- 'DP1.10086.001', ### physical, isotope, and chemical properties
    #dpid2 <- 'DP1.10107.001', ### soil metagenomic sequences
    #dpid3 <- 'DP1.10109.001', ### soil microbial group abundances 
    #dpid4 <- 'DP1.10104.001', ### soil microbial biomass
    #dpid5 <- 'DP1.10108.001', ### soil marker genes
    #dpid6 <- 'DP1.10023.001', ### clip harvest
    #dpid7 <- 'DP1.10033.001', ### litter chem
    #dpid8 <- 'DP1.10067.001',### rootbiomass
    #'DP1.10111.001' ### Site Management data
    #'DP1.10047.001' ### Initial plot characterizations
    #'DP1.00096.001' ### Megapit
    
  )
  
  last_run <- Sys.Date()
  
  ## store date of last run, so local data files can be retrieved
  write.csv(last_run, 'last_run_getdata.csv', row.names=FALSE)
  
  ## retrieve data from API
  data_list <- lapply(soil_products[1:length(soil_products)], FUN = function(x){
    data <- loadByProduct(dpID = x, site = 'all', token=NEON_TOKEN, check.size=FALSE, package = "expanded")
    return(data)
  })
  
  ## write list of data.frames to folder as a compressed .rds file for Cody to use in R
  ### these are like Python "pickle" files (i.e. they preserve the native format of the object, rather than coercing to something like a CSV)
  saveRDS(data_list, paste0('data/input/megapit_', last_run, '.rds'))
  
  ## store CSV files for Pete and non-R users
  ### input object is a list of lists (nested list)...i.e. data_list object
  write_files_from_list(data_list)
  
} else if (get_data == FALSE){
  
  last_run <- read.csv('last_run_getdata.csv')
  last_run_date <- as.character(last_run[1,1])
  
  data_list <- readRDS(paste0('data/input/soil_products_', last_run_date, '.rds'))
  
}


