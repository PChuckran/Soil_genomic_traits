## Downloads and calculates annual moisture averages at the site level

packages <- c('tidyverse', 'jsonlite','neonUtilities', 'httr', 'lubridate')

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


## Going to loop through the collections and get just the samples we need

samples <- read.csv("code/getdata/samples_to_get.csv") %>%
  select(siteID) 

samples$siteID <- as.character(samples$siteID)

# remove duplicates, provides list of site names 

samples <- samples[!duplicated(samples$siteID),]


get_moisture_fun <- function(site_to_get, toke){
  data <- loadByProduct(dpID = 'DP1.00094.001', site = site_to_get, token=toke, check.size=F, package = "basic",
                        startdate = '2018-01-01', enddate = '2019-12-31',
                        timeIndex = 30)
  
  sws_30 <- data$SWS_30_minute
  
  sws_30$endDateTime <- as.Date(sws_30$endDateTime)
  sws_30$doy <- as.character(strftime(sws_30$endDateTime, format = "%j"))
  sws_30$year <- format(as.Date(sws_30$endDateTime, format="%Y-%m-%d"),"%Y")
  
  sws_30_sumDay <- sws_30 %>% 
    select(-c(startDateTime)) %>%
    group_by(verticalPosition, doy, year, endDateTime) %>%
    summarise_all(funs(mean), na.rm = T)
  
  sws_30_Ann <- sws_30_sumDay %>% 
    filter(VSWCMean != "NaN")%>%
    #for each year, gather summary statistics 
    group_by(verticalPosition, year) %>%
    summarise(
      VSWC_ann_mean = mean(VSWCMean, na.rm = T),
      VSWC_ann_var = sd(VSWCMean, na.rm = T),
      VSWC_max = max(VSWCMean, na.rm = T),
      VSWC_min = min(VSWCMean, na.rm = T),
      VSWC_median = median(VSWCMean, na.rm = T),
      VSWC_range = VSWC_max - VSWC_min,
      VSIC_ann_mean = mean(VSICMean, na.rm = T),
      VSIC_ann_var = sd(VSICMean, na.rm = T),
      VSIC_max = max(VSICMean, na.rm = T),
      VSIC_min = min(VSICMean, na.rm = T),
      VSIC_median = median(VSICMean, na.rm = T),
      VSIC_range = VSIC_max - VSIC_min,
      days = n_distinct(doy)
    ) %>% 
    # Average across years
    select(-year)%>%
    group_by(verticalPosition)%>% 
    summarise_all(funs(mean), na.rm = T) %>%
    mutate(siteID = "BARR")
  # for this script, return just the annutal sws_30
  return(sws_30_Ann)
}

moisture_summary <- data.frame()

for(site_rep in samples){
  print(site_rep)
  moisture_summary <- rbind(moisture_summary, get_moisture_fun(site_to_get = site_rep, toke = NEON_TOKEN))
}


# function to convert positions to horizon names
parse_postition <- function(vert){
  if(vert == "501"){
    return("O")
  }
  if(vert == "502"){
    return("M")
  }
  else{
    return("NA")
  }
}

# apply function to convert vertical position to horizon 
moisture_summary$horizon = sapply(moisture_summary$verticalPosition, parse_postition)

  
write.csv(moisture_summary, "data/input/derived/moisture_yearly.csv")
  
