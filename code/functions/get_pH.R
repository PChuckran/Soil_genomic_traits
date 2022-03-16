get_pH <- function(refresh=FALSE){
  if(refresh == TRUE){
    sls_ph <- read.csv('data/input/raw/sls_soilpH.csv')
    sls_ph$collectDate <- as.Date(sls_ph$collectDate)
    plot_ph <- sls_ph %>% 
      mutate(year = as.numeric(str_extract(collectDate, regex("^[0-9]+", TRUE))))%>%
      group_by(plotID, horizon, year, collectDate) %>% 
      summarise(mean_pH_CaCl = mean(soilInCaClpH, na.rm=TRUE),
                mean_pH_water = mean(soilInWaterpH, na.rm=TRUE),
                sd_pH_CaCl = sd(soilInCaClpH, na.rm=TRUE),
                sd_pH_water = sd(soilInWaterpH, na.rm=TRUE))
    
    plot_ph$collectDate <- as.Date(plot_ph$collectDate)
    
    write.csv(plot_ph, file = "data/input/derived/soil_pH.csv")
    return(plot_ph)
  }
  else{
    plot_ph <- read.csv("data/input/derived/soil_pH.csv")%>%
      select(-X)
    return(plot_ph)
  }
}