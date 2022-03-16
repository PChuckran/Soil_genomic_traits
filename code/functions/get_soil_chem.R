get_soil_chem <- function(refresh=FALSE){
  if(refresh == TRUE){
    # Read in soil chemistry data
    soil_chem <- read.csv("data/input/raw/sls_soilChemistry.csv") 
    
    # Make the same ID which has plot and horizon
    soil_chem$soilID <-str_extract(soil_chem$sampleID, regex("[a-z]+_[0-9]+-[a-z]", TRUE))
    
    soil_chem$collectDate <- as.Date(soil_chem$collectDate)
    # Summary of those files
    soil_chem_summary <- soil_chem %>%
      mutate(year = str_extract(collectDate, regex("^[0-9]+", TRUE)))%>%
      group_by(soilID, domainID, plotID, collectDate)%>%
      summarise(mean_cn = mean(CNratio, na.rm = T),
                mean_OrgC = mean(organicCPercent, na.rm = T),
                mean_Org13C = mean(organicd13C, na.rm = T),
                mean_15n = mean(d15N, na.rm = T),
                mean_N = mean(nitrogenPercent, na.rm = T))
    
    #soil_chem_summary$year <- as.numeric(soil_chem_summary$year)
    
    write.csv(soil_chem_summary, file = "data/input/derived/soil_chem.csv")
    return(soil_chem_summary)

  }
  else{
    soil_chem_summary <- read.csv("data/input/derived/soil_chem.csv")%>%
      select(-X)
    return(soil_chem_summary)
  }
}