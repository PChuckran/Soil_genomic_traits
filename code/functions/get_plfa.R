get_plfa <- function(refresh=FALSE){
  if(refresh == TRUE){
    sls_plfa <- read.csv('data/input/derived/plfa.csv')
    sls_plfa$collectDate <- as.Date(as.character(sls_plfa$collectDate))
    sls_plfa$horizon <- str_match(sls_plfa$biomassID, regex("^[A-Z]+_[0-9]+-(\\w)"))[,2]
    sls_plfa <- sls_plfa %>% 
      #mutate(year = as.numeric(str_extract(collectDate, regex("^[0-9]+", TRUE))))%>%
      select(plotID, horizon, lipids_total, lipid_fungi_perc, 
             lipid_gram_pos_perc, lipid_gram_neg_perc, lipid_FtoB)%>%
      group_by(plotID, horizon) %>% 
      summarise_all(funs(mean), na.rm = T)
    
    write.csv(sls_plfa, file = "data/input/derived/soil_plfa.csv")
    return(sls_plfa)
  }
  else{
    soil_plfa.csv <- read.csv("data/input/derived/soil_plfa.csv")%>%
      select(-X)
    return(soil_plfa.csv)
  }
}


