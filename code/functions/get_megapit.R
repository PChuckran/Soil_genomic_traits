get_megapit <- function(refresh=FALSE){
  if(refresh == TRUE){
    
    megapit_bgc <- read.csv("data/input/raw/mgp_perbiogeosample.csv") 
    megapit_bgc <- megapit_bgc %>%
      # Make a new horizon which matches the rest of the data: O/M
      mutate(horizon = ifelse(str_extract(horizonName, "^[A-Z]") == "O", "O", "M")) %>%
      # Cores never go below ~35, so we'll only average values above that
      filter(biogeoTopDepth <25 )
    megapit_bgc_sub <- megapit_bgc %>%
      select(-c(setDate, collectDate, uid, domainID, pitNamedLocation, pitID, horizonID, horizonName,
                biogeoID, biogeoSampleType, biogeoHorizonProportion, laboratoryName, biogeoTopDepth,
                biogeoCenterDepth, biogeoBottomDepth, laboratoryName, labProjID, remarks, publicationDate, resist,
                release)) %>%
      group_by(siteID, horizon) %>%
      summarise_all(funs(mean), na.rm = T)
    
    colnames(megapit_bgc_sub) <- paste("mega", colnames(megapit_bgc_sub), sep = "_")
    megapit_bgc_sub <- rename(megapit_bgc_sub, siteID = mega_siteID)
    megapit_bgc_sub <- rename(megapit_bgc_sub, horizon = mega_horizon)
    megapit_bgc_sub[megapit_bgc_sub == "NaN"] <- NA
    megapit_bgc_sub <- megapit_bgc_sub %>% mutate(mega_cton = mega_carbonTot/mega_nitrogenTot)
    megapit_bgc_sub[megapit_bgc_sub == "Inf"] <- NA
    write.csv(megapit_bgc_sub, file = "data/input/derived/megapit_sum.csv")

    return(megapit_bgc_sub)
  }
  else{
    megapit_bgc_sub <- read.csv("data/input/derived/soil_plfa.csv")%>%
      select(-X)
    return(megapit_bgc_sub)
  }
}
