get_litter <- function(refresh=FALSE){
  if(refresh == TRUE){
    #read in data
    ltr_mass <- read.csv("data/input/raw/ltr_massdata.csv") 
    ltr_chem <- read.csv("data/input/raw/ltr_litterCarbonNitrogen.csv")
    ltr_lignin <- read.csv("data/input/raw/ltr_litterLignin.csv")
    ltr_pertrap <- read.csv("data/input/raw/ltr_pertrap.csv")
    
    #Extract date
    ltr_mass <-ltr_mass %>%
      mutate(year = str_extract(collectDate, regex("^[0-9]+", TRUE)),
             month = str_match(collectDate, "^[0-9]+-([0-9]+)")[,2],
             day = str_match(collectDate, "^[0-9]+-([0-9]+)-([0-9]+)")[,3]
      )
    #Get trap info
    ltr_mass_pertrap <- ltr_pertrap %>%
      select(trapID, trapSize, trapType) %>%
      dplyr::left_join(ltr_mass, .)
    
    
    leaf_plot_sum <- ltr_mass_pertrap %>%
      filter(functionalGroup %in% c('Leaves', 'Needles')) %>% 
      ## Sum up the leaf and needle mass for each trap
      group_by( trapID, domainID, siteID, plotID, year, month, trapSize, collectDate, setDate) %>%
      summarise(leaf_needle_biomass = sum(dryMass, na.rm = TRUE)) %>%
      mutate(numDays = as.numeric(as.Date(collectDate) - as.Date(setDate)),
             g_leaf_per_m_per_d = leaf_needle_biomass / trapSize / numDays) %>%
      group_by(domainID, siteID, plotID)%>%
      summarise(ltr_sum_leafdryMass = sum(leaf_needle_biomass),
                ltr_mean_g_leaf_per_m_per_d = mean(g_leaf_per_m_per_d, na.rm=TRUE),
                ltr_sd_g_leaf_per_m_per_d = sd(g_leaf_per_m_per_d, na.rm=TRUE))
    
    ltr_plot_sum <- ltr_mass_pertrap %>%
      group_by( trapID, domainID, siteID, plotID, year, month, trapSize, collectDate, setDate) %>%
      # sum litter mass for each trap
      summarise(pertrap_biomass = sum(dryMass, na.rm = TRUE))%>%
      # express this as a rate
      mutate(ltr_numDays = as.numeric(as.Date(collectDate) - as.Date(setDate)),
             ltr_g_litt_per_m_per_d = pertrap_biomass / trapSize / ltr_numDays)%>%
      ungroup()%>%
      # Get a mean for the entire plot for a given year 
      group_by(domainID, siteID, plotID)%>%
      # CF: added first collection date and year of collection 
      summarise(
        ltr_sum_g_litt_per_m_per_d = sum(ltr_g_litt_per_m_per_d, na.rm = T),
                ltr_mean_g_litt_per_m_per_d= mean(ltr_g_litt_per_m_per_d, na.rm=TRUE),
                ltr_sd_g_litt_per_m_per_d= sd(ltr_g_litt_per_m_per_d, na.rm=TRUE))

    
    

    # Litter chemistry by plot
    ltr_chem_colllection_plot <- ltr_chem %>%
      mutate(year = str_extract(collectDate, regex("^[0-9]+", TRUE)),
             month = str_match(collectDate, "^[0-9]+-([0-9]+)")[,2],
             day = str_match(collectDate, "^[0-9]+-([0-9]+)-([0-9]+)")[,3]
      )%>%
      group_by(domainID, siteID, plotID) %>% 
      summarise(ltr_mean_carbon = mean(carbonPercent, na.rm=TRUE),
                ltr_mean_nitrogen = mean(nitrogenPercent, na.rm=TRUE),
                ltr_mean_CNratio = mean(CNratio, na.rm=TRUE))
    
    ltr_cn_collection_plot <- dplyr::left_join(ltr_chem_colllection_plot, leaf_plot_sum)%>% 
      mutate(ltr_gCperDay = (ltr_mean_carbon / 100) * ltr_mean_g_leaf_per_m_per_d,
             ltr_gNperDay = (ltr_mean_nitrogen / 100) * ltr_mean_g_leaf_per_m_per_d)
    
    
    ### lignin summary
    ltr_lignin_plot <- ltr_lignin %>% 
      group_by(domainID, plotID) %>% 
      summarise(ltr_mean_lignin = mean(ligninPercent, na.rm=TRUE),
                ltr_mean_cellulose = mean(cellulosePercent, na.rm=TRUE))
    
    ## join litterfall rates with litterfall chemistry
    ltr_plot_join <- dplyr::left_join(ltr_cn_collection_plot, ltr_lignin_plot)%>%
                     dplyr::left_join(ltr_plot_sum, .)
    
    
    ltr_plot_join<-ltr_plot_join%>%
      mutate(ltr_CN_inputs = ltr_gCperDay/ltr_gNperDay)%>%
      group_by(siteID)%>%
      summarise(ltr_site_mean_gCperDay = mean(ltr_gCperDay, na.rm = T),
                ltr_site_se_gCperDay = sd(ltr_gCperDay, na.rm = T)/sqrt(length(ltr_gCperDay)),
                ltr_site_mean_gNperDay = mean(ltr_gNperDay, na.rm = T),
                ltr_site_se_gNperDay = sd(ltr_gNperDay, na.rm = T)/sqrt(length(ltr_gNperDay)),
                ltr_site_mean_CN_inputs = mean(ltr_CN_inputs, na.rm = T),
                ltr_site_se_CN_inputs = sd(ltr_CN_inputs, na.rm = T)/sqrt(length(ltr_CN_inputs)),
                ltr_site_mean_litterRate = mean(ltr_mean_g_litt_per_m_per_d, na.rm = TRUE),
                ltr_site_se_litterRate = sd(ltr_mean_g_litt_per_m_per_d, na.rm = T)/sqrt(length(ltr_mean_g_litt_per_m_per_d))) %>%
      dplyr::left_join(ltr_plot_join, . )
    
    write.csv(ltr_plot_join, "data/input/derived/litter_data.csv")
    return(ltr_plot_join)
  }
  else{
    ltr_plot_join <- read.csv("data/input/derived/litter_data.csv")%>%
      select(-X)
    return(ltr_plot_join)
  }
}
