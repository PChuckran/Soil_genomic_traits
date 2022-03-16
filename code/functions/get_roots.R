get_roots <- function(refresh=FALSE){
  if(refresh == TRUE){
    
    # Read in soil chemistry data
    roots_mass <- read.csv("data/input/raw/bbc_rootmass.csv") 
    root_chem <- read.csv("data/input/raw/bbc_rootChemistry.csv") 
    
    # Make a column for root min and max
    root_chem <- root_chem %>%
      filter(cnPercentQF == "OK")%>%
      select(siteID, plotID, collectDate, d15N, d13C, nitrogenPercent, carbonPercent, CNratio, poolSampleID) %>%
      mutate(year = str_extract(collectDate, regex("^[0-9]+", TRUE)),
             month = str_match(collectDate, "^[0-9]+-([0-9]+)")[,2],
             day = str_match(collectDate, "^[0-9]+-([0-9]+)-([0-9]+)")[,3],
             root_min = str_match(poolSampleID, "([0-9]+)-([0-9]+).POOL$")[,2],
             root_max = str_match(poolSampleID, "([0-9]+)-([0-9]+).POOL$")[,3])
    
    # as numeric 
    root_chem$root_min <- as.numeric(str_replace(root_chem$root_min, "05", ".05"))
    root_chem$root_max <- as.numeric(str_replace(root_chem$root_max, "05", ".05"))
    # add size class, > 2mm is coarse 
    root_chem$rootSizeClass <- ifelse(root_chem$root_min < 2, "Fine", "Coarse")
    # isolate unique coreID to match w mass data
    root_chem$coreID <- as.character(str_match(root_chem$poolSampleID, "BBC.[A-Z]+[0-9]+.[0-9]+"))
    
    ## get mean between  reps
    root_chem <- root_chem %>%
      group_by(siteID, plotID, coreID, collectDate, year, month, day, poolSampleID,
               root_min, root_max, rootSizeClass) %>%
      summarise_all(funs(mean))
    
    # get root min and max, size class, and coreID for root mass data
    roots_mass <- roots_mass %>%
      mutate(root_min = str_extract(sizeCategory, regex("^[0-9]+")),
             root_max = str_extract(sizeCategory, regex("[0-9]+$")))
    roots_mass$root_min <- as.numeric(str_replace(roots_mass$root_min, "05", ".05"))
    roots_mass$root_max <- as.numeric(str_replace(roots_mass$root_max, "05", ".05"))
    roots_mass$rootSizeClass <- ifelse(roots_mass$root_min < 2, "Fine", "Coarse")
    roots_mass$coreID <- as.character(str_match(roots_mass$subsampleID, "BBC.[A-Z]+[0-9]+.[0-9]+"))
    
    # mass sample-columns were pooled for the chem analyses, so I average them here
    root_mass_sub <- roots_mass %>%
      group_by(siteID, plotID, collectDate, rootSizeClass, rootStatus, coreID,
               root_min, root_max) %>%
      summarize(dryMass_mean = mean(dryMass, na.rm = T),
                cores = length(dryMass))
    
    # combine mass and chem
    root_mass_w_chem<- root_chem %>%
      ungroup(.)%>%
      select(-c(collectDate, day, poolSampleID))%>%
      full_join(root_mass_sub, .)
    
    # In 2019 roots were no longer seperated into live/dead. Lucky for us
    ## all our MGs are before 2019 which means we don't have to consider that change
    ## to be safe, Lets remove collections before 2019 to start, then summarize
    
    root_mass_w_chem_wide <- root_mass_w_chem %>%
      filter(year < 2019) %>%
      # also calculate the total root C and N mass
      pivot_wider(., names_from = rootStatus, names_prefix = "mass", values_from = dryMass_mean) %>%
      mutate(liveC = masslive * carbonPercent/100,
             liveN = masslive * nitrogenPercent/100)
    
    # Get a plot mean for each size
    root_mass_w_chem_wide_means <- root_mass_w_chem_wide %>%
      ungroup(.)%>%
      select(-c(coreID, cores))%>%
      group_by( siteID, plotID, collectDate, rootSizeClass, root_min, root_max, year, month) %>%
      summarise_all(funs(mean), na.rm = T)
    
    # For each size class, average the chem data and get sum of mass (+chem) data 
    roots_byClass <- root_mass_w_chem_wide_means %>%
      ungroup(.)%>%
      select(-c(root_min, root_max))%>%
      group_by(siteID, plotID, rootSizeClass, collectDate, year, month)%>%
      summarise(
        root_d15N_m = weighted.mean(d15N, w =masslive, na.rm = T),
        root_d13C_m = weighted.mean(d13C, w =masslive,na.rm = T),
        root_NPerc = weighted.mean(nitrogenPercent, w =masslive, na.rm = T),
        root_CPerc = weighted.mean(carbonPercent, w =masslive, na.rm = T),
        root_CN = weighted.mean(CNratio, w =masslive, na.rm = T),
        root_live_mass = sum(masslive, na.rm = T),
        root_dead_mass = sum(massdead, na.rm = T),
        root_liveC = sum(liveC, na.rm = T),
        root_liveN = sum(liveN, na.rm = T)
      ) %>%
      pivot_wider(., names_from = rootSizeClass, values_from = c(root_d15N_m, root_d13C_m,
                                                                 root_NPerc, root_CPerc,
                                                                 root_CN, root_live_mass, 
                                                                 root_dead_mass,
                                                                 root_liveN, root_liveC))
    
    roots_all <- root_mass_w_chem_wide_means %>%
      ungroup(.)%>%
      select(-c(root_min, root_max, rootSizeClass))%>%
      group_by(siteID, plotID, collectDate, year, month)%>%
      summarise(
        root_d15N_m = weighted.mean(d15N, w =masslive, na.rm = T),
        root_d13C_m = weighted.mean(d13C, w =masslive, na.rm = T),
        root_NPerc = weighted.mean(nitrogenPercent, w =masslive, na.rm = T),
        root_CPerc = weighted.mean(carbonPercent, w =masslive, na.rm = T),
        root_CN = weighted.mean(CNratio, w =masslive, na.rm = T),
        root_live_mass = sum(masslive, na.rm = T),
        root_dead_mass = sum(massdead, na.rm = T),
        root_liveC = sum(liveC, na.rm = T),
        root_liveN = sum(liveN, na.rm = T)
      ) %>%
      left_join(., roots_byClass) %>%
      ungroup(.)%>%
      select(-c(collectDate, month, year))%>%
      group_by(siteID, plotID)%>%
      summarise_all(funs(mean), na.rm = T)
    
    
    write.csv(roots_all, "data/input/derived/root_data.csv")
    return(roots_all)
  }
  else{
    roots_all <- read.csv("data/input/derived/root_data.csv")%>%
      select(-X)
    return(roots_all)
  }
}
