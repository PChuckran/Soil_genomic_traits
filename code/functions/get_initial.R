get_soil_initial <- function(refresh=FALSE){
  if(refresh == TRUE){
    # Read in plot chemistry data
    particle_size <- read.csv("data/input/raw/spc_particlesize.csv") 
    bulk <- read.csv("data/input/raw/spc_bulkdensity.csv") %>%
      select(-c(uid, publicationDate, release, dataQF, laboratoryName, nrcsDescriptionID)) 
    biogeo <- read.csv("data/input/raw/spc_biogeochem.csv") %>%
      select(-c(uid, publicationDate, release, dataQF, laboratoryName, nrcsDescriptionID,
                biogeoIDnrcs)) 
    particle_biogeo_bulk <-full_join(particle_size, biogeo)%>%
      full_join(., bulk)
    # No cores go much beyond 30, so trim all horizons that start after that
    particle_biogeo_bulk <- particle_biogeo_bulk %>% filter(biogeoTopDepth < 30)
    # Get a simplified horizon name (O/M) similar to other files
    particle_biogeo_bulk$horizon <- ifelse(grepl("O", particle_biogeo_bulk$horizonName), "O", "M")
    
    init_summary <- particle_biogeo_bulk %>%
      select(plotID, horizon, sandTotal, siltTotal, clayTotal, bulkDensThirdBar,
             bulkDensOvenDry, waterRetentionThirdBar, caco3Conc, gypsumConc, phCacl2, phH2o, bSatx, brSatx,
             caSatx, clSatx, co3Satx, ecSatp, flSatx, waterSatx, hco3Sx, kSatx, mgSatx,
             naSatx, no2Satx, no3Satx, pSatx, phSp, so4Satx, alKcl, feKcl, mnKcl, alSatCecd33,
             carbonTot, nitrogenTot, ctonRatio, estimatedOC, alOxalate, feOxalate,
             mnOxalate, pOxalate, siOxalate)%>%
      group_by(plotID, horizon) %>%
      summarise_all(funs(mean), na.rm = T)

    write_csv(init_summary, "data/input/derived/init_summary.csv")
    return(init_summary)
    
  }
  else{
    init_summary <- read.csv("data/input/derived/init_summary.csv")
    return(init_summary)
  }
}
