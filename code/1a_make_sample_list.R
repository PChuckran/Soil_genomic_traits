# Generate a list of samples for which metagenomes exist.
## Important for gathering data that is computationally expensive 

packages <- c('tidyverse')

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
file_sources <- list.files(file.path( 'code/functions'), pattern="*.R$",
                           full.names=TRUE, ignore.case=TRUE)

file_sources <- base::tryCatch(invisible(sapply(file_sources, source, .GlobalEnv)),
                               error = function(e) {
                                 base::stop()
                               })

# Get metagenomes
metagenomes <- get_metagenomic()%>%
  mutate(collectDate = as.Date(collectDate))

# Get coor info
coors <- read.csv("data/input/raw/sls_soilCoreCollection.csv") 

# Rename lat-long
coors <- coors%>%
  mutate(lat = decimalLatitude, 
         long = decimalLongitude)

# Merge to create sample list which contains sampleID and description
sample_list <- coors%>%
  select(lat, long, plotID, horizon) %>%
  unique(.)%>%
  left_join(metagenomes, .) %>%
  select(lat, long, siteID, plotID, horizon, dnaSampleID, collectDate)

write.csv(sample_list, file = "code/samples_to_get.csv")
