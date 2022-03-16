# Making the big DF to use in analyses

packages <- c('tidyverse', 'jsonlite','neonUtilities', 'httr', 'ggpmisc', 'GGally',
              'randomForest', 'sp', 'rgeos')

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

# Wrote functions to gather, format, and summarize the data. 
## Metagenomes will be refreshed every time since they are often updated
## but others will only if refresh = TRUE. 
## Otherwise it reads from a csv in the derived folder
metagenomes <- get_metagenomic()%>%
  mutate(collectDate = as.Date(collectDate))
litter <- get_litter(refresh = F)
roots <- get_roots(refresh=F)
pH <- get_pH(refresh=F)%>%
  mutate(collectDate = as.Date(collectDate))
site_env <-get_site_env()
soil_chem <- get_soil_chem(refresh = F)%>%
  mutate(collectDate = as.Date(collectDate))
initial_summary <- get_soil_initial(refresh = F)
initial_summary[initial_summary == "NaN"] <- NA
spei <- read.csv("data/input/derived/spei.csv")%>%
  select(-c(X, collectDate))
moisture <- read.csv("data/input/derived/moisture_yearly.csv")
fungal_to_bacterial <- get_ftob(refresh = F)
plfa <- get_plfa(refresh = T)
mega <- get_megapit(refresh = T)



# I know this isn't elegant but.... meh. I wanted to go through step by step
full_df <- left_join(metagenomes, pH)%>%
  left_join(., soil_chem) %>%
  left_join(., roots)%>%
  left_join(., site_env)%>%
  left_join(., litter, by = c('siteID', 'plotID')) %>%
  left_join(., initial_summary) %>%
  mutate(soil_OrgCtoN = mean_OrgC/mean_N)%>%
  left_join(., spei)%>%
  left_join(., moisture)%>%
  left_join(., fungal_to_bacterial)%>%
  left_join(., plfa) %>%
  left_join(., mega)



# There are a lot of gaps in the plot chemistry data, so I want to fill those
## in with nearby plots. The following finds the nearest plot for each plot
## and creates a key
coors <- read.csv("data/input/raw/sls_soilCoreCollection.csv") %>%
  mutate(long = decimalLongitude,
         lat = decimalLatitude)%>%
  select(plotID, lat, long)

plots <- coors[!duplicated(coors$plotID), ]
sp.plots <- plots
# create spatial object
coordinates(sp.plots) <- ~ long+lat 
# get distances
d <- gDistance(sp.plots, byid=T)
# get the 1st, 2nd, and 3rd closest site 
closest_site <- apply(d, 1, function(x) order(x, decreasing=F)[2])
closest_dist <- apply(d, 1, function(x) sort(x, decreasing=F)[2])
sec_closest_site <- apply(d, 1, function(x) order(x, decreasing=F)[3])
sec_closest_dist <- apply(d, 1, function(x) sort(x, decreasing=F)[3])
third_closest_site <- apply(d, 1, function(x) order(x, decreasing=F)[4])
third_closest_dist <- apply(d, 1, function(x) sort(x, decreasing=F)[4])
fourth_closest_site <- apply(d, 1, function(x) order(x, decreasing=F)[5])
fourth_closest_dist <- apply(d, 1, function(x) sort(x, decreasing=F)[5])

# combine w original data
newdata <- cbind(plots, plots[closest_site,], closest_dist, plots[sec_closest_site,], 
                 sec_closest_dist, plots[third_closest_site,], third_closest_dist, 
                 plots[fourth_closest_site,], fourth_closest_dist)
# rename
colnames(newdata) <- c("refPlotID", "lat", "long", "closest_plot", "closest_lat", "closest_long",
                       "closest_dist", "second_plot", "second_lat", "second_long",
                       "second_dist", "third_plot", "third_lat", "third_long", "third_dist",
                       "fourth_plot", "fourth_lat", "fourth_long", "fourth_dist")
# select relevant info and make into a longer look-up key
site_key <- newdata %>% 
  select(refPlotID, closest_plot, second_plot, third_plot, fourth_plot)%>%
  pivot_longer(., cols = c(closest_plot, second_plot, third_plot), 
               names_to = "plot_dist", values_to = "plotID")

# create an alternate summary df which creates alt values for each reference plot
alt_plot_info <- left_join(site_key, initial_summary)%>%
  filter(is.na(horizon) == F)
alt_plot_info_summary <- alt_plot_info %>%
  group_by(refPlotID, horizon)%>%
  summarise_all(funs(mean), na.rm = T) %>%
  select(-c(plot_dist, plotID))%>%
  mutate(plotID = refPlotID)

alt_plot_info_summary[alt_plot_info_summary == "NaN"] <- NA


filled_df <- coalesce_join(full_df, alt_plot_info_summary, by = c("plotID", "horizon"))

soil_chem_site <- soil_chem %>%
  mutate(siteID = str_extract(plotID, "[A-Z]+"),
         horizon = str_extract(soilID, "[A-Z]$")) %>%
  group_by(siteID, horizon )%>%
  summarise_all(funs(mean), na.rm = T)%>%
  mutate(soil_OrgCtoN = mean_OrgC/mean_N) %>%
  select(siteID, horizon, mean_cn, soil_OrgCtoN)

filled_df <- coalesce_join(filled_df, soil_chem_site, by = c("siteID", "horizon"))

#get some good plot specific coordinates on this bad-bad-boy
filled_df <-left_join(filled_df, plots)%>%
  select(-c(field_longitude, field_latitude))

write_csv(filled_df, "data/input/derived/full_df.csv")
filled_df<-read_csv("data/input/derived/full_df.csv")


  
