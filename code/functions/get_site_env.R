get_site_env <- function(){
  site_meta <- read.csv('data/input/meta/NEON_Field_Site_Metadata_20201204.csv')
  site_env <- site_meta %>%
    rename(siteID = field_site_id) %>%
    select(field_latitude, field_longitude, field_mean_elevation_m, field_mean_annual_precipitation_mm, 
           field_mean_annual_temperature_C, field_avg_number_of_green_days, siteID)
  return(site_env)
  
}