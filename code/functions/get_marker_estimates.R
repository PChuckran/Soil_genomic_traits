get_marker_estimates <- function(refresh=FALSE){
  if(refresh == TRUE){
    # Read in soil chemistry data
    markers <- read.csv("data/input/derived/combined_16s_sizes.csv", header = F)
    colnames(markers) <- c("sample", "AGS_marker", "GC_marker", 
                           "Coding_density_marker", "num_taxa_marker")

    
    markers <- markers %>%
      mutate(perc_identity = str_extract(sample, "[0-9]+$"), 
             soilID = str_extract(sample, "^[A-Z]+_[0-9]+-[A-Z]"),
             siteID = str_extract(sample, "^[A-Z]+"))
    
    markers <- markers %>%
      filter(perc_identity == 99) %>%
      select(soilID, AGS_marker, GC_marker, Coding_density_marker, num_taxa_marker)%>%
      group_by(soilID)%>%
      summarise_all(funs(mean), na.rm = T)
    
    return(markers)
    write.csv(markers, "data/input/derived/rRNA_16S_summary.csv")
  }
  else{
    markers <- read.csv("data/input/derived/rRNA_16S_summary.csv")%>%
      select(-X)
    return(markers)
  }
}