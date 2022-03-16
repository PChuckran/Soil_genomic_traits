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

full <- read.csv("data/input/derived/full_df.csv")
mcc_16s <- read.csv("data/input/Raw_bioinfo/16s/mcc_soilSeqVariantMetadata_16S.csv")

#randomly sample 5
set.seed(316)
mcc_16s_new <- mcc_16s %>%
  filter(plotID %in% full$plotID) %>%
  filter(sequenceCountQF == "OK") %>%
  group_by(plotID)%>%
  sample_n(5, replace = T)
#remove duplicates
mcc_16s_new = mcc_16s_new[!duplicated(mcc_16s_new$downloadFileUrl),]

write.csv(mcc_16s_new, "data/input/Raw_bioinfo/16s/mcc_soilSeqVariantMetadata_16S_subset.csv")
write_delim(as.data.frame(mcc_16s_new$downloadFileUrl), "data/input/Raw_bioinfo/16s/mcc_soilSeqVariantMetadata_16S_urls.txt")


