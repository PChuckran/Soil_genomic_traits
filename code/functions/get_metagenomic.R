## Noticed we were doing these same steps every time, why not a funnction to read them in
get_metagenomic <- function(){
  # Read in traits data
  genomic_traits <- read.csv("data/input/derived/genomic_traits.csv") %>%
    select(-X)
  
  # Red rid of some tails on the labID that could be problematic 
  genomic_traits$ID <- str_replace(genomic_traits$ID, "_mms", "")
  genomic_traits$ID <- str_replace(genomic_traits$ID, "^P", "BMI_P")
  
  # Read in sequencing info
  seq_Info <- read.csv("data/input/raw/mms_metagenomeSequencing.csv") 
  seq_Info <- seq_Info%>%
    rename(ID = internalLabID)
  
  # The files starting with "18S" appeared as "BMI_18S" in the sequence files. This should fix it
  seq_Info$ID <- str_replace(seq_Info$ID, "^18S", "BMI_18S")
  
  # Combine those data by the unique Lab ID
  genomic_traits <- left_join(genomic_traits, seq_Info)
  
  # Make a new ID that captures the plot and horizon
  genomic_traits$soilID <-str_extract(genomic_traits$dnaSampleID, regex("[a-z]+_[0-9]+-[a-z]", TRUE))
  genomic_traits$plotID <-str_extract(genomic_traits$dnaSampleID, regex("[a-z]+_[0-9]+", TRUE))
  genomic_traits$year <- as.character(str_extract(genomic_traits$collectDate, regex("^[0-9]+", TRUE)))
  genomic_traits <- genomic_traits[!duplicated(genomic_traits$ID), ]
  
  coors <- read.csv("data/input/raw/sls_soilCoreCollection.csv") 
  
  coors_light <- coors %>%
    select(plotID, nlcdClass)
  coors_light <- coors_light[!duplicated(coors_light$plotID), ]
  genomic_traits_w_coors <- left_join(genomic_traits, coors_light)%>%
    mutate(horizon = str_extract(soilID, regex("[a-z]$", TRUE))) %>%
    mutate(Nuc_CN = ((GC_bact*(9/8))+((1-GC_bact)*(10/7)))) %>%
    mutate(Genome_CN = AGS*Nuc_CN,
           Genome_C = AGS*((GC_bact*(9))+((1-GC_bact)*(10))),
           Genome_N = AGS*((GC_bact*(8))+((1-GC_bact)*(7))))
  genomic_traits_w_coors$year <- as.numeric(genomic_traits_w_coors$year)
  genomic_traits_w_coors <- genomic_traits_w_coors%>%
    select(-c(uid, domainID, namedLocation, laboratoryName,
              sequencingFacilityID, processedDate, 
              barcodeSequence, instrument_model, sequencingMethod, investigation_type,
              sequencingConcentration, labPrepMethod, sequencingProtocol, illuminaAdapterKit,
              illuminaIndex1, illuminaIndex2, sequencerRunID, qaqcStatus, ncbiProjectID,
              analyzedBy, remarks, dataQF, publicationDate, release))
  
  genomic_traits_w_coors$collectDate <- as.Date(genomic_traits_w_coors$collectDate)
  
  
  return(genomic_traits_w_coors)
}

