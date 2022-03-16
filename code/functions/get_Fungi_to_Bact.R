get_ftob <- function(refresh=FALSE){
  if(refresh == TRUE){
    domain_cts <- read_csv("data/input/derived/counts_by_domain.csv")
    domain_cts$ID <- str_replace(domain_cts$ID, "_mms", "")
    domain_cts$ID <- str_replace(domain_cts$ID, "^P", "BMI_P")
    domain_cts$Domain_gene_counts <- as.numeric(domain_cts$Domain_gene_counts) 
    
    domain_cts_summary <- domain_cts %>% 
      filter(Domain == "Bacteria" | Domain == "Eukaryota") %>%
      select(ID, Domain, Domain_gene_counts)%>%
      group_by(ID, Domain)%>%
      summarise_all(funs(mean), na.rm = T)%>%
      pivot_wider(., names_from = Domain, values_from = Domain_gene_counts)%>%
      mutate(FtoB = Eukaryota/Bacteria)
    
    domain_cts_out <- domain_cts %>%
      filter(Domain == "Bacteria") %>%
      select(ID, GC_adj)%>%
      left_join(domain_cts_summary, .)%>%
      select(-c(Bacteria, Eukaryota))
    
    write.csv(domain_cts_out, file = "data/input/derived/FtoB_summary.csv")
    return(domain_cts_out)
  }
  else{
    domain_cts_out <- read.csv("data/input/derived/FtoB_summary.csv")%>%
      select(-X)
    return(domain_cts_out)
  }
}