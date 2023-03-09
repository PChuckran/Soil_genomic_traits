packages <- c('tidyverse', 'ggpmisc', 'viridis', 'patchwork', 'GGally', 'multcompView', 'ggtext', 'svglite')

## method for loading all packages, and installing ones that are missing if this is launched on a different rig
lapply(packages, function(x){
  #browser()
  ## if the required package is present (TRUE), load it up
  if (require(x, character.only = TRUE)==TRUE){
    library(x, character.only=TRUE)
  } else {
    install.packages(x)
    library(x, character.only=TRUE)
  }
})

full_df <- read.csv("data/input/derived/full_df.csv") %>%
  mutate(Mbp = AGS/1000000)

full_df<- full_df %>%
  #filter(!nlcdClass %in% agg)%>%
  filter(sampleFilteredReadNumber > 2000000)%>%
  mutate(mgmt = ifelse(nlcdClass %in% agg, "Agg", "Wild"),
         GC_100 = GC_bact*100)

outliers <- boxplot(full_df$soil_OrgCtoN, plot=FALSE)$out
full_df<- full_df[-which(full_df$soil_OrgCtoN %in% outliers),]

Fig_2A_CtoN_AA <- full_df %>%
  ggplot(., aes(soil_OrgCtoN, AA_CN))+
  geom_point(alpha = 0.7, size = 1.5, shape =21, fill = "#0C2340")+
  geom_smooth(method = "lm", se = F, color = "#BD3039", formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "bottom",
               color = "#BD3039",
               parse = TRUE,
               size = 2.25)+
  ylab("Amino acid C:N")+
  xlab(expression(Soil~C["extr"]:N["extr"]))+
  theme_pete()+
  theme(text=element_text(family="Arial"))
Fig_2A_CtoN_AA



mod_CN_AA <- lm(AA_CN ~ soil_OrgCtoN, full_df)
summary(mod_CN_AA)



my.formula <- y ~ poly(x, 1, raw=TRUE)

Fig_2B_GC_AA <- full_df %>%
  ggplot(., aes(GC_100, AA_CN))+
  geom_point(alpha = 0.7, size = 1.5, shape =21, fill = "#0C2340")+
  geom_smooth(method = "lm", se = F, color = "#BD3039", formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "right",
               label.y = "top",
               color = "#BD3039",
               parse = TRUE,
               size = 2.25)+
  ylab("Amino acid C:N")+
  xlab("GC (%)")+
  theme_pete()+
  theme(text=element_text(family="Arial"))
Fig_2B_GC_AA

mod_GC_AA <- lm(AA_CN ~ GC_100, full_df)
summary(mod_GC_AA)





codon_aa_key <- read_csv("data/input/derived/AA_w_codon_CN.csv")

codon_aa_key$GC_pairs <- as.character(codon_aa_key$GC)

AA_codon_means <- codon_aa_key %>% 
  group_by(AA)%>%
  summarize(mean_codon_CN = mean(Nuc_CN))

Codon_means <- codon_aa_key %>% 
  group_by(Codon)%>%
  summarize(mean_codon_CN = mean(Nuc_CN))

MG_aa_codon <- read.csv("data/input/derived/codon_usage.csv") %>%
  select(-X)
MG_aa_codon$stop = MG_aa_codon$X.


# Red rid of some tails on the labID that could be problematic 
#MG_aa_codon$ID <- str_replace(MG_aa_codon$ID, "_mms", "")

# Listen, I know this isn't the most elegant solution, but I just took the values from the `AA_codon_means`
# table to get these values. The nice thing is the AA in the df are proportions, so no dividing necessary
MG_aa_codon <- MG_aa_codon %>%
  mutate(theor_codon_cn =
           ((stop*1.363333)+(A*1.175000)+(C*1.280000)+
              (D*1.280000)+(E*1.280000)+(F*1.380000)+
              (G*1.175000)+(H*1.280000)+(I*1.396667)+
              (K*1.380000)+(L*1.313333)+(M*1.330000)+
              (N*1.380000)+(P*1.175000)+(Q*1.280000)+
              (R*1.210000)+(S*1.280000)+(T*1.280000)+
              (V*1.280000)+(W*1.230000)+(Y*1.380000)
           ),
         real_codon_cn = 
           ((AAA*1.43)+(AAC*1.33)+(AAG*1.33)+(AAT*1.43)+(ACA*1.33)+
              (ACC*1.23)+(ACG*1.23)+(ACT*1.33)+(AGA*1.33)+(AGC*1.23)+
              (AGG*1.23)+(AGT*1.33)+(ATA*1.43)+(ATC*1.33)+(ATG*1.33)+
              (ATT*1.43)+(CAA*1.33)+(CAC*1.23)+(CAG*1.23)+(CAT*1.33)+
              (CCA*1.23)+(CCC*1.12)+(CCG*1.12)+(CCT*1.23)+(CGA*1.23)+
              (CGC*1.12)+(CGG*1.12)+(CGT*1.23)+(CTA*1.33)+(CTC*1.23)+
              (CTG*1.23)+(CTT*1.33)+(GAA*1.33)+(GAC*1.23)+(GAG*1.23)+
              (GAT*1.33)+(GCA*1.23)+(GCC*1.12)+(GCG*1.12)+(GCT*1.23)+
              (GGA*1.23)+(GGC*1.12)+(GGG*1.12)+(GGT*1.23)+(GTA*1.33)+
              (GTC*1.23)+(GTG*1.23)+(GTT*1.33)+(TAA*1.43)+(TAC*1.33)+
              (TAG*1.33)+(TAT*1.43)+(TCA*1.33)+(TCC*1.23)+(TCG*1.23)+
              (TCT*1.33)+(TGA*1.33)+(TGC*1.23)+(TGG*1.23)+(TGT*1.33)+
              (TTA*1.43)+(TTC*1.33)+(TTG*1.33)+(TTT*1.43)),
         codon_CN_percOff = (theor_codon_cn-real_codon_cn)/theor_codon_cn
  )




MG_aa_codon$ID <- str_replace(MG_aa_codon$ID, "_mms", "")
MG_aa_codon$ID <- str_replace(MG_aa_codon$ID, "^P", "BMI_P")




genomic_traits <- read.csv("data/input/derived/full_df.csv") %>%
  select(-X)

MG_aa_codon_chem <- left_join(MG_aa_codon, genomic_traits)



MG_aa_codon_chem_sub_long <- MG_aa_codon_chem %>%
  select(AA_CN, mean_cn, mean_OrgC, mean_N,
         theor_codon_cn, real_codon_cn, codon_CN_percOff, horizon)%>%
  pivot_longer(cols=c(theor_codon_cn, real_codon_cn), names_to = "type", values_to = "codon_cn")



MG_codons_long <- MG_aa_codon_chem %>% 
  pivot_longer(c(TTT,TTC,TTA,TTG,CTT,CTC,CTA,
                 CTG,ATT,ATC,ATA,ATG,GTT,GTC,GTA, 
                 GTG,TAT,TAC,TAA,TAG,CAT,CAC,CAA, 
                 CAG,AAT,AAC,AAA,AAG,GAT,GAC,GAA, 
                 GAG,TCT,TCC,TCA,TCG,CCT,CCC,CCA, 
                 CCG,ACT,ACC,ACA,ACG,GCT,GCC,GCA, 
                 GCG,TGT,TGC,TGA,TGG,CGT,CGC,CGA, 
                 CGG,AGT,AGC,AGA,AGG,GGT,GGC,GGA,GGG), names_to = "Codon", values_to = "codon_freq")

MG_codons_long <- codon_aa_key %>%
  select(Codon, AA, GC_content, GC_pairs) %>%
  left_join(MG_codons_long, . )

MG_codons_long <- MG_codons_long %>%
  group_by(AA, ID) %>%
  summarize(AA_freq = sum(codon_freq)) %>%
  left_join(MG_codons_long,.) %>%
  mutate(codon_perc_of_aa = codon_freq/AA_freq)


MG_codons_long$AA <- str_replace(MG_codons_long$AA, pattern = "\\*", replacement = "Stop")

AA_Switch <- Vectorize(function(a) {
  switch(as.character(a),
         'A' = 'Alanine',
         'R' = 'Arginine',
         'N' = 'Asparagine',
         'D' = 'Aspartic Acid',
         'C' = 'Cysteine',
         'E' = 'Glutamic Acid',
         'Q' = 'Glutamine',
         'G' = 'Glycine',
         'H' = 'Histidine',
         'I' = 'Isoleucine',
         'L' = 'Leucine',
         'K' = 'Lysine',
         'M' = 'Methionine',
         'F' = 'Phenylalanine',
         'P' = 'Proline',
         'S' = 'Serine',
         'T' = 'Threonine',
         'W' = 'Tryptophan',
         'Y' = 'Tyrosine',
         'V' = 'Valine',
         'Stop' = 'Stop')
})


my.formula <- y ~ poly(x, 1, raw=TRUE)

Fig_2D_Asp_Glut <- MG_codons_long %>%
  filter(AA == "D" | AA == "E") %>%
  filter(soil_OrgCtoN < 40)%>%
  ggplot(., aes(soil_OrgCtoN, codon_perc_of_aa, group = Codon, color = GC_pairs)) +
  geom_point(alpha = 0.4, size = 2, formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right",
               small.p = T, 
               label.y = "top",
               parse = TRUE,
               size = 2)+
  ylab("Codon frequency for amino acid")+
  xlab(expression(Soil~C["extr"]:N["extr"]))+
  theme_pete()+
  geom_smooth(method = "lm", se = F)+
  #scale_fill_brewer(palette = "Set2")+
  #scale_color_brewer(palette = "Dark2")+
  scale_color_manual(name = "GC base pairs\nin codon", values = c("#007A33", "#BA9653"))+
  facet_wrap(.~AA_Switch(AA), ncol = 1)+
  theme(text=element_text(family="Arial"),
        legend.direction = "horizontal",
        legend.position = c(0,0))

Fig_2D_Asp_Glut


## Just going to look at fourfold substitution sites

# list of amino acids with a fourfold degenerate sites
degen_AA <- c("A", "G", "P", "T", "V")

degenerates <- MG_codons_long %>% filter(AA %in% degen_AA)

degenerates <- degenerates %>%
  mutate(deg_site = str_sub(Codon, start = 3, end = 3))

degenerates$deg_site <- factor(degenerates$deg_site, levels = c("A", "T", "C", "G"))


generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$type=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$type) , ]
  return(Tukey.labels)
}


codon_deg_mod <- aov(codon_perc_of_aa ~ deg_site, degenerates)
codon_deg_anova <- anova(codon_deg_mod)
codon_deg_anova
f2a_tukey <- TukeyHSD(codon_deg_mod,  "deg_site", conf.level=0.95)
f2a_tukey
f2a_labels <- generate_label_df(f2a_tukey, "deg_site")



Fig2C_degen <- degenerates %>%
  ggplot(., aes(deg_site, codon_perc_of_aa))+
  geom_boxplot(alpha = 0.5, color = "black", fill = "grey")+
  ylab("Frequency")+
  xlab("Four-fold degenerate nucleotide")+
  theme_pete()+
  geom_text(data = f2a_labels, aes(x = type, y = .8, label = Letters), inherit.aes = FALSE, color = "darkred")+
  theme(text=element_text(family="Arial"))

Fig2C_degen


output <- matrix(ncol=6, nrow = length(unique(MG_codons_long$Codon)))
i=1
for(codon_to_sub in unique(MG_codons_long$Codon)){
  sub <- MG_codons_long %>% filter(Codon == codon_to_sub)
  GC_to_add <- sub$GC_pairs[1]
  AA_to_add <- sub$AA[1]
  mod <- lm(codon_perc_of_aa~soil_OrgCtoN, sub)
  out <- summary(mod)
  Rsq <- out$r.squared
  slope <- mod$coefficients[2]
  p <- out$coefficients[2,4]
  output[i,] <- c(codon_to_sub, slope, Rsq, p, GC_to_add, AA_to_add)
  i = i +1
}

codon_output <- as.data.frame(output)
colnames(codon_output) <- c("Codon", "Slope", "Rsq", "p", "GC_pairs", "AA")


codon_output <- codon_output %>%
  mutate(firstC = str_sub(Codon, start = 1, end = 1),
         secondC = str_sub(Codon, start = 2, end = 2),
         thirdC = str_sub(Codon, start = 3, end = 3),
         dum = 1,
         Slope = as.numeric(Slope),
         AA_name = AA_Switch(AA),
         for_index = paste(AA, GC_pairs, sep = "_"))

codon_output <- codon_output %>%
  group_by(for_index) %>% 
  mutate(num_unique = length(unique(Codon)),
         splitn = 1:length(unique(Codon)),
         split_dum = splitn/num_unique,
         psig = ifelse(as.numeric(p) < 0.05, "*", " "),
         Codon_label = paste(Codon, psig))


Fig2E <- codon_output %>%
  ggplot(., aes(GC_pairs, dum, fill = Slope, label=Codon_label))+
  geom_col(position = "fill", color = "grey")+
  geom_text(position = "fill", vjust = 1.35, size = 2.5)+
  scale_fill_distiller(palette = "RdBu", name = expression(Slope~with~Soil~C["extr"]:N["extr"]))+
  #"Relationship\nbetween\nCodon Freq.\nand Soil C:N", palette = "RdBu")+
  #scale_fill_distiller(name = "Relationship\nbetween\nCodon Freq. \nand Soil C:N", palette = "RdBu")+
  #theme_tile()+
  facet_grid(AA_name~., switch="both")+
  theme(panel.grid = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "pt"),
        panel.background = element_rect(fill = "#ECE5F1"),
        strip.text.y.left = element_text(size = 10, angle = 0),
        legend.text=element_text(size=5),
        legend.title = element_text(size=7),
        legend.text.align = 1,
        legend.position = "bottom",
        text = element_text(size=8))+
  guides(fill =guide_colourbar(title.position = "top"))+
  ylab("Amino acid")+
  xlab("GC base pairs in codon")

Fig2E



layout <- "
ADE
BDE
C#E
##E
"

Fig2 <- Fig_2A_CtoN_AA+Fig_2B_GC_AA+Fig2C_degen+
  Fig_2D_Asp_Glut+Fig2E+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))
Fig2



ggsave(filename = "paper/figures/Fig2.svg", Fig2, height =9, width = 7)

