
packages <- c('tidyverse', 'ggpmisc', 'viridis', 'patchwork', 'GGally')

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


full_df <- read.csv("data/input/derived/full_df.csv")

agg <- c("cultivatedCrops", "pastureHay")

theme_pete<- function(base_size = 18) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      strip.background =element_blank(),
      strip.placement = "outside",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

# Remove samples with low depth
full_df<- full_df %>%
  #filter(!nlcdClass %in% agg)%>%
  filter(sampleFilteredReadNumber > 2000000)%>%
  mutate(mgmt = ifelse(nlcdClass %in% agg, "Agg", "Wild"),
         GC_100 = GC_bact*100)



outliers <- boxplot(full_df$soil_OrgCtoN, plot=FALSE)$out
full_df<- full_df[-which(full_df$soil_OrgCtoN %in% outliers),]



my.formula <- y ~ poly(x, 3, raw=TRUE)
CtoN_GC <- full_df %>%
  ggplot(., aes(soil_OrgCtoN, GC_100))+
  geom_point(alpha = 0.7, size = 2.5, shape =21, fill = "darkgreen")+
  geom_smooth(method = "lm", se = F, color = "black", formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right",
               label.y = "top",
               parse = TRUE)+
  ylab("GC (%)")+
  xlab("Soil C:N ratio")+
  theme_pete()

CtoN_GC


mod_CN_GC <- lm(GC_100 ~ poly(soil_OrgCtoN, 3, raw = T), full_df)
summary(mod_CN_GC)

my.formula <- y ~ poly(x, 1, raw=TRUE)

GC_AA <- full_df %>%
  ggplot(., aes(GC*100, AA_CN))+
  geom_smooth(method = "lm", se = F, color = "black", formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right",
               label.y = "top",
               parse = TRUE)+
  ylab("Amino Acid C:N")+
  xlab("GC (%)")+
  geom_point(alpha = 0.7, size = 2.5, shape =21, fill = "darkgreen")+
  theme_pete()
GC_AA

mod_GC_AA <- lm(AA_CN ~ GC_100, full_df)
summary(mod_GC_AA)


my.formula <- y ~ poly(x, 1, raw=TRUE)

CtoN_AA <- full_df %>%
  ggplot(., aes(soil_OrgCtoN, AA_CN))+
  geom_smooth(method = "lm", se = F, color = "black", formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right",
               label.y = "bottom",
               parse = TRUE)+
  ylab("Amino Acid C:N")+
  xlab("Soil C:N ratio")+
  geom_point(alpha = 0.7, size = 2.5, shape =21, fill = "darkgreen")+
  theme_pete()
CtoN_AA

mod_CN_AA <- lm(AA_CN ~ soil_OrgCtoN, full_df)
summary(mod_CN_AA)

CtoN_AGS <- full_df %>%
  pivot_longer(., cols = c(AGS, AGS_marker), names_to = "AGS_type", values_to = "AGS_val") %>%
  mutate(Source = ifelse(AGS_type == "AGS", "Metagenome estimated", "Community weighted mean"))%>%
  ggplot(., aes(soil_OrgCtoN, AGS_val, fill = Source))+
  geom_smooth(aes(color = Source), method = "lm", se = F, formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right",
               label.y = "bottom",
               parse = TRUE)+
  ylab("Average Genome Size (bp)")+
  xlab("Soil C:N ratio")+
  scale_fill_manual(values = c("wheat1", "darkgreen"))+
  scale_color_manual(values = c("wheat4", "darkgreen"))+
  geom_point(alpha = 0.7, size = 2.5, shape =21)+
  theme_pete()+
  theme(legend.position = "bottom", 
        legend.direction = "vertical")


CN_AGS_MG <- lm(AGS ~ soil_OrgCtoN, full_df)
summary(CN_AGS_MG)
CN_AGS_m16<- lm(AGS_marker ~ soil_OrgCtoN, full_df)
summary(CN_AGS_m16)

Fig1 <- CtoN_GC + GC_AA + CtoN_AA + CtoN_AGS+
  plot_annotation(tag_levels = 'A')
Fig1


#ggsave(plot = Fig1, filename = "paper/figures/Fig1.png", 
#       height = 8, width = 8)

CtoN_AGS_Gensci <- full_df %>%
  ggplot(., aes(soil_OrgCtoN, AGS))+
  geom_smooth(method = "lm", se = F, color = "black", formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right",
               label.y = "bottom",
               parse = TRUE)+
  ylab("Average\nGenome Size (bp)")+
  xlab("Soil C:N ratio")+
  geom_point(alpha = 0.7, size = 2.5, shape =21, fill = "darkgreen")+
  theme_pete()

CtoN_AGS_Gensci
CtoN_GC + CtoN_AGS_Gensci + CtoN_AA + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 2)

