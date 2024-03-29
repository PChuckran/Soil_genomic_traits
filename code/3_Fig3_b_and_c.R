
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

outliers <- boxplot(full_df$soil_OrgCtoN, plot=FALSE)$out
full_df<- full_df[-which(full_df$soil_OrgCtoN %in% outliers),]

theme_pete<- function(base_size = 12) {
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
  filter(sampleFilteredReadNumber > 2000000)



m1 <- lm(GC_bact ~ mean_pH_water, full_df)
m2 <- lm(GC_bact ~ poly(mean_pH_water, 2, raw=TRUE), full_df)
m3 <- lm(GC_bact ~ poly(mean_pH_water, 3, raw=TRUE), full_df)

AIC(m1, m2, m3)

my.formula <- y ~ poly(x, 1, raw=TRUE)
pH_GC <- full_df %>%
  ggplot(., aes(mean_pH_water, GC_bact*100))+
  geom_point(alpha = 0.7, size = 2.5, shape =21, fill = "black")+
  geom_smooth(method = "lm", se = F, color = "black", 
              formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right",
               label.y = "bottom",
               parse = TRUE)+
  ylab("GC (%)")+
  xlab("Soil pH")+
  theme_pete()

pH_GC

pH_AGS <- full_df %>%
  pivot_longer(., cols = c(AGS, AGS_marker), names_to = "AGS_type", values_to = "AGS_val") %>%
  mutate(Source = ifelse(AGS_type == "AGS", "Metagenome estimated", "Community weighted mean"))%>%
  ggplot(., aes(mean_pH_water, AGS_val, fill = Source))+
  geom_smooth(aes(color = Source), method = "lm", se = F, formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "left",
               label.y = "bottom",
               parse = TRUE)+
  ylab("Average Genome Size (bp)")+
  xlab("Soil pH")+
  scale_fill_manual(values = c("grey", "black"))+
  scale_color_manual(values = c("grey", "black"))+
  geom_point(alpha = 0.7, size = 2.5, shape =21)+
  theme_pete()+
  theme(legend.position = c(.75,.85), 
        legend.direction = "vertical",
        legend.title = element_text(size=rel(0.6)),
        legend.text = element_text(size=rel(0.6)))

pH_AGS

ggsave(plot = pH_GC, filename = "paper/figures/Fig3b.png", 
       height = 4, width = 4.5)

ggsave(plot = pH_AGS, filename = "paper/figures/Fig3c.png", 
       height = 4, width = 4.5)

