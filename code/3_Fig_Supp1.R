
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



my.formula <- y ~ poly(x, 1, raw=TRUE)
precip_AGS <- full_df %>%
  ggplot(., aes(field_mean_annual_precipitation_mm, AGS))+
  geom_smooth(method = "lm", se = F, formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label..,  ..p.value.label.., sep = "~~~")),
               label.x = "right",
               label.y = "bottom",
               parse = TRUE)+
  ylab("Average genome size (bp)")+
  xlab("MAP (mm)")+
  geom_point(alpha = 0.7, size = 2.5, shape =21)+
  theme_classic()

precip_AGS




ggsave(plot = precip_AGS, filename = "paper/figures/Fig_S1.png", 
       height = 4, width = 4)

