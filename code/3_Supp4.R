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

full_df <- read.csv("data/input/derived/full_df.csv")

theme_pete<- function(base_size = 8) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.border = element_rect(size = 1, fill = NA),
      strip.background =element_blank(),
      strip.placement = "outside",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

my.formula <- y ~ poly(x, 1, raw=TRUE)

F4A <- full_df %>%
  ggplot(., aes(GC_16s*100, GC_bact*100))+
  geom_smooth(method = "lm", se = F)+
  xlab("GC content of 16S rRNA genes")+
  ylab("GC content of bacterial contigs")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.25)+
  theme_pete()+
  geom_point()

F4B <- full_df %>%
  ggplot(., aes(GC_16s*100, GC_marker*100))+
  geom_smooth(method = "lm", se = F)+
  xlab("GC content of 16S rRNA genes")+
  ylab(" GC content predicted from
       16S rRNA gene taxonomy")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.25)+
  theme_pete()+
  geom_point()

F4C <- full_df %>%
  ggplot(., aes(GC_bact*100, GC_marker*100))+
  geom_smooth(method = "lm", se = F)+
  xlab("GC content of bacterial contigs")+
  ylab(" GC content predicted from
       16S rRNA gene taxonomy")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.25)+
  theme_pete()+
  geom_point()

SF4 <- F4A+F4B+F4C+plot_annotation(tag_levels = "A")+plot_layout(ncol = 2)
SF4

ggsave(plot = SF4, "paper/figures/Fig_S4.png", width = 7.5, height = 7)

full_df %>%
  ggplot(., aes(AGS, AGS_marker))+
  geom_smooth(method = "lm")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               small.p = T, 
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.25)+
  theme_pete()+
  geom_point()
