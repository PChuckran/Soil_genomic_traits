packages <- c('tidyverse', 'ggpmisc', 'viridis', 'patchwork', 'GGally', 'svglite')

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

theme_vip<- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      strip.background =element_blank(),
      strip.placement = "outside",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

theme_pete<- function(base_size = 10) {
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

vip_data <- read.csv("data/output/gc_variable_importance.csv")

top_vip_data <- vip_data %>%
  top_n(8, wt = importance)

variable_name_to_readable <- Vectorize(function(a) {
  switch(as.character(a),
         "phH2o_poly_1" = "pH",
         "mega_bsesatCecd10" = 'Base saturation',
         "VSIC_median" = "Site salinity (median)",
         "soil_OrgCtoN_poly_1" = "Extractable C:N (plot-level",
         "ctonRatio" = "C:N (plot-level)",
         "VSIC_min" = "Site salinity (min)",
         "mean_cn" = "C:N (collection)",
         "field_avg_number_of_green_days" = "Green days"
  )
  
})

Fig_3A <- top_vip_data%>%
  mutate(variable = variable_name_to_readable(variable))%>%
  mutate(variable = fct_reorder(variable, importance)) %>%
  ggplot(., aes(variable, importance))+
  geom_segment( aes(x=variable, xend=variable, y=0, yend=importance), color="#000E54") +
  geom_point( color="#F76900", size=4, alpha=1) +
  theme_vip()+
  coord_flip()+
  ylab("Variable importance")+
  xlab("")

top_vip_data%>%
  mutate(variable = variable_name_to_readable(variable))



full_df <- read.csv("data/input/derived/full_df.csv")%>%
  mutate(Mbp = AGS/1000000)

outliers <- boxplot(full_df$soil_OrgCtoN, plot=FALSE)$out
full_df<- full_df[-which(full_df$soil_OrgCtoN %in% outliers),]



# Remove samples with low depth
full_df<- full_df %>%
  filter(sampleFilteredReadNumber > 2000000)

my.formula <- y ~ poly(x, 1, raw=TRUE)

Fig_3C_ph_GC <- full_df %>%
  ggplot(., aes(mean_pH_water, GC_bact*100))+
  geom_point(alpha = 0.7, size = 2.5, shape =21, fill = "#F76900")+
  geom_smooth(method = "lm", se = F, color = "#000E54", 
              formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               label.x = "right",
               small.p = T,
               label.y = "bottom",
               parse = TRUE,
               size = 3)+
  ylab("GC (%)")+
  xlab("Soil pH")+
  theme_pete()


Fig_3B_ph_AGS <- full_df %>%
  ggplot(., aes(mean_pH_water, AGS))+
  geom_point(alpha = 0.7, size = 2.5, shape =21, fill = "#F76900")+
  geom_smooth(method = "lm", se = F, color = "#000E54", 
              formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               label.x = "right",
               small.p = T,
               label.y = "top",
               parse = TRUE,
               size = 3)+
  ylab("Average genome size (bp)")+
  xlab("Soil pH")+
  theme_pete()



Fig_3D_ph_cn_traits <- full_df%>%
  ggplot(., aes(mean_pH_water, soil_OrgCtoN))+
  geom_smooth(method = "lm", se = F, color = "grey", formula=my.formula)+
  geom_point(aes(fill = GC_bact*100, size = Mbp), alpha = 0.9, shape =21)+
  ylab(expression(Soil~C["extr"]:N["extr"]))+
  xlab("Soil pH")+
  theme_pete()+
  scale_fill_viridis_c("GC (%)", breaks = c(60, 65, 70))+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               label.x = "left",
               label.y = "bottom",
               small.p = T,
               color = "black",
               parse = TRUE,
               size = 3)+
  scale_radius("Average genome\n size (Mbp)", 
               breaks = c(6, 7, 8),
               labels = c("6", "7", "8"))+
  theme(legend.position = c(.75,.85),
        legend.box = "verticle",
        legend.direction = "horizontal",
        #legend.justification = "right",
        legend.title = element_text(size=rel(0.8)),
        legend.background = element_rect(fill=alpha('white', 0.9)),
        legend.key.size = unit(.25, 'cm'),        
        legend.text=element_text(size=7),
  )



Fig_3 <- ggpubr::ggarrange(Fig_3A, Fig_3B_ph_AGS, Fig_3C_ph_GC, Fig_3D_ph_cn_traits, 
                      labels = c("A","B","C","D"), font.label = list(size = 12, color = "black", face = "plain", family = "Arial"))

Fig_3

ggsave(Fig_3,filename =  "paper/figures/Fig3.svg", height = 8, width = 8)

ggsave(Fig_3,filename =  "paper/figures/Fig3.png", height = 8, width = 8)
