packages <- c('tidyverse', 'patchwork')

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



state <- map_data("world")

full <- read.csv("data/input/derived/full_df.csv") %>%
  mutate(Mbp = AGS/1000000)

theme_pete<- function(base_size = 10, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_minimal(base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
      
    )
}

min_lat = min(full$lat)-2
max_lat = max(full$lat)+2
min_long = min(full$long)-10
max_long = max(full$long)+2

Fig_map <- ggplot() +
  geom_map(
    data = state, map = state,
    aes(long, lat, map_id = region),
    color = "black", fill = "white"
  )+
  geom_jitter(data = full, aes(y = lat, x = long, 
                              fill = GC_bact*100, size = Mbp), 
             alpha = 0.8, shape = 21, width = .5)+
  #scale_fill_manual(values = pal_hmts2)+
  xlab("Longitude") + 
  ylab("Latitude")+
  theme_pete()+
  #scale_shape_manual(values = shape_scale)+
  theme(legend.position = "none",
        legend.box.background = element_rect(colour = "black"),
        legend.text=element_text(size=8))+
  guides(fill=guide_legend(title="System"),
         shape =guide_legend(title="System") )+
  scale_size(name = "Average\nGenome Size")+
  scale_fill_viridis_c("Bacterial\nGC-%")+
  xlim(min_long, max_long)+
  ylim(min_lat, max_lat)

Fig_map


theme_neat<- function(base_size = 10) {
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

pH_cn_gc <- full%>%
  ggplot(., aes(mean_pH_water, soil_OrgCtoN))+
  geom_point(aes(fill = GC_bact*100, size = Mbp), alpha = 0.9, shape =21)+
  ylab("Extractable C:N")+
  xlab("Soil pH (water)")+
  theme_neat()+
  scale_fill_viridis_c("Bacterial\nGC-%")+
  scale_radius("Ave. Genome\n size (Mbp)", 
             breaks = c(5, 6, 7, 8),
             labels = c("5", "6", "7", "8"))+
  theme(legend.position = c(-0.1,-.2),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.justification = c("left", "top"),
        legend.title = element_text(size=rel(0.8)),
        legend.background = element_rect(fill=alpha('white', 0.9))
        )

pH_cn_gc

layout <- c(
  area(t = 1, l = 1, b = 11, r = 8),
  area(t = 6, l = 1, b = 10, r = 3)
)


Fig4 <- Fig_map + pH_cn_gc + 
  plot_layout(design = layout)

Fig4


ggsave(Fig4, filename = "paper/figures/Fig4.png", width = 8, height = 6)
