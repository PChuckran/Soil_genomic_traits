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

theme_map<- function(base_size = 8, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.border = element_rect(size = 1, fill = NA),
      strip.background = element_blank(),
      #panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      #panel.grid.major.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank()
      
    )
}

min_lat = min(full_df$lat)-5
max_lat = max(full_df$lat)+8
min_long = min(full_df$long)-9
max_long = max(full_df$long)+2

state <- map_data("world")

full_df_site <- full_df %>%
  group_by(siteID)%>%
  summarise_all(funs(mean), na.rm = T)

Fig_1A <- ggplot() +
  geom_map(
    data = state, map = state,
    aes(long, lat, map_id = region),
    color = "black", fill = "white"
  )+
  geom_jitter(data = full_df_site, aes(y = lat, x = long, 
                                  fill = GC_bact*100, size = Mbp), 
              alpha = 0.8, shape = 21, width = .5)+
  #scale_fill_manual(values = pal_hmts2)+
  xlab("Longitude") + 
  ylab("Latitude")+
  theme_map()+
  #scale_shape_manual(values = shape_scale)+
  theme(legend.position = c(0.6,0.8),
        legend.box.background = element_rect(colour = "black"),
        legend.text=element_text(size=7),
        legend.box = "vertical",
        legend.direction = "horizontal",
        legend.key.size = unit(.3, 'cm'),
        legend.box.just = "center")+
  #scale_fill_continuous(palette = "Greens", "GC (%)")+
  scale_fill_distiller("Bacterial\nGC (%)", palette = "BuGn", breaks = c(60, 64, 68))+
  xlim(min_long, max_long)+
  ylim(min_lat, max_lat)+
  theme(text=element_text(family="Arial"))+
  scale_size(name = "Average\ngenome size (Mbp)", breaks = c(5.5, 6.5, 7.5))

Fig_1A





agg <- c("cultivatedCrops", "pastureHay")

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

# Remove samples with low depth
full_df<- full_df %>%
  #filter(!nlcdClass %in% agg)%>%
  filter(sampleFilteredReadNumber > 2000000)%>%
  mutate(mgmt = ifelse(nlcdClass %in% agg, "Agg", "Wild"),
         GC_100 = GC_bact*100)



outliers <- boxplot(full_df$soil_OrgCtoN, plot=FALSE)$out
full_df<- full_df[-which(full_df$soil_OrgCtoN %in% outliers),]


my.formula <- y ~ poly(x, 1, raw=TRUE)
Fig_1B_CtoN_GC <- full_df %>%
  ggplot(., aes(soil_OrgCtoN, GC_100))+
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
  ylab("GC (%)")+
  xlab(expression(Soil~C["extr"]:N["extr"]))+
  theme_pete()+
  theme(text=element_text(family="Arial"))

Fig_1B_CtoN_GC

       

mod_CN_GC <- lm(GC_100 ~ poly(soil_OrgCtoN, 3, raw = T), full_df)
summary(mod_CN_GC)



my.formula <- y ~ poly(x, 1, raw=TRUE)

Fig_1C_CtoN_AGS <- full_df %>%
  pivot_longer(., cols = c(AGS, AGS_marker), names_to = "AGS_type", values_to = "AGS_val") %>%
  mutate(Source = ifelse(AGS_type == "AGS", "Metagenome\nestimated", "Community-\nweighted mean"))%>%
  ggplot(., aes(soil_OrgCtoN, AGS_val, fill = Source))+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~"), color = Source),
               small.p = T, 
               label.x = "left",
               label.y = "top",
               parse = TRUE,
               size = 2.25)+
  ylab("Average genome size (bp)")+
  xlab(expression(Soil~C["extr"]:N["extr"]))+
  scale_linetype_manual(values=c("solid", "solid"))+
  scale_fill_manual(values = c("white", "#0C2340"))+
  scale_color_manual(values = c("blue", "#BD3039"))+
  geom_point(alpha = 0.7, size = 1.5, shape =21)+
  geom_smooth(aes(linetype=Source, color = Source), method = "lm", se = F, formula=my.formula)+
  theme_pete()+
  theme(#legend.position = c(.8,.15), 
    legend.position = "bottom",
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        text=element_text(family="Arial"),
        legend.background = element_rect(fill = "transparent"))

Fig_1C_CtoN_AGS


CN_AGS_MG <- lm(AGS ~ soil_OrgCtoN, full_df)
summary(CN_AGS_MG)
CN_AGS_m16<- lm(AGS_marker ~ soil_OrgCtoN, full_df)
summary(CN_AGS_m16)


my.formula <- y ~ poly(x, 1, raw=TRUE)

scg_data <- read.csv("data/input/derived/all_scgs_all.csv")

scg_data <- scg_data %>%
  mutate(ID = str_replace(MGID, "_R_genes", ""))%>%
  mutate(ID = str_replace(ID, "_mms", ""))

full_df <- left_join(full_df, scg_data)

for_fig <- full_df %>%
  filter(contributing_bp > 200000000) %>% 
  filter(! AGS_contig_wt %in% boxplot(full_df$AGS_contig_wt)$out
  )

my.formula <- y ~ poly(x, 1, raw=TRUE)

Fig_1D_AGS_GC <- for_fig %>%
  ggplot(., aes(AGS_contig_wt, GC_bact*100))+
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
  ylab("GC (%) content of Bacteria")+
  xlab("Average Genome Size of Bacteria")+
  theme_pete()+
  theme(text=element_text(family="Arial"))





layout <- "
AABD
AAC#
"


Fig1 <- Fig_1A+Fig_1B_CtoN_GC+Fig_1C_CtoN_AGS+Fig_1D_AGS_GC+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))
#Fig1



ggsave(filename = "paper/figures/Fig1.svg", Fig1, height =4.5, width = 8.5)

