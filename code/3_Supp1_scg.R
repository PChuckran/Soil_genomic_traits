packages <- c('tidyverse', 'car', 'AICcmodavg')

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

full_df <- read.csv("data/input/derived/full_df.csv")%>%
  mutate(Mbp = AGS/1000000)

scg_data <- read.csv("data/input/derived/all_scgs_all.csv")
#scg_data_wt <- read.csv("data/input/derived/all_scgs_wt.csv")
#scg_data_wt <- scg_data_wt %>%
#  rename(AGS_contig_old = AGS_contig_wt)
#scg_data <- scg_data %>%
#  left_join(., scg_data_wt)
meta <- read.csv("data/input/derived/summary_statistic_table.csv")%>%
  select(-X)

meta <-  meta%>%
  mutate(ID = str_replace(Sample_name, "/", ""))%>%
  mutate(ID = str_replace(ID, "_mms", ""))


scg_data <- scg_data %>%
  mutate(ID = str_replace(MGID, "_R_genes", ""))%>%
  mutate(ID = str_replace(ID, "_mms", ""))


test_df <- left_join(scg_data, full_df)

test_df <- left_join(test_df, meta)

another_test <- left_join(full_df, meta)


test_df %>%
  filter(contributing_bp > 200000000) %>%
  ggplot(., aes(AGS, AGS_contig_wt))+
  geom_point()+
  ylab("Average genome size (bp) of Bacteria")+
  #xlab("N50")+
  geom_smooth(method = "lm")

for_fig <- test_df %>%
  filter(contributing_bp > 200000000) %>% 
  filter(! AGS_contig_wt %in% boxplot(test_df$AGS_contig_wt)$out
  )



SFig1C <- for_fig %>%
  ggplot(., aes(mean_pH_water, AGS_contig_wt))+
  ylab("Average genome size (bp) of Bacteria")+
  xlab("pH")+
  geom_point()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label..,  ..p.value.label.., sep = "~~~")),
               label.x = "right",
               small.p = T,
               label.y = "top",
               parse = TRUE)+
  geom_smooth(method = "lm")+
  theme_classic()
SFig1C




n50_pH
summary(lm(AGS_contig_wt~n50, test_df))

ggsave("~/Desktop/n50_pH.png", n50_pH)


k1 <- lm(AGS~AGS_contig_wt, for_fig)
k2 <- lm(AGS ~lipid_FtoB, for_fig)
k3 <- lm(AGS ~lipid_FtoB+AGS_contig_wt, for_fig)

AIC(k1, k2, k3)



my.formula <- y ~ poly(x, 1, raw=TRUE)
lipid_AGS <- full_df %>%
  filter(!is.na(mean_pH_water))%>%
  ggplot(., aes(lipid_FtoB, AGS))+
  geom_smooth(method = "lm", se = F, formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label..,  ..p.value.label.., sep = "~~~")),
               label.x = "right",
               label.y = "bottom",
               small.p = T,
               parse = TRUE)+
  ylab("Average genome size (bp)")+
  xlab("Fungi:Bacteria biomass (PLFA)")+
  geom_point(alpha = 0.9, size = 2.5)+
  theme_classic()

ph_lipid <- full_df %>%
  ggplot(., aes(mean_pH_water, lipid_FtoB))+
  geom_smooth(method = "lm", se = F, formula=my.formula)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label..,  ..p.value.label.., sep = "~~~")),
               label.x = "right",
               small.p = T,
               label.y = "top",
               parse = TRUE)+
  ylab("Fungi:Bacteria biomass (PLFA)")+
  xlab("pH")+
  geom_abline()+
  geom_point(alpha = 0.9, size = 2.5)+
  theme_classic()

FigS1 <- ph_lipid+lipid_AGS+SFig1C+plot_annotation(tag_levels = "A")

FigS1

ggsave(FigS1, filename = "paper/figures/Fig_S1.png")

summary(k2)
AIC(k1, k2)


