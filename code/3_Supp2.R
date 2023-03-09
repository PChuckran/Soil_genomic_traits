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
  ylab("Average genome size (bp)")+
  #xlab("N50")+
  geom_smooth(method = "lm")

for_fig <- test_df %>%
  filter(contributing_bp > 200000000) %>% 
  filter(! AGS_contig_wt %in% boxplot(test_df$AGS_contig_wt)$out
  )



boxplot(for_fig$AGS_contig_wt)$out

ags_pH <- for_fig %>%
  ggplot(., aes(mean_pH_water, AGS_contig_wt))+
  ylab("Average genome size (bp)")+
  xlab("pH")+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()
ags_pH
summary(lm(AGS_contig_wt~mean_pH_water, for_fig))

summary(lm(GC_bact~mean_pH_water, full_df %>%
             filter(horizon == "M")))

summary(lm(GC_bact~soil_OrgCtoN, full_df ))

for_fig %>%
  ggplot(., aes(AGS_contig, GC_bact))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

summary(lm(GC_bact~AGS_contig_wt, for_fig))


ggsave("paper/figures/Fig_S1.png", ags_pH)

test_df %>%
  filter(contributing_bp > 200000000) %>%
  ggplot(., aes(contributing_bp, AGS_contig))+
  #xlim(400, 650)+
  geom_point(data = test_df, aes(contributing_bp, AGS_contig), color = "grey", alpha = 0.5)+
  geom_point()+
  ylab("Average Genome Size")+
  xlab("Read depth (in bp)")+
  theme_linedraw()

n50_pH
summary(lm(AGS_contig_wt~n50, test_df))

ggsave("~/Desktop/n50_pH.png", n50_pH)



for_mod <- test_df %>%
  filter(contributing_bp > 200000000)%>%
  drop_na(mean_pH_water, contributing_bp) %>%
  filter(AGS_contig != Inf)

m0 <- lm(AGS_contig_wt~1, for_mod)
m1 <- lm(AGS_contig_wt~contributing_bp, for_mod)
m2 <- lm(AGS_contig_wt~mean_pH_water, for_mod)
m3 <- lm(AGS_contig_wt~mean_pH_water+contributing_bp, for_mod)

AIC(m0, m1, m2, m3)

model_list = list('m0' =m0, 'm1' =m1, 'm2'=m2, 'm3'=m3)
aictab(model_list)

vif(m3)
confint(m2)


AIC(m1, m2)
summary(m1)
lmer()

n1 <- lm(AGS_contig_wt~mean_pH_water+n50, for_mod)
n2 <- lm(AGS_contig_wt~n50, for_mod)
AIC(n1, n2)
summary(n1)
