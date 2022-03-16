CN_vs_AGS <- combined_gen_chem %>%
  filter(sampleFilteredReadNumber > 2000000)%>%
  ggplot(., aes(mean_OrgC/mean_N, AGS, color = nlcdClass))+
  geom_point(size = 2, alpha = 0.7)+
  ylab("Average Genome Size (bp)")+
  xlab("Soil C:N")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               label.x = "right",
               label.y = "bottom",
               parse = TRUE)+  
  geom_smooth(method = "lm", se = F, color = "black")+
  theme_bw()

ggplotly(CN_vs_AGS)
