library(tidyr)
library(dplyr)
library(ggplot2)
library(extrafont)
font_import(pattern="Montserrat")
loadfonts(device="win")

#experiment 1
data_1<- read.csv('experiment-1/data/generated/for-jags-exp-1.csv')

ggplot(data_1)+
  geom_smooth(aes(x=t,y=rt, col= as.character(is_predictable)))+
  facet_grid(~subject_condition)+
  ylim(0,1000)+
  scale_color_manual(guide=F, values=c("#e41a1c", "#377eb8"))+
  labs(y="Response time (ms)", x="Time (discrete presentations of item)", col = "predictable")+
  theme_bw(base_size = 28, base_family = "Montserrat")

