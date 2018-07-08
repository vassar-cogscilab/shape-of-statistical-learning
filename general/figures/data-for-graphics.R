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


data_1_raw<- read.csv('experiment-1/data/raw/all-data.csv')
block_data<-data_1_raw %>% filter(trial_type == 'letter-matching',triple_type == 'critical',practice == '0',correct == 1, triple_position%in%c(1,3))%>% group_by(block,cond,triple_position) %>% summarise(mean_rt = mean(rt))

ggplot(block_data)+
  geom_point(aes(x=block,y=mean_rt,col=as.character(triple_position)))+
  facet_grid(~cond)




#experiment 2
data_2<- read.csv('experiment-2/data/generated/for-jags-exp-2.csv')


data_2_raw<- read.csv('experiment-2/data/raw/sl_size_data.csv')
block_data<-data_2_raw %>% filter(trial_type == 'letter-matching',triple_type == 'critical',practice == '0',correct == 1, triple_position%in%c(1,3))%>% group_by(block,cond,triple_position) %>% summarise(mean_rt = mean(rt))

ggplot(block_data)+
  geom_point(aes(x=block,y=mean_rt,col=as.character(triple_position)))+
  facet_grid(~cond)

#experiment 3
data_3<- read.csv('experiment-3/data/generated/for-jags-exp-3.csv')



ggplot(simulated.data, aes(x=t, y=rt, color=type)) +
  geom_point(size=5, alpha=0.5)+
  geom_line(data=(plotting.data %>% filter(type %in% c('unpredictable', 'predictable'))), aes(y=y), size=3)+
  ylim(0,1000)+
  scale_color_manual(guide=F, values=c("#e41a1c", "#377eb8"))+
  labs(y="Response time (ms)", x="Time (discrete presentations of item)")+
  theme_bw(base_size = 28, base_family = "Montserrat")
ggsave("model-simulated-data.png", device="png", path="general/", dpi=300, width=15, height=8, units="in")


ggplot(plotting.data %>% filter(type=="relative learning rate"), aes(x=t, y=y)) +
  geom_line(size=3)+
  ylim(0,1)+
  #scale_color_manual(guide=F, values=c("#e41a1c", "#377eb8"))+
  labs(y="Relative learning", x=NULL)+
  theme_bw(base_size = 28, base_family = "Montserrat")
ggsave("model-rlr.png", device="png", path="general/", dpi=300, width=15, height=4, units="in")
