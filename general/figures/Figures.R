source('general/figures/plotData.R')
library(ggplot2)
library(extrafont)
font_import(pattern="Montserrat")


loadfonts(device="win")
##figures

#group level

avg_plot_data_1 <- plot_data_1%>% ungroup()%>%
  group_by(cond,t,predictable,is_predictable)%>% summarise(mean_rt=mean(rt))


ggplot(avg_plot_data_1)+
  #geom_point(aes(x=t, y= mean_rt, col= cond, alpha = is_predictable),size =2.5)+
  geom_point(data = avg_plot_data_1 %>% filter(is_predictable == 0), aes(x=t,y=mean_rt, col = cond),size=2.5, alpha = .2, se =F)+
  geom_point(data = avg_plot_data_1 %>% filter(is_predictable == 1), aes(x=t,y=mean_rt, col = cond),size=2.5, alpha = 1, se =F)+
  guides(alpha = F)+
  ylim(0,1400)+
  scale_color_manual(guide=F, values=c("#e41a1c", "#377eb8", "#38761d"))+
  labs(y="Response time (ms)", x="Time (discrete presentations of item)", col = "predictable")+
  theme_bw(base_size = 28, base_family = "Montserrat")+
  ggtitle('Experiment 1')
ggsave("group-stat-plot-exp-1.png", device="png", path="general/figures/", dpi=300, width=14, height=8, units="in")

avg_plot_data_2 <- plot_data_2 %>% ungroup()%>%
  group_by(set_size,t,predictable, is_predictable)%>% summarise(mean_rt=mean(rt))

ggplot(avg_plot_data_2)+
  #geom_point(aes(x=t, y= mean_rt, col= set_size, alpha = is_predictable), size = 2.5)+
  geom_point(data = avg_plot_data_2 %>% filter(is_predictable == 0), aes(x=t,y=mean_rt, col = set_size),size=2.5, alpha = .2)+
  geom_point(data = avg_plot_data_2 %>% filter(is_predictable == 1), aes(x=t,y=mean_rt, col = set_size),size=2.5, alpha = 1)+
  #geom_smooth(aes(x=t,y=mean_rt, col = set_size), se =F)+
  guides(alpha = F)+
  #facet_grid(~predictable)+
  ylim(0,1400)+
  scale_color_manual(guide=F, values=c("#e41a1c", "#38761d","#377eb8" ))+
  labs(y="Response time (ms)", x="Time (discrete presentations of item)", col = "predictable")+
  theme_bw(base_size = 28, base_family = "Montserrat")+
  ggtitle('Experiment 2')
ggsave("group-stat-plot-exp-2.png", device="png", path="general/figures/", dpi=300, width=14, height=8, units="in")

avg_plot_data_3 <- plot_data_3 %>% ungroup()%>%
  group_by(set_size,t,predictable, is_predictable)%>% summarise(mean_rt=mean(rt))

ggplot(avg_plot_data_3)+
  geom_point(aes(x=t, y= mean_rt, col= set_size, alpha = is_predictable), size = 2.5)+
  geom_point(data = avg_plot_data_3 %>% filter(is_predictable == 0), aes(x=t,y=mean_rt, col = set_size),size=2.5, alpha = .2, se =F)+
  geom_point(data = avg_plot_data_3 %>% filter(is_predictable == 1), aes(x=t,y=mean_rt, col = set_size),size=2.5, alpha = 1, se =F)+
  #geom_smooth(aes(x=t,y=mean_rt, col = set_size), se =F)+
  guides(alpha = F)+
  #facet_grid(~predictable)+
  ylim(0,1400)+
  scale_color_manual(guide=F, values=c("#e41a1c",  "#38761d", "#377eb8"))+
  labs(y="Response time (ms)", x="Time (discrete presentations of item)", col = "predictable")+
  theme_bw(base_size = 28, base_family = "Montserrat")+
  ggtitle('Experiment 3')
ggsave("group-stat-plot-exp-3.png", device="png", path="general/figures/", dpi=300, width=14, height=8, units="in")

