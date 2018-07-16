intercept.adapt <- 770
intercept.learn <- 889

beta.adapt <- 0.2
rate.adapt <- 0.15

beta.learn <- 0.4
rate.learn <- .16
onset.learn <- 24

# plot the two curves ####
layout(1:2)
t <- 0:71
adapt <- (1 + beta.adapt * (exp(-rate.adapt * t) - 1))
learn <- sapply(t, function(x){
  if(x < onset.learn) { return(0); }
  else { return( beta.learn * ( 1 + -exp(-rate.learn * (x - onset.learn)))) }
})
plot(t, intercept.adapt*adapt, col="red", ylim=c(0,2000))
points(t, intercept.learn * adapt * (1 - learn), col="blue")

plot(t, learn, ylim=c(0,1))

# fancier plot with simulated data

library(dplyr)

plotting.data <- data.frame(
  t = rep(0:71,3),
  y = c(intercept.adapt*adapt, intercept.learn * adapt * (1-learn), learn),
  type = c(rep('unpredictable', 72), rep('predictable', 72), rep('relative learning rate', 72))
)

simulated.data <- plotting.data %>% filter(type %in% c('unpredictable', 'predictable')) %>%  mutate(rt = rnorm(n(), y, 100))

library(ggplot2)
library(extrafont)
#font_import(pattern="Montserrat")
#loadfonts(device="win")
library(cowplot)

fit.plot <- ggplot(simulated.data, aes(x=t, y=rt, alpha=type)) +
  geom_point(size=5, color="#e41a1c")+
  geom_line(data=(plotting.data %>% filter(type %in% c('unpredictable', 'predictable'))), aes(y=y), color="#e41a1c", size=3)+
  ylim(0,1000)+
  scale_alpha_discrete(range=c(1,0.2), guide=F)+
  #scale_color_manual(guide=F, values=c("#333333", "#377eb8"))+
  labs(y="Response time (ms)", x=NULL)+
  theme_bw(base_size = 28, base_family = "Montserrat")

rlr.plot <- ggplot(plotting.data %>% filter(type=="relative learning rate"), aes(x=t, y=y)) +
  geom_line(size=3)+
  ylim(0,1)+
  #scale_color_manual(guide=F, values=c("#e41a1c", "#377eb8"))+
  labs(y="Relative learning", x="Time (discrete presentations of item)")+
  theme_bw(base_size = 28, base_family = "Montserrat")

cowplot::plot_grid(fit.plot, rlr.plot, nrow=2, align="V", rel_heights = c(1.8,1))

ggsave("model-simulated-data.png", device="png", path="general/", dpi=300, width=16, height=12, units="in")


ggsave("model-rlr.png", device="png", path="general/", dpi=300, width=15, height=4, units="in")
