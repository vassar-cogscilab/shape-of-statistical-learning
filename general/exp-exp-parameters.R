intercept.adapt <- 1070
intercept.learn <- 1189

beta.adapt <- 0.18
rate.adapt <- 0.28

beta.learn <- 0.4
rate.learn <- .1
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
font_import(pattern="Montserrat")
loadfonts(device="win")

ggplot(simulated.data, aes(x=t, y=rt, color=type)) +
  geom_point(size=5, alpha=0.5)+
  geom_line(data=(plotting.data %>% filter(type %in% c('unpredictable', 'predictable'))), aes(y=y), size=3)+
  ylim(0,1500)+
  scale_color_manual(guide=F, values=c("#e41a1c", "#377eb8"))+
  labs(y="Response time (ms)", x="Time (discrete presentations of item)")+
  theme_bw(base_size = 28, base_family = "Montserrat")
ggsave("model-simulated-data.png", device="png", path="general/", dpi=300, width=10, height=8, units="in")
