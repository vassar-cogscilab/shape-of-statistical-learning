library(runjags)
library(dplyr)
library(tidyr)
library(ggplot2)
library(coda)

dataJags<-function(data){
data.for.jags <-list(rt = data$rt,
  subject_id = data$subject_id,
  is.predictable = data$is_predictable,
  t = data$t,
  N = length(data$rt),
  S = length(unique(data$subject_id)),
  C = length(unique(data$subject_condition)),
  condition = (data%>%group_by(subject_id) %>% summarise(subject_condition = unique(subject_condition)) %>% select(subject_condition) %>% as.matrix)[,1],
  max_t = (data%>%group_by(subject_id)%>% summarise(max_t = max(t)) %>% select(max_t) %>% as.matrix)[,1],
  n.pairs = (data %>% group_by(subject_id) %>% summarise(n.pairs = length(unique(pair))) %>% select(n.pairs) %>% as.matrix)[,1]
)
return(data.for.jags)
}

dat<-read.csv('data/for-jags-exp-1.csv')

data.for.jags<-dataJags(dat)

params.to.monitor <- c('v', 'alpha', 'gamma', 'beta', 'beta.learn','gamma.learn', 'delta')

jags.result <- run.jags('jags-model-expexp-reparam.txt', monitor=params.to.monitor, data=data.for.jags, n.chains=2,
                        burnin=1000, sample=1000, adapt=500)

result <- as.matrix(as.mcmc(jags.result))
save(jags.result, file = 'example_model_output')


plot(jags.result)

summary(jags.result)
