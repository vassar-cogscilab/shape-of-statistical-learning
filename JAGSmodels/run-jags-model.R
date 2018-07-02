library(runjags)
library(dplyr)
library(tidyr)
library(ggplot2)
library(coda)

dataJags<-function(data){
data.for.jags <-list(rt = data$rt,
  subject_id = data$subject_id,
  is.predictable = data$is.predictable,
  t = data$t,
  N = length(data$rt),
  S = length(unique(data$subject_id)),
  C = length(unique(data$subject_condition)),
  condition = data$subject_condition,
  max_t = (data%>%group_by(subject_id)%>% summarise(max_t = max(t)) %>% select(max_t) %>% as.matrix)[,1],
  n.pairs = (data %>% group_by(subject_id) %>% summarise(n.pairs = length(unique(pair))) %>% select(n.pairs) %>% as.matrix)[,1]
)
return(data.for.jags)
}

data.for.jags<-dataJags(sub)

params.to.monitor <- c('v', 'alpha', 'gamma', 'beta', 'beta.learn','gamma.learn', 'delta')

jags.result <- run.jags('jags-model-expexp.txt', monitor=params.to.monitor, data=data.for.jags, n.chains=2,
                        burnin=50, sample=100, adapt=50)

result <- as.matrix(as.mcmc(jags.result))
save(jags.result, file = 'example_model_output')


summary(j)