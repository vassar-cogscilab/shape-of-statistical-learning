#subject_id,rt,pair,target_index,t,is_predictable,subject_condition

library(runjags)
library(dplyr)
library(tidyr)
library(ggplot2)
library(coda)

data.for.jags <-list(rt = sub$rt,
  subject_id = sub$subject_id,
  is.predictable = sub$is.predictable,
  t = sub$t,
  N = length(sub$rt),
  S = length(unique(sub$subject_id)),
  C = length(unique(sub$cond)),
  condition = sub%>%group_by(subject_id) %>% summarise(cond = unique(cond)) %>% select(cond) %>% as.matrix
)
params.to.monitor <- c('sd', 'alpha1', 'alpha2', 'lambda', 'beta', 'gamma', 'delta')

jags.result <- run.jags('jags_model.txt', monitor=params.to.monitor, data=data.for.jags, n.chains=2,
                        burnin=500, sample=100, adapt=1000)

result <- as.matrix(as.mcmc(jags.result))