setwd('~/shape-of-statistical-learning')

library(runjags)
library(parallel)

burnin = 5000
adapt = 5000
sample = 100
thin = 50

load('experiment-3/data/generated/jags-data-exp-3.Rdata')

item.level <- c('alpha', 'item.learned', 'item.beta.learn', 'item.gamma.learn', 'item.delta')
subject.level <- c('sigma', 'subject.adapt.beta', 'subject.adapt.gamma', 'subject.prob.learn', 'subject.beta.learn.mode', 'subject.beta.learn.concentration', 'subject.gamma.learn.mode', 'subject.gamma.learn.sd', 'subject.delta.learn.mode', 'subject.delta.learn.concentration')
condition.level <- c('condition.alpha.mu', 'condition.alpha.sigma', 'condition.adapt.beta.mode', 'condition.adapt.beta.concentration', 'condition.adapt.gamma.mode', 'condition.adapt.gamma.sd', 'condition.prob.learn.mode', 'condition.prob.learn.concentration', 'condition.beta.learn.mode.mode', 'condition.beta.learn.mode.concentration', 'condition.beta.learn.concentration.mode', 'condition.beta.learn.concentration.sd', 'condition.gamma.learn.mode.mode', 'condition.gamma.learn.mode.sd', 'condition.gamma.learn.sd.mode', 'condition.gamma.learn.sd.sd', 'condition.delta.learn.mode.mode', 'condition.delta.learn.mode.concentration', 'condition.delta.learn.concentration.mode', 'condition.delta.learn.concentration.sd') 

params.to.monitor <- c(item.level, subject.level, condition.level)

jags.result <- run.jags('JAGSmodels/jags-model-expexp-latent-learner-pairs.txt', monitor=params.to.monitor, data=data.for.jags, n.chains=detectCores(),
                        burnin=burnin, sample=sample, adapt=adapt, thin=thin, method="parallel")

save(jags.result, file = paste0('experiment-3/data/generated/jags/exp-3-', sample(1:10000000, 1), '.Rdata'))