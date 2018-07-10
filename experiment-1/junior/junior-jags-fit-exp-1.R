setwd('~/shape-of-statistical-learning')

library(runjags)
library(parallel)

burnin = 5000
adapt = 5000
sample = 500
thin = 50

load('experiment-1/data/generated/jags-data-exp-1.Rdata')

params.to.monitor <- c('sigma', 'alpha', 'gamma', 'beta', 'beta.learn', 'gamma.learn', 'delta', 'is.learner')

jags.result <- run.jags('JAGSmodels/jags-model-expexp-latent-learner.txt', monitor=params.to.monitor, data=data.for.jags, n.chains=detectCores(),
                        burnin=burnin, sample=sample, adapt=adapt, thin=thin, method="parallel")

save(jags.result, file = paste0('experiment-1/data/generated/jags/', sample(1:10000000, 1), '.Rdata'))