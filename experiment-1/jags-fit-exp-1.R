library(runjags)
library(dplyr)
library(tidyr)
library(ggplot2)
library(coda)

load('experiment-1/data/generated/jags-data-exp-1.Rdata')

params.to.monitor <- c('sigma', 'alpha', 'gamma', 'beta', 'beta.learn', 'gamma.learn', 'delta', 'is.learner')

jags.result <- run.jags('JAGSmodels/jags-model-expexp-latent-learner.txt', monitor=params.to.monitor, data=data.for.jags, n.chains=2,
                        burnin=1000, sample=5000, adapt=2000, thin=2, method="para")

result <- as.matrix(as.mcmc(jags.result))
save(jags.result, file = 'experiment-1/data/generated/jags/example_model_output.Rdata')