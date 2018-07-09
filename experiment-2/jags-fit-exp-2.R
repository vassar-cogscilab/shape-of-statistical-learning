library(runjags)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(coda)

data <- read_csv('experiment-2/data/generated/for-jags-exp-2.csv')

data <- data %>% filter(subject_id <= 20)

data.for.jags <- list(
  rt = data$rt,
  subject_id = data$subject_id,
  is_predictable = data$is_predictable,
  pair = data$pair,
  t = data$t,
  N = length(data$rt),
  S = length(unique(data$subject_id)),
  C = length(unique(data$subject_condition)),
  condition = data$subject_condition,
  max_t = (data%>%group_by(subject_id)%>% summarise(max_t = max(t)) %>% select(max_t) %>% as.matrix)[,1],
  P = (data %>% group_by(subject_id) %>% summarise(n.pairs = length(unique(pair))) %>% select(n.pairs) %>% as.matrix)[,1]
)

params.to.monitor <- c('sigma', 'alpha', 'subject.adapt.gamma', 'subject.adapt.beta', 'item.learned', 'item.beta.learn', 'item.gamma.learn', 'item.delta', 'subject.prob.learn')

jags.result <- run.jags('JAGSmodels/jags-model-expexp-latent-learner-pairs.txt', monitor=params.to.monitor, data=data.for.jags, n.chains=2,
                        burnin=1000, sample=5000, adapt=2000, thin=2, method="parallel")

result <- as.matrix(as.mcmc(jags.result))
save(jags.result, file = 'experiment-2/data/generated/jags/example_model_output.Rdata')