library(runjags)
library(dplyr)
library(tidyr)
library(ggplot2)
library(coda)

data <- read_csv('experiment-1/data/generated/for-jags-exp-1.csv')

data <- data %>% filter(subject_id <= 20)

data.for.jags <- list(
  rt = data$rt,
  subject_id = data$subject_id,
  is.predictable = data$is_predictable,
  t = data$t,
  N = length(data$rt),
  S = length(unique(data$subject_id)),
  C = length(unique(data$subject_condition)),
  condition = data$subject_condition,
  max_t = (data%>%group_by(subject_id)%>% summarise(max_t = max(t)) %>% select(max_t) %>% as.matrix)[,1],
  n.pairs = (data %>% group_by(subject_id) %>% summarise(n.pairs = length(unique(pair))) %>% select(n.pairs) %>% as.matrix)[,1]
)

params.to.monitor <- c('sigma', 'alpha', 'gamma', 'beta', 'beta.learn', 'gamma.learn', 'delta')

jags.result <- run.jags('JAGSmodels/jags-model-expexp.txt', monitor=params.to.monitor, data=data.for.jags, n.chains=2,
                        burnin=1000, sample=5000, adapt=2000)

result <- as.matrix(as.mcmc(jags.result))
save(jags.result, file = 'experiment-1/data/generated/jags/example_model_output.Rdata')