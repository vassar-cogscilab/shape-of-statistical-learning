library(readr)
library(rstan)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data <- read_csv('experiment-1/data/generated/for-jags-exp-1.csv')

data <- data %>% filter(subject_id <= 20)

model_data <- list(
  n = nrow(data),
  s = max(data$subject_id),
  t = data$t - 1,
  is_predictable = data$is_predictable,
  subject = data$subject_id,
  target_id = data$target_index,
  rt = data$rt
)

fit <- stan(file="experiment-1/stan-model-exp-1.stan", data=model_data, iter=4000, chains=2,
            pars=c('sigma_subject', 'subject_intercept', 'beta_adapt', 'rate_adapt', 'beta_learn', 'rate_learn', 'midpoint_learn'),
            control = list(adapt_delta = 0.90))
pairs(fit, pars=c("sigma_subject"))
