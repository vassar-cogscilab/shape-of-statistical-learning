library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# Read Stan models
exp_log_model_stanc <- stanc(file = "experiment-1/kendal-model/model-comparison/exp_logistic_model.stan",
                             model_name = "exponential_logistic_model")
exp_log_model <- stan_model(stanc_ret = exp_log_model_stanc)
save(exp_log_model, file = "experiment-1/kendal-model/model-comparison/compiled-stan/exp_log_model.Rds",
     compress = "xz", compression_level = 9)
load("experiment-1/kendal-model/model-comparison/compiled-stan/exp_log_model.Rds")

exp_log_model_old_stanc <- stanc(file = "experiment-1/kendal-model/model-comparison/exp_logistic_model_old.stan",
                             model_name = "exponential_logistic_model_old")
exp_log_model_old <- stan_model(stanc_ret = exp_log_model_old_stanc)

exp_log_model_new_stanc <- stanc(file = "experiment-1/kendal-model/model-comparison/exp_logistic_model_new.stan",
                             model_name = "exponential_logistic_model_new")
exp_log_model_new <- stan_model(stanc_ret = exp_log_model_new_stanc)


sh_logn_model_stanc <- stanc(file = "experiment-1/kendal-model/model-comparison/shifted_lognormal_model.stan",
                             model_name = "shifted_lognormal_model")
sh_logn_model <- stan_model(stanc_ret = sh_logn_model_stanc)
save(sh_logn_model, file = "experiment-1/kendal-model/model-comparison/compiled-stan/sh_logn_model.Rds",
     compress = "xz", compression_level = 9)
load("experiment-1/kendal-model/model-comparison/compiled-stan/sh_logn_model.Rds")


exp_stp_model_stanc <- stanc(file = "experiment-1/kendal-model/model-comparison/exp_step_model.stan",
                             model_name = "exponential_step_model")
exp_stp_model <- stan_model(stanc_ret = exp_stp_model_stanc)
save(exp_stp_model, file = "experiment-1/kendal-model/model-comparison/compiled-stan/exp_stp_model.Rds",
     compress = "xz", compression_level = 9)
load("experiment-1/kendal-model/model-comparison/compiled-stan/exp_stp_model.Rds")


exp_exp_model_stanc <- stanc(file = "experiment-1/kendal-model/model-comparison/exp_exp_model.stan",
                             model_name = "exponential_exponential_model")
exp_exp_model <- stan_model(stanc_ret = exp_exp_model_stanc)
save(exp_exp_model, file = "experiment-1/kendal-model/model-comparison/compiled-stan/exp_exp_model.Rds",
     compress = "xz", compression_level = 9)
load("experiment-1/kendal-model/model-comparison/compiled-stan/exp_exp_model.Rds")



# Exponential-Logistic function for data simulation
exp_log <- function(n, V, E, A, D, L, H) {
  if (length(n) == 1) {
    n <- seq_len(n)
  }
  return( (V + E*exp(-A*n)) * (1 - D/(1 + exp(-L*(n-H)))) )
}
sim_exp_log <- function(n, V, E, A, D, L, H, var = 0.15) {
  if (length(n) == 1) {
    n <- seq_len(n)
  }
  sim <- exp_log(n = n, V = V, E = E, A = A, D = D, L = L, H = H) +
         rnorm(length(n), 0, var)
  sim[sim < 0] <- 0
  return(sim)
}



########## Simulated Data
n <- 100
V <- 1
E <- .25
A <- .25
D <- .1
L <- .75
H <- ceiling(n/2)
var <- .15
syn <- sim_exp_log(n = n, V = V, E = E, A = A, D = 0, L = L, H = H, var = var)
syl <- sim_exp_log(n = n, V = V, E = E, A = A, D = D, L = L, H = H, var = var)

sft_exp_log_old <- sampling(exp_log_model_old,
                        data = list('n' = n, 'yn' = syn, 'yl' = syl),
                        refresh = FALSE, chains = 1, iter = 500, seed = 2,
                        control = list(adapt_delta = 0.9, max_treedepth = 10)
                       )

sft_exp_log_new <- sampling(exp_log_model_new,
                        data = list('n' = n, 'yn' = syn, 'yl' = syl),
                        refresh = FALSE, chains = 1, iter = 500, seed = 2,
                        control = list(adapt_delta = 0.9, max_treedepth = 10)
                       )


ylpred_elo <- rstan::extract(sft_exp_log_old)[['ylpred']]
ylpred_eln <- rstan::extract(sft_exp_log_new)[['ylpred']]
plot_pred(ylpred_elo, nr_samples = 0, title = 'Posterior Predictions: Simulated Data', ylim = c(0, 2), lty = 'dashed')
  points(syl, pch = 20)
  par(new = TRUE)
  plot_pred(ylpred_eln, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Exp Logistic Old', 'Exp Logistic New'), lty = c(2, 1), bty = 'n', cex = 1.2)



library("loo")
library("bridgesampling")


# Bayes Factor
bayes_factor(
  bridge_sampler(sft_exp_log_old, silent = TRUE),
  bridge_sampler(sft_exp_log_new, silent = TRUE)
)

loo_exp_log_old <- loo(sft_exp_log_old)
loo_exp_log_new <- loo(sft_exp_log_new)
loo_compare(loo_exp_log_old, loo_exp_log_new)

stacking_weights <- function(fit1, fit2) {
  log_lik_list <- list(
    '1' = extract_log_lik(fit1),
    '2' = extract_log_lik(fit2)
  )

  r_eff_list <- list(
    '1' = relative_eff(extract_log_lik(fit1, merge_chains = FALSE)),
    '2' = relative_eff(extract_log_lik(fit2, merge_chains = FALSE))
  )

  loo_model_weights(log_lik_list, method = 'stacking', r_eff_list = r_eff_list)
}

stacking_weights(sft_exp_log_old, sft_exp_log_new)





########## Real-World Data
load(file = "experiment-1/kendal-model/exp1.Rds")
sub_ids <- unique(exp1$subject_id)[c(78, 7, 10)]
sub_ids <- sub_ids[1]

rK <- length(sub_ids)
gs_idx <- vector()
rNTI <- rep(0, rK)
rY0 <- vector()
rY1 <- vector()
for (i in seq_len(rK)) {
  temp <- exp1[exp1$subject_id == sub_ids[i], ]
  temp0 <- temp[temp$is_predictable == 0, ]$rt/1000
  temp1 <- temp[temp$is_predictable == 1, ]$rt/1000
  mm <- min(length(temp0), length(temp1))
  if (mm >= 50) {
    rNTI[i] <- mm
    rY0 <- c(rY0, temp0[seq_len(mm)])
    rY1 <- c(rY1, temp1[seq_len(mm)])
    gs_idx <- c(gs_idx, i)
  }
}

good_sub_ids <- sub_ids[gs_idx]
rNTI <- rNTI[rNTI > 0]
rNK <- sum(rNTI)
rK <- length(good_sub_ids)


fit_exp_log <- sampling(exp_log_model,
                        data = list('n' = rNTI, 'y' = rY1),
                        refresh = FALSE, chains = 1, iter = 500, seed = 2,
                        control = list(adapt_delta = 0.9, max_treedepth = 10)
                       )

fit_exp_log_old <- sampling(exp_log_model_old,
                        data = list('n' = rNTI, 'yn' = rY0, 'yl' = rY1),
                        refresh = FALSE, chains = 1, iter = 500, seed = 2,
                        control = list(adapt_delta = 0.9, max_treedepth = 10)
                       )

fit_sh_logn <- sampling(sh_logn_model,
                        data = list('n' = rNTI, 'y' = rY1),
                        refresh = FALSE, chains = 1, iter = 500, seed = 2,
                        control = list(adapt_delta = 0.9, max_treedepth = 10)
                       )

fit_exp_stp <- sampling(exp_stp_model,
                        data = list('n' = rNTI, 'y' = rY1),
                        refresh = FALSE, chains = 1, iter = 500, seed = 2,
                        control = list(adapt_delta = 0.9, max_treedepth = 10)
                       )

fit_exp_exp <- sampling(exp_exp_model,
                        data = list('n' = rNTI, 'y' = rY1),
                        refresh = FALSE, chains = 1, iter = 500, seed = 2,
                        control = list(adapt_delta = 0.9, max_treedepth = 10)
                       )





##### Model Comparison


library("loo")
library("bridgesampling")


# Bayes Factor
bayes_factor(
  bridge_sampler(fit_exp_log, silent = TRUE),
  bridge_sampler(fit_sh_logn, silent = TRUE)
)

bayes_factor(
  bridge_sampler(fit_exp_log, silent = TRUE),
  bridge_sampler(fit_exp_stp, silent = TRUE)
)

bayes_factor(
  bridge_sampler(fit_exp_log, silent = TRUE),
  bridge_sampler(fit_exp_exp, silent = TRUE)
)

bayes_factor(
  bridge_sampler(fit_exp_stp, silent = TRUE),
  bridge_sampler(fit_exp_exp, silent = TRUE)
)


# Loo
loo_exp_log <- loo(fit_exp_log)
loo_sh_logn <- loo(fit_sh_logn)
loo_exp_stp <- loo(fit_exp_stp)
loo_exp_exp <- loo(fit_exp_exp)

loo_compare(loo_exp_log, loo_sh_logn)
loo_compare(loo_exp_log, loo_exp_stp)
loo_compare(loo_exp_log, loo_exp_exp)
loo_compare(loo_exp_stp, loo_exp_exp)


# Stacking Weights
stacking_weights <- function(fit1, fit2) {
  log_lik_list <- list(
    '1' = extract_log_lik(fit1),
    '2' = extract_log_lik(fit2)
  )

  r_eff_list <- list(
    '1' = relative_eff(extract_log_lik(fit1, merge_chains = FALSE)),
    '2' = relative_eff(extract_log_lik(fit2, merge_chains = FALSE))
  )

  loo_model_weights(log_lik_list, method = 'stacking', r_eff_list = r_eff_list)
}

stacking_weights(fit_exp_log, fit_sh_logn)
stacking_weights(fit_exp_log, fit_exp_stp)
stacking_weights(fit_exp_log, fit_exp_exp)
stacking_weights(fit_exp_stp, fit_exp_exp)









##########################

plot_pred <- function(ypred, nr_samples = 25, title = '',
                      ylim = c(0, 8), shade_col = 'grey', lty = 'solid') {
  ns <- seq(ncol(ypred))
  sds <- apply(ypred, 2, sd)
  means <- apply(ypred, 2, mean)
  lo <- means - sds
  hi <- means + sds
  plot(means, type = 'l', main = title,
       ylim = ylim, panel.first =
       polygon(c(ns, rev(ns)), c(lo, rev(hi)),
       col = grDevices::adjustcolor(shade_col, alpha = .4), border = NA),
       xlab = 'Trials', ylab = 'Response Time (sec)', axes = FALSE,
       cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex = 1.5,
       lty = lty
  )

  axis(2, las = 1)
  axis(1, at = seq(0, ncol(ypred), 5))

  for (i in sample(nrow(ypred), size = nr_samples)) {
    lines(ypred[i, ], col = 'skyblue')
  }
}


ypred_elo <- rstan::extract(fit_exp_log_old)[['ypred']]
ypred_el <- rstan::extract(fit_exp_log)[['ypred']]
ypred_sl <- rstan::extract(fit_sh_logn)[['ypred']]
ypred_es <- rstan::extract(fit_exp_stp)[['ypred']]
ypred_ee <- rstan::extract(fit_exp_exp)[['ypred']]

plot_pred(ypred_elo, nr_samples = 0, title = 'Posterior Predictions: Real Data', ylim = c(0, 2), lty = 'dashed')
  points(rY1, pch = 20)
  # par(new = TRUE)
  # plot_pred(ypred_el, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Exp Logistic Old', 'Exp Logistic'), lty = c(2, 1), bty = 'n', cex = 1.2)


plot_pred(ypred_sl, nr_samples = 0, title = 'Posterior Predictions: Real Data', ylim = c(0, 2), lty = 'dashed')
  points(rY1, pch = 20)
  par(new = TRUE)
  plot_pred(ypred_el, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Shifted Lognormal', 'Exp Logistic'), lty = c(2, 1), bty = 'n', cex = 1.2)

plot_pred(ypred_el, nr_samples = 0, title = 'Posterior Predictions: Real Data', ylim = c(0, 2), lty = 'dashed')
  points(rY1, pch = 20)
  par(new = TRUE)
  plot_pred(ypred_es, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Exp Logistic', 'Exp Step'), lty = c(2, 1), bty = 'n', cex = 1.2)

plot_pred(ypred_el, nr_samples = 0, title = 'Posterior Predictions: Real Data', ylim = c(0, 2), lty = 'dashed')
  points(rY1, pch = 20)
  par(new = TRUE)
  plot_pred(ypred_ee, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Exp Logistic', 'Exp Exp'), lty = c(2, 1), bty = 'n', cex = 1.2)

plot_pred(ypred_ee, nr_samples = 0, title = 'Posterior Predictions: Real Data', ylim = c(0, 2), lty = 'dashed')
  points(rY1, pch = 20)
  par(new = TRUE)
  plot_pred(ypred_es, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Exp Exp', 'Exp Step'), lty = c(2, 1), bty = 'n', cex = 1.2)
