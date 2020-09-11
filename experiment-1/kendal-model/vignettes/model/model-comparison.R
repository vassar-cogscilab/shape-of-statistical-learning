library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')




########## Read and Load Stan models
os <- ifelse(.Platform$OS.type == "windows", "windows", "unix")

# no learning
no_learning_model_stanc <- stanc(
  file = "stan-models/individual/no_learning_model.stan",
  model_name = "no_learning_model")
no_learning_model <- stan_model(
  stanc_ret = no_learning_model_stanc)
save(no_learning_model,
  file = paste0("stan-models/compiled/", os,
    "/no_learning_model.Rds"),
  compress = "xz", compression_level = 9)

# step learning
step_learning_model_stanc <- stanc(
  file = "stan-models/individual/step_learning_model.stan",
  model_name = "step_learning_model")
step_learning_model <- stan_model(
  stanc_ret = step_learning_model_stanc)
save(step_learning_model,
  file = paste0("stan-models/compiled/", os,
    "/step_learning_model.Rds"),
  compress = "xz", compression_level = 9)

# symmetric logistic learning
symmetric_logistic_learning_model_stanc <- stanc(
  file = "stan-models/individual/symmetric_logistic_learning_model.stan",
  model_name = "symmetric_logistic_learning_model")
symmetric_logistic_learning_model <- stan_model(
  stanc_ret = symmetric_logistic_learning_model_stanc)
save(symmetric_logistic_learning_model,
  file = paste0("stan-models/compiled/", os,
    "/symmetric_logistic_learning_model.Rds"),
  compress = "xz", compression_level = 9)

# asymmetric logistic learning
asymmetric_logistic_learning_model_stanc <- stanc(
  file = "stan-models/individual/asymmetric_logistic_learning_model.stan",
  model_name = "asymmetric_logistic_learning_model")
asymmetric_logistic_learning_model <- stan_model(
  stanc_ret = asymmetric_logistic_learning_model_stanc)
save(asymmetric_logistic_learning_model,
  file = paste0("stan-models/compiled/", os,
    "/asymmetric_logistic_learning_model.Rds"),
  compress = "xz", compression_level = 9)


##### Load Stan Models
load(paste0("stan-models/compiled/", os,
  "/no_learning_model.Rds"))
load(paste0("stan-models/compiled/", os,
  "/step_learning_model.Rds"))
load(paste0("stan-models/compiled/", os,
  "/symmetric_logistic_learning_model.Rds"))
load(paste0("stan-models/compiled/", os,
  "/asymmetric_logistic_learning_model.Rds"))





########## Real-World Data Fits
source("experiment-1/kendal-model/vignettes/model/model-comparison-helper-functions.R")
fit_path <- "experiment-1/kendal-model/vignettes/model/fits/"
library("tictoc")
load(file = "experiment-1/kendal-model/exp1.Rds")

tic()
fit_74 <- data_fit(exp1, sub_ids = c(74), n_chains = 1, n_iterations = 1000, adapt_delta = 0.9, max_treedepth = 10, models = c("step_learning"))
toc()


tic()
fit_171 <- data_fit(exp1, sub_ids = c(171), n_chains = 1, n_iterations = 1000, adapt_delta = 0.9, max_treedepth = 10, models = c("no_learning"))
toc()
save(fit_171, compress = "xz", compression_level = 9,
     file = paste0(fit_path, "pieces/fit_171.Rds"))
load(file = paste0(fit_path, "pieces/fit_171.Rds"))

library("shinystan")
launch_shinystan(fit_171[["171"]][["symmetric_logistic_learning"]])


fit_74[["74"]][["step_learning"]]
fit_171[["171"]][["no_learning"]]

pairs(fit_171[["171"]][["no_learning"]],
  pars = c("V", "E", "A", "S"))
pairs(fit_171[["171"]][["symmetric_logistic_learning"]],
  pars = c("V", "E", "A", "S", "D", "L", "H"))


########## Model Comparison

source("experiment-1/kendal-model/vignettes/model/model-comparison-helper-functions.R")
fit_path <- "experiment-1/kendal-model/vignettes/model/fits/"
load(file = paste0(fit_path, "pieces/fit_171.Rds"))
library("tictoc")


### LOO
tic()
loo_fits_171 <- loo_fit(fit_171, sub_ids = c(171))
toc()
save(loo_fits_171, compress = "xz", compression_level = 9,
     file = paste0(fit_path, "pieces/loo_fits_171.Rds"))
load(file = paste0(fit_path, "pieces/loo_fits_171.Rds"))

loo1no <- loo(fit_171[["171"]][["no_learning"]], save_psis = TRUE)
loo1sym <- loo(fit_171[["171"]][["symmetric_logistic_learning"]], save_psis = TRUE)

loo_compare(loo1no, loo1sym)

plot(loo1no)
plot(loo1sym)

library("bayesplot")
ppc_loo_pit_overlay(
  y = exp1[exp1[["subject_id"]] == 171 & exp1[["is_predictable"]] == 1, "rt"]/1000,
  yrep = rstan::extract(fit_171[["171"]][["no_learning"]])[["ylpred"]],
  lw = weights(loo1no$psis_object)
)
ppc_loo_pit_overlay(
  y = exp1[exp1[["subject_id"]] == 171 & exp1[["is_predictable"]] == 1, "rt"]/1000,
  yrep = rstan::extract(fit_171[["171"]][["symmetric_logistic_learning"]])[["ylpred"]],
  lw = weights(loo1sym$psis_object)
)



fit_74[["74"]][["no_learning"]]
fit_74[["74"]][["step_learning"]]
fit_74[["74"]][["symmetric_logistic_learning"]]
fit_74[["74"]][["asymmetric_logistic_learning"]]

fit_171[["171"]][["no_learning"]]
fit_171[["171"]][["step_learning"]]
fit_171[["171"]][["symmetric_logistic_learning"]]
fit_171[["171"]][["asymmetric_logistic_learning"]]



library("bridgesampling")

log10(bayes_factor(
  bridge_sampler(fit_74[["74"]][["no_learning"]], silent = TRUE),
  bridge_sampler(fit_74[["74"]][["step_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_74[["74"]][["no_learning"]], silent = TRUE),
  bridge_sampler(fit_74[["74"]][["symmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_74[["74"]][["no_learning"]], silent = TRUE),
  bridge_sampler(fit_74[["74"]][["asymmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_74[["74"]][["step_learning"]], silent = TRUE),
  bridge_sampler(fit_74[["74"]][["symmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_74[["74"]][["step_learning"]], silent = TRUE),
  bridge_sampler(fit_74[["74"]][["asymmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_74[["74"]][["symmetric_logistic_learning"]], silent = TRUE),
  bridge_sampler(fit_74[["74"]][["asymmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])



log10(bayes_factor(
  bridge_sampler(fit_171[["171"]][["no_learning"]], silent = TRUE),
  bridge_sampler(fit_171[["171"]][["step_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_171[["171"]][["no_learning"]], silent = TRUE),
  bridge_sampler(fit_171[["171"]][["symmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_171[["171"]][["no_learning"]], silent = TRUE),
  bridge_sampler(fit_171[["171"]][["asymmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_171[["171"]][["step_learning"]], silent = TRUE),
  bridge_sampler(fit_171[["171"]][["symmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_171[["171"]][["step_learning"]], silent = TRUE),
  bridge_sampler(fit_171[["171"]][["asymmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])
log10(bayes_factor(
  bridge_sampler(fit_171[["171"]][["symmetric_logistic_learning"]], silent = TRUE),
  bridge_sampler(fit_171[["171"]][["asymmetric_logistic_learning"]], silent = TRUE)
)[["bf"]])



### WAIC
waic1no <- waic(log_lik_extract(fit_171[["171"]][["no_learning"]]))
waic1sym <- loo(fit_171[["171"]][["symmetric_logistic_learning"]], save_psis = TRUE)
waic_compare(waic1no, loo1sym)



waic_fits_171 <- waic_fit(fit_171, sub_ids = c(171))


plot_loo_fit(loo_fits_171)
plot_waic_fit(waic_fits_171)




########## Visualize Model Fits
load(file = "experiment-1/kendal-model/exp1.Rds")

plot_post_pred(fit_74, exp1, models = "step_learning")
plot_model_est(fit_74, exp1, models = "step_learning")



plot_post_pred(fit_171, exp1)
plot_model_est(fit_171, exp1, models = "no_learning")



###

loo_sig <- loo_significant(loo_fits_171)
loo_sig



########## Plot Prior Estimates
source("experiment-1/kendal-model/vignettes/model/model-comparison-helper-functions.R")
plot_prior_est(ntrials = 100, ndraws = 1000)


y <- seq(0, 1, by = 0.01)
plot(y, dbeta(y, 10, 8))

x <- seq(0, 4, by = 0.01)
plot(x, dgamma(x, 12, 36))
  abline(v = .1, lty = "dashed", col = "blue")
  abline(v = .75, lty = "dashed", col = "blue")

t <- seq_len(70)
mu_sym <- function(t, pars) {
  return(
    (pars[["V"]] + pars[["S"]] + pars[["E"]] * exp(-pars[["A"]] * t)) *
    (1 - pars[["D"]] /
      ( 1 + exp(-pars[["L"]] * (t - pars[["H"]])) ) )
  )
}
pars <- list("V" = 1, "E" = .25, "A" = .25, "S" = 0,
             "D" = .75, "L" = .75, "H" = 40)
plot(t, mu_sym(t, pars))
