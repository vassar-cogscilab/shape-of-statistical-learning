library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')




########## Read and Load Stan models
if (.Platform$OS.type == "windows") { os <- "windows" } else { os <- "unix" }

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
library("tictoc")
load(file = "experiment-1/kendal-model/exp1.Rds")

tic()
fit <- data_fit(exp1)
toc()

fit_path <- "experiment-1/kendal-model/vignettes/model/fits/"
save(fit, compress = "xz", compression_level = 9,
  file = paste0(fit_path, "fit.Rds"))



##### Load Stan Fits
load(file = paste0(fit_path, "fit_no_learning.Rds"))
load(file = paste0(fit_path, "fit_step_learning.Rds"))
load(file = paste0(fit_path, "fit_symmetric_logistic_learning.Rds"))
load(file = paste0(fit_path, "fit_asymmetric_logistic_learning.Rds"))





########## Visualize Model Fits

plot_post_pred(fit_no_learning, nk, yn, yl,
               plot_title = "No Learning Model Data Fit")
plot_post_pred(fit_step_learning, nk, yn, yl,
               plot_title = "Step Learning Model Data Fit")
plot_post_pred(fit_symmetric_logistic_learning, nk, yn, yl,
               plot_title = "Symmetric Logistic Learning Model Data Fit")
plot_post_pred(fit_asymmetric_logistic_learning, nk, yn, yl,
               plot_title = "Asymmetric Logistic Learning Model Data Fit")

plot_model_est(fit_no_learning, nk, yn, yl,
               c("V", "E", "A", "S"),
               plot_title = "No Learning Model Data Fit")
plot_model_est(fit_step_learning, nk, yn, yl,
               c("V", "E", "A", "S", "D", "H"),
               plot_title = "Step Learning Model Data Fit")
plot_model_est(fit_symmetric_logistic_learning, nk, yn, yl,
               c("V", "E", "A", "S", "D", "L", "H"),
               plot_title = "Symmetric Logistic Learning Model Data Fit")
plot_model_est(fit_asymmetric_logistic_learning, nk, yn, yl,
               c("V", "E", "A", "S", "D", "L", "H", "NU", "C", "Q"),
               plot_title = "Asymmetric Logistic Learning Model Data Fit")



########## Model Comparison


### Loo
library("loo")

ll_no_learning <-
  rstan::extract(fit_no_learning)[["log_lik"]]
ll_step_learning <-
  rstan::extract(fit_step_learning)[["log_lik"]]
ll_symmetric_logistic_learning <-
  rstan::extract(fit_symmetric_logistic_learning)[["log_lik"]]
ll_asymmetric_logistic_learning <-
  rstan::extract(fit_asymmetric_logistic_learning)[["log_lik"]]




loo_compare_list <- list()
st <- 0
for (i in 1:k) {
  idx <- (st+1):(st+nk[i])
  loo_no_learning <- loo(ll_no_learning[, idx])
  loo_step_learning <- loo(ll_step_learning[, idx])
  loo_symmetric_logistic_learning <- loo(ll_symmetric_logistic_learning[, idx])
  loo_asymmetric_logistic_learning <- loo(ll_asymmetric_logistic_learning[, idx])

  loo_compare_list[[i]] <- loo_compare(loo_no_learning,
                                       loo_step_learning,
                                       loo_symmetric_logistic_learning,
                                       loo_asymmetric_logistic_learning)
  st = st + nk[i]
}

loo_compare_list[[1]]
loo_compare_list[[2]]
loo_compare_list[[3]]































### Bayes Factor
library("bridgesampling")

bayes_factor(
  bridge_sampler(fit_no_learning, silent = TRUE),
  bridge_sampler(fit_step_learning, silent = TRUE)
)

bayes_factor(
  bridge_sampler(fit_step_learning, silent = TRUE),
  bridge_sampler(fit_symmetric_logistic_learning, silent = TRUE)
)

bayes_factor(
  bridge_sampler(fit_symmetric_logistic_learning, silent = TRUE),
  bridge_sampler(fit_asymmetric_logistic_learning, silent = TRUE)
)


### Stacking Weights
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

stacking_weights(fit_no_learning,
                 fit_step_learning)
stacking_weights(fit_step_learning,
                 fit_symmetric_logistic_learning)
stacking_weights(fit_symmetric_logistic_learning,
                 fit_asymmetric_logistic_learning)
