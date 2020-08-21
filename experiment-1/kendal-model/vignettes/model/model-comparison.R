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

### Load Stan Fits
load(file = paste0(fit_path, "fit.Rds"))



########## Model Comparison

fit_path <- "experiment-1/kendal-model/vignettes/model/fits/"
load(file = paste0(fit_path, "fit_251_258.Rds"))
source("experiment-1/kendal-model/vignettes/model/model-comparison-helper-functions.R")



### Loo
fit_loo_comp <- loo_fit(fit_251_258)
plot_loo_fit(fit_loo_comp)





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
