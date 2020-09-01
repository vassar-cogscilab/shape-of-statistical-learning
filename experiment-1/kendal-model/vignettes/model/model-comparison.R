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
fit_171_74 <- data_fit(exp1, sub_ids = c(171, 74), n_chains = 4, n_iterations = 5000, adapt_delta = 0.9, max_treedepth = 15)
toc()

fit_path <- "experiment-1/kendal-model/vignettes/model/fits/"
save(fit_171_74, compress = "xz", compression_level = 9,
     file = paste0(fit_path, "pieces/fit_171_74.Rds"))

### Load Stan Fits
load(file = paste0(fit_path, "fit.Rds"))


pars <- list(
  "no_learning" = c("V" = 0, "E" = 0, "A" = 0, "S" = 0,
                       "sigma2_n" = 0, "sigma2_l" = 0),
  "step_learning" = c("V" = 0, "E" = 0, "A" = 0, "S" = 0,
                      "D" = 0, "H" = 0,
                      "sigma2_n" = 0, "sigma2_l" = 0),
  "symmetric_logistic_learning" = c("V" = 0, "E" = 0, "A" = 0, "S" = 0,
                                    "D" = 0, "L" = 0, "H" = 0,
                                    "sigma2_n" = 0, "sigma2_l" = 0),
  "asymmetric_logistic_learning" = c("V" = 0, "E" = 0, "A" = 0, "S" = 0,
                                     "D" = 0, "L" = 0, "H" = 0,
                                     "NU" = 0, "C" = 0, "Q" = 0,
                                     "sigma2_n" = 0, "sigma2_l" = 0)
)
npars <- length(pars[[3]])
parm <- pars
parm[["no_learning"]][["sigma2_n"]]
names(parm) <- pars[[3]]
pars[["no_learning"]] <- apply(rstan::extract(fit_251_258[["251"]][["no_learning"]])[[names(pars[["no_learning"]])]], 2, mean)
parm[["no_learning"]] <- lapply(rstan::extract(fit_251_258[["251"]][["no_learning"]])[names(pars[["no_learning"]])], mean)





########## Model Comparison

source("experiment-1/kendal-model/vignettes/model/model-comparison-helper-functions.R")
fit_path <- "experiment-1/kendal-model/vignettes/model/fits/"
load(file = paste0(fit_path, "pieces/fit_171_74.Rds"))


### Loo
tic()
loo_fits_171_74 <- loo_fit(fit_171_74, sub_ids = c(171, 74))
toc()


save(loo_fits_171_74, compress = "xz", compression_level = 9,
     file = paste0(fit_path, "loo_fits_171_74.Rds"))
load(file = paste0(fit_path, "loo_fits.Rds"))

plot_loo_fit(loo_fits_171_74)
             # save_path = "experiment-1/kendal-model/vignettes/model/fits")


########## Visualize Model Fits

plot_post_pred(fit_171_74, exp1)
plot_model_est(fit_171_74, exp1)


###

loo_sig <- loo_significant(loo_fits)
loo_sig
