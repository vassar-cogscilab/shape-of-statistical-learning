library("rstan")
options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')




##### Load Stan Models
cmp_dir <- "home/kendalfoster/sosl/stan-models/cmp/"

load(paste0(cmp_dir, "no_learning_model.Rds"))
load(paste0(cmp_dir, "step_learning_model.Rds"))
load(paste0(cmp_dir, "symmetric_logistic_learning_model.Rds"))
load(paste0(cmp_dir, "asymmetric_logistic_learning_model.Rds"))


##### Read Data
source("home/kendalfoster/sosl/helper-functions/fit-functions.R")
fit_path <- "home/kendalfoster/sosl/exp-1/fit.Rds"
load(file = "home/kendalfoster/sosl/exp-1/exp1.Rds")


##### Run and Save Fits
fit <- data_fit(exp1, n_chains = 4, n_iterations = 10000,
  adapt_delta = 0.9, max_treedepth = 15)
save(fit, compress = "xz", compression_level = 9, file = fit_path)
