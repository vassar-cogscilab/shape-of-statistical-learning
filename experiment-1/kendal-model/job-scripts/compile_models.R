library("rstan")
options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')




########## Read Stan models
raw_dir <- "home/kendalfoster/sosl/stan-models/raw/"
cmp_dir <- "home/kendalfoster/sosl/stan-models/cmp/"

# no learning
model_name <- "no_learning_model"
no_learning_model_stanc <- stanc(
  file = paste0(raw_dir, model_name, ".stan"),
  model_name = model_name)
no_learning_model <- stan_model(
  stanc_ret = no_learning_model_stanc)
save(no_learning_model,
  file = paste0(cmp_dir, model_name, ".Rds"),
  compress = "xz", compression_level = 9)

# step learning
model_name <- "step_learning_model"
step_learning_model_stanc <- stanc(
  file = paste0(raw_dir, model_name, ".stan"),
  model_name = model_name)
step_learning_model <- stan_model(
  stanc_ret = step_learning_model_stanc)
save(step_learning_model,
  file = paste0(cmp_dir, model_name, ".Rds"),
  compress = "xz", compression_level = 9)

# symmetric logistic learning
model_name <- "symmetric_logistic_learning_model"
symmetric_logistic_learning_model_stanc <- stanc(
  file = paste0(raw_dir, model_name, ".stan"),
  model_name = model_name)
symmetric_logistic_learning_model <- stan_model(
  stanc_ret = symmetric_logistic_learning_model_stanc)
save(symmetric_logistic_learning_model,
  file = paste0(cmp_dir, model_name, ".Rds"),
  compress = "xz", compression_level = 9)

# asymmetric logistic learning
model_name <- "asymmetric_logistic_learning_model"
asymmetric_logistic_learning_model_stanc <- stanc(
  file = paste0(raw_dir, model_name, ".stan"),
  model_name = model_name)
asymmetric_logistic_learning_model <- stan_model(
  stanc_ret = asymmetric_logistic_learning_model_stanc)
save(asymmetric_logistic_learning_model,
  file = paste0(cmp_dir, model_name, ".Rds"),
  compress = "xz", compression_level = 9)
