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

plot_model_est(fit_51_100, exp1, sub_ids = c(75, 75, 75))
fit <- fit_51_100
data <- exp1
sub_ids <- c(74, 75)
models = c("no_learning", "step_learning",
           "symmetric_logistic_learning",
           "asymmetric_logistic_learning")

mu_n <- function(nk, pars) {
  return(pars[["V"]] + pars[["E"]] * exp(-pars[["A"]] * 1:nk))
}
mu_nol <- function(nk, pars) {
  return(pars[["V"]] + pars[["S"]] + pars[["E"]] * exp(-pars[["A"]] * 1:nk))
}
mu_stp <- function(nk, pars) {
  return(
    c(pars[["V"]] + pars[["S"]] + pars[["E"]] * exp(-pars[["A"]] *
        1:(round(pars[["H"]])-1)),
      (pars[["V"]] + pars[["S"]] +
        pars[["E"]] * exp(-pars[["A"]] * round(pars[["H"]]):nk)) * (1 - pars[["D"]]))
  )
}
mu_sym <- function(nk, pars) {
  return(
    (pars[["V"]] + pars[["S"]] + pars[["E"]] * exp(-pars[["A"]] * 1:nk)) *
    (1 - pars[["D"]][ind] /
      ( 1 + exp(-pars[["L"]][ind] * (1:nk - pars[["H"]][ind])) ) )
  )
}
mu_asy <- function(nk, pars) {
  return(
    (pars[["V"]] + pars[["S"]] + pars[["E"]] * exp(-pars[["A"]] * 1:nk) ) *
    (1 - pars[["D"]] /
      ( pars[["C"]] + pars[["Q"]] * exp(-pars[["L"]] * (1:nk - pars[["H"]]))
        )^(1/pars[["NU"]]) )
  )
}

ind <- 1

yn <- data[data[["subject_id"]] == sub_ids[ind] &
           data[["is_predictable"]] == 0, "rt"]/1000
yl <- data[data[["subject_id"]] == sub_ids[ind] &
           data[["is_predictable"]] == 1, "rt"]/1000
nk <- min(length(yn), length(yl))
yn <- yn[seq_len(nk)]
yl <- yl[seq_len(nk)]
bad_idx <- which(yn < 0 | yn > 2 | yl < 0 | yl > 2)
n_bad_idx <- length(bad_idx)
if (n_bad_idx > 0) {
  yn <- yn[-bad_idx]
  yl <- yl[-bad_idx]
  nk <- nk - n_bad_idx
}

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

for (m in 1:length(models)) {
  pars[[models[m]]] <- lapply(rstan::extract(
    fit[[as.character(sub_ids[ind])]][[models[m]]])[names(
      pars[[models[m]]])], mean)
}

df <- data.frame(
  x = rep(1:nk, 2),
  y_data = c(yn, yl),
  data_label = c(rep(1, nk), rep(2, nk)),
  model_label = c(rep(3, nk), rep(4, nk)),
  y_nol = c(mu_n(nk, pars[["no_learning"]]),
            mu_nol(nk, pars[["no_learning"]])),
  y_stp = c(mu_n(nk, pars[["step_learning"]]),
            mu_stp(nk, pars[["step_learning"]])),
  y_sym = c(mu_n(nk, pars[["symmetric_logistic_learning"]]),
            mu_sym(nk, pars[["symmetric_logistic_learning"]])),
  y_asy = c(mu_n(nk, pars[["asymmetric_logistic_learning"]]),
            mu_asy(nk, pars[["asymmetric_logistic_learning"]]))
)

df[["sd_nol"]] <- df[["y_nol"]] * sqrt(c(
  rep((pars[["no_learning"]][["sigma2_n"]] - 1) *
      pars[["no_learning"]][["sigma2_n"]], nk),
  rep((pars[["no_learning"]][["sigma2_l"]] - 1) *
      pars[["no_learning"]][["sigma2_l"]], nk)
  ))
df[["sd_stp"]] <- df[["y_stp"]] * sqrt(c(
  rep((pars[["step_learning"]][["sigma2_n"]] - 1) *
      pars[["step_learning"]][["sigma2_n"]], nk),
  rep((pars[["step_learning"]][["sigma2_l"]] - 1) *
      pars[["step_learning"]][["sigma2_l"]], nk)
  ))
df[["sd_sym"]] <- df[["y_sym"]] * sqrt(c(
  rep((pars[["symmetric_logistic_learning"]][["sigma2_n"]] - 1) *
      pars[["symmetric_logistic_learning"]][["sigma2_n"]], nk),
  rep((pars[["symmetric_logistic_learning"]][["sigma2_l"]] - 1) *
      pars[["symmetric_logistic_learning"]][["sigma2_l"]], nk)
  ))
df[["sd_asy"]] <- df[["y_asy"]] * sqrt(c(
  rep((pars[["asymmetric_logistic_learning"]][["sigma2_n"]] - 1) *
      pars[["asymmetric_logistic_learning"]][["sigma2_n"]], nk),
  rep((pars[["asymmetric_logistic_learning"]][["sigma2_l"]] - 1) *
      pars[["asymmetric_logistic_learning"]][["sigma2_l"]], nk)
  ))
df[["dummy"]] <- rep(5, 2*nk)



print(ggplot(df) +
  geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
    aes(x = x, ymin = y_asy-sd_asy, ymax = y_asy+sd_asy,
        fill = factor(model_label, levels = c(3, 4)))) +
  geom_line(size = 1.25, alpha = 0.6,
    aes(x = x, y = y_asy,
        color = factor(model_label, levels = c(3, 4)))) +
  geom_point(shape = 16, size = 1.8, alpha = 0.6,
    aes(x = x, y = y_data,
        fill = factor(data_label, levels = c(1, 2)),
        color = factor(data_label, levels = c(1, 2)))) +
  geom_point(alpha = 0, size = 1, shape = 21,
    aes(x = x, y = y_data,
        fill = factor(dummy, levels = c(5)))) +
  scale_color_manual(values = c("#999999", "#000000", "#80b3ff", "#ff4f4f"),
    guide = FALSE) +
  scale_fill_manual(values = c("#999999", "#000000", "#80b3ff", "#ff4f4f", "#ffffff"),
    labels = c("Non-Learned Data", "Learned Data",
               "Non-Learned Fit \u00B1 1 sd", "Learned Fit \u00B1 1 sd", ""),
               name = NULL) +
  guides(fill = guide_legend(
    override.aes = list(
      size = c(2.5, 2.5, 6, 6, 8),
      shape = c(16, 16, -0x2014L, -0x2014L, -0x2014L),
      color = c("#999999", "#000000", "#80b3ff", "#ff4f4f", "#ffffff"),
      fill = c(NA, NA, "#80b3ff40", "#ff4f4f40", "#ffffff")))) +
  labs(title = "No Learning Model Fit",
       subtitle = paste0("Subject ", sub_ids[ind]),
       x = "Trial Number", y = "Response Time (sec)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    plot.title = element_text(size = 20,
      margin = margin(5, 0, 5, 0, "pt")),
    plot.subtitle = element_text(size = 18,
      margin = margin(5, 0, 20, 0, "pt")),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16,
      margin = margin(10, 0, 0, 0, "pt")),
    axis.title.y = element_text(size = 16,
      margin = margin(0, 10, 0, 0, "pt")),
    legend.position = c(1, 1.1),
    legend.justification = c(1, 1),
    legend.box = "horizontal",
    legend.direction = "vertical",
    legend.background = element_rect(fill = "transparent"),
    legend.text = element_text(size = 12))
)




########## Model Comparison

source("experiment-1/kendal-model/vignettes/model/model-comparison-helper-functions.R")
fit_path <- "experiment-1/kendal-model/vignettes/model/fits/"
load(file = paste0(fit_path, "fit.Rds"))


### Loo
tic()
loo_fits <- loo_fit(fit)
toc()


save(loo_fits, compress = "xz", compression_level = 9,
     file = paste0(fit_path, "loo_fits.Rds"))
load(file = paste0(fit_path, "loo_fits.Rds"))

plot_loo_fit(loo_fits,
             save_path = "experiment-1/kendal-model/vignettes/model/fits")


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

###

loo_sig <- loo_significant(loo_fits)

loo_sig



plot_model_est(fit, exp1, sub_ids = names(loo_sig))
