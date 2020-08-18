library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')





########## Read and Load Stan models
if (.Platform$OS.type == "windows") { os <- "windows" } else { os <- "unix" }

# no learning
no_learning_model_stanc <- stanc(
  file = "stan-models/no_learning_model.stan",
  model_name = "no_learning_model")
no_learning_model <- stan_model(
  stanc_ret = no_learning_model_stanc)
save(no_learning_model,
  file = paste0("stan-models/compiled/", os,
    "/no_learning_model.Rds"),
  compress = "xz", compression_level = 9)

# step learning
step_learning_model_stanc <- stanc(
  file = "stan-models/step_learning_model.stan",
  model_name = "step_learning_model")
step_learning_model <- stan_model(
  stanc_ret = step_learning_model_stanc)
save(step_learning_model,
  file = paste0("stan-models/compiled/", os,
    "/step_learning_model.Rds"),
  compress = "xz", compression_level = 9)

# symmetric logistic learning
symmetric_logistic_learning_model_stanc <- stanc(
  file = "stan-models/symmetric_logistic_learning_model.stan",
  model_name = "symmetric_logistic_learning_model")
symmetric_logistic_learning_model <- stan_model(
  stanc_ret = symmetric_logistic_learning_model_stanc)
save(symmetric_logistic_learning_model,
  file = paste0("stan-models/compiled/", os,
    "/symmetric_logistic_learning_model.Rds"),
  compress = "xz", compression_level = 9)

# # symmetric logistic learning bored
# symmetric_logistic_learning_bored_model_stanc <- stanc(
#   file = "stan-models/symmetric_logistic_learning_bored_model.stan",
#   model_name = "symmetric_logistic_learning_bored_model")
# symmetric_logistic_learning_bored_model <- stan_model(
#   stanc_ret = symmetric_logistic_learning_bored_model_stanc)
# save(symmetric_logistic_learning_bored_model,
#   file = paste0("stan-models/compiled/", os,
#     "/symmetric_logistic_learning_bored_model.Rds"),
#   compress = "xz", compression_level = 9)
#
# # symmetric logistic learning return
# symmetric_logistic_learning_return_model_stanc <- stanc(
#   file = "stan-models/symmetric_logistic_learning_return_model.stan",
#   model_name = "symmetric_logistic_learning_return_model")
# symmetric_logistic_learning_return_model <- stan_model(
#   stanc_ret = symmetric_logistic_learning_return_model_stanc)
# save(symmetric_logistic_learning_return_model,
#   file = paste0("stan-models/compiled/", os,
#     "/symmetric_logistic_learning_return_model.Rds"),
#   compress = "xz", compression_level = 9)

# asymmetric logistic learning
asymmetric_logistic_learning_model_stanc <- stanc(
  file = "stan-models/asymmetric_logistic_learning_model.stan",
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
fit_path <- "experiment-1/kendal-model/vignettes/model/fits/"

load(file = "experiment-1/kendal-model/exp1.Rds")
sub_ids <- unique(exp1$subject_id)[c(78, 7, 10)] # 78, 7, 10

k <- length(sub_ids)
nk <- vector()
yn <- vector()
yl <- vector()
for(i in 1:k) {
  tempn <- exp1[exp1[["subject_id"]] == sub_ids[i] &
                exp1[["is_predictable"]] == 0, "rt"]/1000
  templ <- exp1[exp1[["subject_id"]] == sub_ids[i] &
                exp1[["is_predictable"]] == 1, "rt"]/1000
  mi <- min(length(tempn), length(templ))
  if (mi >= 50) {
    nk <- c(nk, mi)
    yn <- c(yn, tempn[seq_len(mi)])
    yl <- c(yl, templ[seq_len(mi)])
  }
}



fit_no_learning <- sampling(no_learning_model,
  data = list('k' = k, 'nk' = as.array(nk), 'total_length' = sum(nk),
              'yn' = yn, 'yl' = yl),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10))
save(fit_no_learning, compress = "xz", compression_level = 9,
  file = paste0(fit_path, "fit_no_learning.Rds"))

fit_step_learning <- sampling(step_learning_model,
  data = list('k' = k, 'nk' = as.array(nk), 'total_length' = sum(nk),
              'yn' = yn, 'yl' = yl),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10))
save(fit_step_learning, compress = "xz", compression_level = 9,
  file = paste0(fit_path, "fit_step_learning.Rds"))

fit_symmetric_logistic_learning <- sampling(symmetric_logistic_learning_model,
  data = list('k' = k, 'nk' = as.array(nk), 'total_length' = sum(nk),
              'yn' = yn, 'yl' = yl),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10))
save(fit_symmetric_logistic_learning, compress = "xz", compression_level = 9,
  file = paste0(fit_path, "fit_symmetric_logistic_learning.Rds"))

# fit_symmetric_logistic_learning_bored <- sampling(symmetric_logistic_learning_bored_model,
#   data = list('n' = n, 'yn' = yn, 'yl' = yl),
#   refresh = FALSE, chains = 1, iter = 500, seed = 2,
#   control = list(adapt_delta = 0.9, max_treedepth = 10))
# save(fit_symmetric_logistic_learning_bored, compress = "xz", compression_level = 9,
#   file = paste0(fit_path, "fit_symmetric_logistic_learning_bored.Rds"))
#
# fit_symmetric_logistic_learning_return <- sampling(symmetric_logistic_learning_return_model,
#   data = list('n' = n, 'yn' = yn, 'yl' = yl),
#   refresh = FALSE, chains = 1, iter = 500, seed = 2,
#   control = list(adapt_delta = 0.9, max_treedepth = 10))
# save(fit_symmetric_logistic_learning_return, compress = "xz", compression_level = 9,
#   file = paste0(fit_path, "fit_symmetric_logistic_learning_return.Rds"))

fit_asymmetric_logistic_learning <- sampling(asymmetric_logistic_learning_model,
  data = list('k' = k, 'nk' = as.array(nk), 'total_length' = sum(nk),
              'yn' = yn, 'yl' = yl),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10))
save(fit_asymmetric_logistic_learning, compress = "xz", compression_level = 9,
  file = paste0(fit_path, "fit_asymmetric_logistic_learning.Rds"))



##### Load Stan Fits
load(file = paste0(fit_path, "fit_no_learning.Rds"))
load(file = paste0(fit_path, "fit_step_learning.Rds"))
load(file = paste0(fit_path, "fit_symmetric_logistic_learning.Rds"))
load(file = paste0(fit_path, "fit_asymmetric_logistic_learning.Rds"))





########## Visualize Model Fits
library("ggplot2")
plot_post_pred <- function(fit, nk, yn, yl,
                           pred_names = c("ynpred", "ylpred"),
                           plot_title = "Experimental Data Fit") {
  ynpred <- rstan::extract(fit)[[pred_names[1]]]
  ylpred <- rstan::extract(fit)[[pred_names[2]]]

  k <- length(nk)
  st <- 0
  for (i in 1:k) {
    n <- nk[i]
    ynpredi <- ynpred[, (st+1):(st+nk[i])]
    ylpredi <- ylpred[, (st+1):(st+nk[i])]
    yni <- yn[(st+1):(st+nk[i])]
    yli <- yl[(st+1):(st+nk[i])]
    st = st + nk[i]

    df <- data.frame(
      x = rep(seq_len(n), 2),
      y = c(yni, yli),
      data_label = c(rep(1, n), rep(2, n)),
      pred_label = c(rep(3, n), rep(4, n)),
      y_means = c(apply(ynpredi, 2, mean), apply(ylpredi, 2, mean)),
      y_sds = c(apply(ynpredi, 2, sd), apply(ylpredi, 2, sd))
    )
    df$y_hi <- df$y_means + df$y_sds
    df$y_lo <- df$y_means - df$y_sds
    df$dummy <- rep(5, 2*n)

    # factor levels key (since R makes it alphabetical otherwise)
      # 1: Non-Learned Data
      # 2: Learned Data
      # 3: Non-Learned Fit
      # 4: Learned Fit

    print(ggplot(df) +
      geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
        aes(x = x, ymin = y_lo, ymax = y_hi,
            fill = factor(pred_label, levels = c(3, 4)))) +
      geom_line(size = 1.25, alpha = 0.6,
        aes(x = x, y = y_means,
            color = factor(pred_label, levels = c(3, 4)))) +
      geom_point(shape = 16, size = 1.8, alpha = 0.6,
        aes(x = x, y = y,
            fill = factor(data_label, levels = c(1, 2)),
            color = factor(data_label, levels = c(1, 2)))) +
      geom_point(alpha = 0, size = 1, shape = 21,
        aes(x = x, y = y,
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
      labs(title = plot_title,
           x = "Trial Number", y = "Response Time (sec)") +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        plot.title = element_text(size = 20,
          margin = margin(5, 0, 20, 0, "pt")),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16,
          margin = margin(10, 0, 0, 0, "pt")),
        axis.title.y = element_text(size = 16,
          margin = margin(0, 10, 0, 0, "pt")),
        legend.position = c(1, 1.05),
        legend.justification = c(1, 1),
        legend.box = "horizontal",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12))
    )
  }
}


### Print Plots
plot_post_pred(fit_no_learning, nk, yn, yl,
               plot_title = "No Learning Model Data Fit")
plot_post_pred(fit_step_learning, nk, yn, yl,
               plot_title = "Step Learning Model Data Fit")
plot_post_pred(fit_symmetric_logistic_learning, nk, yn, yl,
               plot_title = "Symmetric Logistic Learning Model Data Fit")
# plot_post_pred(fit_symmetric_logistic_learning_bored, yn, yl,
#                plot_title = "Symmetric Logistic Learning Bored Model Data Fit")
# plot_post_pred(fit_symmetric_logistic_learning_return, yn, yl,
#                plot_title = "Symmetric Logistic Learning Return Model Data Fit")
plot_post_pred(fit_asymmetric_logistic_learning, nk, yn, yl,
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
