library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')


########## Read Stan models
# no learning
no_learning_model_stanc <- stanc(
  file = "stan-models/no_learning_model.stan",
  model_name = "no_learning_model")
no_learning_model <- stan_model(
  stanc_ret = no_learning_model_stanc)
save(no_learning_model,
  file = "stan-models/compiled/no_learning_model.Rds",
  compress = "xz", compression_level = 9)

# step learning
step_learning_model_stanc <- stanc(
  file = "stan-models/step_learning_model.stan",
  model_name = "step_learning_model")
step_learning_model <- stan_model(
  stanc_ret = step_learning_model_stanc)
save(step_learning_model,
  file = "stan-models/compiled/step_learning_model.Rds",
  compress = "xz", compression_level = 9)

# symmetric logistic learning
symmetric_logistic_learning_model_stanc <- stanc(
  file = "stan-models/symmetric_logistic_learning_model.stan",
  model_name = "symmetric_logistic_learning_model")
symmetric_logistic_learning_model <- stan_model(
  stanc_ret = symmetric_logistic_learning_model_stanc)
save(symmetric_logistic_learning_model,
  file = "stan-models/compiled/symmetric_logistic_learning_model.Rds",
  compress = "xz", compression_level = 9)

# asymmetric logistic learning
asymmetric_logistic_learning_model_stanc <- stanc(
  file = "stan-models/asymmetric_logistic_learning_model.stan",
  model_name = "asymmetric_logistic_learning_model")
asymmetric_logistic_learning_model <- stan_model(
  stanc_ret = asymmetric_logistic_learning_model_stanc)
save(asymmetric_logistic_learning_model,
  file = "stan-models/compiled/asymmetric_logistic_learning_model.Rds",
  compress = "xz", compression_level = 9)





########## Load Stan Models
load("stan-models/compiled/no_learning_model.Rds")
load("stan-models/compiled/step_learning_model.Rds")
load("stan-models/compiled/symmetric_logistic_learning_model.Rds")
load("stan-models/compiled/asymmetric_logistic_learning_model.Rds")





########## Real-World Data Fits
load(file = "experiment-1/kendal-model/exp1.Rds")
sub_id <- unique(exp1$subject_id)[7] # 78, 7, 10

yn <- exp1[(exp1$subject_id == sub_id & exp1$is_predictable == 0), "rt"]/1000
yl <- exp1[(exp1$subject_id == sub_id & exp1$is_predictable == 1), "rt"]/1000
n <- min(length(ryn), length(ryl))
yn <- ryn[seq_len(rn)]
yl <- ryl[seq_len(rn)]

fit_no_learning <- sampling(no_learning_model,
  data = list('n' = n, 'yn' = yn, 'yl' = yl),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10))
save(fit_no_learning, compress = "xz", compression_level = 9,
  file = "experiment-1/kendal-model/vignettes/model-comparison/fits/fit_no_learning.Rds")

fit_step_learning <- sampling(step_learning_model,
  data = list('n' = n, 'yn' = yn, 'yl' = yl),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10))
save(fit_step_learning, compress = "xz", compression_level = 9,
  file = "experiment-1/kendal-model/vignettes/model-comparison/fits/fit_step_learning.Rds")

fit_symmetric_logistic_learning <- sampling(symmetric_logistic_learning_model,
  data = list('n' = n, 'yn' = yn, 'yl' = yl),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10))
save(fit_symmetric_logistic_learning, compress = "xz", compression_level = 9,
  file = "experiment-1/kendal-model/vignettes/model-comparison/fits/fit_symmetric_logistic_learning.Rds")

fit_asymmetric_logistic_learning <- sampling(asymmetric_logistic_learning_model,
  data = list('n' = n, 'yn' = yn, 'yl' = yl),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10))
save(fit_asymmetric_logistic_learning, compress = "xz", compression_level = 9,
  file = "experiment-1/kendal-model/vignettes/model-comparison/fits/fit_asymmetric_logistic_learning.Rds")





########## Load Stan Fits
load(file = "experiment-1/kendal-model/vignettes/model-comparison/fits/fit_no_learning.Rds")
load(file = "experiment-1/kendal-model/vignettes/model-comparison/fits/fit_step_learning.Rds")
load(file = "experiment-1/kendal-model/vignettes/model-comparison/fits/fit_symmetric_logistic_learning.Rds")
load(file = "experiment-1/kendal-model/vignettes/model-comparison/fits/fit_asymmetric_logistic_learning.Rds")





########## Visualize Model Fits
library("ggplot2")
plot_post_pred <- function(fit, yn, yl,
                           pred_names = c("ynpred", "ylpred"),
                           plot_title = "Experimental Data Fit") {
  n <- min(length(yn), length(yl))
  ynpred <- rstan::extract(fit)[[pred_names[1]]]
  ylpred <- rstan::extract(fit)[[pred_names[2]]]
  df <- data.frame(
    x = rep(seq_len(n), 2),
    y = c(yn, yl),
    data_label = c(rep(1, n), rep(2, n)),
    pred_label = c(rep(3, n), rep(4, n)),
    y_means = c(apply(ynpred, 2, mean), apply(ylpred, 2, mean)),
    y_sds = c(apply(ynpred, 2, sd), apply(ylpred, 2, sd))
  )
  df$y_hi <- df$y_means + df$y_sds
  df$y_lo <- df$y_means - df$y_sds
  df$dummy <- rep(5, 2*n)

  # factor levels key (since R makes it alphabetical)
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
                                 "Non-Learned Fit \u00B1 1 sd",
                                 "Learned Fit \u00B1 1 sd", ""),
                      name = NULL) +
    guides(fill = guide_legend(override.aes = list(size = c(2.5, 2.5, 6, 6, 8),
                                                   shape = c(16, 16, -0x2014L, -0x2014L, -0x2014L),
                                                   color = c("#999999", "#000000", "#80b3ff", "#ff4f4f", "#ffffff"),
                                                   fill = c(NA, NA, "#80b3ff40", "#ff4f4f40", "#ffffff")))) +
    labs(title = plot_title,
         x = "Trial Number", y = "Response Time (sec)") +
    theme_bw() +
    theme(panel.border = element_blank(),
          plot.title = element_text(size = 20,
                                    margin = margin(5, 0, 15, 0, "pt")),
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
}

plot_post_pred(fit_no_learning, yn, yl, plot_title = "No Learning Model Data Fit")

plot(seq_len(10), rnorm(10), pch = "-")

##### Model Comparison


library("loo")
library("bridgesampling")


# Bayes Factor
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
loo_exp_stp <- loo(fit_exp_stp)
loo_exp_exp <- loo(fit_exp_exp)

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


ypred_el <- rstan::extract(fit_exp_log)[['ylpred']]
ypred_es <- rstan::extract(fit_exp_stp)[['ylpred']]
ypred_ee <- rstan::extract(fit_exp_exp)[['ylpred']]


plot_pred(ypred_el, nr_samples = 0, title = 'Posterior Predictions: Real Data', ylim = c(0, 2), lty = 'dashed')
  points(ryl, pch = 20)
  par(new = TRUE)
  plot_pred(ypred_es, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Exp Logistic', 'Exp Step'), lty = c(2, 1), bty = 'n', cex = 1.2)

plot_pred(ypred_el, nr_samples = 0, title = 'Posterior Predictions: Real Data', ylim = c(0, 2), lty = 'dashed')
  points(ryl, pch = 20)
  par(new = TRUE)
  plot_pred(ypred_ee, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Exp Logistic', 'Exp Exp'), lty = c(2, 1), bty = 'n', cex = 1.2)

plot_pred(ypred_ee, nr_samples = 0, title = 'Posterior Predictions: Real Data', ylim = c(0, 2), lty = 'dashed')
  points(ryl, pch = 20)
  par(new = TRUE)
  plot_pred(ypred_es, nr_samples = 0, title = '', ylim = c(0, 2), lty = 'solid', shade_col = 'skyblue')
  legend('topright', legend = c('Exp Exp', 'Exp Step'), lty = c(2, 1), bty = 'n', cex = 1.2)
