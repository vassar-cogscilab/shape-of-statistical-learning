######################### Data Fitting #########################################
data_fit <- function(data, sub_ids = NULL,
                     n_chains = 1, n_iterations = 500,
                     adapt_delta = 0.9, max_treedepth = 10) {
  if(is.null(sub_ids)) {
    sub_ids <- sort(unique(data[["subject_id"]]))
  }
  k <- length(sub_ids)
  out <- list()

  for (i in sub_ids) {
    # read data
    yn <- data[data[["subject_id"]] == i &
               data[["is_predictable"]] == 0, "rt"]/1000
    yl <- data[data[["subject_id"]] == i &
               data[["is_predictable"]] == 1, "rt"]/1000

    # clean data
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

    # run fits
    out[[as.character(i)]] <- list(
      "no_learning" =
        rstan::sampling(no_learning_model,
          data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
          refresh = FALSE, chains = n_chains, iter = n_iterations,
          control = list(adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth),
          init = rep( list( list(
            V = 1, E = 0.25, A = 0.25, S = 0) ),
            n_chains) ),
      "step_learning" =
        rstan::sampling(step_learning_model,
          data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
          refresh = FALSE, chains = n_chains, iter = n_iterations,
          control = list(adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth),
          init = rep( list( list(
            V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, H_raw = 0.5) ),
            n_chains) ),
      "symmetric_logistic_learning" =
        rstan::sampling(symmetric_logistic_learning_model,
          data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
          refresh = FALSE, chains = n_chains, iter = n_iterations,
          control = list(adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth),
          init = rep( list( list(
            V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, L = 6, H_raw = 0.5) ),
            n_chains) ),
      "asymmetric_logistic_learning" =
        rstan::sampling(asymmetric_logistic_learning_model,
          data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
          refresh = FALSE, chains = n_chains, iter = n_iterations,
          control = list(adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth),
          init = rep( list( list(
            V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, L = 6, H_raw = 0.5,
            NU = 1, C = 1, Q = 1) ), n_chains) )
    )
  }
  return(out)
}


# load(file = "experiment-1/kendal-model/exp1.Rds")
# sub_ids <- unique(exp1$subject_id) # 74, 171, 143
#
# k <- length(sub_ids)
# nk <- vector()
# yn <- vector()
# yl <- vector()
# good_ids <- vector()
# for(i in 1:k) {
#   tempn <- exp1[exp1[["subject_id"]] == sub_ids[i] &
#                 exp1[["is_predictable"]] == 0, "rt"]/1000
#   templ <- exp1[exp1[["subject_id"]] == sub_ids[i] &
#                 exp1[["is_predictable"]] == 1, "rt"]/1000
#   mi <- min(length(tempn), length(templ))
#   if (mi >= 50) {
#     tempn <- tempn[seq_len(mi)]
#     templ <- templ[seq_len(mi)]
#     # 3 subjects have is_pred=1 rt > 2 sec (subject_id: 48, 240, 131)
#     bad_idx <- which(tempn < 0 | tempn > 2 | templ < 0 | templ > 2)
#     n_bad_idx <- length(bad_idx)
#     if (n_bad_idx > 0) {
#       tempn <- tempn[-bad_idx]
#       templ <- templ[-bad_idx]
#       mi <- mi - n_bad_idx
#     }
#     nk <- c(nk, mi)
#     yn <- c(yn, tempn)
#     yl <- c(yl, templ)
#     good_ids <- c(good_ids, sub_ids[i])
#   }
# }
# k <- length(good_ids)



######################### Model Comparison #####################################
library("loo")
options(mc.cores = parallel::detectCores())
loo_fit <- function(fit, sub_ids = NULL,
                    models = c("no_learning", "step_learning",
                               "symmetric_logistic_learning",
                               "asymmetric_logistic_learning")) {
  if (is.null(sub_ids)) {
    sub_ids <- names(fit)
  }
  n_sub_ids <- length(sub_ids)
  n_models <- length(models)
  n_row <- n_sub_ids*n_models

  out <- data.frame(subject_id = rep(as.integer(sub_ids), each = n_models),
                    model_label = character(n_row),
                    elpd_diff = double(n_row),
                    se_diff = double(n_row),
                    elpd_loo = double(n_row),
                    se_elpd_loo = double(n_row),
                    p_loo = double(n_row),
                    se_p_loo = double(n_row),
                    looic = double(n_row),
                    se_looic = double(n_row))
  for (i in 1:n_sub_ids) {
    tmp <- fit[[as.character(sub_ids[i])]]
    # tmp <- fit[[sub_ids[i]]]
    loos <- list()
    for (m in 1:n_models) {
      loos[[m]] <- loo::loo(tmp[[models[m]]])
    }
    idx <- 1:n_models + (i-1)*n_models
    looc <- loo::loo_compare(loos)
    out[idx, 2] <- rownames(looc)
    out[idx, -c(1, 2)] <- looc
  }
  return(out)
}

loo_significant <- function(loo_fit) {
  loo_res <- list()
  for (id in sort(unique(loo_fit[["subject_id"]]))) {
    looi <- loo_fit[loo_fit[["subject_id"]] == id, ]
    loo_res[[as.character(id)]] <- looi[["model_label"]][abs(looi[["elpd_diff"]]) <= 5*looi[["se_diff"]]]
  }
  return(loo_res[lengths(loo_res) < 4])
}



######################### Plotting #############################################
library("ggplot2")

plot_loo_fit <- function(loo_fit, sub_ids = NULL, save_path = NULL) {
  if(is.null(sub_ids)) {
    sub_ids <- sort(unique(loo_fit[["subject_id"]]))
  }
  n_sub_ids <- length(sub_ids)
  models <- c("no learning", "step learning",
              "sym logi learning", "asym logi learning")
  n_models <- length(models)
  colors <- c("#b0e7ff", "#b0ffb7", "#ffc8b0", "#f6b0ff")

  loo_fit[["elpd_per"]] <- loo_fit[["elpd_loo"]] /
    rep(loo_fit[["elpd_loo"]][seq(1, by = n_models,
                                  length.out = n_sub_ids)],
        each = n_models)

  loo_fit[["n_pars"]] <- integer(n_sub_ids*n_models)
  loo_fit[["n_pars"]][loo_fit[["model_label"]] == "model1"] <- 6
  loo_fit[["n_pars"]][loo_fit[["model_label"]] == "model2"] <- 8
  loo_fit[["n_pars"]][loo_fit[["model_label"]] == "model3"] <- 9
  loo_fit[["n_pars"]][loo_fit[["model_label"]] == "model4"] <- 12

  print(ggplot(loo_fit) +
    geom_point(size = 3, shape = 16, alpha = 0.5,
               aes(x = subject_id, y = elpd_per, color = model_label)) +
    scale_color_manual(values = colors, labels = models, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = rep(5, n_models)))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Loo Estimates for elpd_loo",
         subtitle = paste("as a percent of the highest",
                          "elpd_loo per individual", sep = "\n"),
         x = "Subject ID Number", y = "Percent of Maximum elpd_loo") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      plot.title = element_text(size = 20,
        margin = margin(5, 0, 5, 0, "pt")),
      plot.subtitle = element_text(size = 16,
        margin = margin(5, 0, 10, 0, "pt")),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16,
        margin = margin(10, 0, 0, 0, "pt")),
      axis.title.y = element_text(size = 16,
        margin = margin(0, 10, 0, 0, "pt")),
      legend.position = c(1, 1),
      legend.justification = c(1, 0),
      legend.box = "horizontal",
      legend.direction = "vertical",
      legend.background = element_rect(fill = "transparent"),
      legend.text = element_text(size = 12))
  )
  if (!is.null(save_path)) {
    ggsave(file = paste0(save_path, "/elpd_percent.png"),
           width = 16, height = 9)
  }

  print(ggplot(loo_fit) +
    geom_pointrange(size = 1, shape = 16, alpha = 0.5,
                    aes(x = subject_id, y = elpd_diff,
                        ymin = elpd_diff - 5*se_diff,
                        ymax = elpd_diff + 5*se_diff,
                        color = model_label)) +
    scale_color_manual(values = colors, labels = models, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = rep(1, n_models)))) +
    coord_cartesian(ylim = c(min(loo_fit[["elpd_diff"]]), 0)) +
    labs(title = "Difference in ELPD Estimates",
         subtitle = "\u00B1 5 se",
         x = "Subject ID Number", y = "elpd_diff") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      plot.title = element_text(size = 20,
        margin = margin(5, 0, 5, 0, "pt")),
      plot.subtitle = element_text(size = 16,
        margin = margin(5, 0, 20, 0, "pt")),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16,
        margin = margin(10, 0, 0, 0, "pt")),
      axis.title.y = element_text(size = 16,
        margin = margin(0, 10, 0, 0, "pt")),
      legend.position = c(1, 1),
      legend.justification = c(1, 0),
      legend.box = "horizontal",
      legend.direction = "vertical",
      legend.background = element_rect(fill = "transparent"),
      legend.text = element_text(size = 12))
  )
  if (!is.null(save_path)) {
    ggsave(file = paste0(save_path, "/elpd_difference.png"),
           width = 16, height = 9)
  }

  print(ggplot(loo_fit) +
    geom_pointrange(size = 1, shape = 16, alpha = 0.5,
                    aes(x = subject_id, y = p_loo,
                        ymin = p_loo - se_p_loo,
                        ymax = p_loo + se_p_loo,
                        color = model_label)) +
    geom_line(size = 1.15, alpha = 1, linetype = "dashed",
              aes(x = subject_id, y = n_pars, color = model_label)) +
    scale_color_manual(values = colors, labels = models, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = rep(1, n_models)))) +
    coord_cartesian(ylim = c(0, max(loo_fit[["p_loo"]], 12))) +
    labs(title = "Effective Number of Parameters",
         subtitle = "\u00B1 1 se with actual number of parameters",
         x = "Subject ID Number", y = "Number of Parameters") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      plot.title = element_text(size = 20,
        margin = margin(5, 0, 5, 0, "pt")),
      plot.subtitle = element_text(size = 16,
        margin = margin(5, 0, 20, 0, "pt")),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16,
        margin = margin(10, 0, 0, 0, "pt")),
      axis.title.y = element_text(size = 16,
        margin = margin(0, 10, 0, 0, "pt")),
      legend.position = c(1, 1),
      legend.justification = c(1, 0),
      legend.box = "horizontal",
      legend.direction = "vertical",
      legend.background = element_rect(fill = "transparent"),
      legend.text = element_text(size = 12))
  )
  if (!is.null(save_path)) {
    ggsave(file = paste0(save_path, "/p_loo.png"),
           width = 16, height = 9)
  }
}

plot_post_pred <- function(fit, data, sub_ids = NULL,
                           pred_names = c("ynpred", "ylpred"),
                           models = c("no_learning", "step_learning",
                                      "symmetric_logistic_learning",
                                      "asymmetric_logistic_learning")) {
  if (is.null(sub_ids)) {
    sub_ids <- names(fit)
  }

  for (ind in 1:length(sub_ids)) {
    # read and clean data
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

    df <- data.frame(
      x = rep(1:nk, 2),
      y_data = c(yn, yl),
      data_label = c(rep(1, nk), rep(2, nk)),
      model_label = c(rep(3, nk), rep(4, nk)),
      y_nol = c(
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[1]]])[[pred_names[1]]], 2, mean),
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[1]]])[[pred_names[2]]], 2, mean)),
      y_stp = c(
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[2]]])[[pred_names[1]]], 2, mean),
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[2]]])[[pred_names[2]]], 2, mean)),
      y_sym = c(
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[3]]])[[pred_names[1]]], 2, mean),
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[3]]])[[pred_names[2]]], 2, mean)),
      y_asy = c(
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[4]]])[[pred_names[1]]], 2, mean),
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[4]]])[[pred_names[2]]], 2, mean)),
      sd_nol = c(
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[1]]])[[pred_names[1]]], 2, sd),
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[1]]])[[pred_names[2]]], 2, sd)),
      sd_stp = c(
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[2]]])[[pred_names[1]]], 2, sd),
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[2]]])[[pred_names[2]]], 2, sd)),
      sd_sym = c(
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[3]]])[[pred_names[1]]], 2, sd),
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[3]]])[[pred_names[2]]], 2, sd)),
      sd_asy = c(
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[4]]])[[pred_names[1]]], 2, sd),
        apply(rstan::extract(fit[[as.character(sub_ids[ind])]][[models[4]]])[[pred_names[2]]], 2, sd))
    )
    df[["dummy"]] <- rep(5, 2*nk)

    # factor levels key (since R makes it alphabetical otherwise)
      # 1: Non-Learned Data
      # 2: Learned Data
      # 3: Non-Learned Fit
      # 4: Learned Fit

    print(ggplot(df) +
      geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
        aes(x = x, ymin = y_nol-sd_nol, ymax = y_nol+sd_nol,
            fill = factor(model_label, levels = c(3, 4)))) +
      geom_line(size = 1.25, alpha = 0.6,
        aes(x = x, y = y_nol,
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

    print(ggplot(df) +
      geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
        aes(x = x, ymin = y_stp-sd_stp, ymax = y_stp+sd_stp,
            fill = factor(model_label, levels = c(3, 4)))) +
      geom_line(size = 1.25, alpha = 0.6,
        aes(x = x, y = y_stp,
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
      labs(title = "Step Learning Model Fit",
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

    print(ggplot(df) +
      geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
        aes(x = x, ymin = y_sym-sd_sym, ymax = y_sym+sd_sym,
            fill = factor(model_label, levels = c(3, 4)))) +
      geom_line(size = 1.25, alpha = 0.6,
        aes(x = x, y = y_sym,
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
      labs(title = "Symmetric Logistic Learning Model Fit",
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
      labs(title = "Asymmetric Logistic Learning Model Fit",
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
  }
}



plot_model_est <- function(fit, data, sub_ids = NULL,
                           models = c("no_learning", "step_learning",
                                      "symmetric_logistic_learning",
                                      "asymmetric_logistic_learning")) {
  if (is.null(sub_ids)) {
    sub_ids <- names(fit)
  }

  # define some functions
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

  for (ind in 1:length(sub_ids)) {
    # read and clean data
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

    # factor levels key (since R makes it alphabetical otherwise)
      # 1: Non-Learned Data
      # 2: Learned Data
      # 3: Non-Learned Fit
      # 4: Learned Fit

    print(ggplot(df) +
      geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
        aes(x = x, ymin = y_nol-sd_nol, ymax = y_nol+sd_nol,
            fill = factor(model_label, levels = c(3, 4)))) +
      geom_line(size = 1.25, alpha = 0.6,
        aes(x = x, y = y_nol,
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

    print(ggplot(df) +
      geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
        aes(x = x, ymin = y_stp-sd_stp, ymax = y_stp+sd_stp,
            fill = factor(model_label, levels = c(3, 4)))) +
      geom_line(size = 1.25, alpha = 0.6,
        aes(x = x, y = y_stp,
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
      labs(title = "Step Learning Model Fit",
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

    print(ggplot(df) +
      geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
        aes(x = x, ymin = y_sym-sd_sym, ymax = y_sym+sd_sym,
            fill = factor(model_label, levels = c(3, 4)))) +
      geom_line(size = 1.25, alpha = 0.6,
        aes(x = x, y = y_sym,
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
      labs(title = "Symmetric Logistic Learning Model Fit",
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
      labs(title = "Asymmetric Logistic Learning Model Fit",
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
  }
}
