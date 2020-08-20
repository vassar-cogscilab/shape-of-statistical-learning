######################### Data Fitting #########################################
data_fit <- function(data, sub_ids = NULL) {
  if(is.null(sub_ids)) {
    sub_ids <- sort(unique(data[["subject_id"]]))
  }
  k <- length(sub_ids)
  out <- list()

  for (i in 1:k) {
    # read data
    yn <- data[data[["subject_id"]] == sub_ids[i] &
               data[["is_predictable"]] == 0, "rt"]/1000
    yl <- data[data[["subject_id"]] == sub_ids[i] &
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
    out[[i]] <- list(
      "no_learning" =
        sampling(no_learning_model,
        data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
        refresh = FALSE, chains = 1, iter = 500,
        control = list(adapt_delta = 0.9, max_treedepth = 10),
        init = list( list(
          V = 1, E = 0.25, A = 0.25, S = 0 )) ),
      "step_learning" =
        sampling(step_learning_model,
        data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
        refresh = FALSE, chains = 1, iter = 500,
        control = list(adapt_delta = 0.9, max_treedepth = 10),
        init = list( list(
          V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, H_raw = 0.5 )) ),
      "symmetric_logistic_learning" =
        sampling(symmetric_logistic_learning_model,
        data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
        refresh = FALSE, chains = 1, iter = 500,
        control = list(adapt_delta = 0.9, max_treedepth = 10),
        init = list( list(
          V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, L = 6, H_raw = 0.5 )) ),
      "asymmetric_logistic_learning" =
        sampling(asymmetric_logistic_learning_model,
        data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
        refresh = FALSE, chains = 1, iter = 500,
        control = list(adapt_delta = 0.9, max_treedepth = 10),
        init = list( list(
          V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, L = 6, H_raw = 0.5,
          NU = 1, C = 1, Q = 1 )) )
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


######################### Plotting #############################################
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
      labs(title = plot_title, subtitle = paste0("Subject ", i),
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



plot_model_est <- function(fit, nk, yn, yl, parnames,
                           plot_title = "Experimental Data Fit") {
  parnames <- c(parnames, "sigma2_n", "sigma2_l")
  n_pars <- length(parnames)
  pars <- list()
  for (i in 1:n_pars) {
    pars[[parnames[i]]] <- apply(rstan::extract(fit)[[parnames[i]]], 2, mean)
  }

  for (ind in 1:k) {
    nki <- nk[ind]
    st <- ifelse(ind > 1, sum(nk[1:(ind-1)]), 0)
    yni <- yn[(st+1):(st+nki)]
    yli <- yl[(st+1):(st+nki)]
    trials <- seq_len(nki)
    sig2_n <- pars[["sigma2_n"]][ind]
    sig2_l <- pars[["sigma2_l"]][ind]
    mu_n <- pars[["V"]][ind] + pars[["E"]][ind] * exp(-pars[["A"]][ind]*trials)
    if (n_pars == 6) { # no learning
      mu_l <- pars[["V"]][ind] + pars[["S"]][ind] +
              pars[["E"]][ind] * exp(-pars[["A"]][ind] * trials)
    } else if (n_pars == 8) { # step learning
      mu_l <- c(pars[["V"]][ind] + pars[["S"]][ind] + pars[["E"]][ind] *
                  exp(-pars[["A"]][ind] * 1:(pars[["H"]][ind]-1)),
                ( pars[["V"]][ind] + pars[["S"]][ind] + pars[["E"]][ind] *
                  exp(-pars[["A"]][ind] * pars[["H"]][ind]:nki) ) *
                  (1 - pars[["D"]][ind]))
    } else if (n_pars == 9) { # symmetric logistic learning
      mu_l <- ( pars[["V"]][ind] + pars[["S"]][ind] +
                pars[["E"]][ind] * exp(-pars[["A"]][ind] * trials) ) *
              ( 1 - pars[["D"]][ind] /
              ( 1 + exp(-pars[["L"]][ind] * (trials - pars[["H"]][ind])) ) )
    } else if (n_pars == 12) { # asymmetric logistic learning
      mu_l <- ( pars[["V"]][ind] + pars[["S"]][ind] +
                pars[["E"]][ind] * exp(-pars[["A"]][ind] * trials) ) *
              ( 1 - pars[["D"]][ind] /
              ( pars[["C"]][ind] + pars[["Q"]] * exp(-pars[["L"]][ind] *
                  (trials - pars[["H"]][ind])) )^(1/pars[["NU"]]) )
    } else {
      stop("unkown number of parameters")
    }

    df <- data.frame(
      x = rep(trials, 2),
      y_data = c(yn[(st+1):(st+nki)], yl[(st+1):(st+nki)]),
      data_label = c(rep(1, nki), rep(2, nki)),
      model_label = c(rep(3, nki), rep(4, nki)),
      y_means = c(mu_n, mu_l),
      y_sds = c((sig2_n - 1) * mu_n*mu_n * sig2_n*sig2_n,
                (sig2_l - 1) * mu_l*mu_l * sig2_l*sig2_l)
    )
    df$y_hi <- df$y_means + df$y_sds
    df$y_lo <- df$y_means - df$y_sds
    df$dummy <- rep(5, 2*nki)

    # factor levels key (since R makes it alphabetical otherwise)
      # 1: Non-Learned Data
      # 2: Learned Data
      # 3: Non-Learned Fit
      # 4: Learned Fit

    print(ggplot(df) +
      geom_ribbon(linetype = "blank", alpha = 0.25, outline.type = NULL,
        aes(x = x, ymin = y_lo, ymax = y_hi,
            fill = factor(model_label, levels = c(3, 4)))) +
      geom_line(size = 1.25, alpha = 0.6,
        aes(x = x, y = y_means,
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
      labs(title = plot_title, subtitle = paste0("Subject ", ind),
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
