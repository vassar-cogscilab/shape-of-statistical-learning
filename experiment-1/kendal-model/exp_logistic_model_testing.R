library("rstan")
library("loo")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

exp_log_model_stanc <- stanc(file = "experiment-1/kendal-model/exp_logistic_model.stan",
                             model_name = "exponential_logistic_model")
exp_log_model <- stan_model(stanc_ret = exp_log_model_stanc)

save(exp_log_model, file = "experiment-1/kendal-model/model-comparison/exp_logistic_model.Rds",
     compress = "xz", compression_level = 9)
load("experiment-1/kendal-model/model-comparison/exp_logistic_model.Rds")



exp_log <- function(trial, V, E, A, D, L, H) {
  return( (V + E*exp(-A*trial)) * (1 - D/(1 + exp(-L*(trial-H)))) )
}
sim_exp_log <- function(trial, V, E, A, D, L, H, var) {
  sim <- exp_log(trial = trial, E = E, A = A, V = V, D = D, L = L, H = H) +
         rnorm(length(trial), 0, var)
  sim[sim < 0] <- 0
  return(sim)
}




########## Plotting Functions
plot_sim_fit_pars <- function(fits, K, V, E, A, P, D, L, H) {
  pars <- data.frame(id = seq_len(K),
                     V = V,
                     E = E,
                     A = A,
                     P = P,
                     D = D,
                     L = L,
                     H = H)
  idx <- c(1:6, 9)
  for (i in seq_len(length(idx))) {
    df <- as.data.frame(fits[seq_len(K) + (idx[i]-1)*K, c(1, 4, 8)])
    colnames(df) <- c("mean", "qlow", "qupp")
    df$id <- seq_len(K)
    df$dumm = c(1, rep(2, K-1))
    par_name <- substring(rownames(df)[1], 1, 1)

    print(ggplot(df) +
      geom_pointrange(alpha = 0.6, shape = 16, color = "#5fc6fa",
                      aes(x = id, y = mean, ymin = qlow, ymax = qupp,
                          size = factor(dumm, levels = c(1, 2)))) +
      geom_point(data = pars, alpha = 0.8, color = "#04547c",
                 shape = 4, size = 6, stroke = 2,
                 aes_string(x = "id", y = par_name)) +
      scale_size_manual(values = c(1.51, 1.5),
                        labels = c("True Value", "Est. Mean"),
                        name = NULL) +
      guides(size = guide_legend(override.aes = list(size = c(1.5, 1.5),
                                                     stroke = c(2, 1),
                                                     shape = c(4, 16),
                                                     lty = c(0, 0),
                                                     color = c("#04547c",
                                                               "#5fc6fa")))) +

      scale_x_discrete(limits = as.character(seq_len(K))) +
      labs(title = paste0("Simulated Parameter Estimates for ", par_name),
           x = "ID of Individual", y = "Estimated Parameter Value",
           subtitle = "Bars represent the 2.5% and 97.5% CI") +
      theme_bw() +
      theme(panel.border = element_blank(),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 16,
                                         margin = margin(0, 5, 15, 5, "pt")),
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
            legend.background = element_rect(fill = "transparent"))
    )
  }
}

plot_sim_fits <- function(fits, K, NTI,
                          V = NA, E = NA, A = NA, D = NA,
                          L = NA, H = NA, Y0 = NA, Y1 = NA,
                          savepath = NULL) {
  st <- 0
  for (i in 1:K) {
    n <- NTI[i]
    Vi <- fits[i+0*K]
    Ei <- fits[i+1*K]
    Ai <- fits[i+2*K]
    Pi <- fits[i+3*K]
    Di <- fits[i+9*K]
    Li <- fits[i+5*K]
    Hi <- fits[i+8*K]
    sigma_2i <- fits[i+7*K]
    sn <- seq_len(n)

    df <- data.frame(x = rep(sn, 2),
                     data = c(Y0[(st+1):(st+n)], Y1[(st+1):(st+n)]),
                     truth = c(exp_log(sn, V = V[i], E = E[i], A = A[i],
                                       D = 0, L = L[i], H = H[i]),
                               exp_log(sn, V = V[i], E = E[i], A = A[i],
                                       D = D[i], L = L[i], H = H[i])),
                     fit = c(exp_log(sn, V = Vi, E = Ei, A = Ai,
                                     D = 0, L = Li, H = Hi),
                             exp_log(sn, V = Vi, E = Ei, A = Ai,
                                     D = Pi * Di, L = Li, H = Hi)),
                     data_label = c(rep(1,n), rep(2,n)),
                     fit_label = c(rep(3,n), rep(4,n)))
    df$fit_min <- df$fit - df$fit*sqrt(sigma_2i*(sigma_2i - 1))
    df$fit_max <- df$fit + df$fit*sqrt(sigma_2i*(sigma_2i - 1))
    df$dumm <- c(rep(1, n), rep(2, n))

    # factor levels key (since R makes it alphabetical)
      # 1: Non-Learned Data/Truth
      # 2: Learned Data/Truth
      # 3: Non-Learned Fit
      # 4: Learned Fit

    print(ggplot(df) +
      geom_ribbon(linetype="blank", alpha=0.25,
                  aes(x = x, ymin = fit_min, ymax = fit_max,
                      fill = factor(fit_label, levels = c(3, 4)))) +
      geom_point(alpha = 0.75, shape = 16,
                 aes(x = x, y = data,
                     color = factor(data_label, levels = c(1, 2)),
                     size = factor(dumm, levels = c(1, 2)))) +
      geom_line(size = 1.25, alpha = 0.8,
                aes(x = x, y = truth,
                    color = factor(data_label, levels = c(1, 2)))) +
      geom_line(size = 1.25, alpha = 0.6,
                aes(x = x, y = fit,
                    color = factor(fit_label, levels = c(3, 4)))) +
      scale_color_manual(values = c("#999999", "#000000",
                                    "#98c1ff", "#ff4f4f"),
                         labels = c("Non-Learned Truth", "Learned Truth",
                                    "Non-Learned Fit", "Learned Fit"),
                         name = NULL) +
      scale_fill_manual(values = c("#98c1ff", "#ff4f4f"),
                        guide = FALSE) +
      scale_size_manual(values = c(1.8, 1.81),
                        labels = c("Non-Learned Data", "Learned Data"),
                        name = NULL) +
      guides(color = guide_legend(order = 2,
                                  override.aes = list(shape = NA,
                                                      size = 1.25,
                                                      lty = rep(1, 4))),
             size = guide_legend(order = 1,
                                 override.aes = list(size = 2,
                                                     shape = c(16, 16),
                                                     color = c("#999999",
                                                               "#000000")))) +
      labs(title = "Simulated Data Fit",
           x = "Trial Number", y = "Response Time (sec)",
           subtitle = paste0("\u00B1 1 standard deviation\n",
                             "Subject ", i, "\n",
                             "P = ", round(Pi, 3))) +
      theme_bw() +
      theme(panel.border = element_blank(),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 16,
                                         margin = margin(0, 5, 15, 5, "pt")),
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
            legend.background = element_rect(fill = "transparent"))
    )

    if (!is.null(savepath)) {
      ggsave(paste0(savepath, "sim_fit_", i, ".png"))
    }

    st = st + n
  }
}

plot_real_fit_pars <- function(fits, K) {
  idx <- c(1:6, 9)
  for (i in seq_len(length(idx))) {
    df <- as.data.frame(fits[seq_len(K) + (idx[i]-1)*K, c(1, 4, 8)])
    colnames(df) <- c("mean", "qlow", "qupp")
    df$id <- seq_len(K)
    # df$dumm = c(1, rep(2, K-1))
    par_name <- substring(rownames(df)[1], 1, 1)

    print(ggplot(df) +
      geom_pointrange(alpha = 0.6, shape = 16, size = 1.5, color = "#5fc6fa",
                      aes(x = id, y = mean, ymin = qlow, ymax = qupp)) +

      scale_x_discrete(limits = as.character(seq_len(K))) +
      labs(title = paste0("Real Parameter Estimates for ", par_name),
           x = "ID of Individual", y = "Estimated Parameter Value",
           subtitle = "Bars represent the 2.5% and 97.5% CI") +
      theme_bw() +
      theme(panel.border = element_blank(),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 16,
                                         margin = margin(0, 5, 15, 5, "pt")),
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
            legend.background = element_rect(fill = "transparent"))
    )
  }
}

plot_real_fits <- function(fits, K, NTI, Y0 = NA, Y1 = NA, savepath = NULL) {
  st <- 0
  for (i in 1:K) {
    n <- NTI[i]
    Vi <- fits[i+0*K]
    Ei <- fits[i+1*K]
    Ai <- fits[i+2*K]
    Pi <- fits[i+3*K]
    Di <- fits[i+9*K]
    Li <- fits[i+5*K]
    Hi <- fits[i+8*K]
    sigma_2i <- fits[i+7*K]
    sn <- seq_len(n)

    df <- data.frame(x = rep(sn, 2),
                     data = c(Y0[(st+1):(st+n)], Y1[(st+1):(st+n)]),
                     fit = c(exp_log(sn, V = Vi, E = Ei, A = Ai,
                                     D = 0, L = Li, H = Hi),
                             exp_log(sn, V = Vi, E = Ei, A = Ai,
                                     D = Pi * Di, L = Li, H = Hi)),
                     data_label = c(rep(1,n), rep(2,n)),
                     fit_label = c(rep(3,n), rep(4,n)))
    df$fit_min <- df$fit - df$fit*sqrt(sigma_2i*(sigma_2i - 1))
    df$fit_max <- df$fit + df$fit*sqrt(sigma_2i*(sigma_2i - 1))

    # factor levels key (since R makes it alphabetical)
      # 1: Non-Learned Data/Truth
      # 2: Learned Data/Truth
      # 3: Non-Learned Fit
      # 4: Learned Fit

    print(ggplot(df) +
      geom_ribbon(linetype="blank", alpha=0.25,
                  aes(x = x, ymin = fit_min, ymax = fit_max,
                      fill = factor(fit_label, levels = c(3, 4)))) +
      geom_point(alpha = 0.75, shape = 16, size = 1.8,
                 aes(x = x, y = data,
                     color = factor(data_label, levels = c(1, 2)))) +
      geom_line(size = 1.25, alpha = 0.6,
                aes(x = x, y = fit,
                    color = factor(fit_label, levels = c(3, 4)))) +
      scale_color_manual(values = c("#999999", "#000000",
                                    "#98c1ff", "#ff4f4f"),
                         labels = c("Non-Learned Data", "Learned Data",
                                    "Non-Learned Fit", "Learned Fit"),
                         name = NULL) +
      scale_fill_manual(values = c("#98c1ff", "#ff4f4f"),
                        guide = FALSE) +
      guides(color = guide_legend(override.aes = list(size = c(2, 2, 1.25, 1.25),
                                                      shape = c(16, 16, NA, NA),
                                                      lty = c(NA, NA, 1, 1)))) +
      labs(title = "Experimental Data Fit",
           x = "Trial Number", y = "Response Time (sec)",
           subtitle = paste0("\u00B1 1 standard deviation\n",
                             "Subject ", i, "\n",
                             "P = ", round(Pi, 3))) +
      theme_bw() +
      theme(panel.border = element_blank(),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 16,
                                         margin = margin(0, 5, 15, 5, "pt")),
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
            legend.background = element_rect(fill = "transparent"))
    )

    if (!is.null(savepath)) {
      ggsave(paste0(savepath, "real_fit_", i, ".png"))
    }

    st = st + n
  }
}



########## Simulated Data
sK <- 1 # number of individuals
sN <- 100 # number of trials for each individual (same for all)
sNTI <- rep(sN, sK) # so Stan can handle multiple individuals
sNK <- sN*sK
sY0 <- rep(NA, sNK) # vector of non-learned response times
sY1 <- rep(NA, sNK) # vector of learned response times
V <- c(1.1, 1.5, .75)
E <- c(.25, .25, .25)
A <- c(.25, .25, .25)
P <- c(1, 0, 1)
D <- c(.1, 0, .1)
L <- c(.4, .4, .4)
H <- c(50, 60, 30)
st <- 0
for (i in 1:sK) {
  x <- seq_len(sNTI[i])
  sY0[(st+1):(st+sNTI[i])] <- sim_exp_log(x, V=V[i], E=E[i], A=A[i], D=0, L=L[i], H=H[i], var=0.15) # generate non-learned response times
  sY1[(st+1):(st+sNTI[i])] <- sim_exp_log(x, V=V[i], E=E[i], A=A[i], D=D[i], L=L[i], H=H[i], var=0.15) # generate learned response times
  st = st + sNTI[i]
}
sim_fit_orig
sim_fit_orig <- summary(sampling(
  exp_log_model,
  data = list('K' = sK, 'NTI' = as.array(sNTI), 'NK' = sNK, 'Y0' = sY0, 'Y1' = sY1),
  refresh = FALSE, chains = 1, iter = 1000, seed = 2,
  control = list(adapt_delta = 0.99, max_treedepth = 10)
))$summary


plot_sim_fit_pars(sim_fit, sK, V, E, A, P, D, L, H)

plot_sim_fits(sim_fit, sK, sNTI, V, E, A, D, L, H, sY0, sY1)

########## Real-World Data
load(file = "experiment-1/kendal-model/exp1.Rds")
sub_ids <- unique(exp1$subject_id)[c(78, 7, 10)]
sub_ids <- sub_ids[3]
rK <- length(sub_ids)
gs_idx <- vector()
rNTI <- rep(0, rK)
rY0 <- vector()
rY1 <- vector()
for (i in seq_len(rK)) {
  temp <- exp1[exp1$subject_id == sub_ids[i], ]
  temp0 <- temp[temp$is_predictable == 0, ]$rt/1000
  temp1 <- temp[temp$is_predictable == 1, ]$rt/1000
  mm <- min(length(temp0), length(temp1))
  if (mm >= 50) {
    rNTI[i] <- mm
    rY0 <- c(rY0, temp0[seq_len(mm)])
    rY1 <- c(rY1, temp1[seq_len(mm)])
    gs_idx <- c(gs_idx, i)
  }
}

good_sub_ids <- sub_ids[gs_idx]
rNTI <- rNTI[rNTI > 0]
rNK <- sum(rNTI)
rK <- length(good_sub_ids)

real_fit <- summary(sampling(
  exp_log_model,
  data = list('K' = rK, 'NTI' = as.array(rNTI), 'NK' = rNK, 'Y0' = rY0, 'Y1' = rY1),
  refresh = FALSE, chains = 1, iter = 500, seed = 2,
  control = list(adapt_delta = 0.9, max_treedepth = 10)
))$summary


plot_real_fit_pars(real_fit, rK)

plot_real_fits(real_fit, rK, rNTI, rY0, rY1)
