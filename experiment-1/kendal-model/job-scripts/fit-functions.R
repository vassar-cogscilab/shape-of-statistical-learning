# add cat statements

##### Data Fitting
data_fit <- function(data, sub_ids = NULL,
                     n_chains = 1, n_iterations = 500,
                     adapt_delta = 0.9, max_treedepth = 10,
                     models = c("no_learning", "step_learning",
                                "symmetric_logistic_learning",
                                "asymmetric_logistic_learning")) {
  if(is.null(sub_ids)) {
    sub_ids <- sort(unique(data[["subject_id"]]))
  }
  k <- length(sub_ids)
  out <- list()

  for (i in sub_ids) {
    cat("fitting subject", i, "out of", number of ids)
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
      "no_learning" = if ("no_learning" %in% models) {
        rstan::sampling(no_learning_model,
          data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
          refresh = FALSE, chains = n_chains, iter = n_iterations,
          control = list(adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth),
          init = rep( list( list(
            V = 1, E = 0.25, A = 0.25, S = 0) ),
            n_chains) )
      },
      "step_learning" = if ("step_learning" %in% models) {
        rstan::sampling(step_learning_model,
          data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
          refresh = FALSE, chains = n_chains, iter = n_iterations,
          control = list(adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth),
          init = rep( list( list(
            V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, H_raw = 0.5) ),
            n_chains) )
      },
      "symmetric_logistic_learning" = if ("symmetric_logistic_learning" %in% models) {
        rstan::sampling(symmetric_logistic_learning_model,
          data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
          refresh = FALSE, chains = n_chains, iter = n_iterations,
          control = list(adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth),
          init = rep( list( list(
            V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, L = 0.4, H_raw = 0.5) ),
            n_chains) )
      },
      "asymmetric_logistic_learning" = if ("asymmetric_logistic_learning" %in% models) {
        rstan::sampling(asymmetric_logistic_learning_model,
          data = list('nk' = nk, 'yn' = yn, 'yl' = yl),
          refresh = FALSE, chains = n_chains, iter = n_iterations,
          control = list(adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth),
          init = rep( list( list(
            V = 1, E = 0.25, A = 0.25, S = 0, D = 0.5, L = 0.4, H_raw = 0.5,
            NU = 1, C = 1, Q = 1) ), n_chains) )
      }
    )
  }
  return(out)
}
