// early exponential adaptation with logistic learning curve
  // learning curve is a symmetric logistic function
  // response times are lognormally distributed
  // fit is_pred=0 and is_pred=1 data simultaneously

data { // these values are input to the sampling function
  int<lower=1> nk; // number of trials performed by the individual
  vector<lower=0, upper=2>[nk] yn; // RT data from the individual for each trial of is_predictable=0
  vector<lower=0, upper=2>[nk] yl; // RT data from the individual for each trial of is_predictable=1
}

parameters { // define parameters (and their bounds) used in the model
  real<lower=0, upper=2> V;  // vertical shift of all repsonse times
  real<lower=0> E;  // overall scale of exponential
  real<lower=0> A;  // adaptation rate
  real S; // shift in median learned response time relative to non-learned

  real<lower=0, upper=1> D;  // scale of learning "drop"
  real<lower=0> L;  // learning rate
  real<lower=0, upper=1> H_raw;  // horizontal shift in onset of learning (raw)

  real<lower=0> sigma2_n; // main component of variance in non-learned response times
  real<lower=0> sigma2_l; // main component of variance in learned response times
}

transformed parameters { // manipulations of the parameters (really, just their bounds)
  real H = H_raw * nk;
}

model {
  // Prior distributions on parameters
  V ~ gamma(2.5, 2.5);
  E ~ gamma(2.5, 10);
  A ~ gamma(2.5, 10);
  S ~ normal(0, 0.5);
  D ~ beta(2.5, 2.5);
  L ~ gamma(1.5, 0.25);
  H_raw ~ beta(1.5, 1.5);

  sigma2_n ~ gamma(3, 2);
  sigma2_l ~ gamma(3, 2);

  for (dummy in 1:1) {
    vector[nk] mu_n;
    vector[nk] mu_l;
    // Likelihood distribution for model
    for (trial in 1:nk) {
      // We log mu so that the paramters directly determine the median so it's
      // easier to interpret; thus we use mu -> log(mu) so that
      // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
      mu_n[trial] = log( V + E*exp(-A*trial) );
      mu_l[trial] = log( V + E*exp(-A*trial) + S ) +
                    log1m( D/( 1 + exp(-L*(trial-H))) );
    }

    yn ~ lognormal( mu_n, 0.5*log(sigma2_n) );
    yl ~ lognormal( mu_l, 0.5*log(sigma2_l) );
  }
}

generated quantities {
  real ynpred[nk];
  real ylpred[nk];
  real log_lik[nk];

  for (dummy in 1:1) {
    vector[nk] mu_n;
    vector[nk] mu_l;
    for (trial in 1:nk) {
      mu_n[trial] = log( V + E*exp(-A*trial) );
      mu_l[trial] = log( V + E*exp(-A*trial) + S ) +
                    log1m( D/( 1 + exp(-L*(trial-H))) );
      log_lik[trial] =
        lognormal_lpdf(yl[trial] | mu_l[trial], 0.5*log(sigma2_l));
    }

    ynpred = lognormal_rng(mu_n, 0.5*log(sigma2_n) );
    ylpred = lognormal_rng(mu_l, 0.5*log(sigma2_l) );
  }
}
