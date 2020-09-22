// early exponential adaptation with step learning curve
  // learning curve is a simple step function
  // response times are lognormally distributed
  // fit is_pred=0 and is_pred=1 data simultaneously

functions { // user-defined functions
  int argmax(vector x) {
    int idx_max = 1;
    int n = num_elements(x);
    for (i in 1:n) {
      if (x[i] > x[idx_max]) idx_max = i;
    }
    return idx_max;
  }
}

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

  real<lower=0, upper=1> D_raw;  // scale of learning "drop" (raw)

  real<lower=0> sigma2_n; // main component of variance in non-learned response times
  real<lower=0> sigma2_l; // main component of variance in learned response times
}

transformed parameters { // manipulations of the parameters (really, just their bounds)
  real D = D_raw * V; // scale of learning "drop" (properly scaled)
  vector[nk] H; // horizontal shift in onset of learning
  vector[nk+1] lp_e;
  vector[nk+1] lp_l;
  lp_e[1] = 0;
  lp_l[1] = 0;
  for (i in 1:nk) {
    lp_e[i+1] = lp_e[i] + lognormal_lpdf(yl[i] | log( V + E*exp(-A*i) + S ), 0.5*log(sigma2_l));
    lp_l[i+1] = lp_l[i] + lognormal_lpdf(yl[i] | log( V + E*exp(-A*i) + S ) + log1m(D), 0.5*log(sigma2_l));
    H[i] = -lbeta(10, 8) + 9*log((i+0.0)/nk) + 7*log(1 - (i+0.0)/nk);
  }
  H += lp_l[nk+1] + head(lp_e, nk) - head(lp_l, nk);
}

model {
  // Prior distributions on parameters
  target += gamma_lpdf(V | 5, 4);
  target += gamma_lpdf(E | 2.5, 10);
  target += gamma_lpdf(A | 2.5, 10);
  target += normal_lpdf(S | 0, 0.1);
  
  target += beta_lpdf(D | 15, 1);

  target += gamma_lpdf(sigma2_n | 3, 2);
  target += gamma_lpdf(sigma2_l | 3, 2);

  {
    vector[nk] mu_n;
    // Likelihood distribution for model
    for (trial in 1:nk) {
      // We log mu so that the paramters directly determine the median so it's
      // easier to interpret; thus we use mu -> log(mu) so that
      // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
      mu_n[trial] = log( V + E*exp(-A*trial) );
    }
    target += lognormal_lpdf(yn | mu_n, 0.5*log(sigma2_n));
    target += log_sum_exp(H);
  }
}

generated quantities {
  real ynpred[nk];
  real ylpred[nk];
  real log_lik[nk];

  {
    vector[nk] mu_n;
    vector[nk] mu_l;
    real h = argmax(H);
    for (trial in 1:nk) {
      mu_n[trial] = log( V + E*exp(-A*trial) );
      mu_l[trial] = log( V + E*exp(-A*trial) + S ) + step(trial - h)*log1m(D);
      log_lik[trial] =
        lognormal_lpdf(yl[trial] | mu_l[trial], 0.5*log(sigma2_l));
    }

    ynpred = lognormal_rng(mu_n, 0.5*log(sigma2_n) );
    ylpred = lognormal_rng(mu_l, 0.5*log(sigma2_l) );
  }
}
