// early exponential adaptation with logistic learning curve
  // learning curve is a symmetric logistic function
  // response times are lognormally distributed
  // fit is_pred=0 and is_pred=1 data simultaneously

data { // these values are input to the sampling function
  int<lower=1> k; // number of individuals
  int<lower=1> nk[k]; // array of number of trials performed by each individual
  int<lower=k> total_length; // total length of all RT data
  vector<lower=0, upper=2>[total_length] yn; // RT data from the individual for each trial of is_predictable=0
  vector<lower=0, upper=2>[total_length] yl; // RT data from the individual for each trial of is_predictable=1
}

transformed data { // manipulations of the data to use in parameter bounds
  vector<lower=1>[k] upperH = to_vector(nk); // (rough) upper bound for H
}

parameters { // define parameters (and their bounds) used in the model
  vector<lower=0, upper=2>[k] V;  // vertical shift of all repsonse times
  vector<lower=0>[k] E;  // overall scale of exponential
  vector<lower=0>[k] A;  // adaptation rate
  vector[k] S; // shift in median learned response time relative to non-learned

  vector<lower=0, upper=1>[k] D;  // scale of learning "drop"
  vector<lower=0>[k] L;  // learning rate
  vector<lower=0, upper=1>[k] H_raw;  // horizontal shift in onset of learning (raw)

  vector<lower=0>[k] sigma2_n; // main component of variance in non-learned response times
  vector<lower=0>[k] sigma2_l; // main component of variance in learned response times
}

transformed parameters { // manipulations of the parameters (really, just their bounds)
  vector[k] H = H_raw .* upperH;
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

  // Loop through each individual
  for (ind in 1:k) {
    int nki = nk[ind];
    vector[nki] mu_n;
    vector[nki] mu_l;
    int st = (ind > 1) ? sum(nk[1:(ind-1)]) : 0;
    // Likelihood distribution for model
    for (trial in 1:nk[ind]) {
      // We log mu so that the paramters directly determine the median so it's
      // easier to interpret; thus we use mu -> log(mu) so that
      // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
      mu_n[trial] = log( V[ind] + E[ind]*exp(-A[ind]*trial) );
      mu_l[trial] = log( V[ind] + E[ind]*exp(-A[ind]*trial) + S[ind] ) +
                    log1m( D[ind]/( 1 + exp(-L[ind]*(trial-H[ind]))) );
    }

    yn[(st+1):(st+nki)] ~ lognormal( mu_n, 0.5*log(sigma2_n[ind]) );
    yl[(st+1):(st+nki)] ~ lognormal( mu_l, 0.5*log(sigma2_l[ind]) );
  }
}

generated quantities {
  real ynpred[total_length];
  real ylpred[total_length];
  real log_lik[total_length];

  // Loop through each individual
  for (ind in 1:k) {
    int nki = nk[ind];
    vector[nki] mu_n;
    vector[nki] mu_l;
    int st = (ind > 1) ? sum(nk[1:(ind-1)]) : 0;
    for (trial in 1:nk[ind]) {
      mu_n[trial] = log( V[ind] + E[ind]*exp(-A[ind]*trial) );
      mu_l[trial] = log( V[ind] + E[ind]*exp(-A[ind]*trial) + S[ind]) +
                    log1m( D[ind]/( 1 + exp(-L[ind]*(trial-H[ind]))) );
      log_lik[st+trial] =
        lognormal_lpdf(yl[st+trial] | mu_l[trial], 0.5*log(sigma2_l[ind]) );
    }

    ynpred[(st+1):(st+nki)] = lognormal_rng(mu_n, 0.5*log(sigma2_n[ind]) );
    ylpred[(st+1):(st+nki)] = lognormal_rng(mu_l, 0.5*log(sigma2_l[ind]) );
  }
}
