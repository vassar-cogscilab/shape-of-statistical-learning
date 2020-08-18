// early exponential adaptation with logistic learning curve
  // learning curve is a symmetric logistic function
  // response times are lognormally distributed
  // fit is_pred=0 and is_pred=1 data simultaneously

functions { // define functions used in model
  real mylognormal_lpdf(real z, real mu, real sigma_2) {
    real lprob = -log(z*sqrt(2*pi()*sigma_2)) - (log(z)-mu)^2 / (2*sigma_2);
    return lprob;
  }
}

data { // these values are input to the sampling function
  int<lower=1> k; // number of individuals
  int<lower=1> nk[k]; // array of number of trials performed by each individual
  int<lower=k> total_length; // total length of all RT data
  vector<lower=0>[total_length] yn; // RT data from the individual for each trial of is_predictable=0
  vector<lower=0>[total_length] yl; // RT data from the individual for each trial of is_predictable=1
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

  vector<lower=0>[k] sigma2_nV;
  vector<lower=0>[k] sigma2_nE;
  vector<lower=0>[k] sigma2_nA;
  vector<lower=0>[k] sigma2_lV;
  vector<lower=0>[k] sigma2_lE;
  vector<lower=0>[k] sigma2_lA;
}

transformed parameters { // manipulations of the parameters (really, just their bounds)
  vector[k] H = H_raw .* upperH;
}

model {
  real mu_n;
  real mu_l;
  real sigma_n;
  real sigma_l;
  int st = 0;

  // Prior distributions on parameters
  V ~ gamma(2.5, 2.5);
  E ~ gamma(2.5, 10);
  A ~ gamma(2.5, 10);
  S ~ normal(0, 0.5);
  D ~ beta(2.5, 2.5);
  L ~ gamma(1.5, 0.25);
  H ~ cauchy(upperH/2, 25);

  sigma2_nV ~ gamma(10, 10);
  sigma2_nE ~ gamma(2, 10);
  sigma2_nA ~ gamma(1, 10);
  sigma2_lV ~ gamma(10, 10);
  sigma2_lE ~ gamma(2, 10);
  sigma2_lA ~ gamma(1, 10);

  // Loop through each individual
  for (ind in 1:k) {
    // Likelihood distribution for model
    for (trial in 1:nk[ind]) { // cannot vectorize because cannot divide by vector
      // We log mu so that the paramters directly determine the median so it's
      // easier to interpret; thus we use mu -> log(mu) so that
      // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
      mu_n = log(V[ind] + E[ind]*exp(-A[ind]*trial));
      mu_l = log(V[ind] + E[ind]*exp(-A[ind]*trial) + S[ind]) +
             log1m(D[ind]/( 1 + exp(-L[ind]*(trial-H[ind]))));
      sigma_n = log(sigma2_nV[ind] + sigma2_nE[ind]*exp(-sigma2_nA[ind]*trial));
      sigma_l = log(sigma2_lV[ind] + sigma2_lE[ind]*exp(-sigma2_lA[ind]*trial));

      target += mylognormal_lpdf(yn[st+trial] | mu_n, sigma_n);
      target += mylognormal_lpdf(yl[st+trial] | mu_l, sigma_l);
    }
    st += nk[ind];
  }
}

generated quantities {
  real mu_n;
  real mu_l;
  real sigma_n;
  real sigma_l;
  real ynpred[total_length];
  real ylpred[total_length];
  real log_lik[total_length];
  int st = 0;

  // Loop through each individual
  for (ind in 1:k) {
    for (trial in 1:nk[ind]) {
      mu_n = log( V[ind] + E[ind]*exp(-A[ind]*trial) );
      mu_l = log( V[ind] + E[ind]*exp(-A[ind]*trial) + S[ind]) +
             log1m(D[ind]/( 1 + exp(-L[ind]*(trial-H[ind]))) );
      sigma_n = log(sigma2_nV[ind] + sigma2_nE[ind]*exp(-sigma2_nA[ind]*trial));
      sigma_l = log(sigma2_lV[ind] + sigma2_lE[ind]*exp(-sigma2_lA[ind]*trial));

      ynpred[st+trial] = lognormal_rng(mu_n, sigma_n);
      ylpred[st+trial] = lognormal_rng(mu_l, sigma_l);
      log_lik[st+trial] = mylognormal_lpdf(yl[st+trial] | mu_l, sigma_l);
    }
    st += nk[ind];
  }
}
