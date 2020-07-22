// early exponential adaptation with no learning curve
  // there is no learning curve
  // response times are lognormally distributed
  // fit is_pred=0 and is_pred=1 data simultaneously

functions { // define functions used in model
  real mylognormal_lpdf(real z, real mu, real sigma_2) {
    real lprob = -log(z*sqrt(2*pi()*sigma_2)) - (log(z)-mu)^2 / (2*sigma_2);
    return lprob;
  }
}

data { // these values are input to the sampling function
  int<lower=1> n; // number of trials performed by the individual
  vector<lower=0>[n] yn; // RT data from the individual for each trial of is_predictable=0
  vector<lower=0>[n] yl; // RT data from the individual for each trial of is_predictable=1
}

parameters { // define parameters (and their bounds) used in the model
  real<lower=0, upper=2> V;  // vertical shift of all repsonse times

  real<lower=0> E;  // overall scale of exponential
  real<lower=0> A;  // adaptation rate

  real S; // shift in median learned response time relative to non-learned

  real<lower=0> sigma_2;  // main component of variance in response times
}

model {
  real mun;
  real mul;
  real sig;

  // Prior distributions on parameters
  V ~ gamma(2.5, 2.5);
  E ~ gamma(2.5, 10);
  A ~ gamma(2.5, 10);
  S ~ normal(0, 0.5);
  sigma_2 ~ gamma(2, 10);
  sig = log( sigma_2 );

  // Likelihood distribution for model
  for (trial in 1:n) { // cannot vectorize because cannot divide by vector
    // We log mu so that the paramters directly determine the median so it's
    // easier to interpret; thus we use mu -> log(mu) so that
    // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
    mun = log( V + E*exp(-A*trial) );
    mul = log( V + E*exp(-A*trial) + S);

    target += mylognormal_lpdf(yn[trial] | mun, sig);
    target += mylognormal_lpdf(yl[trial] | mul, sig);
  }
}

generated quantities {
  real mun;
  real mul;
  real sig;
  real ynpred[n];
  real ylpred[n];
  real log_lik[n];

  sig = log( sigma_2 );

  for (trial in 1:n) {
    mun = log( V + E*exp(-A*trial) );
    mul = log( V + E*exp(-A*trial) + S);
    
    ynpred[trial] = lognormal_rng(mun, sig);
    ylpred[trial] = lognormal_rng(mul, sig);
    log_lik[trial] = mylognormal_lpdf(yl[trial] | mul, sig);
  }
}
