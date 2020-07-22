// early exponential adaptation with logistic learning curve
  // learning curve is a generalized (possibly asymmetric) logistic function
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

  real<lower=0, upper=1> P; // probability of learning
  real<lower=0, upper=1> D;  // scale of learning "drop"
  real<lower=0> L;  // learning rate
  real<lower=0, upper=n> H;  // horizontal shift in onset of learning

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
  P ~ beta(0.01, 0.01);
  D ~ beta(2.5, 2.5);
  L ~ gamma(1.5, 0.25);
  H ~ cauchy(n*1.0/2, 25);
  sigma_2 ~ gamma(2, 10);
  sig = log( sigma_2 );

  // Likelihood distribution for model
  for (trial in 1:n) { // cannot vectorize because cannot divide by vector
    // We log mu so that the paramters directly determine the median so it's
    // easier to interpret; thus we use mu -> log(mu) so that
    // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
    mun = log( V + E*exp(-A*trial) );
    mul = log( V + E*exp(-A*trial) ) + log1m(D/( 1 + exp(-L*(trial-H))) );

    target += mylognormal_lpdf(yn[trial] | mun, sig);
    target += log_sum_exp(
      log(P) + mylognormal_lpdf(yl[trial] | mul, sig),
      log1m(P) + mylognormal_lpdf(yl[trial] | mun, sig)
    );
  }
}

generated quantities {
  real mun;
  real mul;
  real sig;
  real ylpred[n];
  real log_lik[n];

  sig = log( sigma_2 );

  for (trial in 1:n) {
    mun = log( V + E*exp(-A*trial) );
    mul = log( V + E*exp(-A*trial) ) + log1m(D/( 1 + exp(-L*(trial-H))) );

    ylpred[trial] = lognormal_rng(log_sum_exp(log(P) + mul, log1m(P) + mun),
                                  sig);
    log_lik[trial] = log_sum_exp(
      log(P) + mylognormal_lpdf(yl[trial] | mul, sig),
      log1m(P) + mylognormal_lpdf(yl[trial] | mun, sig)
    );
  }
}
