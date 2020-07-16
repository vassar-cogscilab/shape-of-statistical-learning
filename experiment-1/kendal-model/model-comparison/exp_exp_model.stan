// logistic model with early exponential adaptation
  // response times are lognormally distributed
  // fit is_pred=0 and is_pred=1 data simultaneously
  // probability of learning is included
  // fit multiple individuals simultaneously
  // hierarchy of parameters

functions { // define functions used in model
  real mylognormal_lpdf(real z, real mu, real sigma_2) {
    real lprob = -log(z*sqrt(2*pi()*sigma_2)) - (log(z)-mu)^2 / (2*sigma_2);
    return lprob;
  }
}

data { // these values are input to the sampling function
  int<lower=1> n; // vector of number of trials performed by the individual
  vector<lower=0>[n] y; // RT data from the individual for each trial of is_predictable=1
}

parameters { // define parameters (and their bounds) used in the model
  // non-learned parameters
  real<lower=0, upper=2> Vn;  // vertical shift of all repsonse times
  real<lower=0> En;  // overall scale of exponential
  real<lower=0> An;  // adaptation rate

  // learned parameters
  real<lower=0, upper=n> H;  // horizontal shift in onset of learning
  real<lower=0, upper=2> Vl;  // vertical shift of all repsonse times
  real<lower=0> El;  // overall scale of exponential
  real<lower=0> Al;  // adaptation rate

  // sigma squared
  real<lower=0> sigma_2;  // main component of variance in response times
}

model { // define the Bayesian hierarchical model
  // Things we'll use inside the loops
  real mu;
  real sig;

  // Prior distributions on parameters
  Vn ~ gamma(2.5, 2.5);
  En ~ gamma(2.5, 10);
  An ~ gamma(2.5, 10);

  H ~ cauchy(n*1.0/2, 25);
  Vl ~ gamma(2.5, 2.5);
  El ~ gamma(2.5, 10);
  Al ~ gamma(2.5, 10);

  sigma_2 ~ gamma(2, 10);
  sig = log( sigma_2 );

  // Likelihood distribution for model
  for (trial in 1:n) { // cannot vectorize because cannot divide by vector
    // We log mu so that the paramters directly determine the median so it's
    // easier to interpret; thus we use mu -> log(mu) so that
    // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
    if (trial > H) {
      mu = log( Vl + El*exp(-Al*(trial-H)) );
    } else {
      mu = log( Vn + En*exp(-An*trial) );
    }

    target += mylognormal_lpdf(y[trial] | mu, sig);
  }
}

generated quantities {
  real mu;
  real sig;
  real ypred[n];
  real log_lik[n];

  sig = log( sigma_2 );

  for (trial in 1:n) {
    if (trial > H) {
      mu = log( Vl + El*exp(-Al*trial) );
    } else {
      mu = log( Vn + En*exp(-An*trial) );
    }

    ypred[trial] = lognormal_rng(mu, sig);
    log_lik[trial] = mylognormal_lpdf(y[trial] | mu, sig);
  }
}
