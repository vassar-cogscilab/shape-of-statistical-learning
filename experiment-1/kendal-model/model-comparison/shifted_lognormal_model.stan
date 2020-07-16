// In the function block, we define our new probability density function

functions {
  real shiftlognormal_lpdf(real z, real delta, real mu, real sigma) {
    real lprob;

    lprob = (
      -log((z - delta)*sigma*sqrt(2*pi())) -
      (log(z - delta) - mu)^2 / (2*sigma^2)
    );

    return lprob;
  }

  real shiftlognormal_rng(real delta, real mu, real sigma) {
    return delta + lognormal_rng(mu, sigma);
  }
}

// In the data block, we specify everything that is relevant for
// specifying the data. Note that PRIOR_ONLY is a dummy variable used later.

data {
  int<lower=1> n;
  real y[n];
}


// In the parameters block, we specify all parameters we need.
// Although Stan implicitly adds flat prior on the (positive) real line
// we will specify informative priors below.

parameters {
  real<lower=0> r;
  real<lower=0> tau;
  real<lower=0> delta;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> sigma_e;
}

// In the model block, we specify our informative priors and
// the likelihood of the model, unless we want to sample only from
// the prior (i.e., if PRIOR_ONLY == 1)

model {
  real mu;
  // Renormalize Cauchy prior due to truncation to get correct marginal likelihood
  target += cauchy_lpdf(tau | 0, 1) - cauchy_lccdf(0 | 0, 1);
  target += lognormal_lpdf(delta | log(0.50), .5);
  target += lognormal_lpdf(alpha | log(0.50), .5);
  target += lognormal_lpdf(beta | 1, .5);
  target += gamma_lpdf(r | 1, 3);
  target += gamma_lpdf(sigma_e | 0.5, 5);

  for (trial in 1:n) {
    mu = log(alpha + beta * (tau + 1) / (tau + exp(r*trial)));
    target += shiftlognormal_lpdf(y[trial] | delta, mu, sigma_e);
  }
}

// In this block, we make posterior predictions (ypred) and compute
// the log likelihood of each data point given all the others (log_lik)
// which is needed for the computation of loo later

generated quantities {
  real mu;
  real ypred[n];
  real log_lik[n];

  for (trial in 1:n) {
    mu = log(alpha + beta * (tau + 1) / (tau + exp(r*trial)));
    ypred[trial] = shiftlognormal_rng(delta, mu, sigma_e);
    log_lik[trial] = shiftlognormal_lpdf(y[trial] | delta, mu, sigma_e);
  }

}
