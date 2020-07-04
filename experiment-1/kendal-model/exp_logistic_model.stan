// logistic model with early exponential adaptation
  // response times are lognormally distributed
  // fit is_pred=0 and is_pred=1 data simultaneously
  // probability of learning is included
  // fit multiple individuals simultaneously
  // hierarchy of parameters

functions { // define functions used in model
  real mylognormal_lpdf(real z, real mu, real sigma) {
    real lprob = -log(z*sigma*sqrt(2*pi())) - (log(z)-mu)^2 / (2*sigma^2);
    return lprob;
  }
}

data { // these values are input to the sampling function
  int<lower=1> K; // number of individuals
  int<lower=1> NTI[K]; // vector of number of trials performed by each individual
  int<lower=K> NK; // total length of all RT data
  vector<lower=0>[NK] Y0; // RT data from the individual for each trial of is_predictable=0
  vector<lower=0>[NK] Y1; // RT data from the individual for each trial of is_predictable=1
}

transformed data { // manipulations of the data to use in parameter bounds
  vector<lower=1>[K] upperH = to_vector(NTI); // (rough) upper bound for H
}

parameters { // define parameters (and their bounds) used in the model
  vector<lower=0, upper=2>[K] V;  // vertical shift of all repsonse times
  vector<lower=0>[K] v_alpha;
  vector<lower=0>[K] v_beta;

  vector<lower=0>[K] E;  // overall scale of exponential

  vector<lower=0>[K] A;  // adaptation rate

  vector<lower=0, upper=1>[K] P; // probability of learning

  vector<lower=0, upper=1>[K] D;  // scale of learning "drop"
  vector<lower=0>[K] d_alpha;
  vector<lower=0>[K] d_beta;

  vector<lower=0>[K] L;  // learning rate
  vector<lower=0>[K] l_alpha;
  vector<lower=0>[K] l_beta;

  vector<lower=0, upper=1>[K] H_raw;  // horizontal shift in onset of learning
  vector<lower=0, upper=1>[K] h_loc_raw;
  vector<lower=0>[K] h_scale;

  vector<lower=0>[K] sigma;  // sort of standard deviation in response times
}

transformed parameters { // manipulations of the parameters (really, just their bounds)
  vector[K] H  = upperH .* H_raw;
  vector[K] h_loc = upperH .* h_loc_raw;
}

model { // define the Bayesian hierarchical model
  // Things we'll use inside the loops
  real mu0;
  real mu1;
  int N;
  int st = 0;

  for (i in 1:K) { // iterate through each individual
    N = NTI[i];

    // Distributions on the hyperparameters
    v_alpha[i] ~ cauchy(2.5, 1);
    v_beta[i] ~ cauchy(1.75, 1);
    d_alpha[i] ~ cauchy(2.5, 1);
    d_beta[i] ~ cauchy(2.5, 1);
    l_alpha[i] ~ cauchy(4, 1);
    l_beta[i] ~ cauchy(10, 1);
    h_loc[i] ~ cauchy(N*1.0/2, 1);
    h_scale[i] ~ cauchy(25, 1);

    // Prior distributions on parameters
    V[i] ~ gamma(v_alpha[i], v_beta[i]);
    E[i] ~ gamma(2.5, 10);
    A[i] ~ gamma(2.5, 10);
    P[i] ~ beta(0.05, 0.05);
    D[i] ~ beta(d_alpha[i], d_beta[i]);
    L[i] ~ gamma(l_alpha[i], l_beta[i]);
    H[i] ~ cauchy(h_loc[i], h_scale[i]);
    sigma[i] ~ gamma(2, 10);

    // Likelihood distribution for model
    for (trial in 1:N) { // cannot vectorize because cannot divide by vector
      // We log mu so that the paramters directly determine the median so it's
      // easier to interpret; thus we use mu -> log(mu) so that
      // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
      mu0 = log( V[i] + E[i]*exp(-A[i]*trial) );
      mu1 = log( V[i] + E[i]*exp(-A[i]*trial) ) + log1m(D[i]/(1 + exp(-L[i]*(trial-H[i]))) );
      target += mylognormal_lpdf(Y0[st+trial] | mu0, sigma[i]);
      target += log_sum_exp(
        log(P[i]) + mylognormal_lpdf(Y1[st+trial] | mu1, sigma[i]), // learn
        log1m(P[i]) + mylognormal_lpdf(Y1[st+trial] | mu0, sigma[i]) // no learn
      );
    }

    st += N;
  }

}
