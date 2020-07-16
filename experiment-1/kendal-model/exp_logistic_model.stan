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

  vector<lower=0>[K] E;  // overall scale of exponential

  vector<lower=0>[K] A;  // adaptation rate

  vector<lower=0, upper=1>[K] P; // probability of learning

  vector<lower=0, upper=1>[K] D_raw;  // scale of learning "drop"

  vector<lower=0>[K] L;  // learning rate

  vector<lower=0, upper=1>[K] H_raw;  // horizontal shift in onset of learning

  vector<lower=0>[K] sigma_2;  // main component of variance in response times
}

transformed parameters { // manipulations of the parameters (really, just their bounds)
  vector[K] H  = upperH .* H_raw;
  vector[K] D = V .* D_raw;
}

model { // define the Bayesian hierarchical model
  // Things we'll use inside the loops
  real mu0;
  real mu1;
  real sig;
  int N;
  int st = 0;

  for (i in 1:K) { // iterate through each individual
    N = NTI[i];

    // Prior distributions on parameters
    V[i] ~ gamma(2.5, 2.5);
    E[i] ~ gamma(2.5, 10);
    A[i] ~ gamma(2.5, 10);
    P[i] ~ beta(0.01, 0.01);
    D[i] ~ beta(2.5, 2.5);
    L[i] ~ gamma(4, 10);
    H[i] ~ cauchy(N*1.0/2, 25);
    sigma_2[i] ~ gamma(2, 10);

    // Likelihood distribution for model
    for (trial in 1:N) { // cannot vectorize because cannot divide by vector
      // We log mu so that the paramters directly determine the median so it's
      // easier to interpret; thus we use mu -> log(mu) so that
      // median(lognorm(z | log(mu), sigma)) = exp(log(mu)) = mu
      mu0 = log( V[i] + E[i]*exp(-A[i]*trial) );
      mu1 = log( V[i] + E[i]*exp(-A[i]*trial) ) + log1m(1/( (1/D[i]) + exp(-L[i]*(trial-H[i]))) );
      sig = log( sigma_2[i] );
      target += mylognormal_lpdf(Y0[st+trial] | mu0, sig);
      target += log_sum_exp(
        log(P[i]) + mylognormal_lpdf(Y1[st+trial] | mu1, sig), // learn
        log1m(P[i]) + mylognormal_lpdf(Y1[st+trial] | mu0, sig) // no learn
      );
    }

    st += N;
  }

}
