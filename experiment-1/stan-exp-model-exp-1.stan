data {
  int<lower=0> n;
  int<lower=0> s;
  int<lower=0, upper=71> t[n];
  int<lower=0, upper=1> is_predictable[n];
  int<lower=1, upper=s> subject[n];
  int<lower=1, upper=2> target_id[n];
  vector<lower=0>[n] rt;
}

parameters {
  matrix<lower=0>[s,2] subject_intercept;
  vector<lower=0>[s] sigma_subject;
  vector<lower=0, upper=1>[s] beta_adapt;
  vector<lower=0>[s] rate_adapt;
  vector<lower=0, upper=1>[s] beta_learn;
  vector<lower=0>[s] rate_learn;
  vector<lower=0, upper=71>[s] midpoint_learn;
}

transformed parameters {
  vector[n] mu;
  vector[n] intercept;
  vector[n] adapt;
  vector[n] learn;
  vector[n] sigma;
  for(i in 1:n){
    intercept[i] = subject_intercept[subject[i], target_id[i]];
    adapt[i] = (1 + beta_adapt[subject[i]] * (exp(-rate_adapt[subject[i]] * t[i]) - 1));
    learn[i] = beta_learn[subject[i]] / (1 + exp(-rate_learn[subject[i]] * (t[i] - midpoint_learn[subject[i]])));
    sigma[i] = sigma_subject[subject[i]];
    mu[i] = intercept[i] * adapt[i] * (1 - learn[i] * is_predictable[i]);
  }
}

model {
  sigma_subject ~ gamma(1, .01);
  subject_intercept[,1] ~ normal(1000, 250);
  subject_intercept[,2] ~ normal(1000, 250);
  beta_adapt ~ beta(1,1);
  rate_adapt ~ gamma(1,1);
  beta_learn ~ beta(1,1);
  rate_learn ~ gamma(1,1);
  midpoint_learn ~ uniform(0,71);
  rt ~ normal(mu, sigma);
}