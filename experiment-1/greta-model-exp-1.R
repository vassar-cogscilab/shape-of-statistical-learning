library(readr)

data <- read_csv('experiment-1/data/generated/for-jags-exp-1.csv')

n <- 1000#nrow(data)
s <- max(data$subject_id)

##### model ####
library(greta)

# data
rt <- as_data(data$rt)
subject <- data$subject_id
is_predictable <- data$is_predictable
target <- data$target_index
t <- as_data(data$t - 1)

# subject level priors 
beta.adapt <- beta(1,1,dim=s)
rate.adapt <- gamma(1,1, dim=s)

beta.learn <- beta(1,1,dim=s)
rate.learn <- beta(1,1,dim=s)
midpoint.learn <- uniform(0,71, dim=s)

# item level priors
alpha <- normal(1000, 250, dim=c(s,2))

# subject level residuals
sigma.subject <- gamma(1, 0.1, dim=s)

# learning curves
adapt <- greta_array(1, dim=n)
relative.learn <- greta_array(1, dim=n)
intercept <- greta_array(1, dim=n)
mu <- greta_array(1, dim=n)
sigma <- greta_array(1, dim=n)
for(i in 1:n){
  adapt[i] <- (1 + beta.adapt[subject[i]] * (exp(-rate.adapt[subject[i]] * t[i]) - 1))
  relative.learn[i] <- beta.learn[subject[i]] / (1 + exp(-rate.learn[subject[i]] * (t[i] - midpoint.learn[subject[i]])))
  intercept[i] <- alpha[subject[i], target[i]]
  sigma[i] <- sigma.subject[subject[i]]
}

mu <- intercept * adapt * (1 - relative.learn * is_predictable)

# likelihood
distribution(rt) <- normal(mu, sigma)

