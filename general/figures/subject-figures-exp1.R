library(runjags)
library(MCMCpack)
library(mcmc)
library(coda)

## load chains here ??? (need a neat way to do this...)
jags.result4<-as.mcmc(jags.result)

model_mcmc<-combine.mcmc(list(jags.result1,jags.result2, jags.result3,jags.result4))

test<-as.matrix(model_mcmc)
###

#function to take subset of posterior
samplePosterior<-function(mcmc_result, sample_size){
  rows<- sample(1:dim(mcmc_result)[1], sample_size)
  sample_posterior<-mcmc_result[rows,]
  return(sample_posterior)
}

#function to extract collumns with subjects parameters
sampleParams<- function(subject, sample_posterior){
  vars<- colnames(sample_posterior)
  parameters <-c( paste0('alpha[',subject,',1]'), 
                  paste0('alpha[',subject,',2]'), 
                  paste0('beta[',subject,']'),
                  paste0('gamma[',subject,']'),
                  paste0('beta.learn[',subject,']'),
                  paste0('gamma.learn[',subject,']'),
                  paste0('delta[',subject,']'),
                  paste0('is.learner[',subject,']')
                )
  cols<-c()
  for(p in parameters){
    cols<-c(cols,which(vars == p))
  }
  sample_params<- sample_posterior[,cols]
  return(sample_params)
}

sampleParams(1,samplePosterior(test,27))

##function to generate model prediction values based on single set of sample parameters

exponential<-function(alpha,beta,gamma,t){
  rt = alpha*(1+beta*(exp(-gamma*t)-1))
  return(rt)
}

piecewise.exponential<-function(alpha, beta, gamma, beta1, gamma1, delta,learner,t){
  rt = alpha*(1+beta*(exp(-gamma*t)-1))
  if(t>delta){
    rt = rt*(1-(learner *beta1*(1-exp(-gamma1*(t-delta)))))
  }
  return(rt)
}

generateModelRT<-function(params, t){
  alpha_1<-params[1]
  alpha_2<-params[2]
  beta<-params[3]
  gamma<-params[4]
  beta.learn<-params[5]
  gamma.learn<-params[5]
  delta<-params[6]
  is.learner<-params[7]
  
  rtu<- mapply(exponential, t, MoreArgs = list(alpha = alpha_1, beta = beta, gamma = gamma))
  rtp<- mapply(piecewise.exponential, t, MoreArgs = list(alpha = alpha_2, beta= beta, gamma=gamma, beta1= beta.learn, gamma1= gamma.learn, delta=delta,learner=is.learner))
  
  return(list(rtu = rtu, rtp=rtp))
}

##function to plot subject data (returns grob)



##function to add model samples to plot(input grob, output grob)
function(data,subject,mcmc_result, sample_size){
  
  
  
}
