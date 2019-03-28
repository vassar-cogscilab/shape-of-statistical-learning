library(runjags)
library(MCMCpack)
library(mcmc)
library(coda)
library(ggplot2)

#mcmc
load(file = 'general/figures/jags-fits/exp1-result-subset.Rdata')

result.mcmc<-as.mcmc(jags.result)
eff.sample<- effectiveSize(result.mcmc)
eff.sample[which(eff.sample == min(eff.sample))]
autocorr(result.mcmc)
geweke.diag(result.mcmc)

#plotting conditions
result<-as.matrix(result.mcmc)

##extract variables
exp1ExtractPosteriorCondition<-function(result){
  vars= colnames(result)
  condition_posterior<-matrix(ncol = 8)
  for(i in 1:3){
    parameters = c(paste0('rate.learn.mode[',i,']'),
                   paste0('rate.learn.sd[',i,']'),
                   paste0('beta.learn.mode[',i,']'),
                   paste0('beta.learn.concentration[',i,']'),
                   paste0('delta.mode[',i,']'),
                   paste0('delta.concentration[',i,']'),
                   paste0('prob.learn[',i,']')
                   )
    cols <-c()
    for(p in parameters){
      cols = c(cols, which(p == vars))
    }
  
    post<-cbind(result[,cols], rep(i, dim(result[,cols])[1]))
    print(post)
    condition_posterior<- rbind(condition_posterior, post)
  }
  condition_posterior<- as.data.frame(condition_posterior)
  names(condition_posterior)<- c('rate.learn.mode', 'rate.learn.sd', 'beta.learn.mode', 'beta.learn.concentration', 'delta.mode', 'delta.concentration', 'prob.learn', 'condition')
  condition_posterior<-condition_posterior%>%filter(!is.na(condition)) %>% mutate(condition = if_else(condition == 1, 'Known', if_else(condition == 2, 'Novel','Scrambled')))
  
  return(condition_posterior)
}


##plot provability of learning histogram for each condition
conditionPlotProbLearn<-function(condition_posterior){


p <- ggplot(condition_posterior)+
  geom_histogram(aes(x= prob.learn, fill = condition))+
  facet_grid(condition~.)+
  theme_minimal(base_size = 28, base_family = "Montserrat")+
  theme(strip.background =element_rect(fill="#f3f3f3ff"))+
  theme(strip.text = element_text(colour = '#595959'))+
  scale_color_manual(guide=F, values=c("#e41a1c", "#38761d","#377eb8" ))+
  xlab('probability of learning')+
  ylab('')
return(p)
}

conditionPlotProbLearn(condition_posterior) + ggtitle('Experiment 1')
ggsave("condition-prob.learn-exp1.png", device="png", path="general/figures/", dpi=300, width=14, height=8, units="in")

## subject level figures

###find learners
findLearners<-function(model_mcmc,data,exp){
  vars<-colnames(model_mcmc)
  subjects=unique(data$subject_id)
  learner = list()
  for(s in subjects){
    sub_data<- data %>% filter(subject_id == s)
    if(exp == 1){
      par<- paste0('is.learner[',s,']')
      learner[[s]] = mean(model_mcmc[,which(par == vars)])
    }
    if(exp ==2 || exp == 3){
      npairs = length(unique(sub_data$pair_id))
      percent <- numeric(npairs)
    for(i in 1:npairs){
        par<- paste0('item.learned[',s,',',i,']')
      #post[,,i] = model_mcmc[,which(par == vars)]
      #percent[i] = mean(post[,,i])
      percent[i] = mean(model_mcmc[,which(par == vars)])
    }
    learner[[s]]<-percent
    }
  }
  return(learner)
}

exp1<-read.csv('experiment-1/data/generated/simplified-exp-1.csv')
exp2<-read.csv('experiment-2/data/generated/simplified-exp-2.csv')
exp3<-read.csv('experiment-3/data/generated/simplified-exp-3.csv')

test<-findLearners(model_mcmc3,exp3,3)

exp3%>% filter(subject_id,subject_condition,pair_id)

exp3%>% group_by(subject_conditon) %>%

##################3
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


########

##split point
ggplot(condition_posterior)+
  geom_histogram(aes(x= delta.mode, fill = condition))+
  facet_grid(condition~.)+
  theme_minimal(base_size = 28, base_family = "Montserrat")+
  theme(strip.background =element_rect(fill="#f3f3f3ff"))+
  theme(strip.text = element_text(colour = '#595959'))+
  scale_color_manual(guide=F, values=c("#e41a1c", "#38761d","#377eb8" ))


##get a &b to generate condition level distributions
s<-sample(1:dim(result)[1],1000)

getdistvals<-function(params){
  rate.learn.mode<- params[,1]
  rate.learn.sd<- params[,2]
  beta.learn.mode <- params[,3]
  beta.learn.concentration<-params[,4]
  delta.mode <- params[,5]
  delta.concentration<-params[,6]
  prob.learn<-params[,7]
 
  
   
  rate.learn.params<-apply(data.frame(mode = beta.learn.mode, concentration = beta.learn.concentration),MARGIN = 1,getBetaAB)
  return(rate.learn.params)  
}






library(dplyr)
library(tidyr)

temp<-data.frame(t(test[s,]))
temp<-temp%>% gather(key = 'sample', value = 'val')
x<-seq(0.001,.999,by=.001)
temp$x = rep(x, dim(temp)[1])


ggplot(temp)+
  geom_line(aes(x=x,y=val,group = sample))

getBetaAB<-function(params){
  mode = params[1]
  concentration = params[2]
  x<-seq(0.001,.999,by=.001)
  a<-mode * (concentration-2) + 1
  b<-(1 - mode) * (concentration-2) + 1
  
  dist<-dbeta(x,a,b)
  return(dist)
}

