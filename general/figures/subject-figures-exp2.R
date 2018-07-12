library(runjags)
library(MCMCpack)
library(mcmc)
library(coda)
library(ggplot2)
library(extrafont)
font_import(pattern="Montserrat")

loadfonts(device="win")
load(file = 'experiment-2/data/generated/for_jags_exp_2.Rdata')
load(file = 'general/figures/jags-fits/experiment-2-fit.Rdata')
result2<-as.mcmc(jags.result)

model_mcmc2<-combine.mcmc(result2)
  

model_mcmc2<-as.matrix(model_mcmc2)

###

#function to take subset of posterior
samplePosterior<-function(mcmc_result, sample_size){
  rows<- sample(1:dim(mcmc_result)[1], sample_size)
  sample_posterior<-mcmc_result[rows,]
  return(sample_posterior)
}

#function to extract collumns with subjects parameters
sampleParams<- function(subject, npairs, sample_posterior){
  vars<- colnames(sample_posterior)
  sample_params <- list()
  
  for(i in 1:npairs){
    parameters <-c( paste0('alpha[',subject,',',i,',1]'), 
                    paste0('alpha[',subject,',',i,',2]'), 
                    paste0('subject.adapt.beta[',subject,']'),
                    paste0('subject.adapt.gamma[',subject,']'),
                    paste0('item.beta.learn[',subject,',',i,']'),
                    paste0('item.gamma.learn[',subject,',',i,']'),
                    paste0('item.delta[',subject,',',i,']'),
                    paste0('item.learned[',subject,',',i,']')
    )
    cols<-c()
    for(p in parameters){
      cols<-c(cols,which(vars == p))
    }
    sample_params[[i]]<- sample_posterior[,cols]
  }
  return(sample_params)
}

##function to generate model prediction values based on single set of sample parameters

exponential<-function(alpha,beta,gamma,t){
  rt = alpha*(1+beta*(exp(-gamma*t)-1))
  return(rt)
}

piecewise.exponential<-function(alpha, beta, gamma, beta1, gamma1, delta,learner,t){
  rt = alpha*(1+beta*(exp(-gamma*t)-1))
  if(t>delta){
    rlr = (learner *beta1*(1-exp(-gamma1*(t-delta))))
  }
  else{rlr = 0}
  
  return(rt*(1-rlr))
}

generateModelRT<-function(params, t){
  alpha_1<-params[1]
  alpha_2<-params[2]
  beta<-params[3]
  gamma<-params[4]
  beta.learn<-params[5]
  gamma.learn<-params[6]
  delta<-params[7]
  is.learner<-params[8]
  
  rtu<- mapply(exponential, t, MoreArgs = list(alpha = alpha_1, beta = beta, gamma = gamma))
  rtp<- mapply(piecewise.exponential, t, MoreArgs = list(alpha = alpha_2, beta= beta, gamma=gamma, beta1= beta.learn, gamma1= gamma.learn, delta=delta,learner=is.learner))
  
  return(data.frame(rtu = rtu, rtp=rtp, t = t))
}

modelPredict<-function(params, which_pair,subject_data){
  sub_data = subject_data %>% filter(pair == which_pair) 
  model_predict<-generateModelRT(params, sub_data$t)
  model_predict<-inner_join(model_predict, subject_data, by = 't') %>% mutate(rt_predict = if_else(is_predictable == 0, rtu,rtp)) %>% dplyr::select(rt_predict, t, is_predictable, pair)
  
  return(model_predict)
}

##function to 

generatePlotData<-function(subject_data,sample_params){
  model_plot_data_out<- data.frame(rt_predict = c(), t= c(), is_predictable = c(), pair = c(), index = c())
  j=0
  for(p in sample_params){
    j=j+1
    model_plot_list<-apply(p, MARGIN = 1,function(params){return(modelPredict(params,j,subject_data))})
    model_plot_data<- data.frame(rt_predict = c(), t= c(), is_predictable = c(), pair = c(), index = c())
    i=0
    for(m in model_plot_list){
      i=i+1
      m$pair<- rep(j, dim(m)[1])
      m$index<- rep(i, dim(m)[1])
      model_plot_data<-rbind(model_plot_data, data.frame(m))
    }
    model_plot_data_out<- rbind(model_plot_data_out,model_plot_data)
  }
  return(model_plot_data_out) 
}

##function to add model samples to plot(input grob, output grob)
plotSubjectModel<-function(subject_data,sample_posterior){
  subject = unique(subject_data$subject_id)
  npairs = length(unique(subject_data$pair))
  sample_params <- sampleParams(subject,npairs,sample_posterior)
  
  model_plot_data<-generatePlotData(subject_data,sample_params)

  p<-ggplot()+
    geom_line(data=model_plot_data %>% group_by(pair,index) %>% filter(is_predictable == 1),aes(x=t,y=rt_predict,group = index), alpha = .05, col = "#377eb8")+
    geom_line(data=model_plot_data %>% group_by(pair,index) %>% filter(is_predictable == 0),aes(x=t,y=rt_predict,group = index), alpha = .05, col= "#e41a1c")+
    geom_point(data = subject_data, aes(x=t,y=rt,col=as.character(is_predictable)),size =2)+
    facet_wrap(~pair)+
    ggtitle(paste0('Subject ',subject))+
    #ylim(0,1000)+
    scale_color_manual(guide=F, values=c("#e41a1c", "#377eb8"))+
    labs(y="Response time (ms)", x="Time (discrete presentations of item)")+
    theme_bw(base_size = 28, base_family = "Montserrat")

  return(p)  
}

posterior2<-samplePosterior(model_mcmc2,100)

subject_data2<-for_jags_exp_2 %>% filter(subject_id == 42) 

plotSubjectModel(subject_data2,posterior2)

ggsave("subject-42.png", device="png", path="general/figures/", dpi=300, width=12, height=8, units="in")


###find learners
function(model_mcmc,data){
  vars<-colnames(model_mcmc)
  subjects=unique(data$subject_id)
  for(s in subjects){
    sub_data<- data %>% fitler(subject_id == s)
    npairs = length(unique(sub_data$pair))
    post <- array(dim = c(length(subjects), npairs,dim(model_mcmc)[1]))
    percent <- numeric(length(subjects))
    for(i in 1:npairs){
      par<- paste0('item.learned[',subject,',',i,']')
      post[,,i] = model_mcmc[,which(par == var)]
      percent[i] = mean(post[,,i])
    }
  }
}


