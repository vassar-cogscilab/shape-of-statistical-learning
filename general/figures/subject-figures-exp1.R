library(runjags)
library(dplyr)
library(mcmc)
library(coda)
library(ggplot2)
library(extrafont)
font_import(pattern="Montserrat")

loadfonts(device="win")

## load chains
load(file = 'general/figures/jags-fits/exp1-result-subset.Rdata')
simplified_exp_1<- read.csv('experiment-1/data/generated/simplified-exp-1.csv')
result<-as.mcmc(jags.result)

#model_mcmc<-combine.mcmc(jags.result)
model_mcmc<-as.matrix(result)
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

modelPredict<-function(params, subject_data){
  model_predict<-generateModelRT(params, subject_data$t)
  model_predict<-inner_join(model_predict, subject_data, by = 't') %>% mutate(rt_predict = if_else(is_predictable == 0, rtu,rtp)) %>% dplyr::select(rt_predict, t, is_predictable)

  return(model_predict)
}

##function to 

generatePlotData<-function(subject_data,sample_params){
  model_plot_list<-apply(sample_params, MARGIN = 1,function(params){return(modelPredict(params,subject_data))})
  model_plot_data<- data.frame(rt_predict = c(), t= c(), is_predictable = c(), index = c())
  i=0
  for(m in model_plot_list){
    i=i+1
    m$index<- rep(i, dim(m)[1])
    model_plot_data<-rbind(model_plot_data, data.frame(m))
  }
  return(model_plot_data) 
}

##function to add model samples to plot(input grob, output grob)
plotSubjectModel<-function(subject_data,sample_posterior){
  subject = unique(subject_data$subject_id)
  sample_params <- sampleParams(subject,sample_posterior)
  condition = unique(subject_data$subject_condition)

  
  model_plot_data<-generatePlotData(subject_data,sample_params)
  
  if(condition == 1){
    p<-ggplot()+
      geom_line(data=model_plot_data %>% filter(is_predictable == 1),aes(x=t,y=rt_predict,group = index), color ="#e41a1c", alpha = .1)+
      geom_line(data=model_plot_data %>% filter(is_predictable == 0),aes(x=t,y=rt_predict,group = index), color = "#e41a1c", alpha = .05)+
      geom_point(data = subject_data %>% filter(is_predictable == 0), aes(x=t,y=rt), alpha = .2,color="#e41a1c",size =2)+
      geom_point(data = subject_data %>% filter(is_predictable == 1), aes(x=t,y=rt), alpha = 1, color= "#e41a1c",size =2)+  
      ggtitle(paste0('Subject ',subject))+
      scale_color_manual(guide=F, values=c("#e41a1c"))+
      labs(y="Response time (ms)", x="Time (discrete presentations of item)")+
      theme_bw(base_size = 28, base_family = "Montserrat")+
      theme(strip.background =element_rect(fill="#f3f3f3ff"))+
      theme(strip.text = element_text(colour = '#595959'))
    
  }
  if(condition == 2){
    p<-ggplot()+
      geom_line(data=model_plot_data %>% filter(is_predictable == 1),aes(x=t,y=rt_predict,group = index), color ="#38761d", alpha = .1)+
      geom_line(data=model_plot_data %>% filter(is_predictable == 0),aes(x=t,y=rt_predict,group = index), color = "#38761d", alpha = .05)+
      geom_point(data = subject_data %>% filter(is_predictable == 0), aes(x=t,y=rt), alpha = .2,color="#38761d",size =2)+
      geom_point(data = subject_data %>% filter(is_predictable == 1), aes(x=t,y=rt), alpha = 1, color= "#38761d",size =2)+
      ggtitle(paste0('Subject ',subject))+
      scale_color_manual(guide=F, values=c("#38761d"))+
      labs(y="Response time (ms)", x="Time (discrete presentations of item)")+
      theme_bw(base_size = 28, base_family = "Montserrat")+
      theme(strip.background =element_rect(fill="#f3f3f3ff"))+
      theme(strip.text = element_text(colour = '#595959'))
    
  }
  if(condition == 3){
    p<-ggplot()+
      geom_line(data=model_plot_data %>% filter(is_predictable == 1),aes(x=t,y=rt_predict,group = index), color ="#377eb8", alpha = .1)+
      geom_line(data=model_plot_data %>% filter(is_predictable == 0),aes(x=t,y=rt_predict,group = index), color = "#377eb8", alpha = .05)+
      geom_point(data = subject_data %>% filter(is_predictable == 0), aes(x=t,y=rt), alpha = .2,color="#377eb8",size =2)+
      geom_point(data = subject_data %>% filter(is_predictable == 1), aes(x=t,y=rt), alpha = 1, color= "#377eb8",size =2)+      ggtitle(paste0('Subject ',subject))+
      scale_color_manual(guide=F, values=c("#377eb8"))+
      labs(y="Response time (ms)", x="Time (discrete presentations of item)")+
      theme_bw(base_size = 28, base_family = "Montserrat")+
      theme(strip.background =element_rect(fill="#f3f3f3ff"))+
      theme(strip.text = element_text(colour = '#595959'))
  }

  return(p)  

}

unique(simplified_exp_1%>% filter(subject_condition == 2)%>%select(subject_id))

posterior<-samplePosterior(model_mcmc,100)
subject_data<-simplified_exp_1%>% filter(subject_id == 106)
plotSubjectModel(subject_data, posterior)
ggsave("subject-170.png", device="png", path="general/figures/", dpi=300, width=12, height=8, units="in")


