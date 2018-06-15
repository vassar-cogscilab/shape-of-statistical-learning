library(dplyr)
library(tidyr)
library(readr)

data.all <- read.csv('experiment-1/data/raw/all-data.csv')

n_subjects <- length(unique(data.all$subject))

# test data only ####

data.test <- data.all %>% filter(practice==0)

data.test$appearanceCount <- rep(as.vector(sapply(1:72, function(x){return(rep(x,12))})), n_subjects)

# JAGS data only ####

data.for.model <- subset(data.test, correct==1 & triple_type == 'critical' & triple_position%in%c(1,3), c('subject','cond','appearanceCount','triple_position','rt'))

subject.id <- data.frame(original=unique(data.test$subject))
subject.id$original.numeric <- as.numeric(factor(subject.id$original))

data.for.model$subject_id <- as.numeric(factor(data.for.model$subject))
data.for.model$subject_condition <- factor(as.character(data.for.model$cond), levels=c('known','unknown','one'))
data.for.model$subject_condition <- as.numeric(data.for.model$cond)

data.for.model <- data.for.model %>% mutate(is_predictable = as.numeric(factor(triple_position))-1, target_id = is_predictable+1) %>%
  select(subject_id, rt)

