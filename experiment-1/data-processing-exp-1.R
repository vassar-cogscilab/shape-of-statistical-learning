library(dplyr)
library(tidyr)
library(readr)

data.all <- read.csv('experiment-1/data/raw/all-data.csv')

n_subjects <- length(unique(data.all$subject))

# test data only ####

data.test <- data.all %>% filter(practice==0)

data.test$t <- rep(as.vector(sapply(0:71, function(x){return(rep(x,12))})), n_subjects)

# JAGS-relevant data only ####

data.for.model <- subset(data.test, correct==1 & triple_type == 'critical' & triple_position%in%c(1,3), c('subject','cond','t','triple_position','rt'))

subject.id <- data.frame(original=unique(data.test$subject))
subject.id$original.numeric <- as.numeric(factor(subject.id$original))

data.for.model$subject_id <- as.numeric(factor(data.for.model$subject))
data.for.model$subject_condition <- factor(as.character(data.for.model$cond), levels=c('known','unknown','one'))
data.for.model$subject_condition <- as.numeric(data.for.model$cond)

data.for.model <- data.for.model %>% mutate(is_predictable = as.numeric(factor(triple_position))-1, target_index = is_predictable+1, pair=1) %>%
  select(subject_id, rt, pair, target_index, t, is_predictable, subject_condition)

# output csv for plotting

write_csv(data.for.model, path = "experiment-1/data/generated/simplified-exp-1.csv")

# generate list for jags

data.for.jags <- list(
  rt = data.for.model$rt,
  subject_id = data.for.model$subject_id,
  is.predictable = data.for.model$is_predictable,
  t = data.for.model$t,
  N = length(data.for.model$rt),
  S = length(unique(data.for.model$subject_id)),
  C = length(unique(data.for.model$subject_condition)),
  condition = (data.for.model %>% group_by(subject_id) %>% summarise(s_condition = mean(subject_condition)) %>% arrange(subject_id))$s_condition,
  max_t = (data.for.model %>% group_by(subject_id) %>% summarise(max_t = max(t)) %>% select(max_t) %>% as.matrix)[,1],
  n.pairs = (data.for.model %>% group_by(subject_id) %>% summarise(n.pairs = length(unique(pair))) %>% select(n.pairs) %>% as.matrix)[,1]
)

save(data.for.jags, file="experiment-1/data/generated/jags-data-exp-1.Rdata")

