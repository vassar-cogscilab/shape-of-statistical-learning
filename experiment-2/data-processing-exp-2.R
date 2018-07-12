library(dplyr)
library(ggplot2)
library(readr)

all.data <- read_csv2('experiment-2/data/raw/sl_size_data.csv')

# clean columns ####

set.size.condition <- all.data %>% filter(sequence_type == 'predictable') %>% group_by(subject_id, grid) %>% summarize()
set.size.condition$grid.fct <- as.numeric(factor(set.size.condition$grid))
set.size.condition <- set.size.condition %>% mutate(set_size = grid.fct * 2 + 2) %>% select(subject_id, set_size)

all.data <- all.data %>% inner_join(set.size.condition)

# count the number of subjects ####

n.subjects <- length(unique(all.data$subject_id))

# gather data from test trials ####

test.data <- all.data %>% filter(sequence_type == 'unpredictable' | sequence_type == 'predictable') %>% select(subject_id, set_size, block, rt, correct, target, sequence_type)

## add column for which pair, which target

subject.pairs <- test.data %>% group_by(subject_id) %>% filter(row_number() <= 16) %>% select(subject_id, target)
subject.pairs$pair <- rep(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8), n.subjects)
subject.pairs <- subject.pairs %>% group_by(subject_id, target) %>% filter(row_number() <= 1) %>% ungroup() %>% group_by(subject_id) %>% mutate(target_index = 1:n()) %>% ungroup()

test.data <- test.data %>% left_join(subject.pairs, by=c('subject_id', 'target'))

## add column for target repetions

test.data <- test.data %>% group_by(subject_id, target) %>% mutate(t = 1:n()) %>% ungroup() %>% select(-target)

# output data that is JAGS-friendly as CSV

test.data$subject_id <- as.numeric(factor(test.data$subject_id))
test.data$is_predictable <- as.numeric(factor(test.data$sequence_type, levels=c('unpredictable', 'predictable'))) - 1
test.data <- test.data %>% filter(correct == 1) %>% mutate(subject_condition = as.numeric(factor(set_size))) %>% select(-sequence_type, -block, -correct, -set_size)

write_csv(test.data, path="experiment-2/data/generated/simplified-exp-2.csv")

# generate list for jags

data.for.jags <- list(
  rt = test.data$rt,
  subject_id = test.data$subject_id,
  is_predictable = test.data$is_predictable,
  pair = test.data$pair,
  t = test.data$t,
  N = length(test.data$rt),
  S = length(unique(test.data$subject_id)),
  C = length(unique(test.data$subject_condition)),
  condition = test.data$subject_condition,
  max_t = (test.data%>%group_by(subject_id)%>% summarise(max_t = max(t)) %>% select(max_t) %>% as.matrix)[,1],
  P = (test.data %>% group_by(subject_id) %>% summarise(n.pairs = length(unique(pair))) %>% select(n.pairs) %>% as.matrix)[,1]
)

save(data.for.jags, file="experiment-2/data/generated/jags-data-exp-2.Rdata")


