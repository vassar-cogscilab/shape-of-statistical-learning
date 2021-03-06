library(tidyr)
library(dplyr)

#experiment 1
data.all.exp1 <- read.csv('experiment-1/data/raw/all-data.csv')

test.data.exp1<- data.all.exp1 %>% filter(practice == 0) %>%
  mutate(t = rep(as.vector(sapply(0:71, function(x){return(rep(x,12))})), length(unique(subject))))%>%
  filter(correct==1 & triple_type == 'critical' & triple_position%in%c(1,3))%>% 
  mutate(subject_id = as.numeric(factor(subject)),
         subject_condition = factor(as.character(cond), levels=c('known','unknown','one')),
         is_predictable = as.numeric(factor(triple_position))-1)%>%
  group_by(is_predictable)%>%
  mutate(global_index = ((as.numeric(block)+1)*t)/2)

plot_data_1 <- test.data.exp1 %>% mutate(predictable = if_else(is_predictable == 0, 'unpredictable', 'predictable'))

#experiment 2
library(readr)
data.all.exp2 <- read_csv2('experiment-2/data/raw/sl_size_data.csv')


set.size.condition <- data.all.exp2 %>% filter(sequence_type == 'predictable') %>% 
  group_by(subject_id, grid) %>% 
  summarize() %>% ungroup() %>%
  mutate(grid.fct = as.numeric(factor(grid)),
         set_size = grid.fct * 2 + 2) %>% 
  select(subject_id, set_size)

test.data.exp2 <- data.all.exp2 %>% inner_join(set.size.condition)%>%
  filter(sequence_type == 'unpredictable' | sequence_type == 'predictable') %>% 
  select(subject_id, set_size, block, rt, correct, target, sequence_type)

## add column for which pair, which target

subject.pairs <- test.data.exp2 %>% group_by(subject_id) %>% 
  filter(row_number() <= 16) %>% 
  select(subject_id, target) %>%
  mutate(pair = rep(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8), length(unique(subject_id))))%>% 
  group_by(subject_id, target) %>% 
  filter(row_number() <= 1) %>% 
  ungroup() %>% group_by(subject_id) %>% 
  mutate(target_index = 1:n()) %>% 
  ungroup()

test.data.exp2 <- test.data.exp2 %>% left_join(subject.pairs, by=c('subject_id', 'target'))%>% 
  group_by(subject_id, target) %>% mutate(t = 1:n()) %>% ungroup() %>% select(-target) %>%
  mutate(subject_id = as.numeric(factor(subject_id)),
         is_predictable = as.numeric(factor(sequence_type, levels=c('unpredictable', 'predictable'))) - 1) %>%
  filter(correct == 1) %>% mutate(subject_condition = as.numeric(factor(set_size)))%>%
  group_by(is_predictable)%>%
  mutate(global_index = ((as.numeric(block)+1)*t))


plot_data_2 <- test.data.exp2 %>%   mutate(predictable = if_else(is_predictable == 0, 'unpredictable', 'predictable'))%>%
  mutate(set_size = paste0('2 x ', set_size,' grid'))

#experiment 3

data.all.exp3 <- read_csv('experiment-3/data/raw/sl_mouse_click_data.csv')

test.data.exp3 <- data.all.exp3 %>% filter(phase=='rt-test') %>% 
  group_by(prolific_pid, pair, target_type) %>% mutate(t = 1:n()) %>% ungroup() 

target.id <- test.data.exp3 %>% filter(t==1) %>% group_by(prolific_pid) %>% mutate(target_id = 1:n()) %>% ungroup() %>% select(prolific_pid, pair, target_type, target_id)
subject.condition <- target.id %>% group_by(prolific_pid) %>% summarize(max = n()) %>% mutate(subject_condition = as.numeric(factor(max))) %>% select(prolific_pid, subject_condition)

test.data.exp3 = test.data.exp3 %>% 
  inner_join(target.id)%>% 
  inner_join(subject.condition)%>%
  mutate(subject_id = as.numeric(factor(prolific_pid))) %>%
  mutate(pair_id = pair) %>%
  mutate(rt=as.numeric(rt)) %>%
  mutate(is_predictable = as.numeric(factor(target_type, levels=c('unpredictable', 'predictable'))) - 1)%>%
  mutate(global_index = ((as.numeric(block)+1)*t))

plot_data_3 <- test.data.exp3 %>% mutate(predictable = if_else(is_predictable == 0, 'unpredictable', 'predictable'))%>%
  mutate(set_size = paste0(2*subject_condition+2,' pairs'))%>%
  mutate(pair = as.numeric(pair))
