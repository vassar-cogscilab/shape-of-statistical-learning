library(readr)
library(dplyr)

all.data <- read_csv('experiment-3/data/raw/sl_mouse_click_data.csv')

# count number of subjects
n.subjects <- length(unique(all.data$prolific_pid))


# extract testing data

test.data <- all.data %>% filter(phase=='rt-test')

test.data <- test.data %>% group_by(prolific_pid, pair, target_type) %>% mutate(t = 1:n()) %>% ungroup()
target.id <- test.data %>% filter(t==1) %>% group_by(prolific_pid) %>% mutate(target_id = 1:n()) %>% ungroup() %>% select(prolific_pid, pair, target_type, target_id)
test.data <- test.data %>% inner_join(target.id)
subject.condition <- target.id %>% group_by(prolific_pid) %>% summarize(max = n()) %>% mutate(subject_condition = as.numeric(factor(max))) %>% select(prolific_pid, subject_condition)
test.data <- test.data %>% inner_join(subject.condition)

# get model data

model.data <- test.data %>%
  mutate(subject_id = as.numeric(factor(prolific_pid))) %>%
  mutate(pair_id = pair) %>%
  mutate(is_predictable = as.numeric(factor(target_type, levels=c('unpredictable', 'predictable'))) - 1) %>%
  select(subject_id, subject_condition, t, rt, pair_id, target_id, is_predictable)

# output data

write_csv(model.data, path='experiment-3/data/generated/for-jags-exp-3.csv')
