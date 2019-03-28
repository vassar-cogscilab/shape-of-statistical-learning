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
  mutate(rt = as.numeric(rt)) %>%
  mutate(subject_id = as.numeric(factor(prolific_pid))) %>%
  mutate(pair_id = as.numeric(pair)+1) %>%
  mutate(is_predictable = as.numeric(factor(target_type, levels=c('unpredictable', 'predictable'))) - 1) %>%
  select(subject_id, subject_condition, t, rt, pair_id, target_id, is_predictable)

# output data

write_csv(model.data, path='experiment-3/data/generated/simplified-exp-3.csv')

# generate list for jags

data.for.jags <- list(
  rt = model.data$rt,
  subject_id = model.data$subject_id,
  is_predictable = model.data$is_predictable,
  pair = as.numeric(model.data$pair_id),
  t = model.data$t,
  N = length(model.data$rt),
  S = length(unique(model.data$subject_id)),
  C = length(unique(model.data$subject_condition)),
  condition = (model.data %>% group_by(subject_id) %>% summarise(s_condition = mean(subject_condition)) %>% arrange(subject_id))$s_condition,
  max_t = (model.data%>%group_by(subject_id)%>% summarise(max_t = max(t)) %>% select(max_t) %>% as.matrix)[,1],
  P = (model.data %>% group_by(subject_id) %>% summarise(n.pairs = length(unique(pair_id))) %>% select(n.pairs) %>% as.matrix)[,1]
)

save(data.for.jags, file="experiment-3/data/generated/jags-data-exp-3.Rdata")
