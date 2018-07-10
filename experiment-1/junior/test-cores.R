library(parallel)
cores <- detectCores()
print(cores)

writeLines(as.character(cores), paste0('testcoresout-',sample(1:100000,1),'.txt'))