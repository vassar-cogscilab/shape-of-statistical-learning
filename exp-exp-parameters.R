intercept.adapt <- 1070
intercept.learn <- 1189

beta.adapt <- 0.18
rate.adapt <- 0.28

beta.learn <- 0.4
rate.learn <- .1
onset.learn <- 4

# plot the two curves ####
layout(1:2)
t <- 0:71
adapt <- (1 + beta.adapt * (exp(-rate.adapt * t) - 1))
learn <- sapply(t, function(x){
  if(x < onset.learn) { return(0); }
  else { return( beta.learn * ( 1 + -exp(-rate.learn * (x - onset.learn)))) }
})
plot(t, intercept.adapt*adapt, col="red", ylim=c(0,2000))
points(t, intercept.learn * adapt * (1 - learn), col="blue")

plot(t, learn, ylim=c(0,1))