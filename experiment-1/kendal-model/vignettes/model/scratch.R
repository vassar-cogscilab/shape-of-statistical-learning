x <- seq_len(100)
Vn <- 0.71
En <- 0.25
An <- 0.32
H <- 48.12
Vl <- 0.38
El <- 0.17
Al <- 0.33

y <- c( (Vn + En*exp(-An*x[x<=H])), (Vl + El*exp(-Al*(x[x>H]-H))) )
plot(x, y, type = "l")




V <- 0.73
E <- 0.35
A <- 0.29
D <- 0.44
L <- 0.4
H <- 48.70

y1 <- (V + E*exp(-A*x)) * (1 - 1/( (1/D) + exp(-L*(x-H))))
y2 <- (V + E*exp(-A*x)) * (1 - D/( 1 + exp(-L*(x-H))))
plot(x, y1, type = "l", lwd = 4, col = "darkblue")
  lines(x, y2, type = "l", lwd = 4, col = "red")



t <- seq_len(100)
A <- 0
K <- 1
B <- .25
M <- 50
C <- 1.5
Q <- .1
nu <- .1

(K-A)*C^(-1/nu)

logist <- A + (K - A) * ( C + Q * exp( -B*(t-M) ) )^(-1/nu)
logist_sl <- ( (K-A) * B * Q / nu ) * exp(-B*(t-M)) * ( C + Q * exp(-B*(t-M)) )^(-(1+nu)/nu)
plot(t, logist, type = "l", lwd = 4, col = "skyblue", ylim = c(0, (K-A)*C^(-1/nu)))
  abline(h = 0.5*(K-A)*C^(-1/nu), lty = "dashed", lwd = 2, col = "firebrick")
  abline(v = M, lty = "dashed", lwd = 2, col = "firebrick")
  # abline(h = (Q+1)^(1/Q), lty = "dotted", lwd = 2, col = "gold")
  abline(v = t[which.min(abs(logist - 0.5*(K-A)*C^(-1/nu)))], lty = "dotted", lwd = 2, col = "seagreen4")
  lines(t, logist_sl, lwd = 4, col = "gold2")
  abline(v = t[which.max(abs(logist_sl))], lty = "dotted", lwd = 2, col = "gold2")


logist_Y0 <- A + (K - A) * (C + 1)^(-1/nu)




x <- seq(-10, 10, by = 0.1)
y <- 1/x
plot(x, y, type = "l", lwd = 4)




r <- 0.8
tau <- 50
delta <- 0.5
alpha <- 0.19
beta <- 0.56

y_sl <- alpha + beta * ( (tau + 1)/(tau + exp(r * x)) )
plot(x, y_sl, type = "l")







x <- seq(-3, 3, by = 0.01)
plot(x, dnorm(x, 0, .5), type = "l", col = "red")
  lines(x, dnorm(x, 0, 1), type = "l", col = "blue")





x <- seq(0, 3, by = 0.01)
plot(x, dgamma(x, 3, 2), type = "l")
