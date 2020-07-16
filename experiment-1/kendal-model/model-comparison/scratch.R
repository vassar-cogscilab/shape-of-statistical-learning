x <- seq_len(100)
Vn <- 0.71
En <- 0.25
An <- 0.32
H <- 48.12
Vl <- 0.38
El <- 0.17
Al <- 0.33

y <- c( (Vn + En*exp(-An*x[x<=H])), (Vl + El*exp(-Al*x[x>H])) )

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


A <- 1
K <- 2
B <- .25
nu <- 1
Q <- 1/.4
M <- 50
C <- .4

t <- seq_len(100)
logist <- ( (1/C) + 1/C * exp( -B*(t-M) ) )^(-1/nu)
plot(t, logist, type = "l", lwd = 4, col = "skyblue")
  abline(h = 0.5*C^(1/nu), lty = "dashed", lwd = 2, col = "firebrick")
  # abline(h = (Q+1)^(1/Q), lty = "dotted", lwd = 2, col = "gold")
  abline(v = M, lty = "dashed", lwd = 2, col = "firebrick")
  abline(v = t[which.min(abs(logist - 0.5*C^(1/nu)))], lty = "dashed", lwd = 2, col = "seagreen4")









r <- 0.8
tau <- 50
delta <- 0.5
alpha <- 0.19
beta <- 0.56

y_sl <- alpha + beta * ( (tau + 1)/(tau + exp(r * x)) )
plot(x, y_sl, type = "l")
