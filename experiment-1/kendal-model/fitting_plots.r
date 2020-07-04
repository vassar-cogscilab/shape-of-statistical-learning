library("ggplot2")

########## Define Functions ####################################################
exp_log <- function(rt, V, E, A, D, L, H) {
  return( (V + E*exp(-A*rt)) * (1 - D/(1 + exp(-L*(rt-H)))) )
}
sim_exp_log <- function(rt, V, E, A, D, L, H, var) {
  return( exp_log(rt=rt, E=E, A=A, V=V, D=D, L=L, H=H)
          + rnorm(length(rt), 0, var) )
}


########## Plotting Functions ##################################################
plot_fits <- function(k, NTI, y0, y1, fits) {
  st <- 0
  for (i in seq_len(k)) {
    x <- seq_len(NTI[i])
    y0i <- y0[(st+1):(st+NTI[i])]
    y1i <- y1[(st+1):(st+NTI[i])]
    ymin <- min(y0i, y1i) - 0.1
    ymax <- max(y0i, y1i)
    plot(x, y0i, ylim=c(ymin,ymax), pch=20, col='gray60',
         xlab="Trial", ylab="Response Time (sec)",
         main=paste("Subject", i, "Data", sep=" "))
    points(x, y1i, pch=20, col='black')
    lines(x, exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=0, L=fits[[i+9*k]], H=fits[[i+12*k]]), col='#ff87bc')
    lines(x, exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=fits[[i+5*k]]*fits[[i+6*k]], L=fits[[i+9*k]], H=fits[[i+12*k]]), col='#ff5757')
    mtext(paste0('P=', round(fits[[i+5*k]],3), '\n', 'D=', round(fits[[i+6*k]],3), ',  ', 'H=', round(fits[[i+12*k]],3)), side=3, cex=0.75)
    legend("bottomleft", legend=c("No Learn Data", "Learn Data", "No Learn Fit", "Learn Fit"),
           col=c("gray60", "black", "#ff87bc", "#ff5757"),
           pch=c(20,20,NA,NA), lty=c(NA,NA,1,1), cex=0.65)
    st = st + NTI[i]
  }
}

plot_fits_gg <- function(k, NTI, y0, y1, fits, save=FALSE) {
  st <- 0
  for (i in seq_len(k)) {
    n <- NTI[i]
    x = c(seq_len(n), seq_len(n))
    df <- data.frame(x = x,
                     data = c(y0[(st+1):(st+n)], y1[(st+1):(st+n)]),
                     fit = c(exp_log(x[1:n], V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=0, L=fits[[i+9*k]], H=fits[[i+16*k]]),
                             exp_log(x[1:n], V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=fits[[i+5*k]]*fits[[i+6*k]], L=fits[[i+9*k]], H=fits[[i+16*k]])),
                     label_data = c(rep("Non-Learned Data",n), rep("Learned Data",n)),
                     label_fit = c(rep("Non-Learned Fit",n), rep("Learned Fit",n)))
    print(ggplot(df) +
            geom_point(aes(x=x, y=data, color=label_data), shape=20, size=5) +
            geom_line(aes(x=x, y=fit, color=label_fit), size=3) +
            scale_color_manual(values=c("#050d42bf", "#ff1f1f", "#bac3ffbf", "#ffbfbf")) +
            labs(title = paste0("Subject", i, ":  ", "P=", round(fits[[i+5*k]],3)),
                 x = "Trial Number", y = "Response Time (sec)") +
             theme(plot.title = element_text(size=23),
                   axis.text.x = element_text(size = 16),
                   axis.text.y = element_text(size = 16),
                   axis.title.x = element_text(size = 20),
                   axis.title.y = element_text(size = 20),
                   legend.title = element_blank(),
                   legend.position = c(1, 0.9),
                   legend.justification = "right",
                   legend.direction = "vertical",
                   legend.text = element_text(size = 14, angle=0)))
    if (save) {
      ggsave(paste0("Paper/Images/fits_gg_", i, ".png"))
    }
    st = st + n
  }
}

plot_fits_sim <- function(k, NTI, V, E, A, D, L, H, y0, y1, fits) {
  st <- 0
  for (i in seq_len(k)) {
    x <- seq_len(NTI[i])
    y0i <- y0[(st+1):(st+NTI[i])]
    y1i <- y1[(st+1):(st+NTI[i])]
    ymin <- min(y0i, y1i) - 0.1
    ymax <- max(y0i, y1i)
    plot(x, y0i, ylim=c(ymin,ymax), pch=20, col='gray60',
         xlab="Trial", ylab="Response Time (sec)",
         main=paste("Subject", i, "Data", sep=" "))
    points(x, y1i, pch=20, col='black')
    lines(x, exp_log(x, V=V[i], E=E[i], A=A[i], D=0, L=L[i], H=H[i]), col='gray60')
    lines(x, exp_log(x, V=V[i], E=E[i], A=A[i], D=D[i], L=L[i], H=H[i]), col='black')
    lines(x, exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=0, L=fits[[i+9*k]], H=fits[[i+12*k]]), col='#ff87bc')
    lines(x, exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=fits[[i+5*k]]*fits[[i+6*k]], L=fits[[i+9*k]], H=fits[[i+12*k]]), col='#ff5757')
    mtext(paste0('P=', round(fits[[i+5*k]],3), '\n', 'D=', round(fits[[i+6*k]],3), ',  ', 'H=', round(fits[[i+12*k]],3)), side=3, cex=0.75)
    legend("bottomleft", legend=c("No Learn Data", "Learn Data", "No Learn Truth", "Learn Truth", "No Learn Fit", "Learn Fit"),
           col=c("gray60", "black", "gray60", "black", "#ff87bc", "#ff5757"),
           pch=c(20,20,NA,NA,NA,NA), lty=c(NA,NA,1,1,1,1), cex=0.65)
    st = st + NTI[i]
  }
}

plot_fits_sim_gg <- function(k, NTI, V, E, A, D, L, H, y0, y1, fits, save=FALSE) {
  st <- 0
  for (i in seq_len(k)) {
    n <- NTI[i]
    df <- data.frame(x = c(seq_len(n), seq_len(n)),
                     data = c(y0[(st+1):(st+n)], y1[(st+1):(st+n)]),
                     truth = c(exp_log(x, V=V[i], E=E[i], A=A[i], D=0, L=L[i], H=H[i]),
                               exp_log(x, V=V[i], E=E[i], A=A[i], D=D[i], L=L[i], H=H[i])),
                     fit = c(exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=0, L=fits[[i+9*k]], H=fits[[i+16*k]]),
                             exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=fits[[i+5*k]]*fits[[i+6*k]], L=fits[[i+9*k]], H=fits[[i+16*k]])),
                     fit_min = c(exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=0, L=fits[[i+9*k]], H=fits[[i+16*k]]),
                                 exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=fits[[i+5*k]]*fits[[i+6*k]], L=fits[[i+9*k]], H=fits[[i+16*k]]))
                                 - fits[[i+15*k]],
                     fit_max = c(exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=0, L=fits[[i+9*k]], H=fits[[i+16*k]]),
                                 exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=fits[[i+5*k]]*fits[[i+6*k]], L=fits[[i+9*k]], H=fits[[i+16*k]]))
                                 + fits[[i+15*k]],
                     label_data = c(rep("Non-Learned Data",n), rep("Learned Data",n)),
                     label_fit = c(rep("Non-Learned Fit",n), rep("Learned Fit",n)))
    print(ggplot(df) +
            geom_point(aes(x=x, y=data, color=label_data), shape=20) +
            geom_line(aes(x=x, y=truth, color=label_data)) +
            geom_line(aes(x=x, y=fit, color=label_fit)) +
            # geom_ribbon(aes(x=x, ymin=df$fit_min, ymax=df$fit_max, fill=label_fit), linetype="blank", alpha=0.1) +
            scale_color_manual(values=c("black", "#ff1f1f", "gray60", "#ffbfbf")) +
            # scale_fill_manual(values=c("#ff9494", "#ffbfbf")) +
            labs(title = "Simulated Data Fit",
                 x = "Trial Number", y = "Response Time (sec)",
                 subtitle = paste0("Subject", i, "\n",
                                   "P=", round(fits[[i+5*k]],3))) +
             theme(legend.title=element_blank())
          )
    if (save) {
      ggsave(paste0("Paper/Images/sim_gg_", i, ".png"))
    }
    st = st + n
  }
}


plot_fits_gen <- function(k, NTI, y0, y1, fits) {
  st <- 0
  for (i in seq_len(k)) {
    x <- seq_len(NTI[i])
    y0i <- y0[(st+1):(st+NTI[i])]
    y1i <- y1[(st+1):(st+NTI[i])]
    ymin <- min(y0i, y1i) - 0.1
    ymax <- max(y0i, y1i)
    plot(x, y0i, ylim=c(ymin,ymax), pch=20, col='gray60',
         xlab="Trial", ylab="Response Time (sec)",
         main=paste("Subject", i, "Data", sep=" "))
    points(x, y1i, pch=20, col='black')
    lines(x, exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=0, L=fits[[i+9*k]], H=fits[[i+12*k]]), col='#ff87bc')
    lines(x, exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=fits[[i+5*k]]*fits[[i+6*k]], L=fits[[i+9*k]], H=fits[[i+12*k]]), col='#ff5757')
    points(x, sim_exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=0, L=fits[[i+9*k]], H=fits[[i+12*k]], var=fits[[i+15*k]]), pch=18, col='#ff87bc64')
    points(x, sim_exp_log(x, V=fits[[i]], E=fits[[i+3*k]], A=fits[[i+4*k]], D=fits[[i+5*k]]*fits[[i+6*k]], L=fits[[i+9*k]], H=fits[[i+12*k]], var=fits[[i+15*k]]), pch=18, col='#ff575764')
    mtext(paste0('P=', round(fits[[i+5*k]],3), '\n', 'D=', round(fits[[i+6*k]],3)), side=3, cex=0.75)
    legend("bottomleft", legend=c("No Learn Data", "Learn Data", "No Learn Fit", "Learn Fit", "No Learn Gen'd", "Learn Gen'd"),
           col=c("gray60", "black", "#ff87bc", "#ff5757", "#ff87bc", "#ff5757"),
           pch=c(20,20,NA,NA,18,18), lty=c(NA,NA,1,1,NA,NA), cex=0.65)
    st = st + NTI[i]
  }
}
