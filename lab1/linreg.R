setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab1")
library(TMB)

compile("linreg.cpp")

dyn.load(dynlib("linreg"))
# dyn.unload(dynlib("linreg"))


beta0 <- 1
beta1 <- 2
beta2 <- -0.5
sd_err <- 1
n <- 500

set.seed(123)
tmb_data <- list(x1 = sort(rnorm(n) * 2))
tmb_data$x2 <- (tmb_data$x1 - mean(tmb_data$x1))**2

tmb_data$Y <- beta0 + beta1 * tmb_data$x1 + beta2 * tmb_data$x2 + sd_err * rnorm(n)

plot(tmb_data$x1, tmb_data$Y)

parameters <- list(beta0 = 0, beta1 = 1, beta2 = 2, logSigma = 0)

obj <- MakeADFun(tmb_data, parameters, DLL = "linreg")
# obj$hessian <- TRUE

obj$gr(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr)

rep <- obj$report()
# opt$hessian ## <-- FD hessian from optim
# obj$he()    ## <-- Analytical hessian
sdrep <- sdreport(obj)

lines(tmb_data$x1, rep$mu, col = "red", lwd = 2)
lines(tmb_data$x1, rep$mu - qnorm(0.975) * sdrep$sd, col = "red", lwd = 2, lty = 2)
lines(tmb_data$x1, rep$mu + qnorm(0.975) * sdrep$sd, col = "red", lwd = 2, lty = 2)