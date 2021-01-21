setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/lab1")

library(TMB)
compile("linreg.cpp")
dyn.load(dynlib("linreg"))

set.seed(123)
beta_0 <- 1
beta_1 <- 2
beta_2 <- .5
sd_er <- 3
n <- 100

tmb_data <- list(x1 = 1:n, x2 = rpois(n, 10))
tmb_data$y <- beta_0 + beta_1 * tmb_data$x + beta_2 sd_er * rnorm(n)

data <- list(Y = rnorm(10) + 1:10, x = 1:10)
parameters <- list(a = 1, b = 2, logSigma = 3)

obj <- MakeADFun(data, parameters, DLL = "linreg")
obj$hessian <- TRUE
opt <- do.call("optim", obj)

opt
opt$hessian ## <-- FD hessian from optim
obj$he() ## <-- Analytical hessian
sdreport(obj)

lm(tmb_data$y ~ tmb_data$x)

tmb_data$y