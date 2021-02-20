rm(list = ls())

setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/Project1_mess")

load("data\\3LN_redfish.RData")
load("tmb_data.RData")

library(TMB)

compile("fit_mess.cpp", PKG_CXXFLAGS = "-Wno-ignored-attributes")

# compile("fit_origin.cpp")

dyn.load("fit_mess")
dyn.unload("fit_origin")
## set FRV length pattern in catchability;
L50c <- 17
L95c <- 22
b1c <- log(0.95 / 0.05) / (L95c - L50c)
boc <- -b1c * L50c
ql <- 5 * exp(boc + b1c * tmb_data$len_mid) / (1 + exp(boc + b1c * tmb_data$len_mid))
# plot(tmb_data$len_mid,ql,xlab='length',ylab='Catchability Length Pattern',type='l',lwd=2)

A <- tmb_data$A
Y <- tmb_data$Y
parameters <- list(
  # for pla_func()
  log_Linf = log(35),
  log_vbk = log(.35),
  log_len_o = log(.5),
  log_cv_len = log(.1),

  # for computing F
  log_std_log_F = log(.75),
  logit_F_age = log(.5 / .5),
  logit_F_year = log(.5 / .5),
  log_F_main = log(.1),

  # for cohort model
  log_meanR = 3,
  log_std_log_R = log(1),
  logit_log_R = -5,
  log_N0 = rep(0, A - 1),
  log_std_pe = log(.1),

  # for predicted survey index
  log_std_index = log(.5),
  log_std_logq = 3,

  # for predicted catch
  log_std_catch = log(.5),

  # the random effects;
  log_Rec_dev = rep(0, Y),
  log_F_dev = matrix(0, nrow = A - 2, ncol = Y, byrow = T),
  logq = log(ql),
  pe = matrix(0, nrow = A - 1, ncol = Y - 1, byrow = T)
)
parameters$log_N0 <- parameters$log_meanR - 0.2 * (1:(A - 1))

# not estimating correlayion in rec devs or any process error at first
parameters.L <- list(
  # for pla_func()
  log_Linf = log(20),
  log_vbk = log(.01),
  log_len_o = log(1 / 10000),
  log_cv_len = log(.01),

  # for computing F
  log_std_log_F = log(0.01),
  logit_F_age = -5,
  logit_F_year = -5,
  log_F_main = -10,

  # for cohort model
  log_meanR = 1,
  log_std_log_R = log(.01),
  logit_log_R = -5,
  log_N0 = rep(-Inf, A - 1),
  log_std_pe = log(.01),

  # for predicted survey index
  log_std_index = log(.1),
  log_std_logq = log(.01),

  # for predicted catch
  log_std_catch = log(.01)
)

parameters.U <- list(
  # for pla_func()
  log_Linf = log(100),
  log_vbk = log(.8),
  log_len_o = log(1),
  log_cv_len = log(.5),

  # for computing F
  log_std_log_F = log(5),
  logit_F_age = log(.95 / .05),
  logit_F_year = log(.95 / .05),
  log_F_main = log(2),

  # for cohort model
  log_meanR = 100,
  log_std_log_R = log(10),
  logit_log_R = -5,
  log_N0 = rep(Inf, A - 1),
  log_std_pe = log(.5),

  # for predicted survey index
  log_std_index = log(10),
  log_std_logq = log(5),

  # for predicted catch
  log_std_catch = log(2)
)

lower <- unlist(parameters.L)
upper <- unlist(parameters.U)
ql.map <- 1:length(ql)
ql.map[ql.map > 45] <- 45
ql.map[ql.map <= 3] <- NA

pdat <- expand.grid(age = 3:tmb_data$A, year = unique(CLc_vec$Year))
pdat$age[pdat$age > 9] <- 9
F_map <- paste(pdat$year, pdat$age, sep = "_")

log_F_dev_map <- matrix(F_map, nrow = A - 2, ncol = Y, byrow = F)
# log_F_dev_map=matrix(NA,nrow=A-1,ncol=Y,byrow=F)

## first run with log_std_log_F fixed to get starting values;

map <- list(
  log_len_o = factor(NA),
  # log_vbk = factor(NA),
  # log_std_log_F = factor(NA),
  logit_F_age = factor(NA),
  logit_F_year = factor(NA),
  log_std_pe = factor(NA),
  logit_log_R = factor(NA),
  pe = factor(matrix(NA, nrow = A - 1, ncol = Y - 1, byrow = T)),
  logq = factor(ql.map),
  log_F_dev = factor(log_F_dev_map)
)

t.L <- parameters.L
t.L$log_len_o <- NULL
# t.L$log_vbk = NULL
# t.L$log_std_log_F = NULL
t.L$logit_F_age <- NULL
t.L$logit_F_year <- NULL
t.L$log_std_pe <- NULL
t.L$logit_log_R <- NULL
lower <- unlist(t.L)
t.U <- parameters.U
t.U$log_len_o <- NULL
# t.U$log_vbk = NULL
# t.U$log_std_log_F = NULL
t.U$logit_F_age <- NULL
t.U$logit_F_year <- NULL
t.U$log_std_pe <- NULL
t.U$logit_log_R <- NULL
upper <- unlist(t.U)
## random effects;
rname <- c("log_Rec_dev", "log_F_dev", "logq")

obj <- MakeADFun(tmb_data, parameters,
  random = rname, DLL = "fit_mess", map = map,
  # inner.control=list(maxit=500,trace=F),random.start = expression(last.par[random]))
  inner.control = list(maxit = 500, trace = F)
)

length(obj$par)
length(lower)
length(upper)

obj$fn(obj$par)

# obj$gr(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr,
  lower = lower, upper = upper,
  control = list(trace = 0, iter.max = 100, eval.max = 10000)
)

cbind(lower, opt$par, upper)

rep <- obj$report()
sd.rep <- sdreport(obj)

###########  Do the Plotting ##################

source("plots.R")

save.image(file = "fit.RData")

## MISC code #############

ind <- abs(CLc.vec$resid) > 4
CLc.vec[ind, ]

ind1 <- CLc.vec$Year == 2019
CLc.vec[ind1, ]