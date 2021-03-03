rm(list = ls())

setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/Project1_mess")

load("data\\3LN_redfish.RData")
load("tmb_data.RData")

library(TMB)

compile("fit_mess_m.cpp", PKG_CXXFLAGS = "-Wno-ignored-attributes")

# fit_mess_FRV
# fit_mess_m

dyn.load("fit_mess_m")
dyn.unload("fit_mess_m")

## set FRV length pattern in catchability;

FRV <- c(17, 22)
SRV <- c(17, 22)
Spanish_summer <- c(17, 22)
Spanish_3N <- c(17, 22)
qlfix <- 25   # at the length almost 100% fish retained
L50_95 <- list(FRV, SRV, Spanish_summer, Spanish_3N)
b1c <- list()
boc <- list()
ql <- list()

for (i in 1:ns) {
  b1c[[i]] <- log(0.95 / 0.05) / (L50_95[[i]][2] - L50_95[[i]][1])
  boc[[i]] <- -b1c[[i]] * L50_95[[i]][1]
  ql[[i]] <- 5 * exp(boc[[i]] + b1c[[i]] * tmb_data$len_mid) / (1 + exp(boc[[i]] + b1c[[i]] * tmb_data$len_mid))
}

plot(tmb_data$len_mid, 1 / 5 * ql[[1]], xlab = "Length (cm)", ylab = "Catchability length pattern", type = "l", lwd = 2)
abline(v = 22, col = "red", lty = 2)
abline(h = .95, col = "red", lty = 2)
abline(v = 17, col = "red", lty = 2)
abline(h = .50, col = "red", lty = 2)

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
  log_std_index1 = log(.5),
  log_std_index2 = log(.5),
  log_std_index3 = log(.5),
  log_std_index4 = log(.5),

  log_std_logq1 = 3,
  log_std_logq2 = 3,
  log_std_logq3 = 3,
  log_std_logq4 = 3,

  # for predicted catch
  log_std_catch = log(.5),

  # the random effects;
  log_Rec_dev = rep(0, Y),
  log_F_dev = matrix(0, nrow = A - 2, ncol = Y, byrow = TRUE),
  pe = matrix(0, nrow = A - 1, ncol = Y - 1, byrow = TRUE),

  logq1 = log(ql1),
  logq2 = log(ql2),
  logq3 = log(ql3),
  logq4 = log(ql4)
)

parameters$log_N0 <- parameters$log_meanR - 0.2 * (1:(A - 1))
# not estimating correlation in rec_devs or any process error at first

parameters_L <- list(
  # for pla_func()
  log_Linf = log(20),
  log_vbk = log(.01),
  log_len_o = log(1 / 10000),
  log_cv_len = log(.01),

  # for computing F
  log_std_log_F = log(.01),
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
  log_std_index1 = log(.1),
  log_std_index2 = log(.1),
  log_std_index3 = log(.1),
  log_std_index4 = log(.1),

  log_std_logq1 = log(.01),
  log_std_logq2 = log(.01),
  log_std_logq3 = log(.01),
  log_std_logq4 = log(.01),

  # for predicted catch
  log_std_catch = log(.01)
)

parameters_U <- list(
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
  log_std_index1 = log(10),
  log_std_index2 = log(10),
  log_std_index3 = log(10),
  log_std_index4 = log(10),

  log_std_logq1 = log(5),
  log_std_logq2 = log(5),
  log_std_logq3 = log(5),
  log_std_logq4 = log(5),

  # for predicted catch
  log_std_catch = log(2)
)

lower <- unlist(parameters_L)
upper <- unlist(parameters_U)

ql1_map <- 1:length(ql[[1]])
ql1_map[ql1_map > qlfix] <- qlfix
ql1_map[ql1_map <= 3] <- NA

ql2_map <- 1:length(ql[[2]])
ql2_map[ql2_map > qlfix] <- qlfix
ql2_map[ql2_map <= 3] <- NA

ql3_map <- 1:length(ql[[3]])
ql3_map[ql3_map > qlfix] <- qlfix
ql3_map[ql3_map <= 5] <- NA

ql4_map <- 1:length(ql[[4]])
ql4_map[ql4_map > qlfix] <- qlfix
ql4_map[ql4_map <= 5] <- NA

pdat <- expand.grid(age = 3:tmb_data$A, year = unique(CLc_vec$Year))
pdat$age[pdat$age > 9] <- 9
F_map <- paste(pdat$year, pdat$age, sep = "_")

log_F_dev_map <- matrix(F_map, nrow = A - 2, ncol = Y, byrow = FALSE)
# log_F_dev_map=matrix(NA,nrow=A-1,ncol=Y,byrow=F)

### fix some parameters in TMB model ###

## first run with everything fixed to get starting values;
map <- list(
  # log_Linf = factor(NA),
  log_len_o = factor(NA),
  log_vbk = factor(NA),
  log_std_log_F = factor(NA),
  logit_F_age = factor(NA),
  logit_F_year = factor(NA),
  log_std_pe = factor(NA),
  logit_log_R = factor(NA),
  pe = factor(matrix(NA, nrow = A - 1, ncol = Y - 1, byrow = TRUE)),

  logq1 = factor(ql1_map),
  logq2 = factor(ql2_map),
  logq3 = factor(ql3_map),
  logq4 = factor(ql4_map),

  log_F_dev = factor(log_F_dev_map)
)

t_L <- parameters_L
{
  # t_L$log_Linf = NULL
  t_L$log_len_o <- NULL
  t_L$log_vbk <- NULL
  t_L$log_std_log_F <- NULL
  t_L$logit_F_age <- NULL
  t_L$logit_F_year <- NULL
  t_L$log_std_pe <- NULL
  t_L$logit_log_R <- NULL
}
lower <- unlist(t_L)

t_U <- parameters_U
{
  # t_L$log_Linf = NULL
  t_U$log_len_o <- NULL
  t_U$log_vbk <- NULL
  t_U$log_std_log_F <- NULL
  t_U$logit_F_age <- NULL
  t_U$logit_F_year <- NULL
  t_U$log_std_pe <- NULL
  t_U$logit_log_R <- NULL
}
upper <- unlist(t_U)

## random effects;
rname <- c(
  "log_Rec_dev",
  "log_F_dev",
  "logq1",
  "logq2",
  "logq3",
  "logq4",
  "pe"
)

obj <- MakeADFun(tmb_data, parameters,
  random = rname, DLL = "fit_mess_m", map = map,
  # inner.control=list(maxit=500,trace=F),random.start = expression(last.par[random]))
  inner.control = list(maxit = 500, trace = FALSE)
)

# check the number of parameters
length(obj$par)
length(lower)
length(upper)

obj$fn(obj$par)

obj$gr(obj$par) # check gradient

opt <- nlminb(obj$par, obj$fn, obj$gr,
  lower = lower, upper = upper,
  control = list(trace = 0, iter.max = 100, eval.max = 10000)
)

opt$objective

cbind(lower, opt$par, upper)

rep <- obj$report()
sd_rep <- sdreport(obj)

###########  Do the Plotting ##################

source("plots.R")

save.image(file = "fit.RData")

## MISC code #############

ind <- abs(CLc.vec$resid) > 4
CLc.vec[ind, ]

ind1 <- CLc.vec$Year == 2019
CLc.vec[ind1, ]


obj$env$last.par.best