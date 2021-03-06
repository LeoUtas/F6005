setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/project1_mess3")

load("data\\3LN_redfish.RData")
load("tmb_data.RData")

library(TMB)

compile("fit.cpp", PKG_CXXFLAGS = "-Wno-ignored-attributes")

dyn.load("fit")
dyn.unload("fit")
## set FRV length pattern in catchability;
L50c <- 17
L95c <- 22
b1c <- log(0.95 / 0.05) / (L95c - L50c)
boc <- -b1c * L50c
ql <- 5*exp(boc + b1c * tmb_data$len_mid) / (1 + exp(boc + b1c * tmb_data$len_mid))
plot(tmb_data$len_mid, ql, xlab = "length", ylab = "Catchability length pattern", type = "l", lwd = 2)

A <- tmb_data$A
Y <- tmb_data$Y
ns <- length(tmb_data$sf)

parameters <- list(
  log_meanR = 6,
  log_std_log_R = log(1.1),
  log_std_index = rep(log(1.1), ns),
  log_std_catch = log(1.01),
  log_std_logq = rep(log(1.01), ns),
  log_std_pe = log(1.01),
  logit_log_R = 2,
  log_Linf = log(35),
  log_vbk = log(0.35),
  log_len_o = log(.5),
  log_cv_len = log(.2),
  log_std_log_F = 1,
  logit_F_age = log(1.01),
  logit_F_year = 1,
  log_F_main = matrix(log(.05), nrow = A - 2, ncol = Y),
  log_N0 = rep(4, A - 1),
  # the random effects;
  log_Rec_dev = rep(0, Y),
  log_F_dev = matrix(0, nrow = A - 2, ncol = Y),
  logq = log(matrix(ql, nrow = ns, ncol = length(ql), byrow = T)),
  pe = matrix(0, nrow = A - 1, ncol = Y - 1, byrow = T)
)

parameters$log_N0 <- parameters$log_meanR - 0.2 * (1:(A - 1))
parameters$logq[3, ] <- parameters$logq[3, ] - log(10000)

# not estimating correlation in rec_devs or any process error at first
parameters_L <- list(
  log_meanR = 1,
  log_std_log_R = log(1.01),
  log_std_index = rep(log(1.01), ns),
  log_std_catch = log(.1),
  log_std_logq = rep(-Inf, ns),
  log_std_pe = log(.1),
  logit_log_R = -10,
  log_Linf = log(20),
  log_vbk = log(.01),
  log_len_o = log(.01),
  log_cv_len = log(.01),
  log_std_log_F = log(.01),
  logit_F_age = log(.1),
  logit_F_year = log(.1),
  log_F_main = matrix(-Inf, nrow = A - 2, ncol = Y),
  log_N0 = rep(-Inf, A - 1)
)

parameters_U <- list(
  log_meanR = 10,
  log_std_log_R = log(10),
  log_std_index = rep(log(10), ns),
  log_std_catch = log(10),
  log_std_logq = rep(Inf, ns),
  log_std_pe = log(10),
  logit_log_R = 10,
  log_Linf = log(100),
  log_vbk = log(2),
  log_len_o = log(2),
  log_cv_len = log(5),
  log_std_log_F = log(10),
  logit_F_age = log(10),
  logit_F_year = log(10),
  log_F_main = matrix(Inf, nrow = A - 2, ncol = Y),
  log_N0 = rep(Inf, A - 1)
)

lower <- unlist(parameters_L)
upper <- unlist(parameters_U)

ql_map <- 1:length(ql)
ql_map[ql_map > 45] <- 45
ql_map[ql_map <= 3] <- NA
m1 <- matrix(1:ns, nrow = ns, ncol = length(ql), byrow = F)
m2 <- matrix(ql_map, nrow = ns, ncol = length(ql), byrow = T)
qmap_by_survey <- sprintf("%d,%d", m1, m2)
dim(qmap_by_survey) <- dim(m1)
qmap_by_survey[, 1:3] <- NA

pdat <- expand.grid(age = 3:tmb_data$A, year = unique(CLc_vec$Year))
pdat$age[pdat$age > 9] <- 9
F_map <- paste(pdat$year, pdat$age, sep = "_")

log_F_dev_map <- matrix(F_map, nrow = A - 2, ncol = Y, byrow = F)

## code to map mean F's by age and for two time periods
year <- sort(unique(CL_vec$Year))
agec <- c(2, 2, 2, 2, 2, 2, 2, 2)
ind <- year %in% 1995:2009

F_map_main <- matrix(paste("P1_", agec, sep = ""), ncol = Y, nrow = A - 2, byrow = F)
F_map_main[, ind] <- matrix(paste("P2_", agec, sep = ""), ncol = length(year[ind]), nrow = A - 2, byrow = F)

## first run with logit_F_age and logit_F_year fixed to get starting values;

map <- list(
  #log_len_o = factor(NA),
  #log_cv_len = factor(NA),
  # log_Linf = factor(NA),
  # log_vbk = factor(NA),
  log_F_main = factor(F_map_main),
  #logit_F_age = factor(NA),
  #logit_F_year = factor(NA),
  #log_std_pe = factor(NA),
  #logit_log_R = factor(NA),
  #pe = factor(matrix(NA, nrow = A - 1, ncol = Y - 1, byrow = T)),
  logq = factor(qmap_by_survey),
  log_F_dev = factor(log_F_dev_map)
)

t_L <- parameters_L
{
  #t_L$log_len_o <- NULL
  #t_L$log_cv_len <- NULL
  # t_L$log_Linf <- NULL
  # t_L$log_vbk <- NULL
  t_L$log_F_main <- rep(-Inf, length(levels(map$log_F_main)))
  #t_L$logit_F_age <- NULL
  #t_L$logit_F_year <- NULL
  #t_L$log_std_pe <- NULL
  #t_L$logit_log_R <- NULL
}
lower <- unlist(t_L)

t_U <- parameters_U
{
  #t_U$log_len_o <- NULL
  #t_U$log_cv_len <- NULL
  # t_U$log_Linf <- NULL
  # t_U$log_vbk <- NULL
  t_U$log_F_main <- rep(Inf, length(levels(map$log_F_main)))
  #t_U$logit_F_age <- NULL
  #t_U$logit_F_year <- NULL
  #t_U$log_std_pe <- NULL
  #t_U$logit_log_R <- NULL
}
upper <- unlist(t_U)

## random effects;
rname <- c("log_Rec_dev", "log_F_dev", "logq", "pe")

obj <- MakeADFun(tmb_data, parameters,
  random = rname, DLL = "fit", map = map,
  # inner.control=list(maxit=500,trace=F),random.start = expression(last.par[random]))
  inner.control = list(maxit = 500, trace = T)
)

length(obj$par)
length(lower)
length(upper)

obj$fn(obj$par)

obj$gr(obj$par)

# "FallRV" is=0 "SprgRV" is=1 "Spsh3L" is=2 "Spsh3N" is=3  this is the connection between is and sname

opt <- nlminb(obj$par, obj$fn, obj$gr,
  lower = lower, upper = upper,
  control = list(trace = 0, iter.max = 500, eval.max = 10000)
)

rep <- obj$report()

sd_rep <- sdreport(obj)

opt$objective

cbind(lower, opt$par, upper)

save.image(file = "fit28.RData")

###########  Do the Plotting ##################

source("plots.R")




ind <- abs(CLc_vec$resid) > 10
CLc_vec[ind, ]

ind1 <- CLc_vec$Year == 1993
CLc_vec[ind1, ]