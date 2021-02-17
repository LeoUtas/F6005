rm(list = ls())

setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab5_3NO_cod_SURBA")

load("tmb_data.RData")

library(TMB)

compile("fit_mess.cpp", flags = "-Wno-ignored-attributes")

dyn.load("fit_mess")
dyn.unload("fit_mess")

A <- nrow(tmb_data$mat)
Y <- ncol(tmb_data$mat)

parameters <- list(
  log_mean_Rec = 0,
  log_sd_log_Rec = -1,
  log_sd_index = -1,
  log_sd_log_f = -1,
  log_sd_log_s = c(-1, -2),
  log_sd_pe = log(0.1),
  logit_log_Rec = 0,
  log_No = rep(0, A - 1),
  log_Rec_dev = rep(0, Y),
  log_sp = rep(0, A - 2), # age 6 fixed at one so not estimated;
  log_f = rep(-0.5, Y),
  pe = matrix(0, nrow = A, ncol = Y, byrow = TRUE)
)

parameters_L <- list(
  log_mean_Rec = -1,
  log_sd_log_Rec = -5,
  log_sd_index = -5,
  log_sd_log_f = -5,
  log_sd_log_s = c(-5, -5),
  # log_std_log_s = c(-5),
  log_sd_pe = -5,
  logit_log_Rec = -10,
  log_No = rep(-5, A - 1)
)

parameters_U <- list(
  log_mean_Rec = 15,
  log_sd_log_Rec = 0,
  log_sd_index = 0,
  log_sd_log_f = 4,
  log_sd_log_s = c(Inf, Inf),
  # log_std_log_s = c(Inf),
  log_sd_pe = 0,
  logit_log_Rec = 10,
  log_No = rep(25, A - 1)
)

lower <- unlist(parameters_L)
upper <- unlist(parameters_U)
# map=list(log_std_log_s=factor(c('comb','comb')))

## random effects;
rname <- c("log_Rec_dev", "log_sp", "log_f", "pe")

obj <- MakeADFun(tmb_data, parameters,
  random = rname, DLL = "fit_mess",
  inner.control = list(maxit = 100, trace = TRUE)
)

length(lower)
length(upper)
length(obj$par)

obj$gr(obj$par)


opt <- nlminb(obj$par, obj$fn, obj$gr,
  lower = lower, upper = upper,
  control = list(trace = 0, iter.max = 100, eval.max = 10000)
)

opt$message
obj$gr(opt$par)

report <- obj$report()

sd_report <- sdreport(obj)

save.image(file = "fit_mess.RData")
