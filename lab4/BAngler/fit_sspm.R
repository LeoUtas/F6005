setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab4/BAngler")

library(TMB)

compile("sspm_ar1_sd_logC_fixed.cpp", flags = "-Wno-ignored-attributes")

# sspm_ar1_sd_logC_est.cpp
# sspm_ar1_sd_logC_fixed.cpp
# sspm_iid_sd_logC_est.cpp
# sspm_iid_sd_logC_fixed.cpp
# sspm_ar1_sd_logC_fixed_sd_logindex_vec



dyn.load("sspm_ar1_sd_logC_fixed")
# dyn.unload("sspm_ar1_sd_logC_est")

start_q <- aggregate(tmb_data$index, list(iq = tmb_data$iq), mean)

parameters <- list(
  log_r = log(0.2),
  log_K = log(30),
  log_q = log(start_q$x) - log(15),
  log_Po = log(0.6),
  log_Ho = log(0.1),
  log_sd_rw = log(0.2),
  log_sd_log_index = rep(log(0.3), 3),
  # log_sd_log_index = log(0.3),
  log_sd_pe = log(0.1),
  logit_ar_pe = log(0.50 / (1 - 0.50)),
  log_sd_logC = log(0.1),
  log_pe = rep(0, length(tmb.data$C)),
  log_H_dev = rep(0, length(tmb.data$C) - 1)
)

parameters_L <- list(
  log_r = log(0.1),
  log_K = log(10),
  log_q = rep(-Inf, length(unique(tmb.data$iq))),
  log_Po = log(0.0001),
  log_Ho = log(0.0001),
  log_sd_rw = log(0.01),
  log_sd_log_index = rep(log(0.01), length(unique(tmb.data$iq))),
  # log_sd_log_index = log(0.01),
  log_sd_pe = log(0.001),
  logit_ar_pe = log(0.01 / (1 - 0.01)),
  log_sd_logC = log(0.02)
)

parameters_U <- list(
  log_r = log(0.6),
  log_K = log(100),
  log_q = rep(Inf, length(unique(tmb.data$iq))),
  log_Po = log(10),
  log_Ho = log(1),
  log_sd_rw = log(2),
  log_sd_log_index = rep(log(1), length(unique(tmb.data$iq))),
  # log_sd_log_index = log(1),
  log_sd_pe = log(0.35),
  logit_ar_pe = log(0.950 / (1 - 0.950)),
  log_sd_logC = Inf
)

lower <- unlist(parameters_L)
upper <- unlist(parameters_U)
## random effects;
rname <- c("log_pe", "log_H_dev")

sspm <- MakeADFun(tmb_data, parameters,
  random = rname, DLL = "sspm_ar1_sd_logC_fixed",
  inner.control = list(maxit = 100, trace = TRUE)
)

sspm$gr(sspm$par)

sspm_fit <- nlminb(sspm$par, sspm$fn, sspm$gr,
  lower = lower, upper = upper,
  control = list(trace = 0, iter.max = 5000, eval.max = 10000)
)

sspm_fit$message

sspm$gr(sspm_fit$par)

report <- sspm$report()

sd_report <- sdreport(sspm)

save.image(file = "sspm_fit.RData")

sspm_fit

system.time(
  sspm_fit <- nlminb(sspm$par, sspm$fn, sspm$gr,
    lower = lower, upper = upper,
    control = list(trace = 0, iter.max = 5000, eval.max = 10000)
  )
)