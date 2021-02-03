setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab4/BAngler")

# source("xxx",local=F)
rm(list = ls())
library(numDeriv)
library(maxLik)
library(stargazer)

load("tmb.RData")

dat <- list(
  year = tmb.data$year, C = tmb.data$C, log_C = tmb.data$log_C,
  index = tmb.data$index, log_index = tmb.data$log_index,
  iyear = tmb.data$iyear + 1, iq = tmb.data$iq + 1
)
dat$istd <- dat$iq

start.logq <- aggregate(dat$log_index, list(name = dat$iq), mean)$x - log(14)
n.index <- length(unique(dat$iq))

sparm <- unlist(list(
  logq = start.logq,
  logr = log(0.26),
  logK = log(29),
  log_sigma = rep(log(0.2), n.index)
))


spm <- function(r, K, catch) {
  n <- length(catch)
  B <- rep(NA, n)
  B[1] <- 0.6 * K
  for (i in 2:n) {
    B[i] <- B[i - 1] + r * B[i - 1] * (1 - B[i - 1] / K) - catch[i - 1]
  }
  B[B < 0] <- 1e-10
  return(B)
}

spm_pred <- function(x, dat) {
  logq <- x[substr(names(x), 1, 4) == "logq"]
  logr <- x[names(x) == "logr"]
  logK <- x[names(x) == "logK"]

  B <- spm(exp(logr), exp(logK), dat$C)
  log_EIndex <- logq[dat$iq] + log(B[dat$iyear])
  return(log_EIndex)
}

spm_resid <- function(x, dat) {
  log_sigma <- x[substr(names(x), 1, 9) == "log_sigma"]
  resid <- dat$log_index - spm_pred(x, dat)
  resid.std <- resid / exp(log_sigma[dat$istd])
  return(cbind(resid, resid.std))
}

spm_fit <- function(x, dat) {
  log_sigma <- x[substr(names(x), 1, 9) == "log_sigma"]
  resid <- spm_resid(x, dat)[, 1]
  nll <- -sum(dnorm(resid, 0, exp(log_sigma[dat$istd]), TRUE))
  return(nll)
}
spm_logL <- function(x, dat) {
  log_sigma <- x[substr(names(x), 1, 9) == "log_sigma"]
  resid <- spm_resid(x, dat)[, 1]
  nll <- sum(dnorm(resid, 0, exp(log_sigma[dat$istd]), TRUE))
  return(nll)
}

lower <- c(rep(-Inf, n.index), log(0.01), log(1))
upper <- c(rep(Inf, n.index), log(0.4), log(50))

opt.OE <- nlminb(sparm, spm_fit, , lower = lower, upper = upper, dat = dat)
# Hfit = hessian(spm_fit,opt.OE$par,,,dat)

opt.OE <- maxLik(logLik = spm_logL, start = sparm, , , , , dat = dat)
Hfit <- hessian(opt.OE, dat)

summary(opt.OE)

sum(spm_resid(coef(opt.OE), dat)[, 1]**2)


parms <- data.frame(log_est = coef(opt.OE), log_se = stdEr(opt.OE))
parms$est <- exp(parms$log_est)
parms$se <- parms$log_se * parms$est
pname <- row.names(parms)

index_name <- aggregate(tmb.data$name, list(iq = dat$iq), unique)$x
parm_name <- c(paste("q(1000s)", index_name), "r", "K", paste("sigma", index_name))
out <- parms[, 3:4]
out[1:3, ] <- out[1:3, ] * 1000
rownames(out) <- parm_name

stargazer(out,
  digits = 4,
  title = "Parameter Estimates (natural scale)",
  type = "html",
  summary = FALSE,
  out = "out.doc"
)

posv <- rep("topleft", 3)

rep <- list(
  Eindex = exp(spm_pred(coef(opt.OE), dat)),
  std_resid = spm_resid(coef(opt.OE), dat)[, 2]
)

for (i in 1:length(index_name)) {
  u <- index_name[i]
  pos <- posv[i]

  scale <- 1

  ind <- tmb.data$name == u
  iyear <- tmb.data$iyear[ind] + 1
  index <- tmb.data$index[ind] / scale
  Eindex <- rep$Eindex[ind] / scale
  std_resid <- rep$std_resid[ind]

  jpeg(file = paste("OE_fit_", u, ".jpeg", sep = ""), width = 6, height = 6, units = "in", res = 300)

  par(mfcol = c(3, 1), oma = c(1, 1, 0.5, 1), mar = c(3, 4, 2, 0), cex.axis = 1.2, las = 1)

  ylim <- range(index, Eindex)
  xlim <- range(tmb.data$year)

  xi <- tmb.data$year[iyear]
  yi <- rep(NA, length(xi))
  pyi <- yi
  pres <- yi
  ind1 <- xi %in% tmb.data$year
  yi[ind1] <- index
  pyi[ind1] <- Eindex
  pres[ind1] <- std_resid

  plot(xi, yi, type = "b", lwd = 2, pch = 19, xlab = "", ylab = "", ylim = ylim, xlim = xlim)
  lines(xi, pyi, type = "b", lwd = 2)
  mtext(side = 1, line = 2.5, outer = F, "Year")
  mtext(side = 2, line = 3.5, outer = F, "Index", las = 0)
  legend(pos, bty = "n", lwd = 2, pch = c(19, 1), lty = 1, c("Observed", "Predicted"), pt.cex = 1.5)

  plot(xi, pres, type = "b", xlab = "", ylab = "", xlim = xlim)
  abline(h = 0, lty = 2)
  mtext(side = 1, line = 2.5, outer = F, "Year")
  mtext(side = 2, line = 3.5, outer = F, "Standardized residual", las = 0)

  plot(pyi, pres, type = "p", xlab = "", ylab = "")
  abline(h = 0, lty = 2, lwd = 2, col = grey(0.6))
  mtext(side = 1, line = 2.5, outer = F, "Predicted")
  mtext(side = 2, line = 3.5, outer = F, "Standardized residual", las = 0)

  mtext(side = 3, line = -2, outer = T, u, las = 0)

  dev.off()
}

##############  biomass and harvest rel MSY plot ################;

x <- coef(opt.OE)
logq <- x[substr(names(x), 1, 4) == "logq"]
logr <- x[names(x) == "logr"]
logK <- x[names(x) == "logK"]

B <- spm(exp(logr), exp(logK), dat$C)

Bmsy <- exp(logK) / 2
Hmsy <- exp(logr) / 2
biomass <- B
rbiomass <- biomass / Bmsy
exploit <- tmb.data$C / biomass
rexploit <- exploit / Hmsy

ylim <- range(0, rbiomass, rexploit)

par(mar = c(4, 2, 2, 1), mgp = c(2.5, 1, 0))
plot(tmb.data$year, rbiomass, xlab = "Year", ylab = "", type = "n", ylim = ylim)
lines(tmb.data$year, rbiomass, col = "red", type = "b", lwd = 2)
lines(tmb.data$year, rexploit, col = "blue", type = "b", lwd = 2)
mtext(side = 3, line = 0.2, "Relative Biomass and Fishing Mortality")

tmb_data$istd <- tmb_data$iq

library(TMB)

compile("fit_mess.cpp")

dyn.load("fit_mess")
#dyn.unload("fit");

start_q <- aggregate(tmb_data$index, list(iq = tmb_data$iq), mean)

parameters <- list(
  log_r = log(0.2),
  log_K = log(30),
  log_q = log(start_q$x) - log(15),
  log_Po = log(0.6),
  log_Ho = log(0.1),
  log_sd_rw = log(0.2),
  log_sd_log_index = rep(log(0.3), 3),
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
  log_sd_pe = log(0.35),
  logit_ar_pe = log(0.950 / (1 - 0.950)),
  log_sd_logC = Inf
)

lower <- unlist(parameters_L)
upper <- unlist(parameters_U)
## random effects;
rname <- c("log_pe", "log_H_dev")

spm_ar1 <- MakeADFun(tmb_data, parameters,
  random = rname, DLL = "fit_mess",
  inner.control = list(maxit = 100, trace = TRUE)
)

obj$gr(obj$par)

system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr,
    lower = lower, upper = upper,
    control = list(trace = 0, iter.max = 5000, eval.max = 10000)
  )
)

opt$message
obj$gr(opt$par)

rep <- obj$report()

sd.rep <- sdreport(obj)

save.image(file = "fit.RData")