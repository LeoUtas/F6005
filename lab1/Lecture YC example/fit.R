setwd("C:\\home\\CADIGAN\\GradProgram\\2018\\F6005\\Lecture1\\YC")

require(lattice)
library(latticeExtra)
require(reshape2)
library(lme4)
library(cAIC4)
library(xtable)
library(ggplot2)
library(stargazer)
library(TMB)

sdat <- read.table(file = "F3LNO.dat", header = TRUE)
sdat$fYC <- factor(sdat$YC)
sdat$fAge <- factor(sdat$Age)
sdat$fYear <- factor(sdat$Year)

my.padding <- list(
  layout.heights = list(
    top.padding = 0,
    main.key.padding = 0,
    key.axis.padding = 0,
    axis.xlab.padding = 0,
    xlab.key.padding = 0,
    key.sub.padding = 0
  ),
  layout.widths = list(
    left.padding = 1,
    key.ylab.padding = 0,
    ylab.axis.padding = 1,
    axis.key.padding = 0,
    right.padding = 0
  )
)

colv <- c("blue", "cyan", "aquamarine4", "green", "darkgoldenrod4", "red")


mixfit <- lmer(log.index ~ fYC + fAge - 1 + (1 | fYear), REML = FALSE, data = sdat)

summary(mixfit)



################ TMB ############################

compile("YC_YE.cpp")

dyn.load("YC_YE")
# dyn.unload("YC_YE");

tmb.data <- list(
  y = sdat$log.index,
  iAE = as.numeric(sdat$fAge) - 1,
  iCE = as.numeric(sdat$fYC) - 1,
  iYE = as.numeric(sdat$fYear) - 1
)
npar_AE <- length(unique(tmb.data$iAE))
npar_CE <- length(unique(tmb.data$iCE))
npar_YE <- length(unique(tmb.data$iYE))

parameters <- list(
  par_AE = rep(0, npar_AE - 1),
  par_CE = rep(0, npar_CE),
  log_std_index = log(0.3),
  log_std_YE = log(0.1),
  par_YE = rep(0, npar_YE)
)

parameters.L <- list(
  par_age = rep(-Inf, npar_AE - 1),
  par_CE = rep(-Inf, npar_CE),
  log_std_index = -5,
  log_std_YE = -5
)

parameters.U <- list(
  par_age = rep(Inf, npar_AE - 1),
  par_CE = rep(Inf, npar_CE),
  log_std_index = 5,
  log_std_YE = 5
)
lower <- unlist(parameters.L)
upper <- unlist(parameters.U)
obj <- MakeADFun(tmb.data, parameters,
  random = c("par_YE"),
  DLL = "YC_YE", inner.control = list(maxit = 300, trace = F)
)

## check you have the length of bounds right
length(lower)
length(upper)
length(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr,
  lower = lower, upper = upper,
  control = list(trace = 0, iter.max = 2000, eval.max = 10000)
)

rep <- obj$report()
sd.rep <- sdreport(obj)

cbind(ranef(mixfit)$fYear, rep$par_YE)

cbind(sd.rep$par.random, sqrt(sd.rep$diag.cov.random))

cbind(sd.rep$par.fixed, sqrt(diag(sd.rep$cov.fixed)))

uyear <- sort(unique(sdat$Year))
YE <- sd.rep$par.random
sd <- sqrt(sd.rep$diag.cov.random)
YE.L <- sd.rep$par.random + qnorm(0.025) * sd
YE.U <- sd.rep$par.random + qnorm(0.975) * sd
ylim <- range(YE, YE.L, YE.U)

par(mar = c(3, 3, 0.5, 0.5), las = 1, mgp = c(1, 0.5, 0))

plot(uyear, YE, xlab = "", ylab = "", type = "l", lwd = 2, ylim = ylim)

mtext(side = 1, line = 2, "Year", las = 0)
mtext(side = 2, line = 2, "Year Effect", las = 1)
abline(h = 0, lty = 1, lwd = 2)
lines(uyear, YE.L, lwd = 2, col = "grey")
lines(uyear, YE.U, lwd = 2, col = "grey")

ci <- confint(mixfit, parm = names(fixef(mixfit)), method = "Wald")
out <- data.frame(log_est = fixef(mixfit), ci = ci)
mix.rec.dev <- subset(out, substr(rownames(out), 1, 3) == "fYC")
rnames <- rownames(mix.rec.dev)
mix.rec.dev$YC <- as.numeric(substring(rnames, 4, 7))
mix.rec.dev$est <- exp(mix.rec.dev$log_est)
mean.est <- mean(mix.rec.dev$est)
mix.rec.dev$est <- mix.rec.dev$est / mean.est
mix.rec.dev$L95 <- exp(mix.rec.dev$ci.2.5..) / mean.est
mix.rec.dev$U95 <- exp(mix.rec.dev$ci.97.5..) / mean.est

ind <- names(sd.rep$value) == "CE_dev"
mix.rec.dev$TMB.est <- exp(sd.rep$value[ind])
mix.rec.dev$TMB.L95 <- exp(sd.rep$value[ind] + qnorm(0.025) * sd.rep$sd[ind])
mix.rec.dev$TMB.U95 <- exp(sd.rep$value[ind] + qnorm(0.975) * sd.rep$sd[ind])

ylim <- range(mix.rec.dev[, 5:10], na.rm = T)
# ylim[2]=8

gname <- "YC_mixed.jpeg"
jpeg(file = gname, width = 3, height = 3, units = "in", res = 300)
par(mar = c(3, 3.3, 0.2, 1), mgp = c(2, 0.7, 0))
plot(est ~ YC, data = mix.rec.dev, ylab = "", xlab = "", las = 1, type = "l", lwd = 2, ylim = ylim)
mtext(side = 1, line = 1.7, "Cohort")
mtext(side = 2, line = 2.3, "Relative Cohort Strength")

sx <- mix.rec.dev$YC
low <- mix.rec.dev$L95
high <- mix.rec.dev$U95
polygon(c(sx, rev(sx)), c(low, rev(high)), col = rgb(0, 0, 1, alpha = 0.2, maxColorValue = 1), border = NA)
sx <- mix.rec.dev$YC
low <- mix.rec.dev$TMB.L95
high <- mix.rec.dev$TMB.U95
polygon(c(sx, rev(sx)), c(low, rev(high)), col = rgb(1, 0, 0, alpha = 0.2, maxColorValue = 1), border = NA)

vy <- unique(sdat$YC)
vy <- vy[2:(length(vy) - 1)]
#  abline(v=vy,lty=2,col="darkgoldenrod1")
abline(h = 1, lty = 2)
lines(mix.rec.dev$YC, mix.rec.dev$est, lwd = 2, col = "blue")

legend("topleft", col = c("red", "blue"), lty = 1, lwd = 2, legend = c("TMB", "lmer"), bty = "n")

dev.off()