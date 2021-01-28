setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab3")

library(nlme)
library(stargazer)

fname <- c("SR.dat")
sr_data <- read.table(fname, header = TRUE)
sr_data$logrec <- log(sr_data$rec)
sr_data$logssb <- log(sr_data$ssb)

init_Rmax <- 150
init_S50 <- 50
init_alpha <- init_Rmax
init_beta <- init_S50


BH_fit <- nls(logrec ~ log(alpha) + logssb - log(beta + ssb),
  algorithm = "port", lower = c(0, 0), data = sr_data,
  start = list(beta = init_beta, alpha = init_alpha), trace = TRUE
)

BH.arfit <- gnls(logrec ~ log(alpha) + logssb - log(beta + ssb),
  data = sr_data,
  start = list(beta = init_beta, alpha = init_alpha),
  correlation = corAR1(form = ~year)
)

library(TMB)

## Autocorrelated errors;
compile("fit.cpp")
dyn.load("fit")
# dyn.unload("fit")

tmb_data <- list(
  ssb = sr_data$ssb,
  rec = sr_data$rec,
  lssb = log(sr_data$ssb),
  lrec = log(sr_data$rec)
)

parameters <- list(
  lalpha = log(150),
  lbeta = log(50),
  lsd_lrec_me = log(0.1),
  logit_ar_lrec_me = 0
)

obj <- MakeADFun(tmb_data, parameters,
  DLL = "fit",
  inner.control = list(maxit = 500, trace = TRUE)
)

obj$gr(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

obj$gr(opt$par)

rep <- obj$report()
sd_rep <- sdreport(obj)

tab <- cbind(sd_rep$value, sd_rep$sd)
colnames(tab) <- c("Est", "Sd.err")

stargazer(tab,
  type = "html", out = "out.doc", align = TRUE,
  title = "TMB BH AR errors", single.row = TRUE,
  summary = FALSE
)


##########  RW in log alpha ;

n <- nrow(sr_data)

compile("fit_saoRW.cpp")

dyn.load("fit_saoRW")
# dyn.unload("fit_saoRW")

tmb_data <- list(
  ssb = sr_data$ssb,
  rec = sr_data$rec,
  lssb = log(sr_data$ssb),
  lrec = log(sr_data$rec)
)

parameters <- list(
  lsao = log(3),
  lbeta = log(50),
  lsd_lrec_me = log(0.1),
  lsd_lsao_dev = log(0.1),
  lsao_dev = rep(0, n)
)


# map <- list(lsd_lrec_me = factor(NA))

obj <- MakeADFun(tmb_data, parameters,
  DLL = "fit_saoRW",
  random = c("lsao_dev"), inner.control = list(maxit = 500, trace = TRUE)
)

obj$gr(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

obj$gr(opt$par)

rep <- obj$report()
sd_rep <- sdreport(obj)

sname <- names(sd_rep$value)

ind <- sname %in% c("lsao", "beta", "sd_lrec_me", "sd_lsao_dev")

tab <- cbind(sd_rep$value[ind], sd_rep$sd[ind])
colnames(tab) <- c("Est", "Sd.err")

stargazer(tab,
  type = "html", out = "out_saoRW.doc", align = TRUE,
  title = "TMB BH Sao RW", single.row = TRUE,
  summary = FALSE
)

ind <- sname %in% c("lsao_ts")
pdat <- data.frame(
  year = sr_data$year, lsao = sd_rep$value[ind],
  L95 = sd_rep$value[ind] - qnorm(0.975) * sd_rep$sd[ind],
  U95 = sd_rep$value[ind] + qnorm(0.975) * sd_rep$sd[ind]
)

re <- tab[1, 1] + sd_rep$par.random
sdre <- sqrt(sd_rep$diag.cov.random)
pdat$L95.1 <- re - qnorm(0.975) * sdre
pdat$U95.1 <- re + qnorm(0.975) * sdre

gname <- "saoRW.jpeg"
jpeg(file = gname, width = 3, height = 3, units = "in", res = 300)

ylim <- range(pdat$L95, pdat$U95)
par(mar = c(3, 3.3, 0.2, 1), mgp = c(2, 0.7, 0))
plot(lsao ~ year, data = pdat, ylab = "", xlab = "", las = 1, type = "l", lwd = 2, ylim = ylim)
mtext(side = 1, line = 1.7, "Year")
mtext(side = 2, line = 2.3, "ln(Sao)")

sx <- pdat$year
low <- pdat$L95
high <- pdat$U95
polygon(c(sx, rev(sx)), c(low, rev(high)),
  #  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = "yellow", border = NA
)

low <- pdat$L95.1
high <- pdat$U95.1
polygon(c(sx, rev(sx)), c(low, rev(high)),
  #  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = "goldenrod", border = NA
)

lines(pdat$year, pdat$lsao, lty = 1, lwd = 2)

dev.off()

#################  profile rec me ####################;


map <- list(lsd_lrec_me = factor(NA))

psd <- seq(0.05, 0.45, by = 0.01)
nsd <- length(psd)
pfit <- rep(NA, nsd)

for (i in 1:nsd) {
  parameters$lsd_lrec_me <- log(psd[i])

  obj <- MakeADFun(tmb_data, parameters,
    DLL = "fit_saoRW", map = map,
    random = c("lsao_dev"), inner.control = list(maxit = 500, trace = FALSE)
  )

  optp <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 0))
  pfit[i] <- optp$objective
}

pfit <- pfit - opt$objective

gname <- "saoRW_profile.jpeg"
jpeg(file = gname, width = 3, height = 3, units = "in", res = 300)

par(mar = c(3, 3.3, 0.2, 1), mgp = c(2, 0.7, 0))
plot(pfit ~ psd, data = pdat, ylab = "", xlab = "", las = 1, type = "l", lwd = 2)
mtext(side = 1, line = 1.7, substitute(sigma[epsilon]), las = 0)
mtext(side = 2, line = 2.3, substitute(paste(Delta, " nll")), las = 0)
abline(h = qchisq(0.95, 1) / 2, lty = 2)

dev.off()


########## RK  RW in log alpha ;

compile("fit_RK_RW.cpp")

dyn.load("fit_RK_RW")
# dyn.unload("fit_RK_RW")

tmb_data <- list(
  ssb = sr_data$ssb,
  rec = sr_data$rec,
  lssb = log(sr_data$ssb),
  lrec = log(sr_data$rec)
)

parameters <- list(
  lalpha = log(3),
  lbeta = log(0.01),
  lsd_lrec_me = log(0.05),
  lsd_lalpha_dev = log(0.5),
  lalpha_dev = rep(0, n)
)

obj <- MakeADFun(tmb_data, parameters,
  DLL = "fit_RK_RW",
  random = c("lalpha_dev"), inner.control = list(maxit = 500, trace = TRUE)
)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

obj$gr(opt$par)

rep <- obj$report()
sd.rep <- sdreport(obj)

sname <- names(sd.rep$value)

ind <- sname %in% c("alpha_main", "beta", "sd_lrec_me", "sd_lalpha_dev")

tab <- cbind(sd.rep$value[ind], sd.rep$sd[ind])
colnames(tab) <- c("Est", "Sd.err")

stargazer(tab,
  type = "html", out = "out_RKsaoRW.doc", align = TRUE,
  title = "TMB RK Sao RW", single.row = TRUE,
  summary = FALSE
)

ind <- sname %in% c("lsao_ts")
pdat <- data.frame(
  year = sr_data$year, lsao = sd.rep$value[ind],
  L95 = sd.rep$value[ind] - qnorm(0.975) * sd.rep$sd[ind],
  U95 = sd.rep$value[ind] + qnorm(0.975) * sd.rep$sd[ind]
)

re <- opt$par[1] + sd.rep$par.random
sdre <- sqrt(sd.rep$diag.cov.random)
pdat$L95.1 <- re - qnorm(0.975) * sdre
pdat$U95.1 <- re + qnorm(0.975) * sdre

gname <- "RKsaoRW.jpeg"
jpeg(file = gname, width = 3, height = 3, units = "in", res = 300)

ylim <- range(pdat$L95, pdat$U95)
par(mar = c(3, 3.3, 0.2, 1), mgp = c(2, 0.7, 0))
plot(lsao ~ year, data = pdat, ylab = "", xlab = "", las = 1, type = "l", lwd = 2, ylim = ylim)
mtext(side = 1, line = 1.7, "Year")
mtext(side = 2, line = 2.3, "ln(Sao)")

sx <- pdat$year
low <- pdat$L95
high <- pdat$U95
polygon(c(sx, rev(sx)), c(low, rev(high)),
  #  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = "yellow", border = NA
)

lines(pdat$year, pdat$lsao, lty = 1, lwd = 2)

dev.off()

#################  profile rec me ####################;


map <- list(lsd_lrec_me = factor(NA))

psd <- seq(0.01, 0.6, by = 0.01)
nsd <- length(psd)
pfit <- rep(NA, nsd)

for (i in 1:nsd) {
  parameters$lsd_lrec_me <- log(psd[i])

  obj <- MakeADFun(tmb_data, parameters,
    DLL = "fit_RK_RW", map = map,
    random = c("lalpha_dev"), inner.control = list(maxit = 500, trace = FALSE)
  )

  optp <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 0))
  pfit[i] <- optp$objective
}

pfit <- pfit - opt$objective

gname <- "RKsaoRW_profile.jpeg"
jpeg(file = gname, width = 3, height = 3, units = "in", res = 300)

par(mar = c(3, 3.3, 0.2, 1), mgp = c(2, 0.7, 0))
plot(pfit ~ psd, data = pdat, ylab = "", xlab = "", las = 1, type = "l", lwd = 2, ylim = c(0, 2))
mtext(side = 1, line = 1.7, substitute(sigma[epsilon]), las = 0)
mtext(side = 2, line = 2.3, substitute(paste(Delta, " nll")), las = 0)
abline(h = qchisq(0.95, 1) / 2, lty = 2)

dev.off()


## fixed rec ME CV;
map <- list(lsd_lrec_me = factor(NA))

obj <- MakeADFun(tmb_data, parameters,
  DLL = "fit_RK_RW", map = map,
  random = c("lalpha_dev"), inner.control = list(maxit = 500, trace = TRUE)
)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

obj$gr(opt$par)

rep <- obj$report()
sd_rep <- sdreport(obj)

sname <- names(sd_rep$value)

ind <- sname %in% c("alpha_main", "beta", "sd_lrec_me", "sd_lalpha_dev")

tab <- cbind(sd_rep$value[ind], sd_rep$sd[ind])
colnames(tab) <- c("Est", "Sd.err")

stargazer(tab,
  type = "html", out = "out_RKsaoRW_recMEfixed.doc", align = TRUE,
  title = "TMB RK Sao RW", single.row = TRUE,
  summary = FALSE
)

ind <- sname %in% c("lsao_ts")
pdat <- data.frame(
  year = sr_data$year, lsao = sd_rep$value[ind],
  L95 = sd_rep$value[ind] - qnorm(0.975) * sd_rep$sd[ind],
  U95 = sd_rep$value[ind] + qnorm(0.975) * sd_rep$sd[ind]
)

gname <- "RKsaoRW_recMEfixed.jpeg"
jpeg(file = gname, width = 3, height = 3, units = "in", res = 300)

ylim <- range(pdat$L95, pdat$U95)
par(mar = c(3, 3.3, 0.2, 1), mgp = c(2, 0.7, 0))
plot(lsao ~ year, data = pdat, ylab = "", xlab = "", las = 1, type = "l", lwd = 2, ylim = ylim)
mtext(side = 1, line = 1.7, "Year")
mtext(side = 2, line = 2.3, "ln(Sao)")

sx <- pdat$year
low <- pdat$L95
high <- pdat$U95
polygon(c(sx, rev(sx)), c(low, rev(high)),
  #  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = "yellow", border = NA
)

lines(pdat$year, pdat$lsao, lty = 1, lwd = 2)

dev.off()

###############  Added AR(1) alpha process Dec 10 2018;



########## RK  RW in log alpha ;

compile("fit_RK_AR1.cpp")

dyn.load("fit_RK_AR1")
# dyn.unload("fit_RK_AR1")

tmb_data <- list(
  ssb = sr_data$ssb,
  rec = sr_data$rec,
  lssb = log(sr_data$ssb),
  lrec = log(sr_data$rec)
)

parameters <- list(
  lalpha = log(3),
  lbeta = log(0.2),
  lsd_lrec_me = log(0.3),
  lsd_lalpha_dev = log(0.1),
  logit_phi = 0,
  lalpha_dev = rep(0, n)
)

map <- list(lsd_lrec_me = factor(NA))

obj <- MakeADFun(tmb_data, parameters,
  DLL = "fit_RK_AR1", map = map,
  random = c("lalpha_dev"), inner.control = list(maxit = 500, trace = TRUE)
)

obj$gr(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

obj$gr(opt$par)

rep <- obj$report()

plot(tmb_data$rec)
lines(rep$mu)


sd_rep <- sdreport(obj)

sname <- names(sd_rep$value)

ind <- sname %in% c("alpha_main", "beta", "sd_lrec_me", "sd_lalpha_dev")

tab <- cbind(sd_rep$value[ind], sd_rep$sd[ind])
colnames(tab) <- c("Est", "Sd.err")

stargazer(tab,
  type = "html", out = "out_RKsaoRW.doc", align = TRUE,
  title = "TMB RK Sao RW", single.row = TRUE,
  summary = FALSE
)

ind <- sname %in% c("lsao_ts")
pdat <- data.frame(
  year = sr_data$year, lsao = sd_rep$value[ind],
  L95 = sd_rep$value[ind] - qnorm(0.975) * sd_rep$sd[ind],
  U95 = sd_rep$value[ind] + qnorm(0.975) * sd_rep$sd[ind]
)

re <- opt$par[1] + sd_rep$par.random
sdre <- sqrt(sd_rep$diag.cov.random)
pdat$L95.1 <- re - qnorm(0.975) * sdre
pdat$U95.1 <- re + qnorm(0.975) * sdre

gname <- "RKsaoAR.jpeg"
jpeg(file = gname, width = 3, height = 3, units = "in", res = 300)

ylim <- range(pdat$L95, pdat$U95)
par(mar = c(3, 3.3, 0.2, 1), mgp = c(2, 0.7, 0))
plot(lsao ~ year, data = pdat, ylab = "", xlab = "", las = 1, type = "l", lwd = 2, ylim = ylim)
mtext(side = 1, line = 1.7, "Year")
mtext(side = 2, line = 2.3, "ln(Sao)")

sx <- pdat$year
low <- pdat$L95
high <- pdat$U95
polygon(c(sx, rev(sx)), c(low, rev(high)),
  #  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = "yellow", border = NA
)

lines(pdat$year, pdat$lsao, lty = 1, lwd = 2)

dev.off()