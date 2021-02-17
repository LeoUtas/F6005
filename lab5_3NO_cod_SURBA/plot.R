setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab5_3NO_cod_SURBA")

library(lattice)
library(latticeExtra)
library(reshape2)
library(xtable)
library(ggplot2)
library(psych)
library(stargazer)

load(file = "fit_mess.RData")

my_padding <- list(
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

colv <- colors()[seq(26, by = 5, length = 13)]

age_id <- sort(unique(data2$age))
year <- c(1990:2016)
q_age <- exp(tmb_data$log_q)

gname <- "q.jpeg"
jpeg(file = gname, width = 8, height = 8, units = "in", res = 600)
par(mar = c(3, 3, 0.2, 0.2))
plot(age_id, q_age, type = "l", lwd = 2, xlab = "", ylab = "")
mtext(side = 2, outer = FALSE, line = 2, "catchability q")
mtext(side = 1, outer = FALSE, line = 2, "Age")
dev.off()

gname <- "m.jpeg"
jpeg(file = gname, width = 8, height = 8, units = "in", res = 600)
par(mar = c(3, 3, 0.2, 0.2))
plot(age_id, tmb_data$m, type = "l", lwd = 2, xlab = "", ylab = "", ylim = c(0, max(tmb_data$m)))
abline(h = 0.2, lty = 1, lwd = 2, col = "grey")
mtext(side = 2, outer = FALSE, line = 2, "Natural Mortality Rate (m)")
mtext(side = 1, outer = FALSE, line = 2, "Age")
dev.off()

## F comparison plot;

assess_F <- read.table("Fbar.dat",
  header = FALSE, col.names = c("year", "F6_9", "F4_6"), sep = "\t"
)

value_names <- names(sd_report$value)

ind <- value_names == "log_aveF_46"
est_46 <- exp(sd_report$value[ind])
sdi <- sd_report$sd[ind]
low_46 <- exp(sd_report$value[ind] - qnorm(0.975) * sdi)
high_46 <- exp(sd_report$value[ind] + qnorm(0.975) * sdi)
ind <- value_names == "log_aveF_69"
est_69 <- exp(sd_report$value[ind])
sdi <- sd_report$sd[ind]
low_69 <- exp(sd_report$value[ind] - qnorm(0.975) * sdi)
high_69 <- exp(sd_report$value[ind] + qnorm(0.975) * sdi)
ylim <- range(low_46, low_69, high_69, high_46)
ylim[2] <- 2.5

gname <- "aveF.jpeg"
jpeg(file = gname, width = 4, height = 4, units = "in", res = 300)
par(mfcol = c(2, 1), oma = c(0, 1, 0.5, 0.5), mar = c(3, 2.5, 0, 0), las = 1)
plot(year, est_46, type = "n", lwd = 2, xlab = "", ylab = "", xlim = c(1959, 2016), ylim = ylim)

sx <- year
polygon(c(sx, rev(sx)), c(low_46, rev(high_46)), col = "grey", border = NA)
box(lty = 1)
lines(year, est_46, type = "l", lty = 1, lwd = 2)
lines(assess_F$year, assess_F$F4_6, type = "l", lty = 3, lwd = 2)
mtext(side = 2, outer = FALSE, line = 2.5, "Ave F (Ages 4-6)", las = 0)
legend("topleft", lwd = c(2, 2), lty = c(1, 3), bty = "n", legend = c("SURBA++", "Assessment"), cex = 0.75)


plot(year, est_69, type = "n", lwd = 2, xlab = "", ylab = "", xlim = c(1959, 2016), ylim = ylim)
polygon(c(sx, rev(sx)), c(low_69, rev(high_69)), col = "grey", border = NA)
box(lty = 1)
lines(year, est_69, type = "l", lty = 1, lwd = 2)
lines(assess_F$year, assess_F$F6_9, type = "l", lty = 3, lwd = 2)
mtext(side = 2, outer = FALSE, line = 2.5, "Ave F (Ages 6-9)", las = 0)

mtext(side = 1, outer = FALSE, line = 2, "Year")
dev.off()

assmnt_ssb <- c(
  40088, 30523, 14503, 7309, 4905, 5998, 8076, 8496, 8598, 8576, 7577, 7906,
  7780, 9646, 7447, 8110, 7586, 6831, 8224, 8956, 9047, 11594, 17749, 21610, 23304, 24349, 23204
)

rel_ssb <- assmnt_ssb * mean(report$ssb) / mean(assmnt_ssb)

ylim <- range(rel_ssb, report$ssb)

gname <- "ssb.jpeg"
jpeg(file = gname, width = 8, height = 8, units = "in", res = 600)
par(mar = c(3, 3, 0.2, 0.2))
plot(year, report$ssb, type = "l", lwd = 2, xlab = "", ylab = "", ylim = ylim)
lines(year, rel_ssb, col = "red", lwd = 2)
mtext(side = 2, outer = FALSE, line = 2, "Relative SSB")
mtext(side = 1, outer = FALSE, line = 2, "Year")
legend("topright", col = c("black", "red"), lwd = 2, bty = "n", legend = c("SURBA++", "Assessment"))
dev.off()

ssb_matrix <- t(report$SSB_matrix)

gname <- "SSB_3D.jpeg"
jpeg(file = gname, width = 10, height = 10, units = "in", res = 600)

par(mar = c(0, 2, 1, 0), mgp = c(2, 0.7, 0))
persp(year, age_id, ssb_matrix,
  main = "SSB-at-age Surface",
  xlab = "Year", ylab = "Age", zlab = "Weight",
  col = "#4e8fcc",
  theta = 130, phi = 40,
  r = 50,
  d = 0.1,
  expand = 0.5,
  ltheta = 90, lphi = 180,
  shade = 0.25,
  ticktype = "detailed",
  nticks = 5
) # produces the 3-D plot
dev.off()

N <- t(report$N_matrix)

gname <- "N_3D.jpeg"
jpeg(file = gname, width = 10, height = 10, units = "in", res = 600)

par(mar = c(0, 2, 1, 0), mgp = c(2, 0.7, 0))
persp(year, age_id, N,
  main = "Numbers-at-age Surface",
  xlab = "Year", ylab = "Age", zlab = "Number",
  col = "#4e8fcc",
  theta = 130, phi = 40,
  r = 50,
  d = 0.1,
  expand = 0.5,
  ltheta = 90, lphi = 180,
  shade = 0.25,
  ticktype = "detailed",
  nticks = 5
) # produces the 3-D plot

dev.off()

gname <- "Nold_3D.jpeg"
jpeg(file = gname, width = 10, height = 10, units = "in", res = 600)

par(mar = c(0, 2, 1, 0), mgp = c(2, 0.7, 0))
persp(year, age_id[6:13], N[, 6:13],
  main = "Numbers-at-age Surface",
  xlab = "Year", ylab = "Age", zlab = "Number",
  col = "#4e8fcc",
  theta = 130, phi = 40,
  r = 50,
  d = 0.1,
  expand = 0.5,
  ltheta = 90, lphi = 180,
  shade = 0.25,
  ticktype = "detailed",
  nticks = 5
) # produces the 3-D plot
dev.off()

Z <- t(report$Z)

gname <- "Z_3D.jpeg"
jpeg(file = gname, width = 10, height = 10, units = "in", res = 600)

par(mar = c(0, 2, 1, 0), mgp = c(2, 0.7, 0))
persp(year, age_id, Z,
  main = "Z-at-age Surface",
  xlab = "Year", ylab = "Age", zlab = "Z",
  col = "#4e8fcc",
  theta = 130, phi = 40,
  r = 50,
  d = 0.1,
  expand = 0.5,
  ltheta = 90, lphi = 180,
  shade = 0.25,
  ticktype = "detailed",
  nticks = 5
) # produces the 3-D plot
dev.off()

wtm <- t(tmb_data$weight)
pcatch <- wtm * N * (1 - exp(-Z)) * (abs(Z - 0.2)) / Z
rownames(pcatch) <- year
colnames(pcatch) <- age_id
pland <- apply(pcatch, 1, sum)

assmnt_catch <- c(
  28846, 29454, 12752, 10646, 2702, 172,
  174, 383, 547, 919, 1050, 1310, 2194, 4870, 934,
  724, 600, 848, 923, 1083, 946, 867, 734, 1113, 734, 586, 666
)

rel_catch <- assmnt_catch * mean(pland) / mean(assmnt_catch)

ylim <- range(rel_catch, pland)

gname <- "catch_trend.jpeg"
jpeg(file = gname, width = 6, height = 8, units = "in", res = 600)
#  par(mar=c(3,3,0.2,0.2))
par(mfcol = c(2, 1), oma = c(0, 1, 0.5, 0.5), mar = c(3, 2, 0, 0), las = 1)
plot(year, pland, type = "l", lwd = 2, xlab = "", ylab = "", ylim = ylim)
lines(year, rel_catch, col = "red", lwd = 2)
mtext(side = 2, outer = FALSE, line = 2, "Relative Catch", las = 0)
mtext(side = 1, outer = FALSE, line = 2, "Year")
legend("topright", col = c("black", "red"), lwd = 2, bty = "n", legend = c("SURBA++", "Assessment"))


plot(log(rel_catch), log(pland), type = "p", xlab = "", ylab = "")
abline(a = 0, b = 1, col = "red")
mtext(side = 2, outer = FALSE, line = 2, "SURBA++ log relative Catch", las = 0)
mtext(side = 1, outer = FALSE, line = 2, "Assessment log relative Catch")
ltext <- paste(
  "corr=", round(cor(rel_catch, pland), digits = 3), " log corr=",
  round(cor(log(rel_catch), log(pland)), digits = 3)
)
legend("bottomright", legend = ltext, bty = "n", cex = 0.5)

dev.off()

fit <- data2
fit$survey <- "Fall 3LNO"
fit$resid <- report$sd_resid_index
fit$pred <- report$log_index_pred
fit$YC <- fit$year - fit$age

gname <- "sres.jpeg"
jpeg(file = gname, width = 8, height = 10, units = "in", res = 600)
par(mfcol = c(4, 1), oma = c(0, 3, 1, 1), mar = c(3, 2, 0, 2), las = 1)

ind <- !is.na(fit$pred)
x <- fit$year[ind]
y <- fit$resid[ind]

plot(x, y, xlab = "", ylab = "", type = "n")
text(x, y, fit$age[ind])
abline(h = 0, lty = 1)
mres <- tapply(y, x, "mean")
lines(as.numeric(names(mres)), mres, lty = 1, lwd = 1, col = "red")
mtext(side = 4, line = 0.5, outer = FALSE, "Year", las = 0)

ind1 <- ind
x1 <- fit$YC[ind1]
y1 <- fit$resid[ind1]
plot(x1, y1, xlab = "", ylab = "", type = "n")
text(x1, y1, fit$age[ind])
abline(h = 0, lty = 1)
mres <- tapply(y1, x1, "mean")
lines(sort(unique(x1)), mres, lty = 1, lwd = 1, col = "red")
mtext(side = 4, line = 0.5, outer = FALSE, "Cohort", las = 0)

plot(fit$age[ind], y, xlab = "", ylab = "", pch = 3)
abline(h = 0, lty = 2)
mtext(side = 4, line = 0.5, outer = FALSE, "Age", las = 0)
ma <- tapply(y, fit$age[ind], mean)
lines(as.numeric(names(ma)), ma, lty = 1)

plot(fit$pred[ind], y, xlab = "", ylab = "", pch = 3)
abline(h = 0, lty = 1)
mtext(side = 4, line = 0.5, outer = FALSE, "Expected", las = 0)

i <- "Fall 3LNO"
ytext <- paste("Standardized", i, "Residuals")
mtext(side = 2, ytext, line = 1, outer = TRUE, las = 0, cex = 1.2)

dev.off()

fit$colr <- "deepskyblue"
fit$colr[fit$resid > 0] <- "firebrick2"

jDarkGray <- "grey20"
jPch <- 21
resid.plot1 <- xyplot(factor(Age) ~ Year,
  data = fit,
  ylab = "Age",
  scales = list(y = 0:5, x = c(2000, 2005, 2010, 2015)),
  cex = 3 * sqrt(abs(fit$resid) / pi), fill.color = fit$colr,
  col = jDarkGray,
  par.settings = my_padding,
  par.strip.text = list(cex = 0.5),
  panel = function(x, y, ..., cex, fill.color, subscripts) {
    panel.abline(h = 1:6, v = c(2000, 2005, 2010, 2015), col.line = grey(0.9))
    panel.xyplot(x, y,
      cex = cex[subscripts],
      pch = jPch, fill = fill.color[subscripts], ...
    )
  }
)

# print(resid.plot1)

gname <- c("resid_matrix_bubbles.jpeg")
jpeg(file = gname, width = 4, height = 3, units = "in", res = 300)
print(resid.plot1)
dev.off()

vec_func <- function(mdat) {
  vdat <- melt(mdat, id = c("Year"))
  vdat$Age <- as.numeric(substr(vdat$variable, 4, 5))
  vdat$variable <- NULL
  vdat$index <- vdat$value
  vdat$value <- NULL
  return(vdat)
}

pe <- cbind(year, t(rep$pe))
colnames(pe) <- c("Year", paste("Age", age_id, sep = ""))

dat <- vec_func(as.data.frame(pe))
dat$colr <- "deepskyblue"
dat$colr[dat$index > 0] <- "firebrick2"

jDarkGray <- "grey20"
jPch <- 21
resid.plot1 <- xyplot(factor(Age) ~ Year,
  data = dat,
  ylab = "Age",
  scales = list(y = 0:5, x = c(2000, 2005, 2010, 2015)),
  cex = 10 * sqrt(abs(dat$index) / pi), fill.color = dat$colr,
  col = jDarkGray,
  par.settings = my_padding,
  par.strip.text = list(cex = 0.5),
  panel = function(x, y, ..., cex, fill.color, subscripts) {
    panel.abline(h = 1:6, v = c(2000, 2005, 2010, 2015), col.line = grey(0.9))
    panel.xyplot(x, y,
      cex = cex[subscripts],
      pch = jPch, fill = fill.color[subscripts], ...
    )
  }
)

# print(resid.plot1)

gname <- c("pe_matrix_bubbles.jpeg")
jpeg(file = gname, width = 4, height = 3, units = "in", res = 300)
print(resid.plot1)
dev.off()

par_names <- names(opt$par)
ind <- substr(par_names, 1, 6) == "log_sd"
x <- exp(opt$par[ind])
names(x) <- sub("log_sd_", "", names(x))

jpeg(file = "sd_barplot.jpeg", width = 8, height = 8, units = "in", res = 600)

par(mar = c(3, 4, 1, 1), mgp = c(2, 1, 0)) # increase y-axis margin.

barplot(x, names.arg = names(x), horiz = TRUE, xlab = "STD", las = 1, main = "STD")
abline(v = c(0.25, 0.5, 0.75, 1), col = "grey", lty = 2)

dev.off()


ind <- value_names == "log_Rec"
est <- exp(sd_report$value[ind])
sdi <- sd_report$sd[ind]
low <- exp(sd_report$value[ind] - qnorm(0.975) * sdi)
high <- exp(sd_report$value[ind] + qnorm(0.975) * sdi)
ylim <- range(low, high)
ylim[2] <- 1250

jpeg(file = "Rec.jpeg", width = 8, height = 6, units = "in", res = 600)
par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 1, 0))
plot(year, est, type = "l", lwd = 2, xlab = "", ylab = "", ylim = ylim)

sx <- year
polygon(c(sx, rev(sx)), c(low, rev(high)), col = "grey", border = NA)
box(lty = 1)
lines(year, est, type = "l", lty = 1, lwd = 2)
abline(h = mean(est), lty = 2)

mtext(side = 2, line = 2, "Rec. Age 1, MNPT")
mtext(side = 1, outer = FALSE, line = 2, "Year")
dev.off()

jpeg(file = "Sa.jpeg", width = 8, height = 6, units = "in", res = 600)
par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 1, 0))
y <- exp(report$log_s)
plot(age_id, y, type = "b", xlab = "", ylab = "")
abline(h = 1, lty = 2)
mtext(side = 1, line = 2, "Age")
mtext(side = 2, line = 2, "F Age Pattern")
dev.off()


ttl <- list(c("Observed vs predicted Index, Ages 0-12"), cex = 1)
xttl <- list(c("Year"), cex = 1)
yttl <- list(c("Index"), cex = 1)
stripttl <- list(cex = 0.2, font = 1)
ax <- list(cex = 0.5, relation = "free")
my_key <- list(
  space = "top", columns = 2,
  border = FALSE,
  size = 5,
  lines = list(lty = 1, col = 2:1, lwd = 2),
  text = list(c("Observed", "SURBA++"))
)

fit$epred <- exp(fit$pred)

jpeg(file = "obs_pred_byage.jpeg", width = 10, height = 10, units = "in", res = 600)
print(xyplot(index + epred ~ year | factor(age),
  data = fit, type = "l", main = ttl, xlab = xttl, lwd = 2,
  ylab = yttl, strip.text = stripttl, scales = ax, as.table = TRUE, layout = c(3, 5, 1), las = 1,
  key = my_key, col = 2:1, par.strip.text = list(cex = 0.6),
  par.settings = list(layout.heights = list(strip = .8))
))

dev.off()


xttl <- list(c("Year class"), cex = 1)
fit1 <- subset(fit, (YC >= 1982) & (YC <= 2013))
xlim <- range(age_id)

jpeg(file = "obs_pred_byYC.jpeg", width = 10, height = 14, units = "in", res = 600)
print(xyplot(index + epred ~ age | factor(YC),
  data = fit1, xlim = xlim,
  type = "l", main = ttl, xlab = xttl, lwd = 2,
  ylab = yttl, strip.text = stripttl, scales = ax, as.table = TRUE, layout = c(4, 8, 1), las = 1,
  key = my_key, col = 2:1, par.strip.text = list(cex = 0.6),
  par.settings = list(layout.heights = list(strip = .8))
))

dev.off()