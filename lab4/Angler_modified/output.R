setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab4/BAngler")
# source("xxx",local=F)

library("lattice")
library("xtable")

load("tmb_data.RData")
load("sspm_fit.RData")
index_data <- read.table("indices.txt", header = TRUE, fill = FALSE)

################# CPUE plot ############################;

uindex <- unique(index_data$name)
posv <- c("topright", rep("topleft", 7))

for (i in 1:length(uindex)) {
  u <- uindex[i]
  pos <- posv[i]

  scale <- 1
  # index.name = "Index (000s)"
  index_name <- "Index (000s)"
  if (u == "CPUE") {
    scale <- 1
    index_name <- "Index"
  }

  ind <- index_data$name == u
  year <- index_data$year[ind]
  index <- index_data$index[ind] / scale
  index_pred <- report$index_pred[ind] / scale
  std_resid <- report$std_resid[ind]

  jpeg(file = paste("sspm_ar1_sd_logC_fixed", u, ".jpeg", sep = ""), width = 8, height = 8, units = "in", res = 600)

  par(mfcol = c(3, 1), oma = c(1, 1, 1.5, 0.5), mar = c(3, 4, 0, 0), cex.axis = 1.2, las = 1)

  ylim <- range(index, index_pred)
  xlim <- range(index_data$year)

  xi <- tmb_data$year
  yi <- rep(NA, length(xi))
  pyi <- yi
  pres <- yi
  ind1 <- xi %in% year
  yi[ind1] <- index
  pyi[ind1] <- index_pred
  pres[ind1] <- std_resid

  plot(xi, yi, type = "b", lwd = 2, pch = 19, xlab = "", ylab = "", ylim = ylim, xlim = xlim)
  lines(xi, pyi, type = "b", lwd = 2)
  mtext(side = 1, line = 2, outer = FALSE, "Year")
  mtext(side = 2, line = 3.5, outer = FALSE, index_name, las = 0)
  legend(pos, bty = "n", lwd = 2, pch = c(19, 1), lty = 1, c("Observed", "Predicted"), pt.cex = 1.5)

  plot(xi, pres, type = "b", xlab = "", ylab = "", xlim = xlim)
  abline(h = 0, lty = 2)
  mtext(side = 1, line = 2, outer = FALSE, "Year")
  mtext(side = 2, line = 3.5, outer = FALSE, "Std. residual", las = 0)

  plot(pyi, pres, type = "p", xlab = "", ylab = "")
  mtext(side = 1, line = 2.5, outer = FALSE, "Predicted")
  mtext(side = 2, line = 3.5, outer = FALSE, "Std. residual", las = 0)

  mtext(side = 3, line = 0, outer = TRUE, u, las = 0)

  dev.off()
}

##############  biomass and harvest plot ################;

{
  value_names <- names(sd_report$value)

  ind <- value_names == "log_B"
  biomass <- exp(sd_report$value[ind])
  stdi <- sd_report$sd[ind]
  biomass_L95 <- exp(sd_report$value[ind] - qnorm(0.975) * stdi)
  biomass_U95 <- exp(sd_report$value[ind] + qnorm(0.975) * stdi)

  ind <- value_names == "log_H"
  exploit <- exp(sd_report$value[ind])
  stdi <- sd_report$sd[ind]
  exploit_L95 <- exp(sd_report$value[ind] - qnorm(0.975) * stdi)
  exploit_U95 <- exp(sd_report$value[ind] + qnorm(0.975) * stdi)

  jpeg(file = "sspm_ar1_sd_logC_fixed_pop.jpeg", width = 8, height = 10, units = "in", res = 1200)

  par(mfcol = c(2, 1), oma = c(1, 1.5, 0, 0.5), mar = c(2, 2, 1, 0), cex.axis = 1, las = 1)
  xlim <- c(min(tmb_data$year), max(tmb_data$year))
  ylim <- c(min(biomass_L95), max(biomass_U95))
  ylim[2] <- 40

  n <- length(tmb_data$year)

  ytext <- c("Exploitable biomass (Kt)")
  plot(xlim, ylim, xlab = "", ylab = "", type = "n", cex.lab = 1, lwd = 2)

  sx <- tmb_data$year
  low <- biomass_L95
  high <- biomass_U95
  polygon(c(sx, rev(sx)), c(low, rev(high)),
    #  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
    col = "grey", border = NA
  )
  box(lty = 1)
  lines(tmb_data$year, biomass, type = "l", lty = 1, lwd = 2)

  # lines(tmb_data$year,biomass.L95,type='l',lty=2,lwd=2)
  # lines(tmb_data$year,biomass.U95,type='l',lty=2,lwd=2)
  abline(h = report$Bmsy, lwd = 2, col = "red")
  abline(h = 0.5 * report$Bmsy, lwd = 1, col = "red")
  mtext(side = 2, line = 2.5, outer = FALSE, ytext, cex = 1, las = 0)
  p1 <- report$Bmsy + 0.02 * diff(range(ylim))
  text(tmb_data$year[n] - 2, p1, "Bmsy", cex = 0.75)
  p1 <- 0.5 * report$Bmsy + 0.02 * diff(range(ylim))
  text(tmb_data$year[n] - 3, p1, "1/2 Bmsy", cex = 0.75)

  # legend("topleft",bty='n',lwd=2,lty=2:1,c("95% confidence interval","Predicted"),cex=0.75)

  ylim <- c(min(exploit_L95), min(max(exploit_U95)))
  ylim[2] <- 0.4

  ytext <- c("Exploitation rate")
  plot(xlim, ylim, xlab = "", ylab = "", type = "n", cex.lab = 1, lwd = 2)

  sx <- tmb_data$year
  low <- exploit_L95
  high <- exploit_U95
  polygon(c(sx, rev(sx)), c(low, rev(high)),
    #  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
    col = "grey", border = NA
  )
  box(lty = 1)

  lines(tmb_data$year, exploit, type = "l", lty = 1, lwd = 2)

  # lines(tmb_data$year,exploit.L95,type='l',lty=2,lwd=2)
  # lines(tmb_data$year,exploit.U95,type='l',lty=2,lwd=2)
  abline(h = report$Hmsy, lwd = 2, col = "red")
  r.est <- exp(opt$par[1])
  sigma2 <- exp(opt$par[15])**2
  # phi_opt =  r.est/2 - 2*(2-r.est)*sigma2/((4-r.est)**2)
  phi_opt <- report$Hmsy
  abline(h = phi_opt, lwd = 2, col = "red", lty = 2)
  abline(h = 0.667 * phi_opt, lwd = 1, col = "red", lty = 2)
  mtext(side = 2, line = 2.5, outer = FALSE, ytext, cex = 1, las = 0)
  # p1 = rep$Hmsy + 0.07*diff(range(ylim))
  p1 <- report$Hmsy + 0.02 * diff(range(ylim))
  text(tmb_data$year[n] - 2, p1, "Hmsy", cex = 0.75)
  p1 <- 0.667 * report$Hmsy + 0.02 * diff(range(ylim))
  text(tmb_data$year[n] - 3, p1, "2/3 Hmsy", cex = 0.75)

  mtext(side = 1, line = 2, outer = FALSE, "Year", cex = 1.2)

  dev.off()
}

##############  biomass and harvest wrt MSY plot ################;

{
  value_names <- names(sd_report$value)

  ind <- value_names == "log_rB"
  biomass <- exp(sd_report$value[ind])
  stdi <- sd_report$sd[ind]
  biomass_L95 <- exp(sd_report$value[ind] - qnorm(0.975) * stdi)
  biomass_U95 <- exp(sd_report$value[ind] + qnorm(0.975) * stdi)
  ind <- value_names == "log_rH"
  exploit <- exp(sd_report$value[ind])
  stdi <- sd_report$sd[ind]
  exploit_L95 <- exp(sd_report$value[ind] - qnorm(0.975) * stdi)
  exploit_U95 <- exp(sd_report$value[ind] + qnorm(0.975) * stdi)
  jpeg(file = "sspm_ar1_sd_logC_fixed_rpop.jpeg", width = 8, height = 10, units = "in", res = 1200)

  par(mfcol = c(2, 1), oma = c(1, 1.5, 0, 0.5), mar = c(2, 2, 1, 0), cex.axis = 1, las = 1)
  xlim <- c(min(tmb_data$year), max(tmb_data$year))
  ylim <- c(min(biomass_L95), max(biomass_U95))
  ylim[2] <- 2

  n <- length(tmb_data$year)

  ytext <- c("Biomass/Bmsy")
  plot(xlim, ylim, xlab = "", ylab = "", type = "n", cex.lab = 1, lwd = 2)

  sx <- tmb_data$year
  low <- biomass_L95
  high <- biomass_U95
  polygon(c(sx, rev(sx)), c(low, rev(high)),
    col = "grey", border = NA
  )
  box(lty = 1)
  lines(tmb_data$year, biomass, type = "l", lty = 1, lwd = 2)

  abline(h = 1, lwd = 2, col = "red")
  abline(h = 0.5, lwd = 1, col = "red")
  mtext(side = 2, line = 2.5, outer = FALSE, ytext, cex = 1, las = 0)

  ylim <- c(min(exploit_L95), min(max(exploit_U95)))
  ylim[2] <- 2

  ytext <- c("H/Hmsy")
  plot(xlim, ylim, xlab = "", ylab = "", type = "n", cex.lab = 1, lwd = 2)

  sx <- tmb_data$year
  low <- exploit_L95
  high <- exploit_U95
  polygon(c(sx, rev(sx)), c(low, rev(high)),
    col = "grey", border = NA
  )
  box(lty = 1)

  lines(tmb_data$year, exploit, type = "l", lty = 1, lwd = 2)

  abline(h = 1, lwd = 2, col = "red", lty = 1)
  abline(h = 0.667, lwd = 1, col = "red", lty = 1)
  mtext(side = 2, line = 2.5, outer = FALSE, ytext, cex = 1, las = 0)

  mtext(side = 1, line = 2, outer = FALSE, "Year", cex = 1.2)

  dev.off()
}

############  estimated catch ##################;

{
  jpeg("sspm_ar1_sd_logC_fixed_catch.jpeg", width = 10, height = 8, units = "in", res = 1200)

  par(mar = c(4, 4, 1, 1), cex.axis = 1.2, las = 1)

  ylim <- c(min(exp(report$log_C_pred), tmb_data$C), max(exp(report$log_C_pred), tmb_data$C))

  plot(tmb_data$year, exp(report$log_C_pred),
    type = "l", las = 1, cex.lab = 1.2,
    xlab = "Year", ylab = "", lwd = 2, ylim = ylim
  )

  # lines(tmb_data$year,tmb_data$C_L,lwd=2,col='grey')
  # lines(tmb_data$year,tmb_data$C_U,lwd=2,col='grey')
  points(tmb_data$year, tmb_data$C)
  yname <- c("Catch (Kt)")
  mtext(side = 2, line = 2.5, yname, las = 0, cex = 1.2)

  # legend("topright",bty='n',lwd=2,col=c('grey','black','black'),lty=c(1,1,NA),pch=c(NA,NA,1),
  # c("Bounds","Predicted","Reported"))

  dev.off()
}

############  process error ##################;

{
  jpeg("sspm_ar1_sd_logC_fixed_pe.jpeg", width = 8, height = 8, units = "in", res = 1200)

  par(mar = c(3, 4, 0.5, 0.5), cex.axis = 1, las = 1)

  plot(tmb_data$year, report$log_pe, type = "l", las = 1, xlab = "", ylab = "", lwd = 2)
  mtext(side = 2, line = 3, c("Log process error"), las = 0)
  mtext(side = 1, line = 2, c("Year"))
  abline(h = 0, lty = 2)

  dev.off()
}

############ correlation plot ##############;

{
  sd <- sqrt(diag(sd_report$cov.fixed))
  corrm <- diag(1 / sd) %*% sd_report$cov.fixed %*% diag(1 / sd)

  pnames <- names(sd)
  qname <- paste("log_q", levels(tmb_data$name))
  pnames <- replace(pnames, pnames == "log_q", qname)
  stdname <- paste("log_std", levels(tmb_data$name))
  pnames <- replace(pnames, pnames == "log_sd_log_index", stdname)
  rownames(corrm) <- pnames
  colnames(corrm) <- pnames

  require(corrplot)

  jpeg(file = "sspm_ar1_sd_logC_fixed_corr.jpeg", width = 8, height = 8, units = "in", res = 600)

  corrplot.mixed(corrm,
    tl.srt = 45, tl.cex = 0.6, cl.cex = 0.6, tl.col = "black", tl.pos = "lt", diag = "u", number.cex = 0.6, title = "Correlation in parameter estimates",
    mar = c(0, 0, 1, 0)
  )

  dev.off()
}

########### table ########################;

{
  n <- length(tmb_data$year)

  out <- matrix(NA, n, 8)
  value_names <- names(sd_report$value)

  tnames <- c("log_B", "log_H", "log_rB", "log_rH")

  for (i in 1:length(tnames)) {
    ind <- value_names == tnames[i]
    est <- exp(sd_report$value[ind])
    cv <- sd_report$sd[ind]
    out[, 2 * i - 1] <- est
    out[, 2 * i] <- 100 * cv
  }
  # out[,1]=out[,1]
  out[, 3] <- out[, 3] * 100

  xalign <- c("c", rep("r", 8))
  # alignment 'l','c', or 'r';
  xdigits <- c(0, 1, 0, 1, 0, 2, 0, 2, 0)
  colnames(out) <- c("Biomass", "CV(%)", "H (%)", "CV(%)", "B/Bmsy", "CV(%)", "H/Hmsy", "CV(%)")
  rownames(out) <- as.character(tmb_data$year)
  ctext <- c("Model output")
  print(xtable(out, digits = xdigits, caption = ctext, align = xalign),
    type = "html",
    file = "sspm_ar1_sd_logC_fixed_pop.doc", caption.placement = "top"
  )
}

##################### Stats ################################# ;

{
  value_names <- names(sd_report$value)
  par_names <- names(sd_report$par.fixed)
  parms <- sd_report$par.fixed
  parms_sd <- sqrt(diag(sd_report$cov.fixed))

  out <- matrix(NA, 19, 2)

  j <- 1
  ind <- par_names == "log_r"
  est <- exp(parms[ind])
  cv <- parms_sd[ind]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- 2
  ind <- par_names == "log_K"
  est <- exp(parms[ind])
  cv <- parms_sd[ind]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- 3
  ind <- par_names == "log_Po"
  est <- exp(parms[ind])
  cv <- parms_sd[ind]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- 4:6
  ind <- par_names == "log_q"
  est <- exp(parms[ind]) * 1000
  # est[8]=est[8]*1000
  cv <- parms_sd[ind]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- 7:9
  ind <- par_names == "log_sd_log_index"
  est <- exp(parms[ind])
  cv <- parms_sd[ind]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- 10
  ind <- par_names == "log_sd_pe"
  est <- exp(parms[ind])
  cv <- parms_sd[ind]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- j + 1
  ind <- par_names == "log_sd_rw"
  est <- exp(parms[ind])
  cv <- parms_sd[ind]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  value_names <- names(sd_report$value)

  j <- j + 1
  ind <- value_names == "ar_pe"
  est <- sd_report$value[ind]
  cv <- sd_report$sd[ind] / est
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- j + 1
  ind <- value_names == "Hmsy"
  est <- 100 * sd_report$value[ind]
  cv <- 100 * sd_report$sd[ind] / est
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- j + 1
  ind <- value_names == "Bmsy"
  est <- sd_report$value[ind]
  cv <- sd_report$sd[ind] / est
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- j + 1
  ind <- value_names == "MSY"
  est <- sd_report$value[ind]
  cv <- sd_report$sd[ind] / est
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- j + 1
  ind <- value_names == "log_H"
  ind1 <- tmb_data$year == 2015
  est <- 100 * exp(sd_report$value[ind][ind1])
  cv <- sd_report$sd[ind][ind1]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- j + 1
  ind <- value_names == "log_B"
  ind1 <- tmb_data$year == 2015
  est <- exp(sd_report$value[ind][ind1])
  cv <- sd_report$sd[ind][ind1]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- j + 1
  ind <- value_names == "log_rH"
  ind1 <- tmb_data$year == 2015
  est <- exp(sd_report$value[ind][ind1])
  cv <- sd_report$sd[ind][ind1]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  j <- j + 1
  ind <- value_names == "log_rB"
  ind1 <- tmb_data$year == 2015
  est <- exp(sd_report$value[ind][ind1])
  cv <- sd_report$sd[ind][ind1]
  out[j, 1] <- round(est, digits = 3)
  out[j, 2] <- round(cv, digits = 3)

  qname <- paste("q", levels(factor(index_data$name)))
  stdname <- paste("std", levels(factor(index_data$name)))

  rownames(out) <- c(
    "r", "K (000s)", "Po", qname, stdname,
    "sd_pe", "sd_rw", "ar_pe", "Hmsy(%)", "Bmsy (000s)", "MSY (000s)", "H2016(%)", "B2016 (000s)", "H2016/Hmsy", "B2016/Bmsy"
  )
  colnames(out) <- c("est", "CV")

  tname <- paste("Model results, nll = ", round(sspm_fit$objective, digits = 2), sep = "")
  xdigits <- c(1, 3, 3)
  # a digit for rowname, and all rows;
  xalign <- c("c", rep("r", 2))
  # alignment 'l','c', or 'r';

  print(xtable(out, digits = xdigits, caption = tname, align = xalign),
    type = "html",
    file = "sspm_ar1_sd_logC_fixed_stats.doc", caption.placement = "top"
  )
}