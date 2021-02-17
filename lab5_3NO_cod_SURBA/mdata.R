setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab5_3NO_cod_SURBA")

require(reshape2)

vec_func <- function(mdat) {
  vdat <- melt(mdat, id = c("Year"))
  vdat$Age <- as.numeric(substr(vdat$variable, 4, 5))
  vdat$variable <- NULL
  vdat$index <- vdat$value
  vdat$value <- NULL
  return(vdat)
}


age <- 0:12

cnames <- c("Year", paste("Age", age, sep = ""))
mdata <- read.table(file = "FRV.dat", header = FALSE, col.names = cnames)
year <- seq(min(mdata[, 1]), max(mdata[, 1]), by = 1)

## read and make maturities

cnames <- c("Year", paste("Age", 1:12, sep = ""))
all_mat_data <- read.table(file = "mat.txt", header = FALSE, col.names = cnames)
mat_data <- subset(all_mat_data, (Year >= 1990) & (Year <= 2016))

mat <- matrix(0, nrow = length(1990:2016), ncol = length(0:12))
mat[, 2:13] <- as.matrix(mat_data[, 2:13])

## read and make weights

cnames <- c("Year", paste("Age", 3:12, sep = ""))
all_wt_data <- read.table(file = "wt.txt", header = FALSE, col.names = cnames)
wt_data <- subset(all_wt_data, (Year >= 1990) & (Year <= 2016))

vwt_data <- melt(wt_data, id = c("Year"))
vwt_data$Age <- as.numeric(substr(vwt_data$variable, 4, 5)) + 10 / 12
vwt_data$variable <- NULL
vwt_data$wt <- vwt_data$value
vwt_data$value <- NULL
vwt_data$logwt <- log(vwt_data$wt)

vwt_data$map <- as.numeric(as.factor(vwt_data$Year))
uyear <- unique(vwt_data$Year)
nY <- length(unique(vwt_data$Year))

gfit <- function(map, wt, Age, logWinf, logk, logpo) {
  Winf <- exp(logWinf)
  k <- exp(logk)
  po <- exp(logpo)
  pred <- Winf[map] * ((1 - (1 - po) * exp(-k[map] * Age))**3.08)
  return(log(pred))
}

start.parm <- list(logWinf = log(rep(10, nY)), logk = log(rep(0.2, nY)), logpo = log(0.02))
lower.parm <- list(logWinf = log(rep(5, nY)), logk = log(rep(0.01, nY)), logpo = log(0.0002))
upper.parm <- list(logWinf = log(rep(35, nY)), logk = log(rep(0.4, nY)), logpo = log(0.2))
lower <- unlist(lower.parm)
upper <- unlist(upper.parm)

fit <- nls(logwt ~ gfit(map, wt, Age, logWinf, logk, logpo),
  data = vwt_data, start = start.parm,
  algorithm = "port", control = list(maxiter = 5000),
  lower = lower, upper = upper
)

vwt_data$pred_wt <- exp(predict(fit, vwt_data))

pred.dat <- data.frame(Year = rep(uyear, 13), Age = rep(0:12 + 10 / 12, each = nY))
pred.dat$map <- as.numeric(as.factor(pred.dat$Year))
pred.dat$logwt <- predict(fit, newdata = pred.dat)
pred.dat$wt <- exp(pred.dat$logwt)

wtm <- matrix(pred.dat$wt, nrow = nY, ncol = 13, byrow = FALSE)
wtm[, 4:13] <- as.matrix(wt_data[, 2:11])

pred.dat <- data.frame(Year = rep(uyear, 13), Age = rep(0:12 + 0.5, each = nY))
pred.dat$map <- as.numeric(as.factor(pred.dat$Year))
pred.dat$logwt <- predict(fit, newdata = pred.dat)
pred.dat$wt <- exp(pred.dat$logwt)

ave_wt <- aggregate(pred.dat$wt, list(age = pred.dat$Age), mean)
m_age <- 0.2 * ((ave_wt$x / 35)**(-0.305))
m_age <- 0.2 * m_age / m_age[length(m_age)]

FRV.vec <- vec_func(mdata)

temp <- cbind(year, wtm)
colnames(temp) <- colnames(mdata)
wt.vec <- vec_func(as.data.frame(temp))

temp <- cbind(year, mat)
colnames(temp) <- colnames(mdata)
mat.vec <- vec_func(as.data.frame(temp))

data <- data.frame(
  Year = mat.vec$Year,
  Age = mat.vec$Age,
  weight = wt.vec$index,
  mat = mat.vec$index
)
data <- merge(data, FRV.vec, by = c("Year", "Age"), all.x = TRUE)
data$iyear <- as.numeric(factor(data$Year)) - 1
data$iage <- as.numeric(factor(data$Age)) - 1

data <- subset(data, !is.na(data$index))
data <- subset(data, data$index > 0)

data$log_index <- log(data$index)
data$sf <- 10 / 12


A <- length(age)
Y <- length(year)

q_age <- c(0.01, 0.3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
log_q_age <- log(q_age)

tmb.data <- list(
  log_index = data$log_index,
  logq = log_q_age,
  sf = data$sf,
  weight = t(wtm),
  mat = t(mat),
  iyear = data$iyear,
  iage = data$iage,
  m = m_age
)


save(data, tmb.data, file = "tmb.RData")