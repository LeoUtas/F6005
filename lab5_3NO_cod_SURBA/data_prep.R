rm(list = ls())

setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab5_3NO_cod_SURBA")

library(tidyr)

# read and make maturity data

cnames <- c("year", paste("age", 1:12, sep = ""))
all_mat_data <- read.table(file = "mat.txt", header = FALSE, col.names = cnames)
mat_data <- subset(all_mat_data, (year >= 1990) & (year <= 2016))

mat_matrix <- matrix(0, nrow = length(1990:2016), ncol = length(0:12))
mat_matrix[, 2:13] <- as.matrix(mat_data[, 2:13])

## read and make weight data

cnames <- c("year", paste("age", 3:12, sep = ""))
all_wt_data <- read.table(file = "wt.txt", header = FALSE, col.names = cnames)
wt_data <- subset(all_wt_data, (year >= 1990) & (year <= 2016))

vwt_data <- gather(wt_data, "age", "wt", -year)
vwt_data$age <- as.numeric(gsub("age", "", vwt_data$age)) + 10 / 12
vwt_data$log_wt <- log(vwt_data$wt)

vwt_data$map <- as.numeric(as.factor(vwt_data$year))

VonB <- function(map, wt, age, log_Winf, log_k, log_po) {
  Winf <- exp(log_Winf)
  k <- exp(log_k)
  po <- exp(log_po)
  wt_pred <- Winf[map] * ((1 - (1 - po) * exp(-k[map] * age))**3.08)
  return(log(wt_pred))
}

year_id <- unique(vwt_data$year)
year_n <- length(unique(vwt_data$year))

start_parm <- list(log_Winf = log(rep(10, year_n)), log_k = log(rep(0.2, year_n)), log_po = log(0.02))
lower_parm <- list(log_Winf = log(rep(5, year_n)), log_k = log(rep(0.01, year_n)), log_po = log(0.0002))
upper_parm <- list(log_Winf = log(rep(35, year_n)), log_k = log(rep(0.4, year_n)), log_po = log(0.2))
lower <- unlist(lower_parm)
upper <- unlist(upper_parm)

VonB_fit <- nls(log_wt ~ VonB(map, wt, age, log_Winf, log_k, log_po),
  data = vwt_data, start = start_parm,
  algorithm = "port", control = list(maxiter = 5000),
  lower = lower, upper = upper
)

# vwt_data$pred_wt <- exp(predict(fit, vwt_data))

pred_data <- data.frame(year = rep(year_id, 13), age = rep(0:12 + 10 / 12, each = year_n))
pred_data$map <- as.numeric(as.factor(pred_data$year))
pred_data$log_wt_pred <- predict(VonB_fit, newdata = pred_data)
pred_data$wt_pred <- exp(pred_data$log_wt_pred)

wt_matrix <- matrix(pred_data$wt, nrow = year_n, ncol = 13, byrow = FALSE)
wt_matrix[, 4:13] <- as.matrix(wt_data[, 2:11]) # replace predicted weights by observed ones from age 3 to age 12

# make mortality data using Lorenzen method

pred_data <- data.frame(year = rep(year_id, 13), age = rep(0:12 + 0.5, each = year_n))
pred_data$map <- as.numeric(as.factor(pred_data$year))
pred_data$log_wt_pred <- predict(VonB_fit, newdata = pred_data)
pred_data$wt_pred <- exp(pred_data$log_wt_pred)

ave_wt <- aggregate(pred_data$wt_pred, list(age = pred_data$age), mean)
m_age <- 0.2 * ((ave_wt$x / 35)**(-0.305))
m_age <- 0.2 * m_age / m_age[length(m_age)]

# make index data

vec_func <- function(x) {
  vec_data <- gather(x, "age", "index", -year)
  vec_data$age <- as.numeric(gsub("age", "", vec_data$age))
  return(vec_data)
}

cnames <- c("year", paste("age", 1:13, sep = ""))
x <- read.table(file = "FRV.dat", header = FALSE, col.names = cnames)
index_vec <- vec_func(x)

# make tmb_data;

year <- c(1990:2016)

wt_temp <- cbind(year, wt_matrix)
colnames(wt_temp) <- cnames
wt_vec <- vec_func(as.data.frame(wt_temp))

mat_temp <- cbind(year, mat_matrix)
colnames(mat_temp) <- cnames
mat_vec <- vec_func(as.data.frame(mat_temp))

data <- data.frame(
  year = mat_vec$year,
  age = mat_vec$age,
  weight = wt_vec$index,
  mat = mat_vec$index
)

data2 <- merge(data, index_vec, by = c("year", "age"), all.x = TRUE)
data2$iyear <- as.numeric(factor(data2$year)) - 1
data2$iage <- as.numeric(factor(data2$age)) - 1

data2 <- subset(data2, !is.na(data2$index)) # remove NA data
data2 <- subset(data2, data2$index > 0) # remove 0 data

data2$log_index <- log(data2$index)
data2$sf <- 10 / 12

A <- length(1:12)
Y <- length(year)

q_age <- c(0.01, 0.3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
log_q_age <- log(q_age)

tmb_data <- list(
  log_index = data2$log_index,
  log_q = log_q_age,
  sf = data2$sf,
  weight = t(wt_matrix),
  mat = t(mat_matrix),
  iyear = data2$iyear,
  iage = data2$iage,
  m = m_age
)

save(data, tmb_data, file = "tmb_data.RData")

