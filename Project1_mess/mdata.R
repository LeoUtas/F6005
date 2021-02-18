rm(list = ls())

setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/Project1_mess")

library(reshape2)
library(stringr)
library(tidyr)

vec_func <- function(x) {
  vec_data <- gather(x, "Year", "Catch", -Length)
  vec_data$Year <- as.numeric(gsub("X", "", vec_data$Year))
  return(vec_data)
}

load("data\\3LN_redfish.RData")

year <- 1990:2019
len_pop <- 1:54
age_pop <- 1:10

Y <- length(year)
L <- length(len_pop)
A <- length(age_pop)


############### Fall RV indices ###############################;

len_rv <- sort(unique(FRV.vec$Length))

# create RV aggregate (i.e. plus) lower and upper length groups;

Lo <- 10
Hi <- 40

FRV <- FRV[, -c(2:6)] # drop Engels indices;
nY <- ncol(FRV)

# aggregate low length group (i.e., <= 10cm)

ind <- len_rv <= Lo
index_lo <- apply(FRV[ind, 2:nY], 2, sum)
index_len_lo <- min(FRV[ind, 1])

# aggregate high length group (i.e., >= 40cm)

ind <- len_rv >= Hi
index_hi <- apply(FRV[ind, 2:nY], 2, sum)
index_len_hi <- max(FRV[ind, 1])

ind <- (len_rv > Lo) & (len_rv < Hi)

FRVc <- rbind(c(Lo, index_lo), FRV[ind, ], c(Hi, index_hi))
# FRVc[FRVc==0]=1; ## replace one 0 with a 1

FRVc_vec <- vec_func(FRVc)
FRVc_vec <- subset(FRVc_vec, !is.na(FRVc_vec$Catch)) # remove NAs catch
FRVc_vec$Len1 <- FRVc_vec$Length                     # len1 ?
FRVc_vec$Len2 <- FRVc_vec$Length
FRVc_vec$Len1[FRVc_vec$Len1 == Lo] <- index_len_lo
FRVc_vec$Len2[FRVc_vec$Len2 == Hi] <- index_len_hi


############### Campelen indices ###############################;

len_c <- sort(unique(CL.vec$Length))

# aggregate low length group (i.e., <= 20cm)

Lo <- 20
ind <- len_c <= Lo
index_lo <- apply(CL[ind, 2:(Y + 1)], 2, sum)
index_len_lo <- min(CL[ind, 1])

# aggregate high length group (i.e., >= 40cm)

Hi <- 40
ind <- len_c >= Hi
index_hi <- apply(CL[ind, 2:(Y + 1)], 2, sum)
index_len_hi <- max(CL[ind, 1])

ind <- (len_c > Lo) & (len_c < Hi)

CLc <- rbind(c(Lo, index_lo), CL[ind, ], c(Hi, index_hi))
# CLc[CLc==0]=0.1; ## replace one 0 with a 1

CLc_vec <- vec_func(CLc)
CLc_vec <- subset(CLc_vec, Year >= min(year)) # no need
CLc_vec$Len1 <- CLc_vec$Length
CLc_vec$Len2 <- CLc_vec$Length
CLc_vec$Len1[CLc_vec$Len1 == Lo] <- index_len_lo
CLc_vec$Len2[CLc_vec$Len2 == Hi] <- index_len_hi



## ???
## some model dimensions;
A <- length(age_pop)
Y <- length(year)
L <- length(len_pop)
FRV_L <- length(unique(FRV.vec$Length))




## fill a matruity matrix from mat data;

mat_input <- mat_matrix

pmat1 <- rep(0, L)
pmat2 <- rep(0, L)
ind <- len_pop %in% mat_input$Length # ???
pmat1[ind] <- mat_input$P1
pmat2[ind] <- mat_input$P2
ind <- len_pop > max(mat_input$Length)
pmat1[ind] <- 1
pmat2[ind] <- 1

mat_matrix <- cbind(
  matrix(pmat1, nrow = L, ncol = length(1990:1995), byrow = FALSE),
  matrix(pmat2, nrow = L, ncol = length(1996:2019), byrow = FALSE)
)

## Miller and Hyun (2017) natural mortality - this is age x year like fishing mortality
## get rough size at age from Cadigan and Campana

len_age <- 33 * (1 - exp(-0.35 * (age_pop**0.75)))
amean <- mean(LW.parms$a[8:29])
bmean <- mean(LW.parms$b[8:29])
m_age <- amean * (len_age**(-0.305 * bmean))
m_age <- 0.05 * m_age / min(m_age)
par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 1, 0))
plot(age_pop, m_age, xlab = "age", ylab = "M at age", type = "l", lwd = 2)

M <- matrix(m_age, nrow = A, ncol = Y, byrow = FALSE)

## breaks between length bins

len_border <- len_pop[1:(L - 1)] + 0.5

## tmb indices for assessment years, ages and length;

iyear <- as.numeric(factor(year)) - 1
iage <- as.numeric(factor(age_pop)) - 1
ilen <- as.numeric(factor(len_pop)) - 1

# Remove 2 zero's from Fall RV and 6 zeros or really small values from Campelen

FRVc_vec <- subset(FRVc_vec, Catch > 0)
CLc_vec <- subset(CLc_vec, Catch > 0.5)

FRV_iyear <- iyear[match(FRVc_vec$Year, year)]
FRV_ilen1 <- ilen[match(FRVc_vec$Len1, len_pop)]
FRV_ilen2 <- ilen[match(FRVc_vec$Len2, len_pop)]

CL_iyear <- iyear[match(CLc_vec$Year, year)]
CL_ilen1 <- ilen[match(CLc_vec$Len1, len_pop)]
CL_ilen2 <- ilen[match(CLc_vec$Len2, len_pop)]

## fraction of year Fall RV survey occurs;

FRV_sf <- 10.5 / 12

tmb_data <- list(
  A = A,
  Y = Y,
  L = L,
  age = age_pop,
  len_mid = len_pop,
  len_border = len_border,
  log_index = log(FRVc_vec$Catch / 1000),
  sf = FRV_sf,
  weight = wt.mat[1:L, ],
  mat = mat_matrix,
  FRV_iyear = FRV_iyear,
  FRV_ilen1 = FRV_ilen1,
  FRV_ilen2 = FRV_ilen2,
  CL_iyear = CL_iyear,
  CL_ilen1 = CL_ilen1,
  CL_ilen2 = CL_ilen2,
  log_catch = log(CLc_vec$Catch / 1000),
  M = M
)

save(tmb_data, FRVc_vec, CLc_vec, file = "tmb_data.RData")