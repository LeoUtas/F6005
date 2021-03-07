rm(list = ls())

setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/project1_mess3")

vec_func <- function(x) {
  vec_data <- tidyr::gather(x, "Year", "Catch", -Length)
  vec_data$Year <- as.numeric(gsub("X", "", vec_data$Year))
  return(vec_data)
}

load("data\\3LN_redfish.RData")

ns <- length(survey)
year <- 1990:2019
len_pop <- 1:65
age_pop <- 1:10 # consider change age >10

Y <- length(year)
L <- length(len_pop)
A <- length(age_pop)

############### survey index ###############################;

FRV <- c(10, 40)
SRV <- c(10, 40)
Spanish_summer <- c(10, 40)
Spanish_3N <- c(10, 40)
lenbar <- list(FRV, SRV, Spanish_summer, Spanish_3N)

survey_org <- survey

len_survey <- list()
nY <- vector()
RV <- list()
ind <- list()
indlo <- list()
indhi <- list()
index_lo <- list()
index_len_lo <- list()
index_hi <- list()
index_len_hi <- list()
surveyc <- list()
surveyc_vec <- list()
lens <- list()

for (i in 1:ns) {
  len_survey[[i]] <- sort(unique(survey_vec[[i]]$Length))
  nY[i] <- ncol(survey[[i]])

  indlo[[i]] <- len_survey[[i]] <= lenbar[[i]][1]
  index_lo[[i]] <- apply(survey[[i]][indlo[[i]], 2:nY[i]], 2, sum)
  index_len_lo[[i]] <- min(survey[[i]][indlo[[i]], 1])

  indhi[[i]] <- len_survey[[i]] >= lenbar[[i]][2]
  index_hi[[i]] <- apply(survey[[i]][indhi[[i]], 2:nY[i]], 2, sum)
  index_len_hi[[i]] <- max(survey[[i]][indhi[[i]], 1])

  ind[[i]] <- (len_survey[[i]] > lenbar[[i]][1]) & (len_survey[[i]] < lenbar[[i]][2])

  surveyc[[i]] <- rbind(c(lenbar[[i]][1], index_lo[[i]]), survey[[i]][ind[[i]], ], c(lenbar[[i]][2], index_hi[[i]]))
  surveyc_vec[[i]] <- vec_func(surveyc[[i]])
  surveyc_vec[[i]] <- subset(surveyc_vec[[i]], !is.na(surveyc_vec[[i]]$Catch))
  surveyc_vec[[i]] <- subset(surveyc_vec[[i]], Catch > 0)
  surveyc_vec[[i]]$len1 <- surveyc_vec[[i]]$Length
  surveyc_vec[[i]]$len2 <- surveyc_vec[[i]]$Length
  surveyc_vec[[i]]$len1[surveyc_vec[[i]]$len1 == lenbar[[i]][1]] <- index_len_lo[[i]]
  surveyc_vec[[i]]$len2[surveyc_vec[[i]]$len2 == lenbar[[i]][2]] <- index_len_hi[[i]]
  if (i > 2) {
    surveyc_vec[[i]]$len2 <- surveyc_vec[[i]]$len2 + 1
  }  # I am assuming this is how the 2cm bins work;

  lens[[i]] <- sort(unique(survey_vec[[i]]$Length))
}

############### Catch index ###############################;

len_c <- sort(unique(CL_vec$Length))

# aggregate low length group (i.e., <= 20cm)
Lo <- 20
indc <- len_c <= Lo
indexc_lo <- apply(CL[indc, 2:(Y + 1)], 2, sum)
indexc_len_lo <- min(CL[indc, 1])

# aggregate high length group (i.e., >= 40cm)
Hi <- 40
indc <- len_c >= Hi
indexc_hi <- apply(CL[indc, 2:(Y + 1)], 2, sum)
indexc_len_hi <- max(CL[indc, 1])

indc <- (len_c > Lo) & (len_c < Hi)

CLc <- rbind(c(Lo, indexc_lo), CL[indc, ], c(Hi, indexc_hi))
# CLc[CLc==0]=0.1; ## replace one 0 with a 1

CLc_vec <- vec_func(CLc)
CLc_vec <- subset(CLc_vec, Year >= min(year))
CLc_vec <- subset(CLc_vec, Catch > 0.5)
CLc_vec$len1 <- CLc_vec$Length
CLc_vec$len2 <- CLc_vec$Length
CLc_vec$len1[CLc_vec$len1 == Lo] <- indexc_len_lo
CLc_vec$len2[CLc_vec$len2 == Hi] <- indexc_len_hi

############### Maturity ##################################;

mat_input <- mat

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



############### Mortality #################################;

## Miller and Hyun (2017) natural mortality - this is age x year like fishing mortality
## get rough size at age from Cadigan and Campana
len_age <- 33 * (1 - exp(-0.35 * (age_pop**0.75)))
amean <- mean(LW_parms$a[8:29])
bmean <- mean(LW_parms$b[8:29])
m_age <- amean * (len_age**(-0.305 * bmean))
m_age <- 0.05 * m_age / min(m_age)
par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 1, 0))
plot(age_pop, m_age, xlab = "age", ylab = "M at age", type = "l", lwd = 2)

M <- matrix(m_age, nrow = A, ncol = Y, byrow = FALSE)

############### Making tmb_data ###########################;

## breaks between length bins
len_border <- len_pop[1:(L - 1)] + 0.5

## tmb index for assessment years, ages and length;
iyear <- as.numeric(factor(year)) - 1
iage <- as.numeric(factor(age_pop)) - 1
ilen <- as.numeric(factor(len_pop)) - 1

# Remove 2 zero's from Fall RV and 6 zeros or really small values from Campelen

surveyc_vec[[1]]$sname <- "FallRV"
surveyc_vec[[1]]$sf <- 10.5 / 12
surveyc_vec[[2]]$sname <- "SprgRV"
surveyc_vec[[2]]$sf <- 5.5 / 12
surveyc_vec[[3]]$sname <- "Spsh3L"
surveyc_vec[[3]]$sf <- 8 / 12
surveyc_vec[[4]]$sname <- "Spsh3N"
surveyc_vec[[4]]$sf <- 5.5 / 12

allSL_vec <- rbind(surveyc_vec[[1]], surveyc_vec[[2]], surveyc_vec[[3]], surveyc_vec[[4]])
allSL_vec$is <- as.numeric(as.factor(allSL_vec$sname)) - 1
ind1 <- (allSL_vec$is == 0) & (allSL_vec$Year <= 1994) & (allSL_vec$len2 <= 20)
# 1995 FallRV Campelen
ind2 <- (allSL_vec$is == 1) & (allSL_vec$Year <= 1995) & (allSL_vec$len2 <= 20) # 1995 SprgRV Engel
allSL_vec <- subset(allSL_vec, !(ind1 | ind2))

allSL_vec$iyear <- iyear[match(allSL_vec$Year, year)]
allSL_vec$ilen1 <- ilen[match(allSL_vec$len1, len_pop)]
allSL_vec$ilen2 <- ilen[match(allSL_vec$len2, len_pop)]

CL_iyear <- iyear[match(CLc_vec$Year, year)]
CL_ilen1 <- ilen[match(CLc_vec$len1, len_pop)]
CL_ilen2 <- ilen[match(CLc_vec$len2, len_pop)]

sf <- tapply(allSL_vec$sf, allSL_vec$is, unique)

F_ratio <- c(rep(mean(sex_ratio$F), 16), sex_ratio$F)
juv_ratio <- c(rep(mean(sex_ratio$I), 16), sex_ratio$I)

tmb_data <- list(
  A = A,
  Y = Y,
  L = L,
  age = age_pop,
  len_mid = len_pop,
  len_border = len_border,

  log_index = log(allSL_vec$Catch / 1000),
  sf = sf,
  is = allSL_vec$is,

  weight = wt_mat[1:L, ],
  mat = mat_matrix,

  SL_iyear = allSL_vec$iyear,
  SL_ilen1 = allSL_vec$ilen1,
  SL_ilen2 = allSL_vec$ilen2,

  CL_iyear = CL_iyear,
  CL_ilen1 = CL_ilen1,
  CL_ilen2 = CL_ilen2,
  log_catch = log(CLc_vec$Catch / 1000),
  M = M,

  F_ratio = F_ratio,
  juv_ratio = juv_ratio
)

save(tmb_data, surveyc_vec, CLc_vec, allSL_vec, file = "tmb_data.RData")