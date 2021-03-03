rm(list = ls())

setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/Project1_mess")

vec_func <- function(x) {
  vec_data <- tidyr::gather(x, "Year", "Catch", -Length)
  vec_data$Year <- as.numeric(gsub("X", "", vec_data$Year))
  return(vec_data)
}

load("data\\3LN_redfish.RData")

ns <- length(survey)
year <- 1990:2019
len_pop <- 1:65 # 64 or 65
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
for (i in 1:2) {
  survey[[i]] <- survey[[i]][, -c(2:6)] # drop 1990:1995 in the Canadian surveys
}

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

s1_iyear <- iyear[match(surveyc_vec[[1]]$Year, year)]
s1_ilen1 <- ilen[match(surveyc_vec[[1]]$len1, len_pop)]
s1_ilen2 <- ilen[match(surveyc_vec[[1]]$len2, len_pop)]

s2_iyear <- iyear[match(surveyc_vec[[2]]$Year, year)]
s2_ilen1 <- ilen[match(surveyc_vec[[2]]$len1, len_pop)]
s2_ilen2 <- ilen[match(surveyc_vec[[2]]$len2, len_pop)]

s3_iyear <- iyear[match(surveyc_vec[[3]]$Year, year)]
s3_ilen1 <- ilen[match(surveyc_vec[[3]]$len1, len_pop)]
s3_ilen2 <- ilen[match(surveyc_vec[[3]]$len2, len_pop)]

s4_iyear <- iyear[match(surveyc_vec[[4]]$Year, year)]
s4_ilen1 <- ilen[match(surveyc_vec[[4]]$len1, len_pop)]
s4_ilen2 <- ilen[match(surveyc_vec[[4]]$len2, len_pop)]

CL_iyear <- iyear[match(CLc_vec$Year, year)]
CL_ilen1 <- ilen[match(CLc_vec$len1, len_pop)]
CL_ilen2 <- ilen[match(CLc_vec$len2, len_pop)]

## fraction of year when surveys occurs;
sf1 <- 10.5 / 12
sf2 <- 5.5 / 12
sf3 <- 5.5 / 12
sf4 <- 8 / 12

tmb_data <- list(
  A = A,
  Y = Y,
  L = L,
  age = age_pop,
  len_mid = len_pop,
  len_border = len_border,

  lens1 = lens[[1]],
  lens2 = lens[[2]],
  lens3 = lens[[3]],
  lens4 = lens[[4]],

  log_index1 = log(surveyc_vec[[1]]$Catch / 1000),
  log_index2 = log(surveyc_vec[[2]]$Catch / 1000),
  log_index3 = log(surveyc_vec[[3]]$Catch / 1000),
  log_index4 = log(surveyc_vec[[4]]$Catch / 1000),

  sf1 = sf1,
  sf2 = sf2,
  sf3 = sf3,
  sf4 = sf4,

  weight = wt_mat[1:L, ],
  mat = mat_matrix,

  s1_iyear = s1_iyear,
  s1_ilen1 = s1_ilen1,
  s1_ilen2 = s1_ilen2,

  s2_iyear = s2_iyear,
  s2_ilen1 = s2_ilen1,
  s2_ilen2 = s2_ilen2,

  s3_iyear = s3_iyear,
  s3_ilen1 = s3_ilen1,
  s3_ilen2 = s3_ilen2,

  s4_iyear = s4_iyear,
  s4_ilen1 = s4_ilen1,
  s4_ilen2 = s4_ilen2,

  CL_iyear = CL_iyear,
  CL_ilen1 = CL_ilen1,
  CL_ilen2 = CL_ilen2,
  log_catch = log(CLc_vec$Catch / 1000),
  M = M
)

save(tmb_data, surveyc_vec, CLc_vec, file = "tmb_data.RData")