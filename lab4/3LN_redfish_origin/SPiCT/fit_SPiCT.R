#library(devtools)
library(spict)

setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\Lecture4\\3LN_redfish\\data")

catch <- read.table("catch.txt",header=TRUE, fill = FALSE)
catch$catch = catch$catch/1000 ## units in Kilo tonnes

index <- read.table("indices.txt",header=TRUE, fill = FALSE)

setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\Lecture4\\3LN_redfish\\SPiCT")

inp <- list(timeC=catch$year, obsC=catch$catch)
u.index = unique(index$name)

inp$timeI <- list(
index$year[index$name==u.index[1]]+6/12,
index$year[index$name==u.index[2]]+4/12,
index$year[index$name==u.index[3]]+10/12,
index$year[index$name==u.index[4]]+6/12,
index$year[index$name==u.index[5]]+2/12,
index$year[index$name==u.index[6]]+6/12,
index$year[index$name==u.index[7]]+10/12,
index$year[index$name==u.index[8]]+6/12
)


inp$obsI <- list()
inp$obsI[[1]] <- index$index[index$name==u.index[1]] 
inp$obsI[[2]] <- index$index[index$name==u.index[2]]/1000
inp$obsI[[3]] <- index$index[index$name==u.index[3]]/1000
inp$obsI[[4]] <- index$index[index$name==u.index[4]]/1000
inp$obsI[[5]] <- index$index[index$name==u.index[5]]/1000
inp$obsI[[6]] <- index$index[index$name==u.index[6]]/1000
inp$obsI[[7]] <- index$index[index$name==u.index[7]]/1000
inp$obsI[[8]] <- index$index[index$name==u.index[8]]/1000

par(mgp=c(2,1,0),mar=c(2,2,1,0.5))
plotspict.data(inp)

inp$ini$logK <- log(4*20/0.11) 
#inp$ini$logr <- log(0.11)

inp$phases$logn=-1

inp$priors$logm <- c(log(20), 0.25, 1)
inp$mapsdi <- rep(1, 8)

res = fit.spict(inp)

capture.output(summary(res))

plot(res)

plotspict.bbmsy(res)

plotspict.ffmsy(res)

plotspict.biomass(res)

plotspict.catch(res)

plotspict.fb(res)


