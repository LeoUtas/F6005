setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/lab4/3LN_redfish_lab4")

catch <- read.table("catch.txt",header=TRUE, fill = FALSE)
catch$catch = catch$catch/1000 ## units in Kilo tonnes
catch$imap = 1:nrow(catch)-1; ## TMB starts at 0;

index <- read.table("indices.txt",header=TRUE, fill = FALSE)
index$iyear = catch$imap[match(index$year,catch$year)] 
index$iq = as.numeric(as.factor(index$name))-1
## index$iyear matches each survey year with a catch year

tmb.data = list(year = catch$year, C=catch$catch, index = index$index,
 iyear  = index$iyear, iq = index$iq
)

tmb.data$log_C = log(tmb.data$C) 
tmb.data$log_index = log(tmb.data$index)  

tmb.data$E_log_r = log(0.1)
tmb.data$sd_log_r = 0.25 

  
save(tmb.data,file='tmb.RData')

tmb.data
