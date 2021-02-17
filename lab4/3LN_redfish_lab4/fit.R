#source("xxx",local=F)
rm(list = ls())
setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\Lecture4\\3LN_redfish\\data")

load("tmb.RData")


setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\Lecture4\\3LN_redfish")

library(TMB)    

compile("fit.cpp")  

dyn.load("fit")  
#dyn.unload("fit");

parameters <- list(
  log_r = log(0.11), 
  log_K = log(4*20/0.11), 
  log_q = rep(1,length(unique(tmb.data$iq))), 
  log_Po = log(0.5), 
  log_Ho = log(0.1),
  log_sd_rw = log(0.2),
  log_sd_log_index = log(0.3),  
  log_sd_pe = log(0.1),    
  logit_ar_pe = log(0.50/(1-0.50)), 
  log_sd_logC = log(0.1),    
  log_pe = rep(0,length(tmb.data$C)),
  log_H_dev = rep(0,length(tmb.data$C)-1)
)

parameters.L <- list( 
  log_r = log(0.01), 
  log_K = log(4*20/0.4), 
  log_q = rep(-Inf,length(unique(tmb.data$iq))), 
  log_Po = log(0.001),
  log_Ho = log(0.0001),
  log_sd_rw = log(0.01),
  log_sd_log_index = log(0.01),
  log_sd_pe = log(0.05),
  logit_ar_pe = log(0.01/(1-0.01)) 
#  log_sd_logC = log(0.001)
)

parameters.U <- list( 
  log_r = log(0.4), 
  log_K = log(4*20/0.01),  
  log_q = rep(Inf,length(unique(tmb.data$iq))), 
  log_Po = log(10),
  log_Ho = log(1),
  log_sd_rw = log(2),
  log_sd_log_index = log(1),
  log_sd_pe = log(0.35),
  logit_ar_pe = log(0.950/(1-0.950)) 
#  log_sd_logC = log(1)
)

lower = unlist(parameters.L);
upper = unlist(parameters.U);

## random effects;
rname = c("log_pe","log_H_dev")

map = list(log_sd_logC = factor(NA))

obj <- MakeADFun(tmb.data,parameters,random=rname,DLL="fit",
  map=map,inner.control=list(maxit=100,trace=T))  

obj$gr(obj$par)

length(obj$par)
length(lower)
length(upper) 


opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
control = list(trace=0,iter.max=5000,eval.max=10000))

opt$message
obj$gr(opt$par)

rep = obj$report()

sd.rep<-sdreport(obj)

save.image(file="fit.RData")


