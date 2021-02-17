
setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\Lecture5\\3NO_cod_SURBA++")

load("tmb.RData")

library(TMB)    

compile("fit.cpp")  

dyn.load("fit")  
#dyn.unload("fit"); 

A = nrow(tmb.data$mat);
Y = ncol(tmb.data$mat);  

parameters <- list(
  log_meanR = 0, 
  log_std_log_R = -1, 
  log_std_index = -1, 
  log_std_log_f=-1,  
  log_std_log_s = c(-1,-2), 
  log_std_pe=log(0.1),
  logit_log_R = 0,   
  log_No = rep(0,A-1),
  log_Rec_dev = rep(0,Y),
  log_sp = rep(0,A-2), # age 6 fixed at one so not estimated;
  log_f = rep(-0.5,Y),
  pe=matrix(0,nrow=A,ncol=Y,byrow=T) 
)   

parameters.L <- list( 
  log_meanR = -1, 
  log_std_log_R = -5, 
  log_std_index = -5, 
  log_std_log_f=-5,  
  log_std_log_s = c(-5,-5), 
 # log_std_log_s = c(-5), 
  log_std_pe=-5,  
  logit_log_R = -10,      
  log_No = rep(-5,A-1)
)   

parameters.U <- list(
  log_meanR = 15, 
  log_std_log_R = 0, 
  log_std_index = 0, 
  log_std_log_f=4,  
  log_std_log_s = c(Inf,Inf),
  #log_std_log_s = c(Inf),   
  log_std_pe=0,   
  logit_log_R = 10,  
  log_No = rep(25,A-1)   
)

lower = unlist(parameters.L);
upper = unlist(parameters.U);

#map=list(log_std_log_s=factor(c('comb','comb')))

## random effects;
rname = c("log_Rec_dev","log_sp","log_f","pe")

 
obj <- MakeADFun(tmb.data,parameters,random=rname,DLL="fit",
  inner.control=list(maxit=100,trace=T)) 
  
length(lower)
length(upper)
length(obj$par)   

obj$gr(obj$par) 


opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
control = list(trace=0,iter.max=100,eval.max=10000))

opt$message
obj$gr(opt$par)

rep = obj$report()

sd.rep<-sdreport(obj)

save.image(file="fit.RData")
