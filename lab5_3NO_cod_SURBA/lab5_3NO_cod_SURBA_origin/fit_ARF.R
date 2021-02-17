
setwd("C:\\home\\CADIGAN\\GradProgram\\2021\\F6005\\lab5\\3NO_cod_SURBA++")

load("tmb.RData")

library(TMB)    

compile("fit_ARF.cpp")  

dyn.load("fit_ARF")  
#dyn.unload("fit_ARF"); 

A = nrow(tmb.data$mat);
Y = ncol(tmb.data$mat);      

parameters <- list(
  log_meanR = 0, 
  log_std_log_R = -1, 
  log_std_index = -1, 
  log_meanF=log(0.2),  
  log_std_log_F = log(0.1), 
  logit_log_F = c(log(0.9/0.1),log(0.9/0.1)),
  log_std_pe=log(0.1),
  logit_log_R = 0,   
  log_No = rep(0,A-1),
  log_Rec_dev = rep(0,Y),
  log_F_dev=matrix(0,nrow=A-1,ncol=Y,byrow=T),
  pe=matrix(0,nrow=A,ncol=Y,byrow=T) 
)   

parameters.L <- list( 
  log_meanR = -Inf, 
  log_std_log_R = log(0.001), 
  log_std_index = log(0.001), 
  log_meanF=-Inf,  
  log_std_log_F = log(0.001), 
#  logit_log_F = c(-10,-10),   
  logit_log_F = -10,
#  log_std_pe=log(0.1),
  logit_log_R = -10,   
  log_No = rep(-10,A-1)
)   

parameters.U <- list(
  log_meanR = Inf, 
  log_std_log_R = log(2), 
  log_std_index = log(2), 
  log_meanF=Inf,  
  log_std_log_F = log(2), 
#  logit_log_F = c(10,10), 
  logit_log_F = 10,
#  log_std_pe=log(1),
  logit_log_R = 10,   
  log_No = rep(10,A-1)
)

lower = unlist(parameters.L);
upper = unlist(parameters.U);

map=list(
  log_std_pe = factor(NA),
  logit_log_F = factor(c('est',NA)),
  pe=factor(matrix(NA,nrow=A,ncol=Y,byrow=T))
  )

## random effects;
rname = c("log_Rec_dev","log_F_dev")

obj <- MakeADFun(tmb.data,parameters,random=rname,map=map,DLL="fit_ARF",
  inner.control=list(maxit=100,trace=T)) 
  
length(lower)
length(upper)
length(obj$par) 

#obj$env$parList()  

#obj$gr(obj$par) 


opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
control = list(trace=0,iter.max=100,eval.max=10000))

opt$message
obj$gr(opt$par)


cbind(lower,opt$par,upper)


rep = obj$report()

sd.rep<-sdreport(obj)

save.image(file="fit_ARF.RData")
