setwd("C:\\home\\CADIGAN\\GradProgram\\2021\\F6005\\3LN_redfish")

load("data\\3LN_redfish.RData") 
load("tmb_data.RData")

library(TMB)    

compile("fit.cpp")  

dyn.load("fit")  
#dyn.unload("fit"); 

## set FRV length pattern in catchability;
L50c = 17
L95c = 22
b1c = log(0.95/0.05)/(L95c-L50c)
boc = -b1c*L50c
ql = 5*exp(boc + b1c*tmb_data$len_mid)/(1+exp(boc + b1c*tmb_data$len_mid))
#plot(tmb_data$len_mid,ql,xlab='length',ylab='Catchability Length Pattern',type='l',lwd=2)

A = tmb_data$A;
Y = tmb_data$Y;
ns = length(tmb_data$sf)

parameters <- list(
  log_meanR = 3, 
  log_std_log_R = log(2), 
  log_std_index = rep(log(1),ns),
  log_std_catch = log(0.5),
  log_std_logq = rep(0,ns), 
  log_std_pe=log(0.1),
  logit_log_R = -5,  
  log_Linf=log(35),
  log_vbk=log(0.35),  
  log_len_o=log(0.5),
  log_cv_len=log(0.2),
  log_std_log_F = log(1),
  logit_F_age = log(0.25/0.75),
  logit_F_year = log(0.25/0.75),
  log_F_main = log(0.01),        
  log_N0 = rep(0,A-1),    
# the random effects; 
  log_Rec_dev = rep(0,Y),   
  log_F_dev=matrix(0,nrow=A-1,ncol=Y,byrow=T), 
  logq = log(matrix(ql,nrow=ns,ncol=length(ql),byrow=T)),
  pe=matrix(0,nrow=A-1,ncol=Y-1,byrow=T) 
)
parameters$log_N0 = parameters$log_meanR -0.2*(1:(A-1))
parameters$logq[3,]=parameters$logq[3,] - log(10000)   

# not estimating correlayion in rec devs or any process error at first
parameters.L <- list( 
  log_meanR = 1, 
  log_std_log_R = log(0.01), 
  log_std_index = rep(log(0.1),ns), 
  log_std_catch = log(0.01),
  log_std_logq = rep(-Inf,ns), 
  log_std_pe=log(0.01),
  logit_log_R = -5,  
  log_Linf=log(20),
  log_vbk=log(0.01),  
  log_len_o=log(1/10000),
  log_cv_len=log(0.01),
  log_std_log_F = log(0.01),
  logit_F_age = -5,
  logit_F_year = -5,
  log_F_main = -Inf,
  log_N0 = rep(-Inf,A-1) 
)   

parameters.U <- list(
  log_meanR = 100, 
  log_std_log_R = log(10), 
  log_std_index = rep(log(10),ns),
  log_std_catch = log(2),
  log_std_logq = rep(Inf,ns),  
  log_std_pe=log(0.5),
  logit_log_R = -5,  
  log_Linf=log(100),
  log_vbk=log(0.8),  
  log_len_o=log(1),
  log_cv_len=log(0.5),
  log_std_log_F = log(5),
  logit_F_age = log(0.95/0.05),
  logit_F_year = log(0.95/0.05),
  log_F_main = log(2),
  log_N0 = rep(Inf,A-1) 
)

lower = unlist(parameters.L);
upper = unlist(parameters.U);

ql.map = 1:length(ql)
ql.map[ql.map>45]=45 
ql.map[ql.map<=3]=NA
m1 = matrix(1:ns,nrow=ns,ncol=length(ql),byrow=F) 
m2 = matrix(ql.map,nrow=ns,ncol=length(ql),byrow=T)
qmap_by_survey = sprintf('%d,%d', m1, m2)
dim(qmap_by_survey) <- dim(m1)
qmap_by_survey[,1:3]=NA


pdat = expand.grid(age=2:tmb_data$A,year=unique(CLc_vec$Year))
pdat$age[pdat$age>9]=9
F_map = paste(pdat$year,pdat$age,sep="_")

log_F_dev_map=matrix(F_map,nrow=A-1,ncol=Y,byrow=F) 
#log_F_dev_map=matrix(NA,nrow=A-1,ncol=Y,byrow=F) 

## first run with log_std_log_F fixed to get starting values;

map=list(  
log_len_o=factor(NA), 
#log_Linf = factor(NA),  
log_F_main = factor(NA),
#log_vbk = factor(NA),
#log_std_log_F = factor(NA),
#logit_F_age=factor(NA),
#logit_F_year=factor(NA), 
log_std_pe=factor(NA),  
logit_log_R = factor(NA),           
pe=factor(matrix(NA,nrow=A-1,ncol=Y-1,byrow=T)),
logq=factor(qmap_by_survey),
log_F_dev=factor(log_F_dev_map)  
)

t.L = parameters.L
t.L$log_len_o=NULL  
#t.L$log_Linf=NULL     
t.L$log_F_main=NULL    
#t.L$log_vbk = NULL
#t.L$log_std_log_F = NULL
#t.L$logit_F_age=NULL
#t.L$logit_F_year=NULL 
t.L$log_std_pe=NULL  
t.L$logit_log_R = NULL  
lower = unlist(t.L);   

t.U = parameters.U
t.U$log_len_o=NULL  
#t.U$log_Linf=NULL   
t.U$log_F_main=NULL      
#t.U$log_vbk = NULL
#t.U$log_std_log_F = NULL
#t.U$logit_F_age=NULL
#t.U$logit_F_year=NULL 
t.U$log_std_pe=NULL  
t.U$logit_log_R = NULL 
upper = unlist(t.U);

## random effects;
rname = c("log_Rec_dev","log_F_dev","logq")
 
obj <- MakeADFun(tmb_data,parameters,random=rname,DLL="fit",map=map,
  inner.control=list(maxit=500,trace=F),random.start = expression(last.par[random]))
  #inner.control=list(maxit=500,trace=F))
  
length(obj$par)
length(lower)
length(upper)
  
obj$fn(obj$par)
  
#obj$gr(obj$par) 

# "FallRV" is=0 "SprgRV" is=1 "Spsh3L" is=2 "Spsh3N" is=3  this is the connection between is and sname

opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
control = list(trace=0,iter.max=500,eval.max=10000))

opt<-nlminb(opt$par,obj$fn,obj$gr,lower=lower,upper=upper,
control = list(trace=0,iter.max=500,eval.max=10000))

cbind(lower,opt$par,upper)   

rep = obj$report()
   
sd.rep<-sdreport(obj)  

save.image(file="fit.RData")

###########  Do the Plotting ##################

source("plots.R")


## MISC code #############

ind = abs(CLc.vec$resid) > 4
CLc.vec[ind,]

ind1 = CLc.vec$Year==2019    
CLc.vec[ind1,]

