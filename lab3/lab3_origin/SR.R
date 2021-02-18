setwd("C:\\home\\CADIGAN\\GradProgram\\2021\\F6005\\lab3")

library(nlme)  
library(stargazer)

fname = c('SR.dat')
sr.data = read.table(fname,header=T)
sr.data$logrec = log(sr.data$rec)
sr.data$logssb = log(sr.data$ssb)

init.Rmax = 150
init.S50 = 50
init.alpha= init.Rmax
init.beta = init.S50


BH.fit <- nls(logrec ~ log(alpha) + logssb - log(beta + ssb),
           algorithm="port",lower=c(0,0),data=sr.data,
           start = list(beta = init.beta,alpha = init.alpha))

        
BH.arfit <- gnls(logrec ~ log(alpha) + logssb - log(beta + ssb),data=sr.data,
              start = list(beta = init.beta,alpha = init.alpha), 
              correlation = corAR1(form=~year))
              

library(TMB)

## autocorrelated errors;

compile("fit.cpp")  

dyn.load("fit") 
#dyn.unload("fit")

tmb.data = list(
  ssb=sr.data$ssb,    
  rec = sr.data$rec, 
  lssb=log(sr.data$ssb),
  lrec = log(sr.data$rec)    
)

parameters <- list( 
  lalpha = log(150),
  lbeta = log(50),
  lsd_lrec_me = log(0.1),
  logit_ar_lrec_me = 0
) 
 
obj <- MakeADFun(tmb.data,parameters,DLL="fit",
  inner.control=list(maxit=500,trace=T));

obj$gr(obj$par)

opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=1))

obj$gr(opt$par)

rep = obj$report() 
sd.rep<-sdreport(obj) 

tab=cbind(sd.rep$value,sd.rep$sd)
colnames(tab)=c('Est','Sd.err') 

stargazer(tab, type = "html",  out="out.doc",align = TRUE,   
 title="TMB BH AR errors", single.row=TRUE,
 summary=FALSE)
 
 
##########  RW in log alpha ;

n = nrow(sr.data)
 
compile("fit_saoRW.cpp")  

dyn.load("fit_saoRW") 
#dyn.unload("fit_saoRW")

tmb.data = list(
  ssb=sr.data$ssb,    
  rec = sr.data$rec, 
  lssb=log(sr.data$ssb),
  lrec = log(sr.data$rec)    
)  

parameters <- list( 
  lsao = log(3),
  lbeta = log(50),
  lsd_lrec_me = log(0.1), 
  lsd_lsao_dev = log(0.1),
  lsao_dev = rep(0,n)
)


#map <- list(lsd_lrec_me = factor(NA))
 
obj <- MakeADFun(tmb.data,parameters,DLL="fit_saoRW",
  random=c("lsao_dev"),inner.control=list(maxit=500,trace=T));

obj$gr(obj$par)

opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=1))

obj$gr(opt$par)

rep = obj$report() 
sd.rep<-sdreport(obj)

sname=names(sd.rep$value)

ind = sname%in%c("lsao","beta","sd_lrec_me","sd_lsao_dev") 

tab=cbind(sd.rep$value[ind],sd.rep$sd[ind])
colnames(tab)=c('Est','Sd.err') 

stargazer(tab, type = "html",  out="out_saoRW.doc",align = TRUE,   
 title="TMB BH Sao RW", single.row=TRUE,
 summary=FALSE)
 
ind = sname%in%c("lsao_ts")
pdat = data.frame(year=sr.data$year,lsao=sd.rep$value[ind],
L95 = sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind],    
U95 = sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind]) 

re=tab[1,1] + sd.rep$par.random
sdre = sqrt(sd.rep$diag.cov.random)
pdat$L95.1=re-qnorm(0.975)*sdre        
pdat$U95.1=re+qnorm(0.975)*sdre

gname = "saoRW.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)

ylim = range(pdat$L95,pdat$U95)
par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(lsao~year, data=pdat,ylab = '', xlab='',las=1,type='l',lwd=2,ylim=ylim)
  mtext(side=1,line=1.7,"Year") 
  mtext(side=2,line=2.3,"ln(Sao)")
  
  sx = pdat$year
  low = pdat$L95
  high = pdat$U95
  polygon(c(sx,rev(sx)),c(low,rev(high)),
#  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = 'yellow',border = NA)
  
  low = pdat$L95.1
  high = pdat$U95.1
  polygon(c(sx,rev(sx)),c(low,rev(high)),
#  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = 'goldenrod',border = NA)
  
  lines(pdat$year,pdat$lsao,lty=1,lwd=2)
  
dev.off()

#################  profile rec me ####################;


map <- list(lsd_lrec_me = factor(NA))

psd = seq(0.05,0.45,by=0.01)
nsd=length(psd)
pfit = rep(NA,nsd)

for(i in 1:nsd){
  parameters$lsd_lrec_me = log(psd[i])
  
  obj <- MakeADFun(tmb.data,parameters,DLL="fit_saoRW",map=map,
    random=c("lsao_dev"),inner.control=list(maxit=500,trace=F));

  optp<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=0))
  pfit[i] = optp$objective
}

pfit = pfit-opt$objective

gname = "saoRW_profile.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)

  par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(pfit~psd, data=pdat,ylab = '', xlab='',las=1,type='l',lwd=2)
  mtext(side=1,line=1.7,substitute(sigma[epsilon]),las=0)  
  mtext(side=2,line=2.3,substitute(paste(Delta," nll")),las=0) 
  abline(h=qchisq(0.95,1)/2,lty=2)
  
dev.off()


########## RK  RW in log alpha ;

compile("fit_RK_RW.cpp")  

dyn.load("fit_RK_RW") 
#dyn.unload("fit_RK_RW")

tmb.data = list(
  ssb=sr.data$ssb,    
  rec = sr.data$rec, 
  lssb=log(sr.data$ssb),
  lrec = log(sr.data$rec)    
)  

parameters <- list( 
  lalpha = log(3),
  lbeta = log(0.01),
  lsd_lrec_me = log(0.05), 
  lsd_lalpha_dev = log(0.5),
  lalpha_dev = rep(0,n)
) 
 
obj <- MakeADFun(tmb.data,parameters,DLL="fit_RK_RW",
  random=c("lalpha_dev"),inner.control=list(maxit=500,trace=T));

opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=1))

obj$gr(opt$par)

rep = obj$report() 
sd.rep<-sdreport(obj)

sname=names(sd.rep$value)

ind = sname%in%c("alpha_main","beta","sd_lrec_me","sd_lalpha_dev") 

tab=cbind(sd.rep$value[ind],sd.rep$sd[ind])
colnames(tab)=c('Est','Sd.err') 

stargazer(tab, type = "html",  out="out_RKsaoRW.doc",align = TRUE,   
 title="TMB RK Sao RW", single.row=TRUE,
 summary=FALSE)
 
ind = sname%in%c("lsao_ts")
pdat = data.frame(year=sr.data$year,lsao=sd.rep$value[ind],
L95 = sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind],    
U95 = sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind]) 

re=opt$par[1] + sd.rep$par.random
sdre = sqrt(sd.rep$diag.cov.random)
pdat$L95.1=re-qnorm(0.975)*sdre        
pdat$U95.1=re+qnorm(0.975)*sdre

gname = "RKsaoRW.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)

ylim = range(pdat$L95,pdat$U95)
par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(lsao~year, data=pdat,ylab = '', xlab='',las=1,type='l',lwd=2,ylim=ylim)
  mtext(side=1,line=1.7,"Year") 
  mtext(side=2,line=2.3,"ln(Sao)")
  
  sx = pdat$year
  low = pdat$L95
  high = pdat$U95
  polygon(c(sx,rev(sx)),c(low,rev(high)),
#  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = 'yellow',border = NA)
  
  lines(pdat$year,pdat$lsao,lty=1,lwd=2)
  
dev.off()

#################  profile rec me ####################;


map <- list(lsd_lrec_me = factor(NA))

psd = seq(0.01,0.6,by=0.01)
nsd=length(psd)
pfit = rep(NA,nsd)

for(i in 1:nsd){
  parameters$lsd_lrec_me = log(psd[i])
  
  obj <- MakeADFun(tmb.data,parameters,DLL="fit_RK_RW",map=map,
    random=c("lalpha_dev"),inner.control=list(maxit=500,trace=F));

  optp<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=0))
  pfit[i] = optp$objective
}

pfit = pfit-opt$objective

gname = "RKsaoRW_profile.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)

  par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(pfit~psd, data=pdat,ylab = '', xlab='',las=1,type='l',lwd=2,ylim=c(0,2))
  mtext(side=1,line=1.7,substitute(sigma[epsilon]),las=0)  
  mtext(side=2,line=2.3,substitute(paste(Delta," nll")),las=0) 
  abline(h=qchisq(0.95,1)/2,lty=2)
  
dev.off()


## fixed rec ME CV;
map <- list(lsd_lrec_me = factor(NA))

obj <- MakeADFun(tmb.data,parameters,DLL="fit_RK_RW",map=map,
  random=c("lalpha_dev"),inner.control=list(maxit=500,trace=T));

opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=1))

obj$gr(opt$par)

rep = obj$report() 
sd.rep<-sdreport(obj)

sname=names(sd.rep$value)

ind = sname%in%c("alpha_main","beta","sd_lrec_me","sd_lalpha_dev") 

tab=cbind(sd.rep$value[ind],sd.rep$sd[ind])
colnames(tab)=c('Est','Sd.err') 

stargazer(tab, type = "html",  out="out_RKsaoRW_recMEfixed.doc",align = TRUE,   
 title="TMB RK Sao RW", single.row=TRUE,
 summary=FALSE)
 
ind = sname%in%c("lsao_ts")
pdat = data.frame(year=sr.data$year,lsao=sd.rep$value[ind],
L95 = sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind],    
U95 = sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind]) 

gname = "RKsaoRW_recMEfixed.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)

ylim = range(pdat$L95,pdat$U95)
par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(lsao~year, data=pdat,ylab = '', xlab='',las=1,type='l',lwd=2,ylim=ylim)
  mtext(side=1,line=1.7,"Year") 
  mtext(side=2,line=2.3,"ln(Sao)")
  
  sx = pdat$year
  low = pdat$L95
  high = pdat$U95
  polygon(c(sx,rev(sx)),c(low,rev(high)),
#  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = 'yellow',border = NA)
  
  lines(pdat$year,pdat$lsao,lty=1,lwd=2)
  
dev.off()

###############  Added AR(1) alpha process Dec 10 2018;



########## RK  RW in log alpha ;

compile("fit_RK_AR1.cpp")  

dyn.load("fit_RK_AR1") 
#dyn.unload("fit_RK_AR1")

tmb.data = list(
  ssb=sr.data$ssb,    
  rec = sr.data$rec, 
  lssb=log(sr.data$ssb),
  lrec = log(sr.data$rec)    
)  

parameters <- list( 
  lalpha = log(3),
  lbeta = log(0.2),
  lsd_lrec_me = log(0.3), 
  lsd_lalpha_dev = log(0.1),  
  logit_phi = 0,
  lalpha_dev = rep(0,n)
)

map = list(lsd_lrec_me=factor(NA)) 
 
obj <- MakeADFun(tmb.data,parameters,DLL="fit_RK_AR1",map=map,
  random=c("lalpha_dev"),inner.control=list(maxit=500,trace=T));

obj$gr(obj$par)

opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=1))

obj$gr(opt$par)

rep = obj$report() 

plot(tmb.data$rec)
lines(rep$mu)


sd.rep<-sdreport(obj)

sname=names(sd.rep$value)

ind = sname%in%c("alpha_main","beta","sd_lrec_me","sd_lalpha_dev") 

tab=cbind(sd.rep$value[ind],sd.rep$sd[ind])
colnames(tab)=c('Est','Sd.err') 

stargazer(tab, type = "html",  out="out_RKsaoRW.doc",align = TRUE,   
 title="TMB RK Sao RW", single.row=TRUE,
 summary=FALSE)
 
ind = sname%in%c("lsao_ts")
pdat = data.frame(year=sr.data$year,lsao=sd.rep$value[ind],
L95 = sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind],    
U95 = sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind]) 

re=opt$par[1] + sd.rep$par.random
sdre = sqrt(sd.rep$diag.cov.random)
pdat$L95.1=re-qnorm(0.975)*sdre        
pdat$U95.1=re+qnorm(0.975)*sdre

gname = "RKsaoAR.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)

ylim = range(pdat$L95,pdat$U95)
par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(lsao~year, data=pdat,ylab = '', xlab='',las=1,type='l',lwd=2,ylim=ylim)
  mtext(side=1,line=1.7,"Year") 
  mtext(side=2,line=2.3,"ln(Sao)")
  
  sx = pdat$year
  low = pdat$L95
  high = pdat$U95
  polygon(c(sx,rev(sx)),c(low,rev(high)),
#  col = rgb(1,1,0, alpha=0.2, maxColorValue = 1),border = NA)
  col = 'yellow',border = NA)
  
  lines(pdat$year,pdat$lsao,lty=1,lwd=2)
  
dev.off()
  
  
   

