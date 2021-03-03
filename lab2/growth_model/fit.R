
library(TMB)   
library(xtable)    
library(lattice) 
library(latticeExtra) 

setwd("C:\\home\\CADIGAN\\stocks\\models\\redfish\\growth")

load("tmb.RData")
load("data.RData")

compile("fit.cpp")  

dyn.load("fit") 
#dyn.unload("fit")

nstocks = length(unique(tmb.data$istock))   
ngp = length(unique(tmb.data$ui[,1])) 
nobs = length(tmb.data$igroup)           
npred = length(tmb.data$p_igroup)
Linf = 50
k = 0.2
lo = 0.01 
      
parameters <- list( 
  log_Linfp=log(Linf),
  log_kp = log(k),    
  log_deltap = 0, 
  log_Linfp_species=0.35,
  log_kp_species = 0,
  log_Linfp_sex=0,
  log_kp_sex = 0,     
  log_std_log_Linf = log(0.1),  
  log_std_log_k = log(0.1),      
  log_std_log_delta = log(0.1),    
  log_std_me_len = log(0.017), 
  log_std_me_age = log(0.06),     
  log_std_pe = rep(log(0.1),nstocks),
  logit_corr = 0,
  log_alpha = rep(log(7),nstocks),
  log_beta = rep(log(1),nstocks),      
  log_Linf_stock=rep(0,nstocks),
  log_k_stock = rep(0,nstocks),
  log_delta_stock = rep(0,ngp), 
  pe = rep(0,nobs),       
  pe_pred = rep(0,npred),   
  true_age=0.5*tmb.data$age + 0.5*mean(tmb.data$age)       
)           
     
parameters.L <- list(
  log_Linfp = -Inf,
  log_kp = -Inf,  
  log_deltap = -Inf, 
  log_Linfp_species=-Inf,
  log_kp_species = -Inf,
  log_Linfp_sex=-Inf,
  log_kp_sex = -Inf,       
  log_std_log_Linf = log(1e-6),  
  log_std_log_k = log(1e-6),     
  log_std_log_delta = log(1e-6),    
#  log_std_me_len = log(1e-6),      
#  log_std_me_age = log(1e-6),     
  log_std_pe = log(1e-6),
  logit_corr = -Inf,
  log_alpha = -Inf,
  log_beta = -Inf
#  log_alpha = rep(-Inf,nstocks),
#  log_beta = rep(-Inf,nstocks)         
  )    
 
parameters.U <- list( 
  log_Linfp = Inf,
  log_kp = Inf,   
  log_deltap = Inf, 
  log_Linfp_species=Inf,
  log_kp_species = Inf,
  log_Linfp_sex=Inf,
  log_kp_sex = Inf,   
  log_std_log_Linf = log(5),  
  log_std_log_k = log(5),     
  log_std_log_delta = log(5),   
#  log_std_me_len = log(10),      
#  log_std_me_age = log(1),    
  log_std_pe = log(1),
  logit_corr = Inf,  
  log_alpha = Inf,
  log_beta = Inf
#  log_alpha = rep(Inf,nstocks),
#  log_beta = rep(Inf,nstocks)      
)

lower = unlist(parameters.L);
upper = unlist(parameters.U);

map <- list(
log_std_me_len = factor(NA),
log_std_me_age = factor(NA),  
log_alpha=factor(rep('A',nstocks)),  
log_beta=factor(rep('A',nstocks)), 
log_std_pe=factor(rep('A',nstocks)))
                                        
obj <- MakeADFun(tmb.data,parameters,map=map,  
random=c("log_Linf_stock","log_k_stock","log_delta_stock","pe","pe_pred","true_age"),
DLL="fit",inner.control=list(maxit=100,trace=T))

length(lower)
length(upper)
length(obj$par)

opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
control = list(trace=0,iter.max=1000,eval.max=2500))

rep=obj$rep() 

sd.rep<-sdreport(obj,getReportCovariance = FALSE)

save(opt, rep, sd.rep, parameters,parameters.L,parameters.U,map, file = "fit.RData")

temp = t(as.matrix(c(opt$objective,opt$par),nrows=1,byrow=T))
colnames(temp) = c("nll",names(opt$par))

tname = c("Fit")    
xdigits = rep(4,length(temp)+1); #a digit for rowname, and all rows;
xalign = c('c',rep('r',length(temp))); # alignment 'l','c', or 'r';

print(xtable(temp,digits=xdigits,caption=tname,align=xalign),type='html',
file='fit.doc',caption.placement="top") 

par_sd = sqrt(diag(sd.rep$cov.fixed))
corr_matrix = diag(1/par_sd)%*%sd.rep$cov.fixed%*%diag(1/par_sd)

cnames = c("Linf","k","gamma","Linf_sp","k_sp",
"Linf_x","k_x","sigma_inf","sigma_k","sigma_gamma",
"sigma_pe","corr*","alpha","beta") 

rownames(corr_matrix) = cnames     
colnames(corr_matrix) = cnames

require(corrplot)
  

jpeg(file="corr.jpeg",width=5,height=5,units='in',res=300)

corrplot.mixed(corr_matrix, tl.srt=45,tl.cex=0.6,cl.cex=0.6,tl.col='black',tl.pos='lt',diag='u',number.cex=0.6,title="Correlation in log parameter estimates (*logit)",
mar=c(0,0,1,0))

dev.off()


########################################################################;



gp = aggregate(data[,13:16],list(group=data$group),unique)
gp = gp[2:21,]

sex = substr(gp$group,1,1) 
species = substr(gp$group,3,5) 
stock = substr(gp$group,7,11)  
stock.names = as.vector(sort(unique(data$stock)))
#stock.names = gsub('men','nor',stock.names)

value.names = names(sd.rep$value)

ind = value.names=='log_k_stock'
log_kl = rep$log_k_stock + qnorm(0.025)*sd.rep$sd[ind] 
log_ku = rep$log_k_stock + qnorm(0.975)*sd.rep$sd[ind]  
ind = value.names=='log_Linf_stock'
log_Linfl = rep$log_Linf_stock + qnorm(0.025)*sd.rep$sd[ind] 
log_Linfu = rep$log_Linf_stock + qnorm(0.975)*sd.rep$sd[ind]


gname=c('par_vs_k.png')
png(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3,0.5,0.5),mgp=c(2,0.75,0))

pchv = c(1:6,1,2,4,5)
colv = c(rep('red',6),rep('black',4))

plot(rep$log_k_stock,rep$log_Linf_stock,las=1,xlab = "",ylab = '',type = 'p',
pch=pchv,col=colv)

mtext(side=1,line=2,"log(k)")        
mtext(side=2,line=2,expression(paste('log(',L[infinity],')')))

segments(rep$log_k_stock,log_Linfl,rep$log_k_stock,log_Linfu,col=colv)
segments(log_kl,rep$log_Linf_stock,log_ku,rep$log_Linf_stock,col=colv) 

legend("topright",legend=stock.names,bty='n',cex=0.7,pch=pchv,col=colv)

dev.off() 


ox = c(1,11,7,17,2,12,8,18,3,13,4,14,9,19,5,15,10,20,6,16)  
tname = substr(as.character(gp[ox,]$group),3,5) 
#tname = gsub('men','nor',tname)

for(type in 1:3){

if(type==1){ind = value.names=='log_k_gp';gname=c('barplot_k.png');ytext="k"} 
if(type==2){ind = value.names=='log_Linf_gp';gname=c('barplot_Linf.png');ytext=expression(L[infinity])}
if(type==3){ind = value.names=='log_delta_gp';gname=c('barplot_gamma.png');ytext=expression(gamma)}

xl = sd.rep$value[ind][tmb.data$ui[,3]>0]

png(file=gname,width=4,height=3,units='in',res=300)

par(mar=c(3,3,1,0.5),mgp=c(2,0.75,0))

temp = barplot(exp(xl[ox]),names.arg=tname,las=2,col=rep(c('grey','white'),10))
mtext(side=2,line=2,ytext)
abline(v=c(4.9,9.7,12.15,16.95,21.72),lty=1,lwd=4)
x = exp(xl)
if(type==1){ypos = max(x) + 0.06*diff(range(x))} 
if(type==2){ypos = max(x) + 0.07*diff(range(x))}
if(type==3){ypos = max(x) + 0.13*diff(range(x))}

text(temp[1]+1,ypos,'2J3K',xpd=NA)    
text(temp[5]+1,ypos,'3LN',xpd=NA)      
text(temp[9]+0.5,ypos,'3O',xpd=NA)      
text(temp[11]+1.2,ypos,'Unit1',xpd=NA)     
text(temp[15]+1.2,ypos,'Unit2',xpd=NA)     
text(temp[19]+1.2,ypos,'Unit3',xpd=NA)

dev.off()
}


gp = aggregate(data$igroup,list(group=data$group),unique)

res = data.frame(len=data$len,age=data$age,pred=rep$pred,pred_age = rep$true_age,
  resid=rep$resid,std_resid=rep$std_resid,pe=rep$pe,pred_nope=exp(rep$log_pred_nope),
  stock=data$stock,species=data$species,sex=data$sex)
res$group = as.factor(gp$group[data$igroup+1])  
#res$group = gsub('men','nor',res$group)        
#res$stock = gsub('men','nor',res$stock) 
#res$species = gsub('men','nor',res$species)

p_res = data.frame(age=tmb.data$p_age,pred=rep$p_pred)

value.names = names(sd.rep$value)
ind = value.names=='p_log_pred'
log_len_est = rep$p_log_pred
log_len_sd = sd.rep$sd[ind]

p_res$lower = exp(log_len_est + qnorm(0.025)*log_len_sd); 
p_res$upper = exp(log_len_est + qnorm(0.975)*log_len_sd);   
p_res$group = as.factor(gp$group[tmb.data$p_igroup+1]) 
#p_res$group = gsub('men','nor',p_res$group)

ox = order(res$group,res$age)
res=res[ox,]

ylim = range(data$len,p_res$lower,p_res$upper)  
xlim = range(data$age)
xttl <- list(c("Age"), cex=1)
yttl <- list(c("Length (cm)"), cex=1)
stripttl <- list(cex=0.2,font=1) 
#ax <- list(cex=0.75,relation="free",tck=0.3,y=list(rot=0),x=list(draw=FALSE)) 
ax <- list(cex=0.75,alternating=FALSE,tck=0.3,y=list(rot=0))

col.vec = c('black','black') 
my.key <- list(space = "top",columns=2,padding.text=0,height=0.5,
              border = FALSE, size = 1.5,between=0.5,
              lines = list(lty = c(0,1),col=col.vec,lwd=2,type='b',pch=c(3,46)),
              text = list(c("Observed","Predicted")))
theme.novpadding <-
   list(layout.heights =
        list(top.padding = 0,
 	    main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0,
 	    xlab.key.padding = 0.5,
 	    key.sub.padding = 0,
 	    bottom.padding = 0),
        layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 0,
 	    right.padding = 0))


plot1 = xyplot(len~age|group,cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,7,1),las=1,
col=col.vec,par.strip.text = list(cex=0.6),
#par.settings = list(layout.heights = list(strip = .8)))  
par.settings = theme.novpadding)

plot2 = xyplot(pred+lower+upper~age|group,cex=0.5,
data=p_res, type=c("l"),lty=1,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,7,1),las=1,
col=c(grey(0.6),grey(0.6),grey(0.6)),par.strip.text = list(cex=0.6),
#par.settings = list(layout.heights = list(strip = .8)))             
par.settings = theme.novpadding)

gname=c('obs_vs_pred.png')
png(file=gname,width=4,height=5,units='in',res=300)
  plot3 = plot1+plot2
  print(plot3) 
dev.off() 

## stock comparison plot;
theme.novpadding <-
   list(layout.heights =
        list(top.padding =0.1,
 	    main.key.padding = 0,
 	    key.axis.padding = -1,
 	    axis.xlab.padding = 0,
 	    xlab.key.padding = 0.5,
 	    key.sub.padding = 0,
 	    bottom.padding = 0),
        layout.widths =
        list(left.padding = 0.5,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 0,
 	    right.padding = 0))

gp1 = aggregate(data[,13:16],list(group=data$group),unique)
stock1 = substr(gp1$group,7,11)
ss = paste(c('Juvenile','Male','Female')[tmb.data$ui[,3]+1],', ',
c('S. fasciatus','S. mentella')[tmb.data$ui[,2]+1],sep='')

p_res$ss = ss[tmb.data$p_igroup+1]
p_res$stock = stock1[tmb.data$p_igroup+1]

p_res1 = subset(p_res,group!="0_fas_Unit3")
p_res1$fss = factor(p_res1$ss)          
p_res1$fstock = factor(p_res1$stock)

xlim = range(p_res$age)
ylim = range(p_res$pred)
ylim[2]=ylim[2]+1

#col.vec = c('black',grey(0.7),'black',grey(0.7),'black',grey(0.7))
col.vec = c('black','black','black',grey(0.7),grey(0.7),grey(0.7))
my.key <- list(space = "top",columns=3,padding.text=1,height=1,
              border = FALSE, size = 4,between=0.5,cex=0.7,
              lines = list(lty = c(1,2,3,1,2,3),col=col.vec,lwd=2),
              text = list(levels(p_res1$fstock))
)

plot1 = xyplot(pred~age|fss, group=fstock,cex=0.5,
key = my.key,
data=p_res1, type=c("l"),lty=c(1,2,3,1,2,3),xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(2,2,1),las=1,
col=col.vec,par.strip.text = list(cex=0.6),
#par.settings = list(layout.heights = list(strip = .8))) 
par.settings = theme.novpadding,
    panel = function(...) {
      panel.xyplot(...)   
      panel.abline(h=c(20,30,40), lty = 2, col = grey(0.8)) 
})
#  print(plot1)

gname=c('compare_curves.png')
png(file=gname,width=4,height=4,units='in',res=300)
  print(plot1)
 dev.off() 

#################### observed versus predicted ##############;

plot1 = xyplot(age~pred_age|group,cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,7,1),las=1,
col=col.vec,par.strip.text = list(cex=0.6), 
par.settings = theme.novpadding,
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(a=0,b=1, lty = "solid", col = "black") 
})

gname=c('obs_vs_pred_age.png')
png(file=gname,width=4,height=5,units='in',res=300)
  print(plot1) 
dev.off()


gname=c('obs_vs_pred_age_all.png')
png(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3,0.5,0.5),mgp=c(2,0.75,0))

plot(res$pred_age,res$age,las=1,xlab = "",ylab = '',type = 'p',cex=0.5,pch=3)

mtext(side=1,line=2,"Predicted Age")        
mtext(side=2,line=2,"Measured Age")

abline(a=0,b=1,lty=2)

dev.off()


gname=c('density_age_all.png')
png(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3,0.5,0.5),las=1,mgp=c(2,0.5,0))

hist(res$age,breaks=50,prob=T,main='',xlab='Age',ylab='Density')
box(lty=1)
lines(density(res$pred_age))

dev.off()                                      

######### residual plots ####################;
theme.novpadding <-
   list(layout.heights =
        list(top.padding = 0,
 	    main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0,
 	    xlab.key.padding = 0.5,
 	    key.sub.padding = 0,
 	    bottom.padding = 0),
        layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 0,
 	    right.padding = 0))


ylim = range(res$std_resid)  
xlim = range(res$age)
xttl <- list(c("Age"), cex=1)
yttl <- list(c("Standardized Residual"), cex=1)
stripttl <- list(cex=0.2,font=1) 
ax <- list(cex=0.75,relation="free",tck=0.3,y=list(rot=0))

plot1 = xyplot(std_resid~age|group,cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,7,1),las=1,
par.strip.text = list(cex=0.6),
par.settings = theme.novpadding,
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(h=0, lty = "solid", col = "black") 
      panel.loess(x, y,col='black',...)
}) 

gname=c('std_resid.png')
png(file=gname,width=5,height=6,units='in',res=300)
#par(mar=c(2.5,3,0.5,0.5),mgp=c(2,0.75,0))
  print(plot1) 
dev.off() 

res$year=data$year
ylim = range(res$std_resid)  
xlim = range(res$year)
xttl <- list(c("Year"), cex=1)
yttl <- list(c("Standardized Residual"), cex=1)
stripttl <- list(cex=0.2,font=1) 
ax <- list(cex=0.75,relation="free",tck=0.3,y=list(rot=0))

plot1 = xyplot(std_resid~year|group,cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,7,1),las=1,
par.strip.text = list(cex=0.6),
par.settings = theme.novpadding,
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(h=0, lty = "solid", col = "black") 
      #panel.loess(x, y,col='black',...)
})


gname=c('std_resid_year.png')
png(file=gname,width=5,height=6,units='in',res=300)
#par(mar=c(2.5,3,0.5,0.5),mgp=c(2,0.75,0))
  print(plot1) 
dev.off() 

ylim = range(res$std_resid)  
xlim = range(res$age)
xttl <- list(c("Age"), cex=1)
yttl <- list(c("Standardized Residual"), cex=1)
stripttl <- list(cex=0.2,font=1) 
ax <- list(cex=0.75,relation="free",tck=0.3,y=list(rot=0))

plot1 = xyplot(std_resid~age|factor(stock),cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(2,5,1),las=1,
par.strip.text = list(cex=0.6),
par.settings = theme.novpadding,
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(h=0, lty = "solid", col = "black") 
      panel.loess(x, y,col='black',...)
}) 

gname=c('std_resid_sp.png')
png(file=gname,width=5,height=6,units='in',res=300)
#par(mar=c(2.5,3,0.5,0.5),mgp=c(2,0.75,0))
  print(plot1) 
dev.off() 

ylim = range(res$std_resid)  
xlim = range(res$year)
xttl <- list(c("Year"), cex=1)
yttl <- list(c("Standardized Residual"), cex=1)
stripttl <- list(cex=0.2,font=1) 
ax <- list(cex=0.75,relation="free",tck=0.3,y=list(rot=0))

plot1 = xyplot(std_resid~year|factor(stock),cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(2,5,1),las=1,
par.strip.text = list(cex=0.6),
par.settings = theme.novpadding,
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(h=0, lty = "solid", col = "black") 
     # panel.loess(x, y,col='black',...)
}) 

gname=c('std_resid_sp_year.png')
png(file=gname,width=5,height=6,units='in',res=300)
#par(mar=c(2.5,3,0.5,0.5),mgp=c(2,0.75,0))
  print(plot1) 
dev.off() 


ylim = range(res$resid)  
xlim = range(res$age)
xttl <- list(c("Age"), cex=1)
yttl <- list(c("Residual"), cex=1)
stripttl <- list(cex=0.2,font=1) 
ax <- list(cex=0.75,relation="free",tck=0.3,y=list(rot=0))

plot1 = xyplot(resid~age|group,cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,7,1),las=1,
par.strip.text = list(cex=0.6),
par.settings = theme.novpadding,
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(h=0, lty = "solid", col = "black") 
      panel.loess(x, y,col='black',...)
}) 

gname=c('resid.png')
png(file=gname,width=5,height=6,units='in',res=300)
par(mar=c(2.5,3,0.5,0.5),mgp=c(2,0.75,0))
  print(plot1) 
dev.off() 

gname=c('resid_vs_pred.png')
png(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3.5,0.5,0.5),mgp=c(2,0.75,0))

plot(rep$pred,rep$std_resid,las=1,xlab = "",ylab = '',type = 'p',cex=0.5,pch=3)

mtext(side=1,line=2,"Predicted Length (cm)")        
mtext(side=2,line=2.5,"Standardized Residual")

abline(h=0,lty=2,lwd=2,col=grey(0.6))

lo = loess(y~x,data.frame(y=rep$std_resid,x=rep$pred))
pp = data.frame(x = seq(min(rep$pred),max(rep$pred),by=0.1))
lines(pp$x,predict(lo,newdata=pp),lty=1,lwd=2,col=grey(0.6))

dev.off() 

gname=c('resid_vs_pred_nope.png')
png(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3.5,0.5,0.5),mgp=c(2,0.75,0))

plot(res$pred_nope,rep$std_resid,las=1,xlab = "",ylab = '',type = 'p',cex=0.5,pch=3)

mtext(side=1,line=2,"Predicted Length (cm)")        
mtext(side=2,line=2.5,"Standardized Residual")

abline(h=0,lty=2,lwd=2,col=grey(0.6))

lo = loess(y~x,data.frame(y=rep$std_resid,x=res$pred_nope))
pp = data.frame(x = seq(min(res$pred_nope),max(res$pred_nope),by=0.1))
lines(pp$x,predict(lo,newdata=pp),lty=1,lwd=2,col=grey(0.6))

dev.off()


######### process error plots ####################;

ylim = range(res$pe)  
xlim = range(res$age)
xttl <- list(c("Age"), cex=1)
#yttl <- list(c("Process Error"), cex=1)
yttl <- list(expression(paste(epsilon,' Individual Variation')), cex=1)


stripttl <- list(cex=0.2,font=1) 
ax <- list(cex=0.75,alternating=FALSE,tck=0.3,y=list(rot=0))

plot1 = xyplot(pe~age|group,cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,7,1),las=1,
par.strip.text = list(cex=0.6),
par.settings = theme.novpadding,
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(h=0, lty = "solid", col = "black") 
      panel.loess(x, y,col='black',...)
}) 

gname=c('pe_vs_age_by_group.png')
png(file=gname,width=7,height=10,units='in',res=300)
#par(mar=c(2.5,3,0.5,0.5),mgp=c(2,0.75,0))
  print(plot1) 
dev.off() 

ylim = range(res$pe)  
xlim = range(res$age)
xttl <- list(c("Age"), cex=1)
#yttl <- list(c("Process Error"), cex=1)
stripttl <- list(cex=0.2,font=1) 
ax <- list(cex=0.75,alternating=FALSE,tck=0.3,y=list(rot=0))

plot1 = xyplot(pe~age|factor(stock),cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(2,5,1),las=1,
par.strip.text = list(cex=0.6),
par.settings = theme.novpadding,
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(h=0, lty = "solid", col = "black") 
      panel.loess(x, y,col='black',...)
}) 

gname=c('pe_vs_age_by_stock.png')
png(file=gname,width=4,height=5,units='in',res=300)
par(mar=c(2.5,3,0.5,0.5),mgp=c(2,0.75,0))
  print(plot1) 
dev.off()            


ylim = range(res$pe)  
xlim = range(res$pred)
xttl <- list(c("Predicted Length (cm)"), cex=1)
yttl <- list(c("Process Error"), cex=1)
stripttl <- list(cex=0.2,font=1) 
ax <- list(cex=0.75,relation="free",tck=0.3,y=list(rot=0))

plot1 = xyplot(pe~pred|group,cex=0.5,
data=res, type=c("p"),pch=3,xlab=xttl,lwd=2,ylim=ylim,xlim=xlim,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,7,1),las=1,
par.strip.text = list(cex=0.6),
par.settings = list(layout.heights = list(strip = .8)),
    panel = function(x, y,col,...) {
      panel.xyplot(x, y,col=grey(0.6),...)   
      panel.abline(h=0, lty = "solid", col = "black") 
      panel.loess(x, y,col='black',...)
}) 

gname=c('pe_vs_pred_group.png')
png(file=gname,width=7,height=10,units='in',res=300)
par(mar=c(2.5,3,0.5,0.5),mgp=c(2,0.75,0))
  print(plot1) 
dev.off() 


gname=c('pe_vs_pred.png')
png(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3.5,0.5,0.5),mgp=c(2,0.75,0))

pe = rep$pe
plot(rep$pred,pe,las=1,xlab = "",ylab = '',type = 'p',cex=0.5,pch=3)

mtext(side=1,line=2,"Predicted Length (cm)")        
mtext(side=2,line=2.5,expression(paste(epsilon,' Individual Variation')))

abline(h=0,lty=2,lwd=2,col=grey(0.6))

lo = loess(y~x,data.frame(y=pe,x=rep$pred))
pp = data.frame(x = seq(min(rep$pred),max(rep$pred),by=0.1))
lines(pp$x,predict(lo,newdata=pp),lty=1,lwd=2,col=grey(0.6))

dev.off() 

########## The Keys #################;
## stock comparison plot;
theme.novpadding <-
   list(layout.heights =
        list(top.padding =0.1,
 	    main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0.5,
 	    xlab.key.padding = 0,
 	    key.sub.padding = 0,
 	    bottom.padding = 0.5),
        layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 1,
 	    axis.key.padding = -0.5,
 	    right.padding = 0))

value.names = names(sd.rep$value)

ind = value.names=='p_log_pred'
log_len_est = rep$p_log_pred
log_len_sd = sd.rep$sd[ind]

#gp$group = gsub('men','nor',gp$group)

for(i in 1:20){
  
  ind = tmb.data$p_igroup==i
  lli = log_len_est[ind]
  sdi = log_len_sd[ind]
  agei = tmb.data$p_age[ind]
  nagei=length(agei)
  lc = 0:50
  lc1 = c(-Inf,log(seq(0.5,49.5,by=1)),Inf)
  nlc1 = length(lc1)

  ki = matrix(NA,length(agei),length(lc))
  
  for(j in 1:length(agei)){
    ml = lli[j]
    sl = sdi[j]
#    if(sl>1){sl = sdi[order(sdi)][nagei-1]}
    ki[j,] = pnorm(lc1[2:nlc1],ml,sl) - pnorm(lc1[1:(nlc1-1)],ml,sl)
  }
  
  rownames(ki) = agei
  colnames(ki) = lc
#  round(ki,digits=2)

  p_resi = subset(p_res,group==gp$group[i+1])  


gname=paste('keys\\',as.character(gp$group[i+1]),'.png',sep='')
png(file=gname,width=5,height=5,units='in',res=300)

#par(mar=c(3,3,1,0.5),mgp=c(2,0.75,0))

ct = c(1,0.75,seq(0.5,0,by=-0.01))

plot1 = levelplot(ki, scales=list(tck=0, x=list(rot=90),cex=0.5), col.regions=grey(ct),
par.settings = theme.novpadding,
at=ct, main=as.character(gp$group[i+1]), xlab='Age', ylab='Length (cm)')
plot1 = plot1 + layer(llines(1:length(p_resi$age),p_resi$pred+0.5,col='black',lwd=3))

print(plot1)



dev.off()   
    
ctext = as.character(gp$group[i+1])
xdigits = c(1,rep(2,nagei))
xalign = c('c',rep('r',nagei)); 
gname=paste('keys\\',as.character(gp$group[i+1]),'.doc',sep='')
  
  print(xtable(t(ki),digits=xdigits,caption=ctext,align=xalign),type='html',
  file=gname,caption.placement="top") 
  
}   
    
  