#source("xxx",local=F)

setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\Lecture4\\3LN_redfish")  

require("lattice")    
require("xtable")

load("fit.RData")
setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\Lecture4\\3LN_redfish\\data")
load("tmb.RData") 

index.dat <- read.table("indices.txt",header=TRUE, fill = FALSE)


setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\Lecture4\\3LN_redfish\\figs")


################# CPUE plot ############################;

uindex = unique(index.dat$name)
posv = c('topright',rep("topleft",7))

for(i in 1: length(uindex)){

u = uindex[i]
pos = posv[i]

scale=1000
index.name = "Index (000s)"
#index.name = "Index"
if(u=='CPUE'){scale=1;index.name="Index"}

ind = index.dat$name==u
year = index.dat$year[ind]
index = index.dat$index[ind]/scale
Eindex = rep$Eindex[ind]/scale
std.resid = rep$std_resid[ind]

jpeg(file=paste('fit_',u,'.jpeg',sep=''),width=4,height=4,units='in',res=300)

par(mfcol=c(3,1),oma=c(1,1,1.5,0.5),mar=c(3,4,0,0),cex.axis=1.2,las=1)

ylim=range(index,Eindex)
xlim=range(index.dat$year)

xi = tmb.data$year
yi=rep(NA,length(xi));pyi=yi;pres=yi
ind1 = xi%in%year
yi[ind1]=index 
pyi[ind1]=Eindex 
pres[ind1]=std.resid

plot(xi,yi,type='b',lwd=2,pch=19,xlab='',ylab='',ylim=ylim,xlim=xlim)
lines(xi,pyi,type='l',lwd=2,col='red')
mtext(side=1,line=2,outer=F,"Year")            
mtext(side=2,line=3.5,outer=F,index.name,las=0)
legend(pos,bty='n',lwd=2,pch=c(19,NA),lty=1,c("Observed","Predicted"),pt.cex=1.5,col=c('black','red'))

plot(xi,pres,type='b',xlab='',ylab='',xlim=xlim)
abline(h=0,lty=2)
mtext(side=1,line=2,outer=F,"Year")            
mtext(side=2,line=3.5,outer=F,"Std. residual",las=0)

plot(pyi,pres,type='p',xlab='',ylab='')
mtext(side=1,line=2.5,outer=F,"Predicted")            
mtext(side=2,line=3.5,outer=F,"Std. residual",las=0)

mtext(side=3,line=0,outer=T,u,las=0)

dev.off() 
}

##############  biomass and harvest plot ################;
 
value.names = names(sd.rep$value)

ind = value.names == "log_B";
biomass = exp(sd.rep$value[ind]);

stdi =  sd.rep$sd[ind]
biomass.L95 = exp(sd.rep$value[ind] - qnorm(0.975)*stdi); 
biomass.U95 = exp(sd.rep$value[ind] + qnorm(0.975)*stdi);

ind = value.names == "log_H";
exploit = exp(sd.rep$value[ind]);

stdi =  sd.rep$sd[ind]
exploit.L95 = exp(sd.rep$value[ind] - qnorm(0.975)*stdi); 
exploit.U95 = exp(sd.rep$value[ind] + qnorm(0.975)*stdi);


jpeg(file="pop.jpeg",width=4,height=5,units='in',res=1200)

par(mfcol=c(2,1),oma=c(1,1.5,0,0.5),mar=c(2,2,1,0),cex.axis=1,las=1)
xlim <- c(min(tmb.data$year),max(tmb.data$year))
ylim <- c(min(biomass.L95),max(biomass.U95))
ylim[2]=1500

n=length(tmb.data$year)

ytext <- c("Exploitable biomass (Kt)")
plot(xlim,ylim,xlab='',ylab='',type='n',cex.lab=1,lwd=2)

sx = tmb.data$year
low = biomass.L95
high = biomass.U95
polygon(c(sx,rev(sx)),c(low,rev(high)),
  col = 'grey',border = NA)       
  box(lty=1)   
lines(tmb.data$year,biomass,type='l',lty=1,lwd=2)  
  
abline(h=rep$Bmsy,lwd=2,col='red')    
abline(h=0.5*rep$Bmsy,lwd=1,col='red')    
mtext(side=2,line=2.5,outer=F,ytext,cex=1,las=0)
p1 = rep$Bmsy + 0.02*diff(range(ylim))
text(tmb.data$year[n]-2,p1,"Bmsy",cex=0.75)    
p1 = 0.5*rep$Bmsy + 0.02*diff(range(ylim))
text(tmb.data$year[n]-3,p1,"1/2 Bmsy",cex=0.75) 

ylim <- c(min(exploit.L95),min(max(exploit.U95)))
ylim[2]=0.5

ytext <- c("Exploitation rate") 
plot(xlim,ylim,xlab='',ylab='',type='n',cex.lab=1,lwd=2)

sx = tmb.data$year
low = exploit.L95
high = exploit.U95
polygon(c(sx,rev(sx)),c(low,rev(high)),
  col = 'grey',border = NA)
  box(lty=1)   

lines(tmb.data$year,exploit,type='l',lty=1,lwd=2)  

abline(h=rep$Hmsy,lwd=2,col='red')
r.est = exp(opt$par[1])
sigma2 = exp(opt$par[15])**2 
phi_opt =  rep$Hmsy
abline(h=phi_opt,lwd=2,col='red',lty=2) 
abline(h=0.667*phi_opt,lwd=1,col='red',lty=2)     
mtext(side=2,line=2.5,outer=F,ytext,cex=1,las=0)  
#p1 = rep$Hmsy + 0.07*diff(range(ylim)) 
p1 = rep$Hmsy + 0.02*diff(range(ylim))
text(tmb.data$year[n]-2,p1,"Hmsy",cex=0.75)   
p1 = 0.667*rep$Hmsy + 0.02*diff(range(ylim))
text(tmb.data$year[n]-3,p1,"2/3 Hmsy",cex=0.75)

mtext(side=1,line=2,outer=F,"Year",cex=1.2)

dev.off()
 

##############  biomass and harvest wrt MSY plot ################;

value.names = names(sd.rep$value)

ind = value.names == "log_rB";
biomass = exp(sd.rep$value[ind]);

stdi =  sd.rep$sd[ind]
biomass.L95 = exp(sd.rep$value[ind] - qnorm(0.975)*stdi); 
biomass.U95 = exp(sd.rep$value[ind] + qnorm(0.975)*stdi);

ind = value.names == "log_rH";
exploit = exp(sd.rep$value[ind]);

stdi =  sd.rep$sd[ind]
exploit.L95 = exp(sd.rep$value[ind] - qnorm(0.975)*stdi); 
exploit.U95 = exp(sd.rep$value[ind] + qnorm(0.975)*stdi);


jpeg(file="rpop.jpeg",width=4,height=5,units='in',res=1200)

par(mfcol=c(2,1),oma=c(1,1.5,0,0.5),mar=c(2,2,1,0),cex.axis=1,las=1)
xlim <- c(min(tmb.data$year),max(tmb.data$year))
ylim <- c(min(biomass.L95),max(biomass.U95))
ylim[2]=3

n=length(tmb.data$year)

ytext <- c("Biomass/Bmsy")
plot(xlim,ylim,xlab='',ylab='',type='n',cex.lab=1,lwd=2)

sx = tmb.data$year
low = biomass.L95
high = biomass.U95
polygon(c(sx,rev(sx)),c(low,rev(high)),
  col = 'grey',border = NA)       
  box(lty=1)   
lines(tmb.data$year,biomass,type='l',lty=1,lwd=2)  
  
abline(h=1,lwd=2,col='red')    
abline(h=0.5,lwd=1,col='red')    
mtext(side=2,line=2.5,outer=F,ytext,cex=1,las=0)  

ylim <- c(min(exploit.L95),min(max(exploit.U95)))
ylim[2]=6

ytext <- c("H/Hmsy") 
plot(xlim,ylim,xlab='',ylab='',type='n',cex.lab=1,lwd=2)

sx = tmb.data$year
low = exploit.L95
high = exploit.U95
polygon(c(sx,rev(sx)),c(low,rev(high)),
  col = 'grey',border = NA)
  box(lty=1)   

lines(tmb.data$year,exploit,type='l',lty=1,lwd=2)  

abline(h=1,lwd=2,col='red',lty=1) 
abline(h=0.667,lwd=1,col='red',lty=1)     
mtext(side=2,line=2.5,outer=F,ytext,cex=1,las=0)  

mtext(side=1,line=2,outer=F,"Year",cex=1.2)

dev.off()

############  estimated catch ##################;

jpeg("catch.jpeg",width=5,height=4,units='in',res=1200)

par(mar=c(4,4,1,1),cex.axis=1.2,las=1)

ylim = c(min(exp(rep$log_EC),tmb.data$C),max(exp(rep$log_EC),tmb.data$C))

plot(tmb.data$year,exp(rep$log_EC),type='l',las=1,cex.lab=1.2,
xlab='Year',ylab='',lwd=2,ylim=ylim)
 
#lines(tmb.data$year,tmb.data$C_L,lwd=2,col='grey')
#lines(tmb.data$year,tmb.data$C_U,lwd=2,col='grey')
points(tmb.data$year,tmb.data$C)
yname = c("Catch (Kt)")
mtext(side=2,line=2.5,yname,las=0,cex=1.2)

dev.off()     

############  process error ##################;

jpeg("pe.jpeg",width=3,height=3,units='in',res=1200)

par(mar=c(3,4,0.5,0.5),cex.axis=1,las=1)

plot(tmb.data$year,rep$log_pe,type='l',las=1,xlab='',ylab='',lwd=2) 
mtext(side=2,line=3,c('Log process error'),las=0)     
mtext(side=1,line=2,c('Year'))
abline(h=0,lty=2)   

dev.off()

############ correlation plot ##############;

sd = sqrt(diag(sd.rep$cov.fixed))
corrm = diag(1/sd)%*%sd.rep$cov.fixed%*%diag(1/sd)

pnames = names(sd) 
qname = paste('log_q',levels(index.dat$name))
pnames = replace(pnames,pnames=="log_q",qname)
#stdname = paste('log_std',levels(tmb.data$name))   
#pnames = replace(pnames,pnames=="log_sd_log_index",stdname)
rownames(corrm)=pnames
colnames(corrm)=pnames

require(corrplot)  

jpeg(file="corr.jpeg",width=5,height=5,units='in',res=300)

corrplot.mixed(corrm, tl.srt=45,tl.cex=0.6,cl.cex=0.6,tl.col='black',tl.pos='lt',diag='u',number.cex=0.6,title="Correlation in parameter estimates",
mar=c(0,0,1,0))

dev.off()


########### table ########################;
n=length(tmb.data$year)

out = matrix(NA,n,8);

value.names = names(sd.rep$value)

tnames=c("log_B","log_H","log_rB","log_rH")

for(i in 1:length(tnames)){

  ind = value.names == tnames[i];
  est = exp(sd.rep$value[ind]);
  cv=sd.rep$sd[ind]
  out[,2*i-1] = est 
  out[,2*i] = 100*cv
}        
#out[,1]=out[,1]
out[,3]=out[,3]*100

xalign = c('c',rep('r',8)); # alignment 'l','c', or 'r';
xdigits=c(0,1,0,1,0,2,0,2,0)
colnames(out) = c('Biomass','CV(%)','H (%)','CV(%)','B/Bmsy','CV(%)','H/Hmsy','CV(%)' )
rownames(out) = as.character(tmb.data$year)
ctext = c('Model output')
print(xtable(out,digits=xdigits,caption=ctext,align=xalign),type='html',
file='pop.doc',caption.placement="top")
######################################################;

value.names = names(sd.rep$value)
par.names = names(sd.rep$par.fixed)
parms = sd.rep$par.fixed
parms.sd = sqrt(diag(sd.rep$cov.fixed))

out = matrix(NA,22,2)

j=1
ind = par.names == "log_r";
est = exp(parms[ind]);
cv=parms.sd[ind]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)

j=2   
ind = par.names == "log_K";
est = exp(parms[ind]);
cv=parms.sd[ind]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3) 

j=3   
ind = par.names == "log_Po";
est = exp(parms[ind]);
cv=parms.sd[ind]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)  

j=4:11   
ind = par.names == "log_q";
est = exp(parms[ind])
cv=parms.sd[ind]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)

j=12   
ind = par.names == "log_sd_log_index";
est = exp(parms[ind]);
cv=parms.sd[ind]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)   

j=j+1   
ind = par.names == "log_sd_pe";
est = exp(parms[ind]);
cv=parms.sd[ind]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)   

j=j+1   
ind = par.names == "log_sd_rw";
est = exp(parms[ind]);
cv=parms.sd[ind]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3) 

value.names = names(sd.rep$value)

j=j+1
ind = value.names == "ar_pe";
est = sd.rep$value[ind];
cv=sd.rep$sd[ind]/est
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3) 

j=j+1
ind = value.names == "Hmsy";
est = 100*sd.rep$value[ind];
cv=100*sd.rep$sd[ind]/est
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)   

j=j+1
ind = value.names == "Bmsy";
est = sd.rep$value[ind];
cv=sd.rep$sd[ind]/est
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)    

j=j+1
ind = value.names == "MSY";
est = sd.rep$value[ind];
cv=sd.rep$sd[ind]/est  
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)   

j=j+1
ind = value.names == "log_H"
ind1 = tmb.data$year==2013;
est = 100*exp(sd.rep$value[ind][ind1]);
cv=sd.rep$sd[ind][ind1]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)     

j=j+1
ind = value.names == "log_B"
ind1 = tmb.data$year==2013;
est = exp(sd.rep$value[ind][ind1]);
cv=sd.rep$sd[ind][ind1]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)  

j=j+1
ind = value.names == "log_rH"
ind1 = tmb.data$year==2013;
est = exp(sd.rep$value[ind][ind1]);
cv=sd.rep$sd[ind][ind1]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)     

j=j+1
ind = value.names == "log_rB"
ind1 = tmb.data$year==2013;
est = exp(sd.rep$value[ind][ind1]);
cv=sd.rep$sd[ind][ind1]
out[j,1] = round(est,digits=3)
out[j,2] = round(cv,digits=3)

qname = paste('q',levels(index.dat$name))
#stdname = paste('std',levels(tmb.data$name))
stdname = 'std'

rownames(out) = c("r","K (000s)","Po",qname,stdname,
"sd_pe","sd_rw","ar_pe","Hmsy(%)","Bmsy (000s)","MSY (000s)","H2016(%)","B2013 (000s)","H2013/Hmsy","B2013/Bmsy")
colnames(out) = c("est","CV")

taic = 2*opt$objective + 2*length(opt$par)

tname = paste("Model results, nll = ",round(opt$objective,digits=2),
" AIC = ",round(taic,digits=2),sep="")    
xdigits = c(1,3,3); #a digit for rowname, and all rows;
xalign = c('c',rep('r',2)); # alignment 'l','c', or 'r';

print(xtable(out,digits=xdigits,caption=tname,align=xalign),type='html',
  file='stats.doc',caption.placement="top")

