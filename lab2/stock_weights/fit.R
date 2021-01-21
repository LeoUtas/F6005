setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/F6005/lab2/stock_weights")

require(reshape2) 
library(stringr) 
library(TMB)   
require(lattice)  
library(latticeExtra)
library('ggplot2') 
                
vec_func = function(mdat){
  vdat = melt(mdat,id=c("year"))
  vdat$age = vdat$variable
  vdat$variable=NULL
  vdat$index=vdat$value
  vdat$value=NULL
  return(vdat)
}
 
RV_len.dat = read.csv('RV_length.csv',header=TRUE,na.strings = ".") 
RV_len.dat = subset(RV_len.dat,!is.na(RV_len.dat$mean))   
RV_len.dat = subset(RV_len.dat,RV_len.dat$age<=21)

# just use these to comapre with model predictions
RV_wt.dat = read.csv('RV_weight.csv',header=TRUE,na.strings = ".") 
RV_wt.dat = subset(RV_wt.dat,!is.na(RV_wt.dat$mean))   
RV_wt.dat = subset(RV_wt.dat,RV_wt.dat$age<=21)
  

RV_len.dat$cohort = RV_len.dat$year-RV_len.dat$age

RV_len.dat$CV = 0.3/sqrt(RV_len.dat$n_otolith)

year = sort(unique(RV_len.dat$year)) 
age = sort(unique(RV_len.dat$age))

pred.dat = expand.grid(year=min(year):max(year),age=min(age):max(age))
pred.dat$cohort = pred.dat$year-pred.dat$age
## first two cohorts in pred.dat not in sub.catch
#pred.dat$cohort[pred.dat$cohort<=1949]=1949

year = sort(unique(pred.dat$year)) 
age = sort(unique(pred.dat$age))
cohort = sort(unique(pred.dat$cohort))

Y=length(year)
A=length(age)

tmb.data = list(
  se = RV_len.dat$CV,
  x = log(RV_len.dat$mean),
  ia = as.numeric(factor(RV_len.dat$age,levels=age))-1,
  iy = as.numeric(factor(RV_len.dat$year,levels=year))-1,
  ic = as.numeric(factor(RV_len.dat$cohort,levels=cohort))-1, 
  iap = as.numeric(factor(pred.dat$age,levels=age))-1,
  iyp = as.numeric(factor(pred.dat$year,levels=year))-1,
  icp = as.numeric(factor(pred.dat$cohort,levels=cohort))-1
)


#dyn.unload("fit"); 

compile("fit.cpp")  

dyn.load(dynlib("fit"))

parameters <- list( 
  age_eff=tapply(log(RV_len.dat$mean),RV_len.dat$age,mean),
  year_eff = rep(0,length(unique(tmb.data$iyp))),       
  cohort_eff = rep(0,length(unique(tmb.data$icp))),
  log_std=rep(log(0.1),3),                 
  logit_ar=rep(3,4),      
  dev=matrix(0,nrow=Y,ncol=A,byrow=T)
)  

parameters.L <- list(
  age_eff=rep(-Inf,length(parameters$age_eff)),
  log_std=rep(-Inf,3),                 
  logit_ar=rep(-10,4)  
)

parameters.U <- list(
  age_eff=rep(Inf,length(parameters$age_eff)),
  log_std=rep(2,3),                 
  logit_ar=rep(10,4)  
)

lower = unlist(parameters.L);
upper = unlist(parameters.U);

obj <- MakeADFun(tmb.data,parameters,
                 random=c("year_eff","cohort_eff","dev"), 
                 DLL = "fit",                 
                 inner.control=list(maxit=500,trace=F))

length(lower)
length(upper)
length(obj$par)

obj$gr(obj$par)

opt <- nlminb(obj$par,obj$fn,obj$gr,upper=upper,lower=lower,
                  control = list(trace=10,eval.max=2000,iter.max=1000))
                  
rep = obj$report()

sd.rep<-sdreport(obj)

### Do Some plotting, my code for this is not always good, lattice is clunky compared to ggplot;

my.padding <- list(layout.heights = list( 
                        top.padding = 0, 
                        main.key.padding = 0, 
                        key.axis.padding = 0, 
                        axis.xlab.padding = 0, 
                        xlab.key.padding = 0, 
                        key.sub.padding = 0, 
                        strip=0.8), 
                layout.widths = list( 
                        left.padding = 1,
                        key.ylab.padding = 0, 
                        ylab.axis.padding = 1, 
                        axis.key.padding = 0, 
                        right.padding = 1) 
                ) 

dat=RV_len.dat
dat$std.resid = rep$std_resid  
dat$clr = sign(dat$std.resid)
dat$pred = exp(rep$Ex)

png(file='resid_bubbles.png',width=5,height=5,units='in',res=300)

ggplot(dat, aes(x=year, y=age, color=factor(clr))) + 
    geom_point(aes(size = 1.5*abs(std.resid)),alpha=0.6, show.legend = FALSE) +
    scale_size_continuous(range = c(0,10)) + 
    theme(plot.margin = margin(2, 2, 2, 2)) +
    ggtitle("RV Length Standardized Residuals") +
    scale_color_manual(values=c('blue','red'))+ 
    scale_y_discrete(name ="Age",limits=3:21)+ 
    scale_x_discrete(name ="year",limits=seq(1960,2020,by=5))
    
dev.off()

tdat=pred.dat
tdat$pred=exp(sd.rep$value)
tdat$pred_L = exp(sd.rep$value - qnorm(0.975)*sd.rep$sd) 
tdat$pred_U = exp(sd.rep$value + qnorm(0.975)*sd.rep$sd)
tdat$group='all'

yr=cbind(tapply(tdat$pred_L,tdat$age,min),tapply(tdat$pred_U,tdat$age,max)) 
yr1=cbind(tapply(dat$mean,dat$age,min),tapply(dat$mean,dat$age,max))
ind = yr1[,1]<yr[,1]
yr[ind,1]=yr1[ind,1]     
ind = yr1[,2]>yr[,2]
yr[ind,2]=yr1[ind,2]

holdRange <- vector('list', length(age)) 
for(i in 1:length(holdRange)){ 
     holdRange[[i]] <-yr[i,] 
     }

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Stock Length (cm)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free"),limits=holdRange))
nr = ceiling(length(2:22)/3) 
  my.key <- list(space = "top",columns=2,
                border = FALSE,
                size = 3,between=0.5,
                lines = list(lty = c(1,1),col=c('white','red'),lwd=c(2,2)),
                points = list(pch = c(1,NA),col=c('black','red')),
                text = list(c("Observed","Predicted")))

P1 = xyplot(mean~year|as.factor(age),
data=dat,col=c('black'), type=c("p"),cex=0.5,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,
as.table=TRUE, layout=c(3,nr,1),
key=my.key,par.settings = my.padding,
par.strip.text = list(cex=0.6),
)

my.panel.bands <- function(x, y, upper, lower, fill, col,
 subscripts, ..., font, fontface)
 {
 upper <- upper[subscripts]
 lower <- lower[subscripts]
 panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
 col = 'pink', border = FALSE,
 ...)
 }
 
 P2 = xyplot(pred~year|as.factor(age),groups=group,  
 upper = tdat$pred_U, lower = tdat$pred_L,
data=tdat, xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,
as.table=TRUE, layout=c(3,nr,1),alternating=FALSE,
key=my.key,par.settings = my.padding,
par.strip.text = list(cex=0.6), 
 panel = function(x, y, ...){
 panel.superpose(x, y, panel.groups = my.panel.bands, type='l',...)
 panel.xyplot(x, y, type='l', cex=0.6, lty=1,lwd=2,col='red',...)
 }
)

add_axes <- function() {
  library(grid)
  library(lattice)
  l <- trellis.currentLayout()
  pan <- which(l[nrow(l), ]==0)
  if(length(pan) > 0) {
    g <- grid.ls(print=FALSE)
    # use an existing panel as a template for ticks 
    ticks <- grid.get(g$name[grep("ticks.bottom.panel", g$name)][[1]])
    # use an existing panel as a template for labels
    labels <- grid.get(g$name[grep("ticklabels.bottom.panel", g$name)][[1]])
    ax <- grobTree(ticks, labels)
    invisible(sapply(pan, function(x) {
      trellis.focus("panel", x, nrow(l)-1, clip.off=TRUE)
      grid.draw(ax)
      trellis.unfocus()  
    }))
  }
}

#P2 + as.layer(P1)

gname = "obs_pred.jpeg"
jpeg(file=gname,width=5,height=6,units='in',res=300)       
  #print(P3)
  P2 + as.layer(P1)
  add_axes()
dev.off()

gname = paste('model_effects1.jpeg',sep='')
jpeg(file=gname,width=4,height=4,units='in',res=300)

par(mfcol=c(3,1),oma=c(1,3,0.5,1),mar=c(2,2,0.5,2),cex.axis=1,cex.lab=1,las=1,mgp=c(1,1,0))

ylim = range(rep$dev,rep$year_eff,rep$cohort_eff)

x <- age
y <- rep$age_eff
plot(x,y,xlab="",ylab="",type='l') 
abline(h=0,lty=1)
mtext(side=4,line=0.5,outer=F,'Age',las=0)

x <- year
y <- rep$year_eff
plot(x,y,xlab="",ylab="",type='l',ylim=ylim) 
abline(h=0,lty=1)
mtext(side=4,line=0.5,outer=F,'Year',las=0)

x <- sort(unique(dat$cohort))
y <- rep$cohort_eff
plot(x,y,xlab="",ylab="",type='l',ylim=ylim) 
abline(h=0,lty=1)
mtext(side=4,line=0.5,outer=F,'cohort',las=0)

ytext <- c("Model Effects")
mtext(side=2,ytext,line=1,outer=T,las=0,cex=1.2)

dev.off()

colnames(rep$dev)=age
rownames(rep$dev)=year  
                  
tdat = vec_func(t(rep$dev)) 
names(tdat) = c('age','year','dev')
tdat$clr = sign(tdat$dev)

png(file='model_effects2.png',width=5,height=5,units='in',res=300)

ggplot(tdat, aes(x=year, y=age, color=factor(clr))) + 
    geom_point(aes(size = abs(dev)),alpha=0.6, show.legend = FALSE) +
    scale_size_continuous(range = c(0,5)) + 
    theme(plot.margin = margin(2, 2, 2, 2)) +
    ggtitle("Model Effects: Year x Age Deviations") +
    scale_color_manual(values=c('blue','red'))+ 
    scale_y_discrete(name ="Age",limits=1:21)+ 
    scale_x_discrete(name ="year",limits=seq(1975,2018,by=5))
    
dev.off() 


gname = paste('resid_4P.jpeg',sep='')
jpeg(file=gname,width=4,height=5,units='in',res=300)  
 
par(mfcol=c(4,1),oma=c(1,3,1,1),mar=c(3,2,0,2),cex.axis=1.2,cex.lab=1.2,las=1)

x <- dat$year
y <- dat$std.resid

uage = sort(unique(dat$age))

plot(x,y,xlab="",ylab="",type='n')
abline(v=seq(1960,2020,by=5),lty=1,col=grey(0.9))
text(x,y,dat$age,cex=0.5)
mres <- tapply(y,x,"mean")
lines(sort(unique(x)),mres,lty=1,lwd=3,col='red') 
abline(h=0,lty=1)
mtext(side=4,line=0.5,outer=F,'Year',las=0)

x=dat$cohort
plot(x,y,xlab="",ylab="",pch=3)
abline(h=0,lty=1)     
mres <- tapply(y,x,"mean")
lines(sort(unique(x)),mres,lty=1,lwd=3,col='red')
mtext(side=4,line=0.5,outer=F,'Cohort',las=0)

plot(dat$age,y,xlab="",ylab="",pch=3)
abline(h=0,lty=2)
mtext(side=4,line=0.5,outer=F,'Age',las=0)
lines(uage,tapply(y,dat$age,mean),lty=1)

plot(dat$pred,y,xlab="",ylab="",pch=3)
abline(h=0,lty=1)
mtext(side=4,line=0.5,outer=F,'Expected',las=0)

ytext <- c("Standardized residuals")
mtext(side=2,ytext,line=1,outer=T,las=0,cex=1.2)

dev.off() 


#predicted CW for assessment;

# these values come from the weight~length segmented regression
a = -11.5789
b1 = 2.9583
b2 = 3.3048
bp = 3.7874

tdat=pred.dat
tdat$pllen=sd.rep$value
tdat$pred = exp(a + b1*tdat$pllen)
ind = tdat$pllen>=bp  
tdat$pred[ind] = exp(a + (b1-b2)*bp + b2*tdat$pllen[ind])
tdat$plen = exp(sd.rep$value)

SW.pred=dcast(tdat[,c(1,2,5)],year~age)
year = SW.pred[,1]

ind = 15:ncol(SW.pred)  
mwt = exp(-0.4*(ind-14))
mwt=mwt/sum(mwt)

SW.pred = cbind(SW.pred[,2:14],apply(SW.pred[,ind],1,weighted.mean,w=mwt))
colnames(SW.pred) = 1:14

year.old = 1959:(min(year)-1)
year.all = c(year.old,year)

mean.old = apply(SW.pred[1:5,],2,mean)
SW.pred.old = matrix(mean.old,nrow=length(year.old),ncol=ncol(SW.pred),byrow=T)

rownames(SW.pred)=year
rownames(SW.pred.old)=year.old 
colnames(SW.pred.old)=1:14

SW.pred.all = rbind(SW.pred.old,SW.pred) 
SW.pred.all = matrix(unlist(SW.pred.all),nrow=length(year.all),ncol=ncol(SW.pred))

colnames(SW.pred)=1:14

save(SW.pred.all,file='SW.pred.RData')

year=year.all 

mdat = data.frame(cbind(year,SW.pred.all))
names(mdat) = c('year',1:14)
SW.vec = vec_func(mdat)

WM=dcast(RV_wt.dat[,c(1,2,4)],year~age)
year1 = WM[,1]
mdat = data.frame(cbind(year1,WM[,2:14]))
names(mdat) = c('year',1:13)

Wo.vec = vec_func(mdat)
levels(Wo.vec$age)=1:14
Wo.vec = rbind(Wo.vec,c(2018,14,NA)) 


xttl <- list(c("Year"), cex=1)
yttl <- list(c("Predicted Stock Weights (kg)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
nr = ceiling(length(2:14)/3)


yr=cbind(tapply(SW.vec$index,SW.vec$age,min),tapply(SW.vec$index,SW.vec$age,max)) 
yr1=cbind(tapply(Wo.vec$index,Wo.vec$age,min,na.rm=T),tapply(Wo.vec$index,Wo.vec$age,max,na.rm=T))

ind = yr1[,1]<yr[,1]
yr[ind,1]=yr1[ind,1]     
ind = yr1[,2]>yr[,2]
yr[ind,2]=yr1[ind,2]

age=1:14

holdRange <- vector('list', length(age)) 
for(i in 1:length(holdRange)){ 
     holdRange[[i]] <-yr[i,] 
     }

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Stock Weight (kg)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free"),limits=holdRange))
nr = ceiling(length(1:14)/3) 
  my.key <- list(space = "top",columns=2,
                border = FALSE,
                size = 3,between=0.5,
                lines = list(lty = c(1,1),col=c('white','red'),lwd=c(2,2)),
                points = list(pch = c(1,NA),col=c('black','red')),
                text = list(c("Observed","Predicted")))


P1 = xyplot(index~year|as.factor(age),
data=Wo.vec,col=c('black'), type=c("p"),cex=0.5,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,
as.table=TRUE, layout=c(3,nr,1),
key=my.key,par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
 
P2 = xyplot(index~year|as.factor(age),
data=SW.vec, xlab=xttl, scales=ax,col='red',type='l',lwd=2,
ylab=yttl, strip.text=stripttl, las=1,
as.table=TRUE, layout=c(3,nr,1),alternating=FALSE,
key=my.key,par.settings = my.padding,
par.strip.text = list(cex=0.6), 
)

#P2 + as.layer(P1)

gname = "SW_obs_pred.jpeg"
jpeg(file=gname,width=5,height=6,units='in',res=300)       
  #print(P3)
  P2 + as.layer(P1)
  add_axes()
dev.off()

## Compare CW and SW ###;


setwd("C:\\home\\CADIGAN\\stocks\\models\\3Ps_cod\\catch_weights\\3+")

load('catch_weights.RData')

         
setwd("C:\\home\\CADIGAN\\stocks\\models\\3Ps_cod\\stock_weights")  

yearm1 = year[1:(length(year)-1)]
catch.wt = cbind(matrix(0,nrow=length(yearm1),ncol=2),catch.wt[,1:ncol(catch.wt)])

mdat = cbind(yearm1,catch.wt)
names(mdat)=c('year',names(catch.wt))
W.vec = vec_func(mdat)                
names(W.vec)=c("year","age","CW")                  
W.vec$SW = SW.vec$index[SW.vec$year < max(year)] 

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Weight (kg)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
nr = ceiling(length(1:14)/3) 
  my.key <- list(space = "top",columns=2,
                border = FALSE,
                size = 3,between=0.5,
                lines = list(lty = c(1,1),col=c('blue','red'),lwd=c(2,2)),
                text = list(c("Stock","Catch")))


P1 = xyplot(SW+CW~year|as.factor(age),
data=W.vec,col=c('blue','red'), type='l',lwd=2,lty=1,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,
as.table=TRUE, layout=c(3,nr,1),
key=my.key,par.settings = my.padding,
par.strip.text = list(cex=0.6),
)

gname = "compare_wts.jpeg"
jpeg(file=gname,width=5,height=6,units='in',res=300)
  P1
  add_axes()
dev.off()                 

