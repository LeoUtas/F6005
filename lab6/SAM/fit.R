 # SAM AT https://github.com/fishfollower/SAM 
library(stockassessment)  
library(ggplot2)  
library(reshape2)   

## read NAFO assessment results

assess.ssb = read.table("C:\\home\\CADIGAN\\stocks\\models\\SAM\\3LNO_Amplaice\\assmt_SSB.txt",
header=FALSE,col.names=c('year','ssb'),sep = "\t")
assess.ssb$ssb = assess.ssb$ssb/1000

assess.F = read.table("C:\\home\\CADIGAN\\stocks\\models\\SAM\\3LNO_Amplaice\\assmt_F.txt",
header=FALSE,col.names=c('year','est'),sep = "\t")

nafo.ssb = read.table("C:\\home\\CADIGAN\\stocks\\models\\SAM\\3LNO_Amplaice\\NAFOJ_ssb.txt",
header=FALSE,col.names=c('year','est','low','high'),sep = "\t")

nafo.F = read.table("C:\\home\\CADIGAN\\stocks\\models\\SAM\\3LNO_Amplaice\\NAFOJ_F.txt",
header=FALSE,col.names=c('year','est','low','high'),sep = "\t")

 
setwd("C:\\home\\CADIGAN\\stocks\\models\\SAM\\3LNO_Amplaice\\data")
 
cn <- read.ices("cn.dat")/1000
cw <- read.ices("cw.dat")
dw <- read.ices("dw.dat")
lf <- read.ices("lf.dat")
lw <- read.ices("lw.dat")
mo <- read.ices("mo.dat")
nm <- read.ices("nm.dat")
pf <- read.ices("pf.dat")
pm <- read.ices("pm.dat")
sw <- read.ices("sw.dat")
surveys <- read.ices("survey.dat")

dat <- setup.sam.data(
  surveys=surveys,
  residual.fleet=cn, 
  prop.mature=mo, 
  stock.mean.weight=sw, 
  catch.mean.weight=cw, 
  dis.mean.weight=dw, 
  land.mean.weight=lw,
  prop.f=pf, 
  prop.m=pm, 
  natural.mortality=nm, 
  land.frac=lf)
                      
conf <- defcon(dat)

conf$fbarRange <- c(9,14)

par <- defpar(dat,conf)

conf$keyLogFpar <- matrix(
         c(rep(-1,15),
            c(0:6,rep(6,8)),
            c(7:13,rep(13,8)),
            c(14:20,rep(20,8))
           ), nrow=4, byrow=TRUE)

conf$keyVarObs <- matrix(
         c( -1,-1,-1,-1,0,  1,  1,  1,  1,1,1,1,1,1,1,
            2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,
            5,  6,  6,  6, 6,7,7,7,7,7,7,7,8,8,8,
            9,10,10,10,10,11,11,11,11,11,11,11,11,11,11
           ), nrow=4, byrow=TRUE)

par <- defpar(dat,conf)
par$logFpar =  rep(0,length(par$logFpar))

fit <- sam.fit(dat,conf,par)

dataplot(fit)

tble = summary(fit,digits=3)

write.table(ssbtable(fit),file='..//ssb.dat') 
write.table(fbartable(fit),file='..//fbar.dat') 

par(mar=c(3,3,1,1))
ssbplot(fit)
lines(assess.ssb$year,assess.ssb$ssb,lwd=2,lty=1,col='red') 
lines(nafo.ssb$year,nafo.ssb$est,lwd=2,lty=1,col='blue')
legend("topright",c('SAM','ADAPT','NAFOJ'),col=c('black','red','blue'),lty=1,bty='n',lwd=2)

cv=fit$sdrep$sd[names(fit$sdrep$value)=="logssb"]
plot(assess.F$year,cv)

par(mar=c(4,4,0.5,0.5))
fbarplot(fit,partial=TRUE,page=c(9,10,13,14))   
lines(assess.F$year,assess.F$est,lwd=2,lty=1,col='red') 
lines(nafo.F$year,nafo.F$est,lwd=2,lty=1,col='blue')
legend("topright",c('SAM','ADAPT','NAFOJ'),col=c('black','red','blue'),lty=1,bty='n',lwd=2)


recplot(fit)
 
catchplot(fit)

F = faytable(fit)
year=as.numeric(rownames(F)) 
uage=as.numeric(colnames(F))

temp = expand.grid(Year=year, Age=uage)
temp$ssb = as.vector(F)
            
  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
  persp(year,uage,F, 
      main="F-at-age Surface",
      xlab='Year',ylab='Age',zlab='F',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
      

par(mar=c(2,2,0,0),mgp=c(2,1,0))
fitplot(fit)


res <- residuals(fit)     
par(mar=c(4,4,0.5,0.5),mgp=c(2,1,0))
plot(res)
#empirobscorrplot(res)

pres=procres(fit) 
pres$age[pres$fleet==2] = pres$age[pres$fleet==2]+5   
par(mar=c(4,4,0.5,0.5),mgp=c(2,1,0))
plot(pres)
 
retro <- retro(fit,year=4)   
rho <- mohn(retro, lag=0)

 par(mar=c(3,3,0.5,0.5),mgp=c(2,1,0))
  ssbplot(retro,drop=0,xlim=c(1999,2017),las=0) 
  legend("bottomright", legend=bquote(rho == .(rho[2])), bty="n")
  
  par(mar=c(3,3,0.5,0.5),mgp=c(2,1,0))
  fbarplot(retro,drop=0,xlim=c(1999,2017),las=0)
  legend("bottomright", legend=bquote(rho == .(rho[3])), bty="n")
  
  sim=simstudy(fit, nsim=10)


qest = data.frame(survey=rep(c("FRV","SRV","JRV"),each=15),age=rep(1:15,3))

temp = partable(fit)[1:21,]
t2 = conf$keyLogFpar
temp1 = rbind(temp[t2[2,]+1,],temp[t2[3,]+1,],temp[t2[4,]+1,])
qest = cbind(qest,temp1)[,c(1,2,5,6,7)]
colnames(qest) = c('survey','age','est','Low','High')

#jpeg(file='q.jpeg',width=4,height=4,units='in',res=300)

p1=ggplot(qest, aes(x=age, y=est,fill=survey,color=survey)) +
    geom_line(aes(y = est, color=survey),size=1) + 
    labs(x = "Age",y = "Index Q")+
    geom_smooth(aes(ymin = Low, ymax = High,fill=survey,color=survey),stat = "identity",alpha=0.2)
print(p1)
  
#dev.off()
 
 
 
 