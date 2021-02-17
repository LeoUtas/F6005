require(reshape2)

vec_func = function(mdat){
  vdat = melt(mdat,id=c("Year"))
  vdat$Age = as.numeric(substr(vdat$variable,4,5))
  vdat$variable=NULL
  vdat$index=vdat$value
  vdat$value=NULL
  return(vdat)
} 


setwd("C:\\home\\CADIGAN\\GradProgram\\2018\\F6004\\week8\\F3LNO_SURBA")

age=0:12

cnames = c('Year',paste('Age',age,sep=""))
mdat = read.table(file='FRV.dat',header=FALSE,col.names=cnames)
year=seq(min(mdat[,1]),max(mdat[,1]),by=1)

##read and make maturities

cnames = c('Year',paste('Age',1:14,sep=""))
allmat.dat = read.table(file='mat.txt',header=FALSE,col.names=cnames)
mat.dat = subset(allmat.dat,(Year>=1990)&(Year<=2016))

mat = matrix(0,nrow=length(1990:2016),ncol=length(0:12))
mat[,2:13] = as.matrix(mat.dat[,2:13])

## read and make weights

cnames = c('Year',paste('Age',3:12,sep=""))
allwt.dat = read.table(file='wt.txt',header=FALSE,col.names=cnames) 
wt.dat = subset(allwt.dat,(Year>=1990)&(Year<=2016))

vwt.dat = melt(wt.dat,id=c("Year")) 
vwt.dat$Age = as.numeric(substr(vwt.dat$variable,4,5)) + 10/12
vwt.dat$variable=NULL
vwt.dat$wt=vwt.dat$value
vwt.dat$value=NULL
vwt.dat$logwt = log(vwt.dat$wt)

vwt.dat$map = as.numeric(as.factor(vwt.dat$Year))
uyear=unique(vwt.dat$Year)
nY = length(unique(vwt.dat$Year))

gfit = function(map,wt,Age,logWinf,logk,logpo){
  Winf=exp(logWinf);k=exp(logk);po=exp(logpo)
  pred = Winf[map]*((1 - (1-po)*exp(-k[map]*Age))**3.08)
return(log(pred))}

start.parm = list(logWinf=log(rep(10,nY)),logk=log(rep(0.2,nY)),logpo=log(0.02)) 
lower.parm = list(logWinf=log(rep(5,nY)),logk=log(rep(0.01,nY)),logpo=log(0.0002)) 
upper.parm = list(logWinf=log(rep(35,nY)),logk=log(rep(0.4,nY)),logpo=log(0.2))
lower=unlist(lower.parm)
upper=unlist(upper.parm)

fit <- nls(logwt ~ gfit(map,wt,Age,logWinf,logk,logpo),
  data=vwt.dat,start = start.parm,
  algorithm="port",control=list(maxiter=5000),
  lower=lower,upper=upper)
  
vwt.dat$pred_wt = exp(predict(fit,vwt.dat))
  
pred.dat = data.frame(Year=rep(uyear,13),Age = rep(0:12+10/12,each=nY))  
pred.dat$map = as.numeric(as.factor(pred.dat$Year))   
pred.dat$logwt = predict(fit,newdata=pred.dat)
pred.dat$wt = exp(pred.dat$logwt)

wtm = matrix(pred.dat$wt,nrow=nY,ncol=13,byrow=FALSE)
wtm[,4:13] = as.matrix(wt.dat[,2:11])

pred.dat = data.frame(Year=rep(uyear,13),Age = rep(0:12+0.5,each=nY))
pred.dat$map = as.numeric(as.factor(pred.dat$Year))   
pred.dat$logwt = predict(fit,newdata=pred.dat)
pred.dat$wt = exp(pred.dat$logwt)

ave_wt = aggregate(pred.dat$wt,list(age = pred.dat$Age),mean)
m_age = 0.2*((ave_wt$x/35)**(-0.305))
m_age = 0.2*m_age/m_age[length(m_age)]


FRV.vec = vec_func(mdat)

temp = cbind(year,wtm)
colnames(temp) = colnames(mdat)
wt.vec = vec_func(as.data.frame(temp))  

temp = cbind(year,mat)
colnames(temp) = colnames(mdat)
mat.vec = vec_func(as.data.frame(temp)) 

data=data.frame(
  Year=mat.vec$Year, 
  Age=mat.vec$Age,   
  weight = wt.vec$index,
  mat = mat.vec$index
)
data=merge(data,FRV.vec,by=c("Year","Age"),all.x=TRUE)
data$iyear = as.numeric(factor(data$Year))-1          
data$iage = as.numeric(factor(data$Age))-1

data=subset(data,!is.na(data$index))
data=subset(data,data$index>0)

data$log_index = log(data$index)
data$sf=10/12


setwd("C:\\home\\CADIGAN\\GradProgram\\2018\\F6005\\Lecture5\\3NO_cod_SURBA++")

A = length(age)
Y = length(year)        

q_age = c(0.01,0.3,1,1,1,1,1,1,1,1,1,1,1)
log_q_age = log(q_age)

tmb.data = list(
  log_index = data$log_index,
  logq=log_q_age,
  sf=data$sf,
  weight = t(wtm),
  mat = t(mat),
  iyear = data$iyear,
  iage = data$iage,
  m=m_age
)


save(data,tmb.data,file='tmb.RData')
