setwd("C:\\home\\CADIGAN\\GradProgram\\2021\\F6005\\3LN_redfish") 
require(reshape2)
library(stringr)    

                
vec_func = function(mdat){
  vdat = melt(mdat,id=c("Length"))
  vdat$Year = as.numeric(str_replace_all(vdat$variable,"X",""))
  vdat$variable=NULL
  vdat$Catch=vdat$value
  vdat$value=NULL
  return(vdat)
}

load("data\\3LN_redfish.RData") 

year = 1990:2019 
len.rv = sort(unique(FRV.vec$Length))
len.pop=1:54
age.pop=1:10

Y = length(year)
L = length(len.pop)
A = length(age.pop)

#########  create RV aggregate (i.e. plus) lower and upper length groups; 
Lo = 10
Hi = 40

FRV.orig = FRV
# drop Engels indices;
FRV = FRV.orig[,-c(2:6)]
nY = ncol(FRV)

ind = len.rv <= Lo
index.lo = apply(FRV[ind,2:nY],2,sum)  
index.len.lo = min(FRV[ind,1]) 
      
ind = len.rv >= Hi
index.hi = apply(FRV[ind,2:nY],2,sum) 
index.len.hi = max(FRV[ind,1])

ind = (len.rv > Lo)&(len.rv < Hi)

FRVc = rbind(c(Lo,index.lo),FRV[ind,],c(Hi,index.hi))  
#FRVc[FRVc==0]=1; ## replace one 0 with a 1

FRVc.vec = vec_func(FRVc)
FRVc.vec = subset(FRVc.vec,!is.na(FRVc.vec$Catch))
FRVc.vec$Len1 = FRVc.vec$Length 
FRVc.vec$Len2 = FRVc.vec$Length
FRVc.vec$Len1[FRVc.vec$Len1==Lo]=index.len.lo 
FRVc.vec$Len2[FRVc.vec$Len2==Hi]=index.len.hi


#########  create catch plus length group;
len.c = sort(unique(CL.vec$Length)) 

Lo=20
ind = len.c <= Lo
index.lo = apply(CL[ind,2:(Y+1)],2,sum)  
index.len.lo = min(CL[ind,1]) 
 
Hi = 40
      
ind = len.c >= Hi
index.hi = apply(CL[ind,2:(Y+1)],2,sum) 
index.len.hi = max(CL[ind,1])

ind = (len.c > Lo)&(len.c < Hi)
 
CLc = rbind(c(Lo,index.lo),CL[ind,],c(Hi,index.hi)) 
#CLc[CLc==0]=0.1; ## replace one 0 with a 1

CLc.vec = vec_func(CLc) 
CLc.vec = subset(CLc.vec,Year>=min(year))
CLc.vec$Len1 = CLc.vec$Length 
CLc.vec$Len2 = CLc.vec$Length  
CLc.vec$Len1[CLc.vec$Len1==Lo]=index.len.lo 
CLc.vec$Len2[CLc.vec$Len2==Hi]=index.len.hi

## some model dimensions;
A = length(age.pop)  
Y = length(year)
L = length(len.pop)
FRV_L = length(unique(FRV.vec$Length))

## fill a matruity matrix from mat data;

mat.input=mat

pmat1 = rep(0,L) 
pmat2 = rep(0,L)
ind = len.pop %in% mat.input$Length
pmat1[ind]=mat.input$P1             
pmat2[ind]=mat.input$P2    
ind = len.pop > max(mat.input$Length) 
pmat1[ind]=1            
pmat2[ind]=1    

mat = cbind(matrix(pmat1,nrow=L,ncol=length(1990:1995),byrow=F),
matrix(pmat2,nrow=L,ncol=length(1996:2019),byrow=F))

## Miller and Hyun (2017) natural mortality - this is age x year like fishing mortality
## get rough size at age from Cadigan and Campana

len_age = 33*(1-exp(-0.35*(age.pop**0.75)))
amean = mean(LW.parms$a[8:29])      
bmean = mean(LW.parms$b[8:29])
m_age = amean*(len_age**(-0.305*bmean))
m_age = 0.05*m_age/min(m_age)
par(mar=c(3,3,0.5,0.5),mgp=c(2,1,0))        
plot(age.pop,m_age,xlab='age',ylab='M at age',type='l',lwd=2)
   
M = matrix(m_age,nrow=A,ncol=Y,byrow=F)

## breaks between length bins
len_border = len.pop[1:(L-1)]+0.5

## tmb indices for assessment years and ages
iyear = as.numeric(factor(year))-1          
iage = as.numeric(factor(age.pop))-1
ilen = as.numeric(factor(len.pop))-1

#Remove 2 zero's from FRV and 6 zeros or really small values from catch 

FRVc.vec = subset(FRVc.vec,Catch>0) 
CLc.vec = subset(CLc.vec,Catch>0.5)
 
FRV_iyear = iyear[match(FRVc.vec$Year,year)] 
FRV_ilen1 = ilen[match(FRVc.vec$Len1,len.pop)] 
FRV_ilen2 = ilen[match(FRVc.vec$Len2,len.pop)] 

## fraction of year Fall Survey occurs;         
FRV_sf=10.5/12

CL_iyear = iyear[match(CLc.vec$Year,year)] 
CL_ilen1 = ilen[match(CLc.vec$Len1,len.pop)] 
CL_ilen2 = ilen[match(CLc.vec$Len2,len.pop)] 


tmb.data = list(
  A=A,
  Y=Y,
  L=L,  
  age=age.pop,
  len_mid = len.pop, 
  len_border = len_border,
  log_index = log(FRVc.vec$Catch/1000),   
  sf=FRV_sf,
  weight = wt.mat[1:L,],
  mat = mat,
  FRV_iyear = FRV_iyear,
  FRV_ilen1 = FRV_ilen1,   
  FRV_ilen2 = FRV_ilen2, 
  CL_iyear = CL_iyear,
  CL_ilen1 = CL_ilen1,   
  CL_ilen2 = CL_ilen2, 
  log_catch = log(CLc.vec$Catch/1000),
  M=M
)


save(tmb.data,FRVc.vec,CLc.vec,file='tmb.RData')
