library('ggplot2') 
require(reshape2) 
library(stringr)   
library('ggridges')  
library('gridExtra')    
library(patchwork)
library(corrplot)    
library(xtable)

library(viridis)
library(hrbrthemes) 

len.pop=1:65
age.pop=1:10     
age.Fdev=3:10
year = 1990:2019
value.names = names(sd.rep$value)

sname=tapply(allSL_vec$sname,allSL_vec$is,unique) 
ns = length(sname)
#table(allSL_vec$sname,allSL_vec$is)

######### plot log_SDs ###############;

par.names = names(opt$par)
ind = substr(par.names,1,7) == 'log_std'
x = exp(opt$par[ind])

cx = rep('darkorange3',length=length(x))
ind1 = names(x) == 'log_std_index'   
cx[ind1]='grey' 
ind1 = names(x) == 'log_std_logq'   
cx[ind1]='green'

names(x) = sub("log_std_","",names(x))  
xn = names(x) 

xn[xn=='index']=sname 
xn[xn=='logq']=paste('logq',sname,sep="_")

jpeg(file='figs//sd_barplot.jpeg',width=3,height=4,units='in',res=300)

par(mar=c(3,7,0,0.5),mgp=c(1.5,0.5,0)) # increase y-axis margin.

barplot(x,names.arg=xn,horiz=TRUE,xlab="SD's",las=1,col=cx)
abline(v=c(0.25,0.5,0.75,1),col='grey',lty=2)  

dev.off()
 

## corr plot #######

dsd = sqrt(diag(sd.rep$cov.fixed))
corr.matrix = diag(1/dsd)%*%sd.rep$cov.fixed%*%diag(1/dsd)
rownames(corr.matrix) = names(dsd)                
colnames(corr.matrix) = names(dsd)

 
jpeg(file='figs//corr_parm.jpeg',width=7,height=7,units='in',res=300)
par(mar=c(0,5,0,1))

corrplot.mixed(corr.matrix, tl.pos='lt',tl.srt=45,tl.cex=0.6,cl.cex=0.75,tl.col='black',number.cex=0.5)

dev.off() 


## plot of proportion at length for each age;
#jpeg(file='figs//PLA.jpeg',width=5,height=5,units='in',res=300)
#  image(age.pop,len.pop,t(rep$pla),xlab="Age",ylab="Length (cm)",las=1)
#  lines(age.pop,rep$mean_len)
#dev.off()

## plot of proportion at length for each age;
jpeg(file='figs//PLA.jpeg',width=5,height=5,units='in',res=300) 
  par(mfrow=c(3,2),oma=c(1,1,0.5,0),mar=c(2,2,0,0),mgp=c(1,0.5,0),cex=1)
  image(age.pop,len.pop,t(rep$pla),xlab="",ylab="",las=1)
  lines(age.pop,rep$mean_len)
  legend("topleft","Jan 1",bty='n')
      
  image(age.pop,len.pop,t(rep$plac),xlab="",ylab="",las=1)
  lines(age.pop,rep$mean_len) 
  legend("topleft","Catch",bty='n')
  for(i in 1:ns){
    image(age.pop,len.pop,t(rep$plas[i,,]),xlab="",ylab="",las=1)
    lines(age.pop,rep$mean_len)  
  legend("topleft",sname[i],bty='n')
  }  
  mtext(side=1,line=0,outer=T,"Age") 
  mtext(side=2,line=0,outer=T,"Length (cm)")
  
dev.off()

# plot Catch at lengths, observed and model predicted;
CLc_vec$ECatch = 1000*exp(rep$Elog_catch)

jpeg(file='figs//Catch_ridges.jpeg',width=7,height=8,units='in',res=300)
p1=ggplot(CLc_vec)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = Catch, fill=as.factor(Year)),
                      stat = "identity", scale = 3,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1990:2019, expand=c(0,0))+
  scale_x_continuous(limits=c(20,40), breaks=c(20,25,30,35,40,45), expand=c(0,0))+
  xlab('Length (mm)')+
  ylab('Year')+    
  ggtitle('Observed')+ 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=14))

p2=ggplot(CLc_vec)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = ECatch, fill=as.factor(Year)),                                         
                      stat = "identity", scale = 3,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1990:2019, expand=c(0,0))+ 
  scale_x_continuous(limits=c(20,40), breaks=c(20,25,30,35,40,45), expand=c(0,0))+
  xlab('Length (mm)')+ 
  ylab('')+     
  ggtitle('Predicted')+ 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10), 
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5,size=14))

print(p1|p2) 
dev.off()  

####### CL at length ##########

gname <- 'figs//log_catch_at_len.jpeg'
jpeg(file=gname,width=6,height=8,units='in',res=300)    

p1=ggplot(CLc_vec) + geom_line(aes(x = Year, y = log(ECatch/1000)),size=1) +
    geom_point(aes(x=Year, y = log(Catch/1000))) +
    facet_wrap(~Length, ncol = 3,scales = "free_y") + geom_hline(yintercept = 0, linetype="dashed") + 
    xlab('Year')+
    ylab('Log Catch')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
print(p1 + theme(axis.text.x = element_text(angle = 45)))

dev.off()      

gname <- 'figs//catch_at_len.jpeg'
jpeg(file=gname,width=6,height=8,units='in',res=300)    

p1=ggplot(CLc_vec) + geom_line(aes(x = Year, y = ECatch/1000),size=1) +
    geom_point(aes(x=Year, y = Catch/1000)) +
    facet_wrap(~Length, ncol = 3,scales = "free_y") + geom_hline(yintercept = 0, linetype="dashed") + 
    xlab('Year')+
    ylab('Catch (000s)')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
print(p1 + theme(axis.text.x = element_text(angle = 45)))

dev.off()


CLc.vec = CLc_vec
CLc.vec$resid = rep$std_resid_catch
CLc.vec$size = abs(CLc.vec$resid)
CLc.vec$clr = '+'
CLc.vec$clr[rep$std_resid_catch<0] = '-'

jpeg(file='figs//Catch_resid_bubble.jpeg',width=6,height=4,units='in',res=300)
p1=ggplot(CLc.vec,aes(x=Year, y=Length, size=size, color=clr)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 6), name="Std. Resid") + 
    scale_color_manual(values=c("red", "blue"))

print(p1) 
dev.off()

# 3 panel plot

jpeg(file='figs//Catch_resid_3P.jpeg',width=5,height=7,units='in',res=300)
p1 = ggplot(CLc.vec,aes(x=Year, y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') +
     stat_summary(fun=mean, aes(group=1), geom="line",size=1) +
     xlab('Year')+ ylab('Std. Residual')    
p2 = ggplot(CLc.vec,aes(x=Length, y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') +
     stat_summary(fun=mean, aes(group=1), geom="line",size=1) +
     xlab('Length')+ ylab('Std. Residual')     
p3 = ggplot(CLc.vec,aes(x=log(ECatch), y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') +
     xlab('Log Expected Catch')+ ylab('Std. Residual')
print(p1 / p2 / p3)
dev.off()

####### F at age ##########

vec_func = function(mdat){
  vdat = melt(mdat,id=c("Age"))
  vdat$Year = as.numeric(str_replace_all(vdat$variable,"X",""))
  vdat$variable=NULL
  return(vdat)
}

Fmat = data.frame(cbind(age.pop,rep$F))
colnames(Fmat) = c('Age',year)

pdat=  vec_func(Fmat)
pdat = subset(pdat,Age<=20)

gname <- 'figs//F_at_age.jpeg'
jpeg(file=gname,width=6,height=5,units='in',res=300)    

p1=ggplot(subset(pdat,Age>2)) + geom_line(aes(x = Year, y = value),size=1) +
    facet_wrap(~Age, ncol = 2) +          
    xlab('Year')+
    ylab('F')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
#print(p1 + theme(axis.text.x = element_text(angle = 45)))
print(p1)

dev.off()     


gname <- 'figs//F_at_age1.jpeg'
jpeg(file=gname,width=6,height=5,units='in',res=300)    

p1=ggplot(subset(pdat,Age>2)) + geom_line(aes(x = Year, y = value),size=1) +
    facet_wrap(~Age, ncol = 2,scales = "free_y") +          
    xlab('Year')+
    ylab('F')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
#print(p1 + theme(axis.text.x = element_text(angle = 45)))  
print(p1)

dev.off() 


Fdev = data.frame(cbind(age.Fdev,rep$log_F_dev))
colnames(Fdev) = c('Age',year) 

pdat= vec_func(Fdev)  
pdat$size = abs(pdat$value)
pdat$clr = '+'
pdat$clr[pdat$value<0] = '-'

jpeg(file='figs//Fdev_bubble.jpeg',width=6,height=4,units='in',res=300)
p1=ggplot(pdat,aes(x=Year, y=Age, size=size, color=clr)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 6), name="Std. Resid") + 
    scale_color_manual(values=c("red", "blue"))

print(p1) 
dev.off()

## plot of pop N;
jpeg(file='figs//Pop.jpeg',width=7,height=7,units='in',res=300)

par(mfrow=c(2,2),oma=c(3,1,1,1),mar=c(0.5,4,1,1),mgp=c(2,1,0),cex=1)

  plot(year,exp(rep$log_Rec),type='l',xlab="",ylab="",lwd=2,las=1,xaxt='n') 
  mtext(side=3,line=0,outer=F,c("Recruitment")) 
  plot(year,rep$ssb,type='l',xlab="",ylab="",lwd=2,las=1,xaxt='n')   
  mtext(side=3,line=0,outer=F,c("SSB(Kt)")) 

  image(year,age.pop,t(rep$N_matrix),xlab="",ylab="Age",las=1)
  mtext(side=3,line=0,outer=F,"Numbers at age") 
  image(year,len.pop,t(rep$NL),xlab="",ylab="Length (cm)",las=1)
  mtext(side=3,line=0,outer=F,"Numbers at length") 
  
  mtext(side=1,line=1,outer=T,"Year") 
  
dev.off()


### Biomass plot with CIs ###

ind = value.names=='log_ssb'
pdat = data.frame(year=year,est=exp(sd.rep$value[ind]))
pdat$L = exp(sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind]) 
pdat$U = exp(sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind])
pdat$type='SSB'

ind = value.names=='log_biomass'
pdat1 = data.frame(year=year,est=exp(sd.rep$value[ind]))
pdat1$L = exp(sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind]) 
pdat1$U = exp(sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind])
pdat1$type='TSB'

pdat = rbind(pdat,pdat1)

ylim = range(pdat$L,pdat$U)
ylim[2] = min(ylim[2],max(pdat$est)*1.5)

jpeg(file='figs//biomass.jpeg',width=5,height=4,units='in',res=300)

p1=ggplot(pdat, aes(x=year, y=est,fill = type,color=type)) +
    geom_line(aes(y = est, group=type, color=type),size=1) + 
    labs(x = "Year",y = "Biomass (Kt)")+
    geom_smooth(aes(ymin = L, ymax = U,fill=type,color=type),stat = "identity",alpha=0.2) +
    coord_cartesian(ylim = ylim)
print(p1)
    
dev.off()       
       
###  Harvest rate plots with CIs ###

ind = value.names=='log_harvest_rate'
pdat = data.frame(year=year,est=100*exp(sd.rep$value[ind]))
pdat$L = 100*exp(sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind]) 
pdat$U = 100*exp(sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind])

ylim = range(pdat$L,pdat$U)
ylim[2] = min(ylim[2],max(pdat$est)*1.5)

jpeg(file='figs//harvest_rate.jpeg',width=5,height=4,units='in',res=300)

p1=ggplot(pdat, aes(x=year, y=est)) +
    geom_line(aes(y = est),size=1) + 
    labs(x = "Year",y = "Hrvest Rate (%)")+
    geom_smooth(aes(ymin = L, ymax = U,color='grey'),stat = "identity",alpha=0.2) +
    coord_cartesian(ylim = ylim)  + guides(color = FALSE)
print(p1)
    
dev.off()  


## plot of Catch N;
jpeg(file='figs//Catch.jpeg',width=7,height=7,units='in',res=300)

par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(3,4,1,1),mgp=c(2,1,0),cex=1) 

  image(year,age.pop,t(rep$F),xlab="",ylab="Age",las=1)
  mtext(side=3,line=0.5,outer=F,"F at age")          
  image(year,age.pop,t(rep$Z),xlab="",ylab="Age",las=1)
  mtext(side=3,line=0.5,outer=F,"Z at age") 
  image(year,age.pop,t(rep$CNA),xlab="",ylab="Age",las=1)
  mtext(side=3,line=0.5,outer=F,"Catch at age")  
  image(year,len.pop,t(rep$CL),xlab="",ylab="Length (cm)",las=1)
  mtext(side=3,line=0.5,outer=F,"Catch at Length")   
  mtext(side=1,line=0,outer=T,"Year") 
  
dev.off()

pdat = data.frame(sname=rep(sname[1],tmb_data$L),qest = exp(rep$logq[1,]),len=len.pop)
for(i in 2:4){pdat = rbind(pdat,data.frame(sname=rep(sname[i],tmb_data$L),qest = exp(rep$logq[i,]),len=len.pop))}

jpeg(file='figs//Q.jpeg',width=5,height=5,units='in',res=300)    

p1=ggplot(pdat) + geom_line(aes(x = len, y = qest),size=1) +
    facet_wrap(~sname, ncol = 2,scales = "free_y") +          
    xlab('Length')+
    ylab('Q')+       
   # coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
#print(p1 + theme(axis.text.x = element_text(angle = 45)))  
print(p1)
dev.off()

###########  Survey Index Fit and Residual Plots ###########;
allSL_vec$ECatch = 1000*exp(rep$Elog_index) 
allSL_vec$resid = rep$std_resid_index

allSL_vec$size = abs(allSL_vec$resid)
allSL_vec$clr = '+'
allSL_vec$clr[allSL_vec$resid<0] = '-'

gname = paste('figs//index_resid_bubble.jpeg',sep='')
jpeg(file=gname,width=6,height=5,units='in',res=300)
p1=ggplot(allSL_vec,aes(x=Year, y=Length, size=size, color=clr)) +
    geom_point(alpha=0.5) +               
    facet_wrap(~sname, ncol = 2,scales = "free_y") +         
    scale_size_continuous(range = c(.1, 4), name="Std. Resid")  +   
    theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")))+
    scale_color_manual(values=c("red", "blue"))  

print(p1)

dev.off()


gname = paste('figs//index_resid_3P.jpeg',sep='')
jpeg(file=gname,width=7,height=8,units='in',res=300)

p1 = ggplot(allSL_vec,aes(x=Year, y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') +
     stat_summary(fun=mean, aes(group=1), geom="line",size=1, color='red') + 
    facet_wrap(~sname, nrow = ns,scales = "free_y") + 
    theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")))+
     xlab('Year')+ ylab('Std. Residual')    
p2 = ggplot(allSL_vec,aes(x=Length, y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') +
     stat_summary(fun=mean, aes(group=1), geom="line",size=1) + 
    facet_wrap(~sname, nrow = ns,scales = "free_y") +   
    theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")))+
     xlab('Length')+ ylab('Std. Residual')     
p3 = ggplot(allSL_vec,aes(x=log(ECatch), y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') + 
    facet_wrap(~sname, nrow = ns,scales = "free_y") +   
    theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")))+
     xlab('Log Expected Catch')+ ylab('Std. Residual')
print(p1 | p2 | p3 )

dev.off()


for(i in 1:ns){

Sdat = allSL_vec[allSL_vec$sname==sname[i],]
# plot survey index at lengths, observed and model predicted;

gname = paste('figs//details//',sname[i],'_ridges.jpeg',sep='')
jpeg(file=gname,width=6,height=8,units='in',res=300)
p1=ggplot(Sdat)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = Catch, fill=as.factor(Year)),
                      stat = "identity", scale = 5,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1991:2019, expand=c(0,0))+
  scale_x_continuous(limits=c(10,40), breaks=c(10,15,20,25,30,35,40,45), expand=c(0,0))+ 
  xlab('Length (mm)')+
  ylab('Year')+    
  ggtitle(paste(sname[i],'observed Index'))+ 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=14))

p2=ggplot(Sdat)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = ECatch, fill=as.factor(Year)),
                      stat = "identity", scale = 5,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1991:2019, expand=c(0,0))+
  scale_x_continuous(limits=c(10,40), breaks=c(10,15,20,25,30,35,40,45), expand=c(0,0))+ 
  xlab('Length (mm)')+
  ylab('')+    
  ggtitle(paste(sname[i],'predicted Index'))+ 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),  
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5,size=14))

print(p1|p2)

dev.off()    

####### Index at length ##########

gname = paste('figs//details//log_',sname[i],'_at_len.jpeg',sep='')
jpeg(file=gname,width=6,height=8,units='in',res=300)    

p1=ggplot(Sdat) + geom_line(aes(x = Year, y = log(ECatch/1000)),size=1) +
    geom_point(aes(x=Year, y = log(Catch/1000))) +
    facet_wrap(~Length, ncol = 4,scales = "free_y") + geom_hline(yintercept = 0, linetype="dashed") +        
    xlab('Year')+
    ylab('Log Index')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
print(p1 + theme(axis.text.x = element_text(angle = 45)))

dev.off()     

gname = paste('figs//details//',sname[i],'_at_len.jpeg',sep='')
jpeg(file=gname,width=6,height=8,units='in',res=300)    

p1=ggplot(Sdat) + geom_line(aes(x = Year, y = ECatch/1000),size=1) +
    geom_point(aes(x=Year, y = Catch/1000)) +
    facet_wrap(~Length, ncol = 4,scales = "free_y") + geom_hline(yintercept = 0, linetype="dashed") +        
    xlab('Year')+
    ylab('Index (000s)')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
print(p1 + theme(axis.text.x = element_text(angle = 45)))

dev.off()   

}


####################  Table #################################;

np = length(opt$par)   
rnames = names(opt$par) 
ind = substr(rnames,1,5) == 'logit'           

out.tab = matrix(NA,np,2)
out.tab[,1] = exp(opt$par) 
out.tab[ind,1] = opt$par[ind] 
out.tab[,2] = sqrt(diag(sd.rep$cov.fixed))
rnames = sub("log_","",rnames) 

ind = substr(rnames,1,9) == 'std_index'
rnames[ind] = paste(rnames[ind],sname,sep="_")    
ind = substr(rnames,1,8) == 'std_logq'
rnames[ind] = paste(rnames[ind],sname,sep="_")    
ind = substr(rnames,1,2) == 'N0'
rnames[ind] = paste(rnames[ind],age.Fdev,sep="_")

rownames(out.tab)=rnames

t1 = 2*opt$objective+ 2*length(opt$par)
t2 = 2*opt$objective + log(length(rep$resid_index) + length(rep$resid_catch))*length(opt$par)
tname = paste("Model results, nll = ",round(opt$objective,digits=2),
', AIC = ',round(t1,2),", BIC = ",round(t2,2),sep="") 

tf = function(x){return(x)}

print(xtable(out.tab,digits=c(0,2,1),caption=tname,
align=c('l','r','r')),type='html',
  file='figs//model_output.doc',caption.placement="top",
  sanitize.rownames.function=tf)
    
vnames=c("log_biomass","log_ssb","log_Rec","log_harvest_rate")  

out.tab = matrix(NA,tmb_data$Y,8)
for(i in 1:4){

  ind = names(sd.rep$value) %in% vnames[i]

  out.tab[,2*(i-1) + 1] = exp(sd.rep$value[ind]) 
  out.tab[,2*(i-1) + 2] = sd.rep$sd[ind]
}
out.tab[,7:8]= 100*out.tab[,7:8]

cnames = c("biomass",'cv',"ssb",'cv',"Rec",'cv',"H(%)",'cv')
colnames(out.tab)=cnames
rownames(out.tab) = year

tname = "Stock Results"

tf = function(x){return(x)}

print(xtable(out.tab,digits=c(0,1,2,1,2,1,2,2,2),caption=tname,
align=c('l',rep('r',8))),type='html',
  file='figs//stock_table.doc',caption.placement="top",
  sanitize.rownames.function=tf)


