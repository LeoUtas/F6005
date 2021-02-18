library('ggplot2') 
require(reshape2) 
library(stringr)   
library('ggridges')  
library('gridExtra')    
library(patchwork)

library(viridis)
library(hrbrthemes) 



len.pop=1:54
age.pop=1:10
year = 1990:2019

## plot of proportion at length for each age;
jpeg(file='figs//PLA.jpeg',width=5,height=5,units='in',res=300)
  image(age.pop,len.pop,t(rep$pla),xlab="Age",ylab="Length (cm)",las=1)
  lines(age.pop,rep$mean_len)
dev.off()

# plot FRV index at lengths, observed and model predicted;
FRVc.vec$ECatch = 1000*exp(rep$Elog_index)

jpeg(file='figs//FLA_obs_ridges.jpeg',width=6,height=8,units='in',res=300)
p1=ggplot(FRVc.vec)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = Catch, fill=as.factor(Year)),
                      stat = "identity", scale = 5,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1991:2019, expand=c(0,0))+
  scale_x_continuous(limits=c(10,40), breaks=c(10,15,20,25,30,35,40,45), expand=c(0,0))+
  xlab('Length (mm)')+
  ylab('Year')+       
  ggtitle('Observed Index')+  
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        plot.margin=unit(c(0.25,1,0.5,0.25),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75))

print(p1)

dev.off()


jpeg(file='figs//FLA_pred_ridges.jpeg',width=6,height=8,units='in',res=300)
p1=ggplot(FRVc.vec)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = ECatch, fill=as.factor(Year)),
                      stat = "identity", scale = 5,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1991:2019, expand=c(0,0))+
  scale_x_continuous(limits=c(10,40), breaks=c(10,15,20,25,30,35,40,45), expand=c(0,0))+
  xlab('Length (mm)')+
  ylab('Year')+       
  ggtitle('Predicted Index')+  
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        plot.margin=unit(c(0.25,1,0.5,0.25),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75)) 

print(p1)

dev.off()

# plot Catch at lengths, observed and model predicted;
CLc.vec$ECatch = 1000*exp(rep$Elog_catch)

jpeg(file='figs//Catch_obs_ridges.jpeg',width=6,height=8,units='in',res=300)
p1=ggplot(CLc.vec)  +
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
        plot.title = element_text(hjust = 0.5,size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        plot.margin=unit(c(0.25,1,0.5,0.25),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75))  

print(p1)
dev.off()            

jpeg(file='figs//Catch_pred_ridges.jpeg',width=6,height=8,units='in',res=300)
p1=ggplot(CLc.vec)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = ECatch, fill=as.factor(Year)),
                      stat = "identity", scale = 3,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1990:2019, expand=c(0,0))+
  scale_x_continuous(limits=c(20,40), breaks=c(20,25,30,35,40,45), expand=c(0,0))+
  xlab('Length (mm)')+
  ylab('Year')+    
  ggtitle('Predicted')+ 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        plot.margin=unit(c(0.25,1,0.5,0.25),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75)) 

print(p1)
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

####### FRV at length ##########

gname <- 'figs//log_FRV_at_len.jpeg'
jpeg(file=gname,width=6,height=8,units='in',res=300)    

p1=ggplot(FRVc.vec) + geom_line(aes(x = Year, y = log(ECatch/1000)),size=1) +
    geom_point(aes(x=Year, y = log(Catch/1000))) +
    facet_wrap(~Length, ncol = 4,scales = "free_y") + geom_hline(yintercept = 0, linetype="dashed") +        
    xlab('Year')+
    ylab('Log Index')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
print(p1 + theme(axis.text.x = element_text(angle = 45)))

dev.off()     

gname <- 'figs//FRV_at_len.jpeg'
jpeg(file=gname,width=6,height=8,units='in',res=300)    

p1=ggplot(FRVc.vec) + geom_line(aes(x = Year, y = ECatch/1000),size=1) +
    geom_point(aes(x=Year, y = Catch/1000)) +
    facet_wrap(~Length, ncol = 4,scales = "free_y") + geom_hline(yintercept = 0, linetype="dashed") +        
    xlab('Year')+
    ylab('Index (000s)')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
print(p1 + theme(axis.text.x = element_text(angle = 45)))

dev.off()     

####### CL at length ##########

gname <- 'figs//log_catch_at_len.jpeg'
jpeg(file=gname,width=6,height=8,units='in',res=300)    

p1=ggplot(CLc.vec) + geom_line(aes(x = Year, y = log(ECatch/1000)),size=1) +
    geom_point(aes(x=Year, y = log(Catch/1000))) +
    facet_wrap(~Length, ncol = 4,scales = "free_y") + geom_hline(yintercept = 0, linetype="dashed") + 
    xlab('Year')+
    ylab('Log Catch')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
print(p1 + theme(axis.text.x = element_text(angle = 45)))

dev.off()      

gname <- 'figs//catch_at_len.jpeg'
jpeg(file=gname,width=6,height=8,units='in',res=300)    

p1=ggplot(CLc.vec) + geom_line(aes(x = Year, y = ECatch/1000),size=1) +
    geom_point(aes(x=Year, y = Catch/1000)) +
    facet_wrap(~Length, ncol = 4,scales = "free_y") + geom_hline(yintercept = 0, linetype="dashed") + 
    xlab('Year')+
    ylab('Catch (000s)')+       
    coord_cartesian(xlim = c(1991,2018))+
    theme(strip.text = element_text(size = 8, margin = margin()))
  
print(p1 + theme(axis.text.x = element_text(angle = 45)))

dev.off()



## plot of pop N;
jpeg(file='figs//Pop.jpeg',width=7,height=7,units='in',res=300)

par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2,1,0),cex=1)

  plot(year,exp(rep$log_Rec),type='l',xlab="Year",ylab="",lwd=2,las=1) 
  mtext(side=2,line=3,outer=F,c("Recruitment")) 
  plot(year,rep$ssb,type='l',xlab="Year",ylab="",lwd=2,las=1)   
  mtext(side=2,line=3,outer=F,c("SSB(Kt)")) 

  image(year,age.pop,t(rep$N_matrix),xlab="Year",ylab="Age",las=1)
  mtext(side=3,line=0.5,outer=F,"Numbers at age") 
  image(year,len.pop,t(rep$NLs),xlab="Year",ylab="Length (cm)",las=1)
  mtext(side=3,line=0.5,outer=F,"Numbers at length") 
  
dev.off()

value.names = names(sd.rep$value)

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

jpeg(file='figs//biomass.jpeg',width=5,height=4,units='in',res=300)

p1=ggplot(pdat, aes(x=year, y=est,fill = type,color=type)) +
    geom_line(aes(y = est, group=type, color=type),size=1) + 
    labs(x = "Year",y = "Biomass (Kt)")+
    geom_smooth(aes(ymin = L, ymax = U,fill=type,color=type),stat = "identity",alpha=0.2) +
    coord_cartesian(ylim = c(0,5000))
print(p1)
    
dev.off()       
       
###  Harvest rate plots with CIs ###

ind = value.names=='log_harvest_rate'
pdat = data.frame(year=year,est=exp(sd.rep$value[ind]))
pdat$L = exp(sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind]) 
pdat$U = exp(sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind])

jpeg(file='figs//harvest_rate.jpeg',width=5,height=4,units='in',res=300)

p1=ggplot(pdat, aes(x=year, y=est)) +
    geom_line(aes(y = est),size=1) + 
    labs(x = "Year",y = "Hrvest Rate")+
    geom_smooth(aes(ymin = L, ymax = U,color='grey'),stat = "identity",alpha=0.2) +
    coord_cartesian(ylim = c(0,0.03))  + guides(color = FALSE)
print(p1)
    
dev.off()  


## plot of Catch N;
jpeg(file='figs//Catch.jpeg',width=7,height=7,units='in',res=300)

par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2,1,0),cex=1) 

  image(year,age.pop,t(rep$F),xlab="",ylab="Age",las=1)
  mtext(side=3,line=0.5,outer=F,"F at age")          
  image(year,age.pop,t(rep$Z),xlab="",ylab="Age",las=1)
  mtext(side=3,line=0.5,outer=F,"Z at age") 
  image(year,age.pop,t(rep$CNA),xlab="Year",ylab="Age",las=1)
  mtext(side=3,line=0.5,outer=F,"Catch at age")  
  image(year,len.pop,t(rep$CL),xlab="Year",ylab="Length (cm)",las=1)
  mtext(side=3,line=0.5,outer=F,"Catch at Length")  
  
dev.off()


jpeg(file='figs//Q.jpeg',width=3,height=3,units='in',res=300)
  par(mar=c(4,4,1,1),mgp=c(2,1,0))
  plot(len.pop,exp(rep$logq),type='l',xlab="Length",ylab="",lwd=2,las=1)
  mtext(side=2,line=3,outer=F,"Catchability")   
dev.off()


CLc.vec$resid = rep$std_resid_catch
CLc.vec$size = abs(CLc.vec$resid)
CLc.vec$clr = '+'
CLc.vec$clr[rep$std_resid_catch<0] = '-'

jpeg(file='figs//Catch_resid_bubble.jpeg',width=6,height=4,units='in',res=300)
p1=ggplot(CLc.vec,aes(x=Year, y=Length, size=size, color=clr)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 4), name="Std. Resid") + 
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


FRVc.vec$resid = rep$std_resid_index
FRVc.vec$size = abs(FRVc.vec$resid)
FRVc.vec$clr = '+'
FRVc.vec$clr[rep$std_resid_index<0] = '-'

jpeg(file='figs//FRV_resid_bubble.jpeg',width=6,height=4,units='in',res=300)
p1=ggplot(FRVc.vec,aes(x=Year, y=Length, size=size, color=clr)) +
    geom_point(alpha=0.5) +
    scale_size_continuous(range = c(.1, 4), name="Std. Resid")  + 
    scale_color_manual(values=c("red", "blue"))  

print(p1)
dev.off()

jpeg(file='figs//FRV_resid_3P.jpeg',width=5,height=7,units='in',res=300)
p1 = ggplot(FRVc.vec,aes(x=Year, y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') +
     stat_summary(fun=mean, aes(group=1), geom="line",size=1) +
     xlab('Year')+ ylab('Std. Residual')    
p2 = ggplot(FRVc.vec,aes(x=Length, y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') +
     stat_summary(fun=mean, aes(group=1), geom="line",size=1) +
     xlab('Length')+ ylab('Std. Residual')     
p3 = ggplot(FRVc.vec,aes(x=log(ECatch), y=resid)) + geom_point(alpha=0.5) +
     geom_hline(aes(yintercept = 0), linetype='dashed') +
     xlab('Log Expected Catch')+ ylab('Std. Residual')
print(p1 / p2 / p3 )
dev.off()
