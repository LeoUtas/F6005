setwd("C:\\home\\CADIGAN\\GradProgram\\2021\\F6005\\3LN_redfish\\data")

library('ggplot2') 
require(reshape2) 
library(stringr)   
library('ggridges')  
library('gridExtra') 
                
vec_func = function(mdat){
  vdat = melt(mdat,id=c("Length"))
  vdat$Year = as.numeric(str_replace_all(vdat$variable,"X",""))
  vdat$variable=NULL
  vdat$Catch=vdat$value
  vdat$value=NULL
  return(vdat)
}

#landings;

landings = read.table('landings.txt',header=TRUE)
landings$Catch = landings$Catch/1000        
landings$TAC = landings$TAC/1000

gname <- 'landings_ts.jpeg'
jpeg(file=gname,width=4,height=4,units='in',res=300)

p1 = ggplot(landings)+
  geom_line(aes(x = Year, y = Catch),size=1)+ 
  geom_point(aes(x = Year, y = TAC),color='red')+
  labs(x="Year", y = "Landings (Kt)")
    
print(p1)
  
dev.off() 

################## Catch @ Length ####################################################; 

CL = read.table('C@L.txt',header=TRUE)
CL.vec = vec_func(CL)
CL.vec$symbol=19
CL.vec$symbol[CL.vec$Catch==0]=8
CL.vec$size = sqrt(CL.vec$Catch) 
CL.vec$size[CL.vec$Catch==0]=5 

jpeg(file='catch_bubbles.jpeg',width=7,height=7,units='in',res=300)

 ggplot(CL.vec, aes(x=Year, y=Length)) + 
    geom_point(aes(size = size, shape=factor(symbol)),alpha=0.6, show.legend = FALSE) +
    scale_shape_manual(values=c(8,19))+
    scale_size_continuous(range = c(0,8))  
  #  scale_y_continuous(breaks=1:15)              
  #  scale_x_continuous(breaks=seq(1990,2015,by=5))
       
dev.off()

jpeg(file='catch_ridges.jpeg',width=7,height=7,units='in',res=300)

ggplot(CL.vec)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = Catch, fill=as.factor(Year)),
                      stat = "identity", scale = 3,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1991:2019, expand=c(0,0))+
  scale_x_continuous(limits=c(10,45), breaks=c(10,15,20,25,30,35,40,45), expand=c(0,0))+
  xlab('Length (cm)')+
  ylab('Year')+  
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        plot.margin=unit(c(0.25,1,0.5,0.25),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75))
 
dev.off()


################## FALL RV @ Length ####################################################; 

FRV = read.table('Fall_RV.txt',header=TRUE)
FRV.vec = vec_func(FRV)
FRV.vec$symbol=19
FRV.vec$symbol[FRV.vec$Catch==0]=8
FRV.vec$size = sqrt(FRV.vec$Catch) 
FRV.vec$size[FRV.vec$Catch==0]=25 

jpeg(file='FRV_bubbles.jpeg',width=7,height=7,units='in',res=300)

ggplot(FRV.vec, aes(x=Year, y=Length)) + 
    geom_point(aes(size = size, shape=factor(symbol)),alpha=0.6, show.legend = FALSE) +
    scale_shape_manual(values=c(8,19))+
    scale_size_continuous(range = c(0,10))
    
dev.off()

jpeg(file='FRV_ridges.jpeg',width=7,height=7,units='in',res=300)

ggplot(FRV.vec)  +
  theme_ridges(grid=T,center_axis_labels=T,font_size = 10)+
  geom_density_ridges(aes(x = Length, y = Year, height = Catch, fill=as.factor(Year)),
                      stat = "identity", scale = 5,alpha= .8, show.legend = F)+    
  scale_y_continuous(breaks=1991:2019, expand=c(0,0))+
  scale_x_continuous(limits=c(4,45), breaks=c(10,15,20,25,30,35,40,45), expand=c(0,0))+
  xlab('Length (mm)')+
  ylab('Year')+  
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(hjust = 0.5,size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        plot.margin=unit(c(0.25,1,0.5,0.25),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75))
 
dev.off()

############### maturity ogives #######

mat = read.table('maturity.txt',header=TRUE)


gname <- 'maturity_ogives.jpeg'
jpeg(file=gname,width=5,height=4,units='in',res=300)

p1 = ggplot(mat)+
  geom_line(aes(x = Length, y = P1,color='black'),size=1)+ 
  geom_line(aes(x = Length, y = P2,color='red'),size=1)+ 
  labs(x="Length", y = "Female Proportion Mature") +
  scale_colour_manual(name = 'Time Periods', 
         values =c('black'='black','red'='red'), labels = c('1991-1995','1996+'))
    
print(p1)
  
dev.off() 

############### Weight ~ length #######

year = 1990:2019 
len = 1:55
pdat = expand.grid(year=year,length = len)

LW.parms = read.table('LW_parms.txt',header=TRUE)

pdat$a = rep(LW.parms$a,times=length(len)) 
pdat$b = rep(LW.parms$b,times=length(len))

pdat$weight = pdat$a*(pdat$length**pdat$b)/1000

library(viridis)
library(hrbrthemes)

personal_theme = theme(plot.title = element_text(hjust = 0.5),
  axis.title.x = element_text(hjust = 0.5),
  axis.title.y = element_text(hjust = 0.5))   
 
p1 = ggplot(pdat, aes(x = length, y = weight, color=as.factor(year))) +
    geom_line() +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Weight ~ Length") +
    theme_ipsum() +  
    labs(x="Length(cm)", y = "Weight (Kg)")

gname <- 'W-L.jpeg'
jpeg(file=gname,width=6,height=4,units='in',res=300)    
p1 + personal_theme
dev.off()

wt.mat = matrix(pdat$weight,
  nrow = length(len),
  ncol = length(year),
  byrow=TRUE)
rownames(wt.mat) = len
colnames(wt.mat)=year

wt.vec = pdat[,names(pdat) %in% c("year","length","weight")]


save(landings,CL,CL.vec,FRV,FRV.vec,mat,wt.mat,wt.vec,LW.parms,file='3LN_redfish.RData')



