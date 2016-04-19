## ---- echo=FALSE---------------------------------------------------------
library(knitr)

opts_chunk$set(comment=NA, 
               fig.width =8, 
               fig.height=4,
               fig.path  ="../tex/",
               warning=FALSE, 
               message=FALSE, 
               error  =FALSE, 
               echo   =FALSE, 
               eval   =TRUE,
               cache  =TRUE)

## ---- echo=FALSE---------------------------------------------------------
library(FLCore)
library(ggplotFL)
library(plyr)
library(reshape)
library(diags)
library(R4MFCL)
library(kobe)

mf2FLPar=FLife:::mf2FLPar

theme_set(theme_bw(10))

dirMy ="/home/laurie/MEGAsync/mse/albn"
dirDat=paste(dirMy,"/data",sep="")
dirInp="/home/laurie/Desktop/rfmos/iccat/kobe/Inputs/albn/2013/mfcl"

#lfile=paste(dirInp,"/alt",1:7,"/length.fit",sep="")
#ffile=paste(dirInp,"/alt",1:7,"/albN.frq",  sep="")
#rfile=paste(dirInp,"/alt",1:7,"/plot.rep",  sep="")
load(paste(dirDat,"/om.RData",sep=""))

i=0

## ----eval=FALSE----------------------------------------------------------
## purl("/home/laurie/Desktop/scrs-2016/papers/scrs-2016-024/R/scrs-2016-024-figs.Rmd",
##      "/home/laurie/Desktop/scrs-2016/papers/scrs-2016-024/R/scrs-2016-024-figs.R")

## ----dgs-----------------------------------------------------------------
dgs=mdply(file.path(dirInp,tolower(names(om)),"plot.rep"),
          function(x) diags(x,"mfcl"))

flt=c("ESP_BBrec","EsFr_TR","EsFr_BBear","PRT_BB",
      "JPN_LLtrg","JPN_LLtra","JPN_LLbyc",
      "TAI_LL1","TAI_LL2","TAI_LL3",
      "KrPaCu_LL","Other_SU")

nms=c("Spain BB","Spain France\nTroll","Spain France\nEarly BB",
      "Portugal\nBB","Japan\nLL Target","Japan\nLL tra","Japan\nLL Bycatch",
      "Chinese-Taipei\nLL Early","Chinese-Taipei\nLL Mid","Chinese-Taipei\nLL Late",
      "KoreaPaCu LL","Other surface")

nm2=paste("flt",1:12,sep="")

dgs=transform(dgs,fleet=flt[name])

uCV=daply(subset(dgs,year>=1975),.(fleet),
      with, var(effDev.value))^.5

# rep=mlply(paste(dirInp,"/alt",1:7,"/plot.rep",sep=""),
#           function(x) read.rep(x))
# 
# caa=mlply(paste(dirInp,"/alt",1,sep=""),
#           function(x) {
#             rep=read.rep(paste(x,"plot.rep",sep="/"))
#             read.ests(rep,paste(x,"ests.rep",sep="/"))})
            
rep=read.rep("/home/laurie/Desktop/rfmos/iccat/kobe/Inputs/albn/2013/mfcl/alt1/plot.rep")
paa=read.ests(rep,"/home/laurie/Desktop/rfmos/iccat/kobe/Inputs/albn/2013/mfcl/alt1/ests.rep")[[1]]

paa=sweep(paa,1:2,apply(paa,1:2,sum),"/")

dat=ddply(melt(apply(paa[61:81,,],2:3,mean)), .(X2), 
             transform, sel=value/max(value))
dat=transform(dat,fleet=flt[X2])

## ----fig.height=4,fig.width=6--------------------------------------------
ggplot(dat)+
  geom_line(aes(X1,sel,col=fleet),size=1.5)+
  kobe:::theme_ms(12,legend.position="bottom")

paa=cast(dat[,c(1,4,5)],X1~fleet,value="sel")

save(paa,uCV,file=paste(dirDat,"/oem.RData",sep="/"))

## ----fig.height=7.25,fig.width=6-----------------------------------------
ggplot(dat)+
  geom_line(aes(X1,sel,col=fleet),size=1.5)+
  facet_wrap(~fleet,ncol=2)+
  kobe:::theme_ms(12,legend.position="none")

## ----ursd,fig.height=8.5,fig.width=6.5-----------------------------------
dgs=transform(dgs,Scenario=names(om)[X1])

ggplot(dgs)+
  geom_point(aes(year,residual,fill=Scenario),
             col="black",shape=21,size=1.25,lwd=0.25)  +
  geom_errorbar(aes(year,ymin=0,ymax=residual,fill=Scenario),size=0)+
  geom_smooth(aes(year,residual,col=Scenario),span=1,se=F)+
  facet_wrap(~name,ncol=2,scale="free")+
  theme_bw()+theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=2))

## ----uhat,fig.height=9,fig.width=7---------------------------------------
ggplot(dgs,aes(hat,residual,fill=Scenario))+
  geom_point(col="black",shape=21,size=1.25,lwd=0.5)  +
  geom_hline(aes(yintercept=0))+
  geom_smooth(aes(hat,residual,col=Scenario),span=1,se=F)+
  facet_wrap(~name,ncol=2,scale="free")+
  theme_bw()+theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=2))

## ----uar, fig.width=7,fig.height=8---------------------------------------
ggplot(dgs,aes(residual,residualLag,fill=Scenario))+
  geom_hline(aes(yintercept=0))+
  geom_point(col="black",shape=21,size=1.25,lwd=0.5)  +
  geom_smooth(aes(residual,residualLag,col=Scenario),method="lm",se=F)+
  facet_wrap(~name,ncol=3,scale="free")+
  theme_bw()+theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=2))

## ----qq,fig.height=7.25,fig.width=8--------------------------------------
qq=ddply(dgs,.(Scenario,name), with, 
         data.frame(qqnorm(stdz(residual),plot=F)))

ggplot(qq)                                        +
  geom_point( aes(x,y,  fill=Scenario),
              col="black",shape=21,size=1.25,lwd=0.5)            +
  stat_smooth(aes(x,y,group=Scenario),method="lm",
              se=T,fill="blue", alpha=0.1)        +
  facet_wrap(~name,scale="free",ncol=3)           +
  theme(legend.position="bottom")                 +
  theme_bw()  

## ----oem,fig.height=7.25,fig.width=6-------------------------------------
library(FLife)
source('~/Desktop/flr/git/biodyn/R/biodyn-oem.R')

stk=as.data.frame(window(stock(om[[1]]),start=1985)/
                    mean(window(stock(om[[1]]),start=1985)))

object=window(om[[1]],start=1985)
set.seed(7889)
cv    =exp(noise(100,FLQuant(0,dimnames=dimnames(m(object)[1,])),.3,0))
set.seed(7889)
ar    =exp(noise( 100,FLQuant(0,dimnames=dimnames(m(object)[1,])),.3,.75))
vr    =(cv-1)%*%FLQuant(seq(1.25,0.75,length.out=dim(cv)[2]),
                    dimnames=dimnames(cv[,,,,,1]))+1
trend =FLQuant(cumprod(1+rep(0.02,dim(fbar(object))[2])),
                             dimnames=dimnames(fbar(object)))
hyp   =FLQuant((stock(object)%/%stock(object)[,"1990"])^-0.9)
u     =FLQuants("Unbiased"      =oem(object,cv),
                "Hyperstability"=oem(object,cv,q=hyp),
                "Trend"         =oem(object,cv,q=trend),
                "AR"            =oem(object,ar),
                "Variable"      =oem(object,vr),
                "Juvenile"      =oem(object,cv,sel=mat(object)[,1]),
                "Mature"        =oem(object,cv,sel=1-mat(object)[,1]),
                "Numbers"       =oem(object,cv,mass=FALSE))

u=FLQuants(llply(u,function(x) x/mean(x)))

plot(u)+
  geom_line(aes(year,data,col=iter),
            data=subset(as.data.frame(u),iter%in%c(2,11)))+
  geom_line(aes(year,data),data=stk,col="blue")+
  facet_wrap(~qname,ncol=2)+
  theme_bw()+theme(legend.position="none")

## ----sel-----------------------------------------------------------------
sel=mf2FLPar(as.data.frame(paa))

save(sel,file=paste(dirDat,"/sel.RData",sep="/"))

