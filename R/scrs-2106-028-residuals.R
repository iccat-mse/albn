## ----knitr, echo=FALSE---------------------------------------------------
library(knitr)

opts_chunk$set(comment=NA, 
               fig.width =6, 
               fig.height=6,
               fig.path  ="../tex/",
               warning=FALSE, 
               message=FALSE, 
               error  =FALSE, 
               echo   =FALSE, 
               eval   =TRUE,
               cache  =TRUE)

iFig=0

## ----init----------------------------------------------------------------
library(reshape)
library(plyr)

library(FLCore)
library(ggplotFL)
library(diags)
library(kobe)

library(mpb)

library(corrplot)
library(ggbiplot)
library(plotrix)

options(digits=4)

#source('~/Desktop/flr/git/kobe/R/kobe-quads.R')

theme_set(theme_bw(10))

dirMy  ="/home/laurie/Desktop/scrs-2016/papers/scrs-2016-028-albn-prelim"
dirDat =file.path(dirMy,"data")
dirInp =file.path(dirMy,"inputs")

load(file.path(dirDat,"bds.RData"))

dgs=ldply(bds, function(x) x@diags)
names(nms)=nm2
dgs=transform(dgs,name=nms[name])
names(dgs)[1]="Scenario"

## -----rsdl1,fig.height=4-------------------------------------------------
ggplot(ddply(dgs, .(Scenario,name), with, data.frame(obs=stdz(obs),hat=stdz(hat))))+
    geom_abline(aes(intercept=0,slope=1),col="orange")        +
    geom_point( aes(obs,hat),fill="salmon",col="grey",pch=21) +
    stat_smooth(aes(obs,hat),method="lm", se=F)               +
    theme_ms(14,legend.position="bottom")                     +
    xlab("Fitted") + ylab("Observed")                         +
    facet_grid(Scenario~name,scale="free")                    +
      theme_bw()

## -----rsdl2,fig.height=4,fig.width=6-------------------------------------
dat=ddply(subset(dgs,!is.na(residual)),.(Scenario,name), 
          transform, residual=stdz(residual,na.rm=T))

ggplot(aes(year,residual),data=subset(dat,!is.na(residual)))  +
  geom_hline(aes(yintercept=0),col="orange")                  +
  geom_point(fill="salmon",col="grey",pch=21)                 +
  stat_smooth(method="loess",se=F)                            +
  theme_ms(10,legend.position="bottom")                       +
  facet_grid(Scenario~name,scale="free",space="free")

## -----rsdl3,fig.height=4-------------------------------------------------
ggplot(dgs)                                           +
  geom_point( aes(qqx,qqy),fill="salmon",col="grey",pch=21)+
  stat_smooth(aes(qqx,qqHat),
              method="lm",se=T,fill="blue", alpha=0.1) +
  theme_ms(14,legend.position="bottom")                +
  facet_grid(Scenario~name)+
  xlab(expression(QQ[x]))+ylab(expression(QQ[y]))

## -----rsdl4,fig.height=5,fig.width=8-------------------------------------
ggplot(aes(hat, residual),
       data=subset(dgs,!is.na(hat) & !is.na(residual)))    +
  stat_quantile(quantiles = c(0.05, 0.5, 0.95),col="grey") +
  geom_hline(aes(yintercept=0),col="orange")               +
  geom_point(fill="salmon",col="grey",pch=21)              +
  facet_grid(Scenario~name,scale="free")                   +
  stat_smooth(method="loess",se=F,span=1)                  +
             theme_ms(14,legend.position="bottom")


## -----rsdl5,fig.height=4,fig.width=7-------------------------------------
ggplot(dgs)                                               +
  geom_point( aes(residual,residualLag),fill="salmon",col="grey",pch=21) +
  stat_smooth(aes(residual,residualLag),method="lm",se=F) +
  geom_hline(aes(yintercept=0),col="orange")              +
  xlab(expression(Residual[t]))                           + 
  ylab(expression(Residual[t+1]))                         +
  theme_ms(14,legend.position="bottom")                   +
  facet_grid(Scenario~name,scale="free")                  +
  theme_bw()  

## -----rsdl6,fig.height=4,fig.width=6-------------------------------------
ggplot(aes(year,obs),
       data=subset(ddply(dgs,.(Scenario,name), 
                         transform,obs=obs/mean(hat),
                                   hat=hat/mean(hat)),
                   !is.na(obs))) +
  geom_point(fill="salmon",col="grey",pch=21)+
  stat_smooth(method="loess",se=F)                +
  facet_grid(Scenario~name,scale="free",space="free")          +
  geom_line(aes(year,hat),col="black")   +
  scale_x_continuous(breaks=seq(0,3000,10))     +
  theme_bw()+
  theme_ms(10,legend.position="bottom")

## -----hat1,fig.height=4,fig.width=6--------------------------------------
ggplot(aes(year,obs),
       data=subset(ddply(dgs,.(Scenario,name), 
                         transform,obs=obs/mean(hat),
                                   hat=hat/mean(hat)),
                   !is.na(obs)))                      +
  geom_point(fill="salmon",col="grey",pch=21)         +
  geom_point(aes(year,obs))                           +
  stat_smooth(aes(year,obs),method="loess",se=F)      +
  geom_line(aes(year,hat),col="orange")               +
  facet_grid(Scenario~name,scale="free",space="free") +
  theme_bw()+
  theme_ms(10)+xlab("Year")+ylab("Stock")

## -----hat2,fig.height=5--------------------------------------------------

p=ldply(bds,function(x) model.frame(FLQuants(
      stock=FLQuant(seq(0,1,length.out=101)*params(x)["k"]), 
      yield=computePrd(x,FLQuant(seq(0,1,length.out=101)*params(x)["k"]))),drop=TRUE))
s=ldply(bds,function(x) model.frame(mcf(FLQuants(
       stock=stock(x),yield=catch(x)))))
  
ggplot(p)+
  geom_path(aes( stock,yield,col=.id),size=1.25)+
  geom_path(aes( stock,yield),data=s,col="grey")+
  geom_point(aes(stock,yield,fill=.id),data=s,shape=21,size=1.0)+
  theme_bw()+
  scale_x_continuous(limits=c(0,2.0e6))+
  scale_y_continuous(limits=c(0,6.5e4))+
  theme(legend.position="none")

## ------------------------------------------------------------------------
load("/home/laurie/Desktop/scrs-2016/papers/scrs-2016-028-albn-prelim/data/bds.RData")
iJk   =FLQuants(llply(cpue,jackknife))
jks   =rbind(
          expand.grid(run=1,idx=2:3),
          expand.grid(run=2,idx=4:5),
          expand.grid(run=3,idx=6))
key   =list(2:3,4:5,6)

jkSmry=NULL
jkTS=mdply(jks,function(run,idx){
  
  x  =bds[[run]]
  u  =cpue[key[[run]]]
  bd =fit(x,u)
  uJk=mpb:::FLQuantJKs(u)
  
  i       =seq(length(key[[run]]))[key[[run]]==idx]
  uJk[[i]]=iJk[[idx]]
  jk      =fit(bd,uJk)

  bias=jackSummary(params(bd)[c("r","k")],params(jk)[c("r","k")])
  bias=bias[[4]]/bias[[1]]
  jkSmry<<-rbind(jkSmry,cbind(run=run,cpue=idx,
                         as.data.frame(bias,drop=T)))
  
  ts=data.frame(model.frame(params(jk)[c("r","k")],drop=TRUE))
  ts=transform(ts,year=
            as.numeric(dimnames(u[[i]][!is.na(u[[i]])])$year))[,-3]
  
  return(ts)})

save(jkSmry,file=file.path(dirDat,"jkSmry.RData"))
save(jkTS,  file=file.path(dirDat,"jkTS.RData"))

## ----jk-r, fig.width=7,fig.height=4--------------------------------------
jkTS=transform(jkTS,name=nms[idx],Scenario=unique(dgs$Scenario)[run])

ggplot(ddply(jkTS,.(Scenario),transform, r=(r-mean(r,na.rm=T))/mean(r,na.rm=T)))+
  geom_hline(aes(yintercept=0),col="red")+
  geom_point(aes(year,r),fill="salmon",col="grey",pch=21)+
  facet_grid(Scenario~name,scale="free",space="free")+
  scale_y_continuous(labels=percent)+
  theme_bw()

## ----jk-k, fig.width=7,fig.height=4--------------------------------------
ggplot(ddply(jkTS,.(Scenario),
             transform, k=(k-mean(k,na.rm=T))/mean(k,na.rm=T)))+
  geom_hline(aes(yintercept=0),col="red")                      +
  geom_point(aes(year,k),fill="salmon",col="grey",pch=21)      +
  facet_grid(Scenario~name,scale="free",space="free")          +
  scale_y_continuous(labels=percent)                           +
  theme_bw()

## ------------------------------------------------------------------------
xv=function(assessment,bd,u){
  maxYr=range(bd)["maxyear"]
  ft   =fit(window(bd,end=assessment),u)
  if (assessment<maxYr)
   bd=mpb:::fwd(window(ft,end=maxYr),catch=window(catch(bd),start=assessment))
  
  model.frame(mcf(FLQuants(bd,
                  stock  =stock,
                  harvest=mpb:::harvest,
                  Bmsy   =function(x) stock(  x)%/%refpts(x)["bmsy"],
                  Fmsy   =function(x) mpb:::harvest(x)%/%refpts(x)["fmsy"])),
    drop=TRUE)}

xvs=rbind(cbind(Scenario=1,mdply(data.frame(assessment=2003:2011),xv,
                                 bd=bds[[1]],u=cpue[2:3])),
          cbind(Scenario=2,mdply(data.frame(assessment=2003:2011),xv,
                                 bd=bds[[2]],u=cpue[4:5])),
          cbind(Scenario=3,mdply(data.frame(assessment=2003:2011),xv,
                                 bd=bds[[3]],u=cpue[6])))

## ----xv-b----------------------------------------------------------------
xvs$Scenario=unique(dgs$Scenario)[xvs$Scenario]
ggplot(xvs)+
  geom_path(aes(year,stock,col=factor(assessment)))+
  facet_grid(Scenario~.,scale="free")+
  theme_bw()+
  theme(legend.position="none")+
  xlab("Year")+ylab("Stock")

## ----x-bmsy--------------------------------------------------------------

ggplot(xvs)+
  geom_path(aes(year,Bmsy,col=factor(assessment)))+
  facet_grid(Scenario~.,scale="free")+
  theme_bw()+
  theme(legend.position="none")+
  xlab("Year")+ylab(expression(B/B[msy]))

## ----x-f-----------------------------------------------------------------

ggplot(xvs)+
  geom_path(aes(year,harvest,col=factor(assessment)))+
  facet_grid(Scenario~.,scale="free")+
  theme_bw()+
  theme(legend.position="none")+
  xlab("Year")+ylab("Harvest Rate")

## ----xv-fmsy-------------------------------------------------------------

ggplot(xvs)+
  geom_path(aes(year,Fmsy,col=factor(assessment)))+
  facet_grid(Scenario~.,scale="free")+
  theme_bw()+
  theme(legend.position="none")+
  xlab("Year")+ylab(expression(F/F[msy]))

save(xvs,file=file.path(dirDat,"xvs.RData"))

## ----kobe, fig.height=8,fig.height=7-------------------------------------
library(kobe)
dat=subset(xvs,year>2002)

kobePhase(dat,xlim=c(0,2),ylim=c(0,2))+
  geom_path(aes(Bmsy,Fmsy,col=factor(assessment)))+
  facet_wrap(~Scenario,ncol=2)+
  theme(legend.position="bottom")

