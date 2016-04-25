## ----knitr, echo=FALSE---------------------------------------------------
library(knitr)

opts_chunk$set(comment=NA, 
               fig.width =6, 
               fig.height=5,
               fig.path  ="../tex/",
               warning=FALSE, 
               message=FALSE, 
               error  =FALSE, 
               echo   =FALSE, 
               eval   =TRUE,
               cache  =TRUE)

iFig=4

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

##Please change the dirs to those of your own choice
dirMy  ="/home/laurie/Desktop/scrs-2016/papers/scrs-2016-028-albn-prelim"
dirDat =file.path(dirMy,"data")
dirInp =file.path(dirMy,"inputs")

## Prettify names
nms=c("Chinese-Tai\n All","Chinese-Tai\n Early","Chinese-Tai\n Late",
                          "Japanese LL\n Early","Japanese LL\n Late", 
      "Troll Combined")
nm2=c("Chinese-Tai Early","Chinese-Tai Late","Chinese-Tai All",
      "Japanese LL Early","Japanese LL Late", 
      "Troll Combined")

load(file.path(dirDat,"fit.RData"))

## ------------------------------------------------------------------------
## K Profile
#source('~/Desktop/flr/git/mp/R/aspic-profile.R')
setControl(bds[[1]],min=0.25, max=2.25)=params(bds[[1]])
#bds[[1]]=fit(bds[[1]],cpue[2:3])

pK1=profile(bds[[1]],cpue[2:3],which="k",range=seq(0.70,1.3, length.out=201),
            run=FALSE)
pK1=fit(pK1,cpue[2:3])

## ------------------------------------------------------------------------
setControl(bds[[2]],min=0.5, max=2.0)=params(bds[[2]])

pK2=profile(bds[[2]],cpue[4:5],which="k",range=seq(0.75,1.25, length.out=201),
            run=FALSE)
pK2=fit(pK2,cpue[4:5])

## ------------------------------------------------------------------------
setControl(bds[[3]],min=0.25, max=2.0)=params(bds[[3]])

pK3=profile(bds[[3]],cpue[6],  which="k",range=seq(0.90,1.15, length.out=201),
            run=FALSE)
pK3=fit(pK3,cpue[6])

## ------------------------------------------------------------------------
pK=rbind(
data.frame(Scenario="1",CPUE=nm2[1],ss=c(mpb::::calcObjFn(pK1,cpue[2:3])[[1]]["ss"]),
           K=c(params(pK1)["k"]),r=c(params(pK1)["r"])),
data.frame(Scenario="1",CPUE=nm2[2],ss=c(mpb::::calcObjFn(pK1,cpue[2:3])[[2]]["ss"]),
           K=c(params(pK1)["k"]),r=c(params(pK1)["r"])),
data.frame(Scenario="2",CPUE=nm2[4],ss=c(mpb::::calcObjFn(pK2,cpue[4:5])[[1]]["ss"]),
           K=c(params(pK2)["k"]),r=c(params(pK2)["r"])),
data.frame(Scenario="2",CPUE=nm2[5],ss=c(mpb::::calcObjFn(pK2,cpue[4:5])[[2]]["ss"]),
           K=c(params(pK2)["k"]),r=c(params(pK2)["r"])),
data.frame(Scenario="3",CPUE=nm2[6],ss=c(mpb::::calcObjFn(pK3,cpue[  6])[[1]]["ss"]),
           K=c(params(pK3)["k"]),r=c(params(pK3)["r"])))

options(digits=3)
ggplot(pK)+
  geom_line(aes(K,ss,col=CPUE))+
  facet_wrap(~Scenario,scale="free_y",ncol=1)+
  theme_bw()

## ------------------------------------------------------------------------
ggplot(subset(pK))+
  geom_point(aes(K,r,col=CPUE))+
  #facet_wrap(~Scenario,scale="free_y",ncol=1)+
  scale_x_log10()+scale_y_log10()

## ---- pk1----------------------------------------------------------------
library(gam)
dat=transform(pK,k=log(K),r=log(r))
par=lm(r~k,data=dat)
dat=subset(pK,(residuals(par)^2)<0.001)
smh=data.frame(k=seq(11.8,13.75,length.out=201))
smh=cbind(r=predict(par,newdata=smh,type="response"),smh)

ggplot(smh)+
  geom_point(aes(k,r),data=transform(pK,k=log(K),r=log(r)),col="red")+
  geom_line(aes(k,r))
  
pK1=bds[[1]]

control(pK1)=propagate(control(pK1), dim(smh)[1])
params( pK1)=propagate(params( pK1), dim(smh)[1])
pK1@stock   =propagate(stock(  pK1), dim(smh)[1])

params(pK1)["r"]=exp(smh[,"r"])
params(pK1)["k"]=exp(smh[,"k"])

control(pK1)[,"val"]=params(pK1)
control(pK1)[,"min"]=params(pK1)*0.25
control(pK1)[,"max"]=params(pK1)*1.75
control(pK1)["r","phase"]= 2
control(pK1)["k","phase"]=-1

pK1=fwd(pK1,catch=catch(pK1),starvationRations=2.5)
pK1=fit(pK1,cpue[2:3])

## ---- pk2----------------------------------------------------------------
smh=data.frame(k=seq(11.5,12.75,length.out=201))
smh=cbind(r=predict(par,newdata=smh,type="response"),smh)

pK2=iter(bds[[2]],1)

control(pK2)=propagate(control(pK2), dim(smh)[1])
params( pK2)=propagate(params( pK2), dim(smh)[1])
pK2@stock   =propagate(stock(  pK2), dim(smh)[1])

params(pK2)["r"]=exp(smh[,"r"])
params(pK2)["k"]=exp(smh[,"k"])

control(pK2)[,"val"]=params(pK2)
control(pK2)[,"min"]=params(pK2)*0.5
control(pK2)[,"max"]=params(pK2)*2.9
control(pK2)["r","phase"]= 2
control(pK2)["k","phase"]=-1

pK2=mpb::::fwd(pK2,catch=catch(pK2),starvationRations=2.5)
pK2=fit(pK2,cpue[4:5])

## ------------------------------------------------------------------------
pK=rbind(
data.frame(Scenario="1",CPUE=nm2[1],ss=c(mpb::::calcObjFn(pK1,cpue[2:3])[[1]]["ss"]),
           K=c(params(pK1)["k"]),r=c(params(pK1)["r"])),
data.frame(Scenario="1",CPUE=nm2[2],ss=c(mpb::::calcObjFn(pK1,cpue[2:3])[[2]]["ss"]),
           K=c(params(pK1)["k"]),r=c(params(pK1)["r"])),
data.frame(Scenario="2",CPUE=nm2[4],ss=c(mpb::::calcObjFn(pK2,cpue[4:5])[[1]]["ss"]),
           K=c(params(pK2)["k"]),r=c(params(pK2)["r"])),
data.frame(Scenario="2",CPUE=nm2[5],ss=c(mpb::::calcObjFn(pK2,cpue[4:5])[[2]]["ss"]),
           K=c(params(pK2)["k"]),r=c(params(pK2)["r"])),
data.frame(Scenario="3",CPUE=nm2[6],ss=c(mpb::::calcObjFn(pK3,cpue[  6])[[1]]["ss"]),
           K=c(params(pK3)["k"]),r=c(params(pK3)["r"])))

options(digits=3)
ggplot(pK)+#subset(pK,ss<10))+
  geom_line(aes(K,ss,col=CPUE))+
  facet_wrap(~Scenario,scale="free_y",ncol=1)+
  guides(col = guide_legend(ncol=4))+
  theme_bw()+
  theme(legend.position="bottom")

## ------------------------------------------------------------------------
dat=transform(subset(pK,Scenario==2&ss<1.5),k=log(K),r=log(r))
par=lm(r~k,data=dat)
dat=subset(dat,(residuals(par)^2)<0.001)
smh=data.frame(k=seq(11.5,13.0,length.out=201))
smh=cbind(r=predict(par,newdata=smh,type="response"),smh)

ggplot()+
  geom_point(aes(k,r),data=transform(dat,k=k,r=r),col="red")+
  geom_line(aes(k,r),data=smh)

pK2=bds[[2]]

control(pK2)=propagate(control(pK2), dim(smh)[1])
params( pK2)=propagate(params( pK2), dim(smh)[1])
pK2@stock   =propagate(stock(  pK2), dim(smh)[1])

params(pK2)["r"]=exp(smh[,"r"])
params(pK2)["k"]=exp(smh[,"k"])

control(pK2)[,"val"]=params(pK2)
control(pK2)[,"min"]=params(pK2)*0.25
control(pK2)[,"max"]=params(pK2)*1.75
control(pK2)["r","phase"]= 2
control(pK2)["k","phase"]=-1

pK2=fwd(pK2,catch=catch(pK2),starvationRations=2.5)
pK2=fit(pK2,cpue[4:5])

## ------------------------------------------------------------------------
pK=rbind(
data.frame(Scenario="1",CPUE=nm2[1],ss=c(mpb::::calcObjFn(pK1,cpue[2:3])[[1]]["ss"]),
           K=c(params(pK1)["k"]),r=c(params(pK1)["r"])),
data.frame(Scenario="1",CPUE=nm2[2],ss=c(mpb::::calcObjFn(pK1,cpue[2:3])[[2]]["ss"]),
           K=c(params(pK1)["k"]),r=c(params(pK1)["r"])),
data.frame(Scenario="2",CPUE=nm2[4],ss=c(mpb::::calcObjFn(pK2,cpue[4:5])[[1]]["ss"]),
           K=c(params(pK2)["k"]),r=c(params(pK2)["r"])),
data.frame(Scenario="2",CPUE=nm2[5],ss=c(mpb::::calcObjFn(pK2,cpue[4:5])[[2]]["ss"]),
           K=c(params(pK2)["k"]),r=c(params(pK2)["r"])),
data.frame(Scenario="3",CPUE=nm2[6],ss=c(mpb::::calcObjFn(pK3,cpue[  6])[[1]]["ss"]),
           K=c(params(pK3)["k"]),r=c(params(pK3)["r"])))

options(digits=3)
ggplot(pK)+#subset(pK,ss<10))+
  geom_line(aes(K,ss,col=CPUE))+
  facet_wrap(~Scenario,scale="free_y",ncol=1)+
  guides(col = guide_legend(ncol=4))+
  theme_bw()+
  theme(legend.position="bottom")

## ------------------------------------------------------------------------
library(gam)
dat=transform(pK,k=log(K),r=log(r))
par=lm(r~k,data=dat)
dat=subset(pK,(residuals(par)^2)<0.001)
smh=data.frame(k=seq(12.0,14.0,length.out=201))
smh=cbind(r=predict(par,newdata=smh,type="response"),smh)

ggplot(smh)+
  geom_point(aes(k,r),data=transform(pK,k=log(K),r=log(r)),col="red")+
  geom_line(aes(k,r))

## ------------------------------------------------------------------------
ss=mpb::::calcObjFn(pK1,cpue[4:5])[[1]]["ss"]+mpb::::calcObjFn(pK1,cpue[4:5])[[2]]["ss"]
params(bds[[1]])[]=params(pK1)[,seq(dims(params(pK1))$iter)[ss==min(ss)]]
setControl(bds[[1]],min=0.75,max=1.25)=params(bds[[1]])
bds[[1]]=fit(bds[[1]],cpue[2:3])

ss=mpb::::calcObjFn(pK2,cpue[4:5])[[1]]["ss"]+mpb::::calcObjFn(pK2,cpue[4:5])[[2]]["ss"]
params(bds[[2]])[]=params(pK2)[,seq(dims(params(pK2))$iter)[ss==min(ss)]]
setControl(bds[[2]],min=0.75,max=1.25)=params(bds[[2]])
bds[[2]]=fit(bds[[2]],cpue[4:5])


save(bds,cpue,nm2,nms,file=file.path(dirDat,"bds.RData"))

