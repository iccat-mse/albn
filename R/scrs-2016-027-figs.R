## ----knitr_init, echo=FALSE, cache=FALSE---------------------------------
library(knitr)
library(rmdformats)

## Global options
options(max.print="75")
opts_chunk$set(echo       =FALSE,
               cache      =TRUE,
               cache.path ="../cache/validation",
               fig.path   ="../tex/",
               prompt     =FALSE,
               tidy       =TRUE,
               comment    =NA,
               message    =FALSE,
               warnings   =FALSE,
               fig.width  =7,
               fig.height =7)
opts_knit$set(width=75)

iFig=0

## ----eval=FALSE----------------------------------------------------------
## purl("/home/laurie/Desktop/scrs-2016/papers/scrs-2016-027/R/scrs-2016-027-figs.Rmd",
##      "/home/laurie/Desktop/scrs-2016/papers/scrs-2016-027/R/scrs-2016-027-figs.R")

## ----init----------------------------------------------------------------
library(stringr)
library(plyr)
library(reshape)
library(mpb)
library(ggplotFL)

library(scales)

##Please change the dirs to those of your own choice
dirMy="/home/laurie/Desktop/scrs-2016/papers/scrs-2016-027"
dirInp=file.path(dirMy,"inputs")
dirDat=file.path(dirMy,"data")
dirTex=file.path(dirMy,"tex")

# dirKobe="http://rscloud.iccat.int/kobe/Inputs"
# dirKobe="/home/laurie/Desktop/rfmos/iccat/kobe/Inputs"
dirKobe ="/home/laurie/MEGAsync/mse/albn/inputs/aspic"

#source('~/Desktop/flr/git/mp/R/aspic-constructors.R')

assessment=data.frame(stk   =c("albn","albs","yft"),
                      yr    =c(2013,2011,2011),
                      stringsAsFactors=F)
scen=rbind.fill(
   expand.grid(stk="albn",method="aspic",  run=paste("run",c(1:7),   sep=""),
                              stringsAsFactors=F)
  ,expand.grid(stk="albs",method="aspic",  run=paste("run",c(2,7),   sep=""),
                              stringsAsFactors=F) 
  ,expand.grid(stk="yft", method="aspic",  run=paste("run",c(9:12)[c(1,3)],  sep=""),
                              stringsAsFactors=F))

scen=merge(assessment,scen)
scen=arrange(scen,stk,method,run,yr)
scen=mdply(scen,function(stk,yr,method,run) 
  data.frame(dir=paste(dirKobe,stk,yr,method,run,sep="/")))

rm(assessment)

asp=aspics(file.path(scen$dir[1:7],"aspic.inp"))
est=cast(read.csv(file.path(dirMy,"inputs/estimates.tex"),sep=""),
           scen~stat,value="data")

library(doParallel)
library(foreach)

cl=makeCluster(4)
registerDoParallel(cl)

asp=fit(asp)

est=cbind(est[1:7,],
          ldply(asp[1:7],function(x)
              cbind(model.frame(params(x)[1:3]),
                    model.frame(refpts(x)),
                    model.frame(x@objFn)[,-3])))[,-c(6:7,10,14)]
names(est)[c(5:7,10)]=c("wg","msy2","k2","rerun")

save(asp,est,file=file.path(dirDat,"aspic.RData"),compress="xz")

# cpue
cpue=llply(asp[1:7],index,FALSE)
names(cpue)=names(asp[1:7])

# biodyn
bds=mpb:::aspics2biodyns(asp)

## ----bet-cpue,fig.height=8,fig.width=7-----------------------------------
library(gam)

mpb:::plotIndex(asp[1:7])

## ----bet-cor-------------------------------------------------------------
library(corrplot)

u=ldply(asp[1:7],index)
u=u[!duplicated(u[,c("name","year")]),]

cr=cor(cast(u,year~name,value="index")[,-1],
       use="pairwise.complete.obs")
dimnames(cr)=list(unique(u$name), unique(u$name))
cr[is.na(cr)]=0
corrplot(cr,diag=F,order="hclust",addrect=2)  +          
             theme(legend.position="bottom")  

## ----aspicAlbn-----------------------------------------------------------
plot(asp[1:7])+
  theme_bw()

## ----pK1-----------------------------------------------------------------
i=1
asp[[i]]=fit(asp[[i]])
pr1=profile(asp[[i]],which="k",
                range=seq(0.90,1.10,length.out=51))           

prK=cbind(scen="albn1",pr1)
p1=ggplot(pr1)+geom_line(aes(k,rss))+
      geom_point(aes(k,wg),   data=subset(est,scen=="albn1"),col="blue",size=4)+ 
      geom_point(aes(k,rerun),data=subset(est,scen=="albn1"),col="red", size=2)+
  theme_bw()+xlab("K")+ylab("Sum of Squares")
#p1

## ----pK2-----------------------------------------------------------------
i=2
params(asp[[i]])[c("msy","k")]=unlist(c(subset(est,scen=="albn2")[c("msy","k")]))
setControl(asp[[i]])=params(asp[[i]])
asp[[i]]=fit(asp[[i]])

pr2=profile(asp[[i]],which="k",
                range=seq(0.50,1.50,length.out=101))           

prK=rbind.fill(prK,cbind(scen="albn2",pr2))
p2=ggplot(subset(pr2,rss<5.5))+
    geom_line(aes(k,rss))+
    geom_point(aes(k,wg),   data=subset(est,scen=="albn2"),col="blue",size=4)+ 
    geom_point(aes(k,rerun),data=subset(est,scen=="albn2"),col="red", size=2)+
  theme_bw()+xlab("K")+ylab("Sum of Squares")
#p2

## ----pK3-----------------------------------------------------------------
i=3
params(asp[[i]])[c("msy","k")]=unlist(c(subset(est,scen=="albn3")[c("msy","k")]))
setControl(asp[[i]])=params(asp[[i]])
asp[[i]]=fit(asp[[i]])
pr3=profile(asp[[i]],which="k",
                range=seq(0.90,1.50,length.out=51))           

prK=rbind.fill(prK,cbind(scen="albn3",pr3))
p3=ggplot(pr3)+geom_line(aes(k,rss))+
    geom_point(aes(k,wg),   data=subset(est,scen=="albn3"),col="blue",size=4)+ 
    geom_point(aes(k,rerun),data=subset(est,scen=="albn3"),col="red", size=2)+
  theme_bw()+xlab("K")+ylab("Sum of Squares")
#p3

## ----pK4-----------------------------------------------------------------
i=4
params(asp[[i]])[c("msy","k")]=as.numeric(subset(est,scen=="albn4",c("msy","k")))
setControl(asp[[i]])=params(asp[[i]])
asp[[i]]=fit(asp[[i]])
pr4=profile(asp[[i]],which="k",
                range=seq(0.90,1.20,length.out=101),min=0.9,max=1.1)           

p4=ggplot(pr4)+
    geom_line(aes(k,rss),col="grey")+
    geom_line(aes(k,rss,col=Level),
           data=transform(pr4,Level=ifelse(rss<6,"Low","High")))+
    geom_point(aes(k,wg),data=subset(est,scen=="albn4"),col="blue",size=4)+ 
    geom_point(aes(k,rss),col="red",
                       data=cbind(model.frame(asp[[i]]@objFn["rss"]),
                                             model.frame(params(asp[[i]])["k"])))+
  theme_bw()+xlab("K")+ylab("Sum of Squares")
#p4

## ----pK5-----------------------------------------------------------------
i=5
params(asp[[i]])[c("msy","k")]=unlist(c(subset(est,scen=="albn5")[c("msy","k")]))
setControl(asp[[i]])=params(asp[[i]])
asp[[i]]=fit(asp[[i]])
pr5=profile(asp[[i]],which="k",
                range=seq(0.75,1.25,length.out=51),min=0.75,max=1.5)           

prK=rbind.fill(prK,cbind(scen="albn5",pr5))
p5=ggplot(pr5)+
    geom_line(aes(k,rss))+
    geom_point(aes(k,wg),   data=subset(est,scen=="albn5"),col="blue",size=4)+ 
    geom_point(aes(k,rerun),data=subset(est,scen=="albn5"),col="red", size=2)+
  theme_bw()+xlab("K")+ylab("Sum of Squares")
#p5

## ----pK6-----------------------------------------------------------------
i=6
params(asp[[i]])[c("msy","k")]=unlist(c(subset(est,scen=="albn6")[c("msy","k")]))
setControl(asp[[i]])=params(asp[[i]])
asp[[i]]=fit(asp[[i]])
pr6=profile(asp[[i]],which="k",
                range=seq(0.8,1.2,length.out=51),min=0.75,max=1.5)           

prK=rbind.fill(prK,cbind(scen="albn6",pr6))
p6=ggplot(pr6)+
    geom_line(aes(k,rss))+
    geom_point(aes(k,wg),   data=subset(est,scen=="albn6"),col="blue",size=4)+ 
    geom_point(aes(k,rerun),data=subset(est,scen=="albn6"),col="red", size=2)+
  theme_bw()+xlab("K")+ylab("Sum of Squares")
#p6

## ----pK7-----------------------------------------------------------------
i=7
params(asp[[i]])[c("msy","k")]=unlist(c(subset(est,scen=="albn7")[c("msy","k")]))
setControl(asp[[i]])=params(asp[[i]])
asp[[i]]=fit(asp[[i]])
pr7=profile(asp[[i]],which="k",
                range=seq(0.9,1.2,length.out=51),min=0.75,max=1.5)           

prK=rbind.fill(prK,cbind(scen="albn7",pr7))
p7=ggplot(pr7)+
    geom_line(aes(k,rss))+
    geom_point(aes(k,wg),   data=subset(est,scen=="albn7"),col="blue",size=4)+ 
    geom_point(aes(k,rerun),data=subset(est,scen=="albn7"),col="red", size=2)+
  theme_bw()+xlab("K")+ylab("Sum of Squares")
#p7

## ----pAlbn, fig.height=8,fig.width=6-------------------------------------
library(kobe)

kobe:::multiplot(p1+ylab("LAV")+
                   scale_x_continuous(limits=c(5e5,15e5),breaks=c(5e5,10e5,15e5)),
                 p2+ylab("LAV")+
                   scale_x_continuous(limits=c(5e5,15e5),breaks=c(5e5,10e5,15e5)),
                 p3+ylab("LAV")+
                   scale_x_continuous(limits=c(5e5,15e5),breaks=c(5e5,10e5,15e5)),
                 p4+ylab("LAV")+scale_x_continuous(limits=c(3.25e5,4.25e5))+
                        theme(legend.position="none"),
                 p5+ylab("LAV")+
                   scale_x_continuous(limits=c(5e5,15e5),breaks=c(5e5,10e5,15e5)),
                 
                 p6+ylab("LAV")+
                   scale_x_continuous(limits=c(5e5,15e5),breaks=c(5e5,10e5,15e5)),
                 p7+ylab("LAV")+
                   scale_x_continuous(limits=c(5e5,15e5),breaks=c(5e5,10e5,15e5)),
                 p4+ylab("LAV")+scale_x_continuous(limits=c(3.25e5,4.25e5))+
                               coord_cartesian(ylim=c(4.6,5.0))+
                        theme(legend.position="none"),
                 cols=2)

## ----fwd-----------------------------------------------------------------
cf1=ldply(asp,function(x){
  x  =fit(x)
  res=mpb:::aspic2biodyn(x)
  control(res)["p",       "phase"]=-1
  res=fwd(res,catch=catch(res))
  
  ##time series
  cbind(Method="biodyn",model.frame(mcf(FLQuants(res,"stock"))))
  })

## ----fwd2----------------------------------------------------------------
#source('~/Desktop/flr/git/mp/R/biodyn-F.R')
cf2=ldply(asp,function(x){
  x  =fit(x)
  res=mpb:::aspic2biodyn(x)
  
  stock=FLQuant(c(params(res)["k"]*params(res)["b0"]),
                  dimnames=dimnames(catch(x)))
  res@stock=mpb:::nr(catch(x),stock,
               params(res)["r"],
               params(res)["k"],
               params(res)["b0"],
               tolVal=1e-6,niter=20,yieldFlag=TRUE)[["B"]]

  ##time series
  cbind(Method="biodyn",model.frame(mcf(FLQuants(res,"stock"))))
  })

## ----fit-----------------------------------------------------------------
cf3=ldply(asp,function(x){
  x  =fit(x)
  res=mpb:::aspic2biodyn(x,
                  phase=c("b0"=-1,"r"=2,"k"=-1,"p"=-1,"q"=2,"sigma"=1),
                  min=0.5,max=2)
  res=fit(res,index(x,FALSE))
  
  ##time series
  cbind(Method="biodyn",model.frame(mcf(FLQuants(res,"stock"))))
  })

## ----fit2----------------------------------------------------------------
cf4=ldply(asp,function(x){
  x  =fit(x)
  res=mpb:::aspic2biodyn(x)
  control(res)["p", "phase"]=-2
  
  control(res)["r", "phase"]=2
  control(res)["k", "phase"]=2
  control(res)["k",c("min","max")]=control(res)["k","val"]*c(0.9,1.25)
  #control(res)["r",c("min","max")]=control(res)["r","val"]*c(0.8,2.5)
  res=fit(res,index(x,FALSE))
  
  ##time series
  cbind(Method="biodyn",model.frame(mcf(FLQuants(res,"stock"))))
  })

## ----cf------------------------------------------------------------------
cfAsp=ldply(asp,function(x){
  x  =fit(x)
  
  ##time series
  cbind(Method="aspic", model.frame(mcf(FLQuants(x,  "stock"))))})

cf=rbind.fill(cbind("Type"="Projection","Form"="Rate",    cf1),
              cbind("Type"="Projection","Form"="Logistic",cf2),
              cbind("Type"="Fit",       "Form"="Rate",    cf3),
              cbind("Type"="Fit",       "Form"="Logistic",cf4))
rm(cf1,cf2,cf3,cf4)

## ----chk-----------------------------------------------------------------
ggplot(subset(cf,Method=="biodyn"&Type!="Fit"))+
  geom_path(aes(year,stock),data=cfAsp)+
  geom_path(aes(year,stock),colour="red")+
  facet_grid(.id~Form)+xlab("Year")+ylab("Stock Biomass")+
  theme_bw()

## ----chk2----------------------------------------------------------------
ggplot(subset(cf,Method=="biodyn"&Type=="Fit"))+
  geom_path(aes(year,stock),data=cfAsp)+
  geom_path(aes(year,stock),colour="red")+
  facet_grid(.id~Form)+xlab("Year")+ylab("Stock Biomass")+
  theme_bw()

## ----lav-----------------------------------------------------------------
lav=ldply(asp, function(x){
  x  =fit(x)
  res=mpb:::aspic2biodyn(x)
  res=fit(res,index(x,FALSE))
  
  cbind(aspic =sum(ldply(mpb:::calcObjFn(x), function(x) x["lav"])),
        biodyn=sum(ldply(mpb:::calcObjFn(res,index(x,FALSE)), function(x) x["lav"])))

  #ddply(x@diags,  .(name),with,sum(abs(residual), na.rm=TRUE))
  #ddply(res@diags,.(name),with,sum(abs(residual), na.rm=TRUE))
  })

## ------------------------------------------------------------------------
names(lav)[1]="Run"

library(xtable)

scen=ldply(asp[1:7],function(run) data.frame(cpue=ac(unique(index(run)$name))))
names(scen)[1]="Scenario"
#dimnames(scen)[[1]]=NULL
print(
  xtable(scen,
         caption="2013 Assessment Scenarios"),
  file=file.path(dirTex,"table1.tex"))

#dimnames(lav)[[1]]=NULL
print(
  xtable(lav,
         caption="Objective function values (least absolute value)"),
  file=file.path(dirTex,"table2.tex"))

## ---- eval=FALSE---------------------------------------------------------
## fn=function(x) cbind(model.frame(params(x)),
##                      ll     =model.frame(x@ll),
##                      model.frame(refpts(x))[,-4],
##                      stock  =c(stock(  x)[,ac(range(x)['maxyear'])]%/%bmsy(x)),
##                      harvest=c(harvest(x)[,ac(range(x)['maxyear'])]%/%fmsy(x))),
## 
## pK=profile(res,index(x,FALSE),which="k",range=seq(0.9,1.1,length.out=51),run=FALSE)
## pK=fit(pK,index(x,FALSE))

## ------------------------------------------------------------------------
## Assessment Scenarios
scen=ldply(asp[1:7],function(run) data.frame(cpue=ac(unique(index(run)$name))))
names(scen)[1]="run"
scen$cpue=ac(scen$cpue)

## Indices
idx=ldply(asp[1:7],index)
idx=idx[!duplicated(idx[,c("name","year")]),c("name","year","index")]
idx=FLQuants(dlply(idx,.(name), with,as.FLQuant(data.frame(year=year,data=index))))
iJk=FLQuants(llply(idx,jackknife))

save(cpue,asp,scen,idx,file="/home/laurie/Desktop/temp/tmp.RData")

## ------------------------------------------------------------------------
#load("/home/laurie/Desktop/temp/tmp.RData")
names(cpue)=unique(scen$cpue)

jkSmry<-NULL
jkTS=mdply(scen,function(run,cpue){
  #run =scen[7,1];cpue=scen[7,2]

  x      =asp[[run]]
  iU     =unique(index(x)$name)
  u      =idx[iU]
  names(u)=iU

  bd     =mpb:::aspic2biodyn(x,"biodyn")
  bd     =fit(bd,u)

  uJk     =mpb:::FLQuantJKs(u)
  names(uJk)=names(u)
  
  uJk[[cpue]]=iJk[[cpue]]
  names(uJk)=iU
  jk     =fit(bd,uJk)

#   nms=c("r","k",
#         "cnow","bnow","fnow",
#         "bthen","fthen",
#         "bnowthen","fnowthen",
#         "msy","bmsy","fmsy","cmsy",
#         "bbmsy","ffmsy",
#         "bk","fr",
#         "bratio","fratio",
#         "slopeb","slopef")  
#   
#   bias=jackSummary(bd@mng[nms,1,],jk@mng[nms,1,])
#   bias=bias[[4]]/bias[[1]]
   bias=jackSummary(params(bd)[c("r","k")],params(jk)[c("r","k")])
   bias=bias[[4]]/bias[[1]]
   jkSmry<<-rbind(jkSmry,cbind(run=run,cpue=cpue,
                         as.data.frame(bias,drop=T)))
  
  ts=data.frame(model.frame(params(jk)[c("r","k")],drop=TRUE))
  ts=transform(ts,year=
            as.numeric(dimnames(u[[ac(cpue)]][!is.na(u[[ac(cpue)]])])$year))[,-3]
  
  return(ts)})

save(jkSmry,file=file.path(dirDat,"jkSmry.RData"))
save(jkTS,  file=file.path(dirDat,"jkTS.RData"))

## ----jk-r, fig.width=12,fig.height=6-------------------------------------
ggplot(ddply(jkTS,.(run),transform, r=(r-mean(r,na.rm=T))/mean(r,na.rm=T)))+
  geom_hline(aes(yintercept=0),col="red")+
  geom_point(aes(year,r))+
  facet_grid(run~cpue,scale="free",space="free")+
  scale_y_continuous(labels=percent)+
  theme_bw()

## ----jk-k, fig.width=12,fig.height=6-------------------------------------
ggplot(ddply(jkTS,.(run),transform, k=(k-mean(k,na.rm=T))/mean(k,na.rm=T)))+
  geom_hline(aes(yintercept=0),col="red")+
  geom_point(aes(year,k))+
  facet_grid(run~cpue,scale="free",space="free")+
  scale_y_continuous(labels=percent)+
  theme_bw()

## ------------------------------------------------------------------------
xv=ldply(asp[1:7],function(x) mdply(data.frame(assessment=2003:2011), 
                    function(assessment,bd,cpue){
  bd. =window(bd,end=assessment)
  
  bd. =fit(bd.,cpue)
  if (range(bd)["maxyear"]>assessment)
   bd. =fwd(bd.,catch=window(catch(bd),start=assessment))
  
  model.frame(mcf(FLQuants(bd.,
                  stock  =stock,
                  harvest=harvest,
                  Bmsy   =function(x) stock(  x)%/%refpts(x)["bmsy"],
                  Fmsy   =function(x) harvest(x)%/%refpts(x)["fmsy"])),
    drop=TRUE)
  
  },
    bd=mpb:::aspic2biodyn(x,"biodyn"),cpue=index(x,FALSE)))

## ----xv-b----------------------------------------------------------------
ggplot(xv)+
  geom_path(aes(year,stock,col=factor(assessment)))+
  facet_wrap(~.id,ncol=2)+
  theme_bw()+theme(legend.position="none")+
  xlab("Year")+ylab("Stock")

## ----x-bmsy--------------------------------------------------------------
ggplot(xv)+
  geom_path(aes(year,Bmsy,col=factor(assessment)))+
  facet_wrap(~.id,ncol=2)+
  theme_bw()+theme(legend.position="none")+
  xlab("Year")+ylab(expression(B/B[msy]))

## ----x-f-----------------------------------------------------------------
ggplot(xv)+
  geom_path(aes(year,harvest,col=factor(assessment)))+
  facet_wrap(~.id,ncol=2,scale="free")+
  theme_bw()+theme(legend.position="none")+
  xlab("Year")+ylab("Harvest Rate")

## ----xv-fmsy-------------------------------------------------------------
ggplot(xv)+
  geom_path(aes(year,Fmsy,col=factor(assessment)))+
  facet_wrap(~.id,ncol=2)+
  theme_bw()+theme(legend.position="none")+
  xlab("Year")+ylab(expression(F/F[msy]))

## ----kobe, fig.height=9,fig.height=8-------------------------------------
library(kobe)
dat=subset(xv,year>2002)

kobePhase(dat,xlim=c(0,2),ylim=c(0,2))+
  geom_path(aes(Bmsy,Fmsy,col=factor(assessment)))+
  facet_wrap(~.id,ncol=3)+
  theme(legend.position="bottom")

## ------------------------------------------------------------------------
names(cpue)=unique(scen$cpue)

## Assessment Scenarios
scen=ldply(asp[1:7],function(run) data.frame(cpue=ac(unique(index(run)$name))))
names(scen)[1]="run"
scen$cpue=ac(scen$cpue)

## Indices
idx=ldply(asp[1:7],index)
idx=idx[!duplicated(idx[,c("name","year")]),c("name","year","index")]
idx=FLQuants(dlply(idx,.(name), with,as.FLQuant(data.frame(year=year,data=index))))
iJk=FLQuants(llply(idx,jackknife))

jkSmry<-NULL
jkTS=mdply(scen,function(run,cpue){
  #run =scen[7,1];cpue=scen[7,2]

  print(paste(run,cpue))
  x      =asp[[run]]
  iU     =unique(index(x)$name)
  u      =idx[iU]
  names(u)=iU

  bd     =mpb:::aspic2biodyn(x)
  bd     =fit(bd,u,lav=TRUE)

  uJk     =mpb:::FLQuantJKs(u)
  names(uJk)=names(u)
  
  uJk[[cpue]]=iJk[[cpue]]
  names(uJk)=iU
  jk     =fit(bd,uJk,lav=TRUE)

#   nms=c("r","k",
#         "cnow","bnow","fnow",
#         "bthen","fthen",
#         "bnowthen","fnowthen",
#         "msy","bmsy","fmsy","cmsy",
#         "bbmsy","ffmsy",
#         "bk","fr",
#         "bratio","fratio",
#         "slopeb","slopef")  
#   
#   bias=jackSummary(bd@mng[nms,1,],jk@mng[nms,1,])
#   bias=bias[[4]]/bias[[1]]
   bias=jackSummary(params(bd)[c("r","k")],params(jk)[c("r","k")])
   bias=bias[[4]]/bias[[1]]
   jkSmry<<-rbind(jkSmry,cbind(run=run,cpue=cpue,
                         as.data.frame(bias,drop=T)))
  
  ts=data.frame(model.frame(params(jk)[c("r","k")],drop=TRUE))
  ts=transform(ts,year=
            as.numeric(dimnames(u[[ac(cpue)]][!is.na(u[[ac(cpue)]])])$year))[,-3]
  
  return(ts)})

save(jkSmry,file=file.path(dirDat,"jkSmry.RData"))
save(jkTS,  file=file.path(dirDat,"jkTS.RData"))

## ----jk-k-lav, fig.width=12,fig.height=6---------------------------------
ggplot(ddply(jkTS,.(run),transform, k=(k-mean(k,na.rm=T))/mean(k,na.rm=T)))+
  geom_hline(aes(yintercept=0),col="red")+
  geom_point(aes(year,k))+
  facet_grid(run~cpue,scale="free",space="free")+
  scale_y_continuous(labels=percent)+
  theme_bw()

## ------------------------------------------------------------------------
xv=ldply(asp[1:7],function(x) mdply(data.frame(assessment=2003:2011), 
                    function(assessment,bd,cpue){
  bd. =window(bd,end=assessment)
  
  bd. =fit(bd.,cpue,lav=TRUE)
  if (range(bd)["maxyear"]>assessment)
   bd. =fwd(bd.,catch=window(catch(bd),start=assessment))
  
  model.frame(mcf(FLQuants(bd.,
                  stock  =stock,
                  harvest=harvest,
                  Bmsy   =function(x) stock(  x)%/%refpts(x)["bmsy"],
                  Fmsy   =function(x) harvest(x)%/%refpts(x)["fmsy"])),
    drop=TRUE)
  
  },
    bd=mpb:::aspic2biodyn(x),cpue=index(x,FALSE)))

## ----xv-b-lav------------------------------------------------------------
ggplot(xv)+
  geom_path(aes(year,stock,col=factor(assessment)))+
  facet_wrap(~.id,ncol=2)+
  theme_bw()+theme(legend.position="none")+
  xlab("Year")+ylab("Stock")

## ----xv-bmsy-lav---------------------------------------------------------
ggplot(xv)+
  geom_path(aes(year,Bmsy,col=factor(assessment)))+
  facet_wrap(~.id,ncol=2)+
  theme_bw()+theme(legend.position="none")+
  xlab("Year")+ylab(expression(B/B[msy]))

## ----xv-f-lav------------------------------------------------------------
ggplot(xv)+
  geom_path(aes(year,harvest,col=factor(assessment)))+
  facet_wrap(~.id,ncol=2,scale="free")+
  xlab("Year")+ylab("Harvest Rate")+
  theme_bw()+theme(legend.position="none")

## ----xv-fmsy-lav---------------------------------------------------------
ggplot(xv)+
  geom_path(aes(year,Fmsy,col=factor(assessment)))+
  facet_wrap(~.id,ncol=2)+
  theme_bw()+theme(legend.position="none")+
  xlab("Year")+ylab(expression(F/F[msy]))

## ----xv-kobe-lav, fig.height=9,fig.height=8------------------------------
library(kobe)
dat=subset(xv,year>2002)

kobePhase(dat,xlim=c(0,2),ylim=c(0,2))+
  geom_path(aes(Bmsy,Fmsy,col=factor(assessment)))+
  facet_wrap(~.id,ncol=3)+
  theme(legend.position="bottom")

