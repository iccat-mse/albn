## ----knitr, echo=FALSE---------------------------------------------------
library(knitr)

opts_chunk$set(comment=NA, 
               fig.width =6, 
               fig.height=6,
               fig.path  ="../tex/",
               warnings=FALSE, 
               message =FALSE, 
               error   =FALSE, 
               echo    =FALSE, 
               eval    =TRUE,
               cache   =TRUE)

## ----init, echo=FALSE----------------------------------------------------
library(FLCore)
library(FLBRP)
library(mpb)
library(popbio)
library(gam)

##Please change the dirs to those of your own choice
dirAlb="/home/laurie/MEGAsync/mse/albn"
dirTex="/home/laurie/Desktop/scrs-2016/papers/scrs-2016-026-albn-sa/tex"
dirDat="/home/laurie/Desktop/scrs-2016/papers/scrs-2016-026-albn-sa/data"

load(paste(dirAlb,"data","om.RData",sep="/"))

#source('~/Desktop/flr/git/mp/R/biodyn-oem.R')

options(digits=3)

nits=101

iFig=0

## ------------------------------------------------------------------------
##with fit
chkFit<-function(i,bd,u,nits=101){
  pK=profile(bd[[i]],u[[i]],which="k",range=seq(0.5,2.0,length.out=nits),run=FALSE)
  pK=fit(pK,u[[i]])

  dat=data.frame(ss=c(mpb:::calcObjFn(pK,u[[i]])[[1]]["ss"]),
                 k =c(params(pK)["k"]),
                 r =c(params(pK)["r"]))

  dat}

##without fit
chk<-function(bd,u){
  data.frame(ss=c(mpb:::calcObjFn(bd,u)[[1]]["ss"]),
             k =c(params(bd)["k"]),
             r =c(params(bd)["r"]))}

##refit
refit<-function(b,u,d,fit=FALSE){
  if (dim(params( b))[2]==1)
    params( b)=propagate(params( b),nits)
  if (dim(control( b))[3]==1)
    control(b)=propagate(control(b),nits)
  b@stock   =propagate(stock(  b),nits)
  b@catch   =propagate(catch(  b),nits)
  #b         =propagate(b,nits)
  
  control(b)["sigma1","phase"]=2
  control(b)["r","phase"]=3
  control(b)["r","val"]=d$r
  control(b)["r","min"]=control(b)["r","val"]*0.75
  control(b)["r","max"]=control(b)["r","val"]*1.75
  control(b)["k","phase"]=-1
  control(b)["k","val"]=d$k
  control(b)["k","min"]=control(b)["k","val"]*1.0
  control(b)["k","max"]=control(b)["k","max"]*1.1
  
  params(b)[]=control(b)[,"val"]
  b=mpb:::fwd(b,catch=catch(b))
  
  if (fit)
    b=fit(b,u)
  b}
  
gamRK<-function(r,k,ss) {
        dat=data.frame(r=log(r),k=log(k))[ss<5,]
        mdl=gam(r~lo(k),dat=dat)
        new=expand.grid(k=seq(min(dat$k),max(dat$k),length.out=nits))
        dat=cbind(r=predict(mdl,newdata=new),new)
        dat=transform(dat,k=exp(k),r=exp(r))
    
        dat}

## ----fig.height=10,fig.width=7-------------------------------------------
cpue=FLQuants(llply(om,stock))
juve=FLQuants(llply(om,function(x) 
              oem(x,sel=1-mat(x),cv=0)))
mat =FLQuants(llply(om,function(x) 
              oem(x,sel=mat(x),  cv=0)))

ggplot()+
  geom_line(aes(year,data),col="red",  
            data=as.data.frame(cpue,drop=TRUE))+
  geom_line(aes(year,data),col="blue", 
            data=as.data.frame(juve,drop=TRUE))+
  geom_line(aes(year,data),col="green",
            data=as.data.frame(mat,drop=TRUE))+
  facet_wrap(~qname,ncol=2)+
  theme_bw()+
  scale_y_continuous(breaks=c(0,2.5e8,5e8))

## ------------------------------------------------------------------------
## fit to all data and logistic production function
mult=FLPar(rep(c(0.5,1,2.0),each=4),
                  dimnames=dimnames(biodyn()@control[,-1]))
mult[c("r","k"),"val"]=1.15
r=c(1,0.6,1,1,1,1)
k=c(1,1,1,1,1,1)

bd=mpb:::biodyns(mlply(1:7,function(i,om,u,mult){
  res=mpb:::FLStock2biodyn(om[[i]],mult=mult)

  setParams( res)     =u[[i]]
  params(    res)["r"]=params( res)["r"]*r[i]
  params(    res)["k"]=params( res)["k"]*k[i]
  setControl(res)     =params(res)
  res                 =fit(res,u[[i]])
  res},om=om,u=cpue,mult=mult))

params(bd[[7]])    =params(bd[[1]])
setControl(bd[[7]])=params(bd[[7]])
bd[[7]]            =fit(bd[[7]],cpue[[7]])

params(bd[[5]])    =params(bd[[1]])
setControl(bd[[5]])=params(bd[[5]])
bd[[5]]            =fit(bd[[5]],cpue[[5]])

p1=plot(bd,fn=list("Stock"=stock))+
  theme_bw()+theme(legend.position="none")

## ------------------------------------------------------------------------
p2=ggplot(ldply(bd,function(x) x@diags))+
  geom_point(aes(year,residual,col=.id))+
  theme_bw()+theme(legend.position="none")

## ------------------------------------------------------------------------
pK=mdply(1:7, function(i) chkFit(i,bd,cpue))

p3=ggplot(subset(pK,ss<6))+
  geom_point(aes(log(k),log(r),col=X1))+
  theme_bw()+theme(legend.position="none")

## ------------------------------------------------------------------------
dat=rbind(cbind(X1=factor(2),with(subset(pK,X1==2),       gamRK(r,k,ss))),
          cbind(X1=factor(3),with(subset(pK,X1==3),       gamRK(r,k,ss))),
          cbind(X1=factor(1),with(subset(pK,!(X1%in%2:3)),gamRK(r,k,ss))),
          cbind(X1=factor(4),with(subset(pK,!(X1%in%2:3)),gamRK(r,k,ss))),
          cbind(X1=factor(5),with(subset(pK,!(X1%in%2:3)),gamRK(r,k,ss))),
          cbind(X1=factor(6),with(subset(pK,!(X1%in%2:3)),gamRK(r,k,ss))),
          cbind(X1=factor(7),with(subset(pK,!(X1%in%2:3)),gamRK(r,k,ss))))

p4=ggplot(dat)+
  geom_point(aes(k,r,col=X1))+
  theme_bw()+theme(legend.position="none")

## ------------------------------------------------------------------------
save(bd,cpue,dat,file="/home/laurie/Desktop/temp/tmp.RData",compress="xz")

## ------------------------------------------------------------------------
# load("/home/laurie/Desktop/temp/tmp.RData")
# 
# b=bd[[1]]
# u=cpue[[1]]
# d=subset(dat,X1==1)

reBd=mpb:::biodyns(mlply(1:7,function(x,b,u,d) 
  refit(b[[x]],u[[x]],subset(d,X1==x),fit=TRUE),b=bd,u=cpue,d=dat))

## ------------------------------------------------------------------------
prf=mdply(1:7, function(x,bd,u) chk(bd[[x]],u[[x]]),bd=reBd,u=cpue)
p5=ggplot(subset(prf,ss<10.0))+
  geom_line(aes(k,ss,col=X1))+
  theme_bw()+theme(legend.position="none")

## ------------------------------------------------------------------------
min=mdply(1:7,function(i) {
  res=model.frame(mpb:::calcObjFn(reBd[[i]],cpue[[i]])[[1]])
  subset(res,min(res$ss)==res$ss)})

## ------------------------------------------------------------------------
bd=mpb:::biodyns(mlply(1:7, function(x) iter(reBd[[x]],as.numeric(min[1,"iter"]))))
p6=plot(bd,fn=list("Stock"=stock))+
  theme_bw()+theme(legend.position="none")

## ------------------------------------------------------------------------
kobe:::multiplot(p1,p3,p5,p2,p4,p6,cols=2)

cf=cbind(Index="Total",ldply(bd,function(x) as.data.frame(stock(x),drop=T)))

## ------------------------------------------------------------------------
## fit to all data and logistic production function
mult=FLPar(rep(c(0.5,1,2.0),each=4),
                  dimnames=dimnames(biodyn()@control[,-1]))
mult[c("r","k"),"val"]=1.15
r=c(1,1,1,1,1,1)
k=c(1,1,1,1,1,1)

bdJ=mpb:::biodyns(mlply(1:7,function(i,b,u){
  
  res           =biodyn(catch=catch(b[[i]]))
  #bug
  params(res)[] =params(b[[i]])[c("r","k","p","b0")]
  res           =fwd(res,catch=catch(res))
  setParams( res)=u[[i]]
  setControl(res)=params(res)
  res=fit(res,u[[i]])
  
  res},b=bd,u=juve))

params(bdJ[[2]])    =params(bdJ[[1]])
setControl(bdJ[[2]])=params(bdJ[[1]])
bdJ[[2]]            =fit(bdJ[[2]],juve[[2]])

params(bdJ[[3]])    =params(bdJ[[1]])
setControl(bdJ[[3]])=params(bdJ[[1]])
bdJ[[3]]            =fit(bdJ[[3]],juve[[3]])

params(bdJ[[6]])    =params(bdJ[[1]])
setControl(bdJ[[6]])=params(bdJ[[1]])
bdJ[[6]]            =fit(bdJ[[6]],juve[[6]])

p1=plot(bdJ,fn=list("Stock"=stock))+
  theme_bw()+theme(legend.position="none")

p2=ggplot(ldply(bdJ,function(x) x@diags))+
  geom_point(aes(year,residual,col=.id))+
  theme_bw()+theme(legend.position="none")

pK=mdply(1:7, function(i) chkFit(i,bdJ,cpue))

p3=ggplot(subset(pK,ss<8&log(k)<22.2))+
  geom_point(aes(log(k),log(r),col=X1))+
  theme_bw()+theme(legend.position="none")
  
dat=ddply(subset(pK,ss<8&log(k)<22.3),.(X1),with,gamRK(r,k,ss))

options(digits=2)
p4=ggplot(dat)+
  geom_point(aes(k,r,col=X1))+
  theme_bw()+theme(legend.position="none")

reBd=mpb:::biodyns(mlply(1:7,function(x,b,u,d) 
  refit(b[[x]],u[[x]],subset(dat,X1==x)),b=bdJ,u=juve,d=dat))

prf=mdply(1:7, function(x,b,u) chk(b[[x]],u[[x]]),b=reBd,u=juve)
p5=ggplot(subset(prf,ss<10.0))+
  geom_line(aes(k,ss,col=X1))+
  theme_bw()+theme(legend.position="none")

min=mdply(1:7,function(i) {
  res=model.frame(mpb:::calcObjFn(reBd[[i]],cpue[[i]])[[1]])
  subset(res,min(res$ss)==res$ss)})

bd=mpb:::biodyns(mlply(1:7, function(x) iter(reBd[[x]],as.numeric(min[1,"iter"]))))
p6=plot(bd,fn=list("Stock"=stock))+
  theme_bw()+theme(legend.position="none")

kobe:::multiplot(p1,p3,p5,p2,p4,p6,cols=2)

cf=rbind(cf,
         cbind(Index="Juvenile",ldply(bd,function(x) as.data.frame(stock(x),drop=T))))        

## ------------------------------------------------------------------------
of=ldply(om,function(x) as.data.frame(stock(x),drop=T))      
of$.id=as.numeric(as.factor(of$.id))

cf=merge(cf,of,by=c("year",".id"))
cf.=ddply(cf,.(Index,.id),transform,om=data.y/mean(data.y),
                                    mp=data.x/mean(data.x))
ggplot(cf.)+geom_line(aes(year,mp,group=.id))+
            geom_line(aes(year,om,group=.id),col="red")+
  facet_wrap(~Index,ncol=1)+
  theme_bw()+theme(legend.position="none")

