## ---- echo=FALSE---------------------------------------------------------
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

## ---- echo=FALSE---------------------------------------------------------
library(ggplot2)
library(scales)
library(plyr)
library(reshape)
library(FLCore)
library(FLash)
library(mpb)
library(FLBRP)
library(ggplotFL)
library(kobe)
library(FLife)

dirMy ="/home/laurie/MEGAsync/mse/albn"
dirDat=file.path(dirMy,"data")

theme_set(theme_bw(12))

load(paste(dirDat,"om.RData",sep="/"))

i=0

## ------------------------------------------------------------------------
#source("/home/laurie/Desktop/flr/git/biodyn/R/biodyn-production.R")

rodFn=FLife:::rodFn

hrate=function(x) catch(x)/stock(x)

FLBRP2biodyn=function(from){
  
  ## age based regference points
  msy =c(from@refpts["msy",   "yield"])
  bmsy=c(from@refpts["msy",   "biomass"])
  k   =c(from@refpts["virgin","biomass"])
  
  # biomass based reference points
  p   =optimise(function(p,bmsy,k) 
    (bmsy-k*(1/(1+p))^(1/p))^2, c(0.001,5), bmsy=bmsy,k=k)$minimum
  k=bmsy/((1/(1+p))^(1/p))
  r=msy/(k*(1/(1+p))^(1/p+1))
  b0=mean(stock(from)[,1:5]/k)
  
  bd=biodyn()
  
  bd@params=FLPar(c(r=r,k=k,p=p,b0=b0))
  
  bd@catch=catch.obs(from)
  bd@stock=from@stock.obs
  
  range(bd)["minyear"]=dims(bd@catch)$minyear
  range(bd)["maxyear"]=dims(bd@catch)$maxyear
  
  bd}

productionFn=function(object,stk="missing",
                      slots=c("landings.sel","discards.sel",
                              "stock.wt","landings.wt","discards.wt",
                              "m","mat",
                              "harvest.spwn","m.spwn")){
  
  nyr=dim(biomass.obs(object))[2]
  
  prd=biomass.obs(object)[,-nyr]+
    catch.obs(  object)[,-nyr]-
    biomass.obs(object)[,-1]
  
  ref               =FLBRP:::refpts(object)[1,]
  dimnames(ref)[[1]]="biomass"
  ref               =propagate(ref,nyr-1)
  ref[]             =NA
  ref[,"biomass"]   =biomass.obs(object)[,-nyr]
  refpts(object)    =ref
  
  if (!missing(stk)){
    for (i in slots[!(slots%in%c("landings.sel","discards.sel"))]){
      slot(object,i)=propagate(slot(object,i),dims(ref)$iter)
      slot(object,i)[]=slot(stk,i)[,-nyr]
    }
    
    if (any(c("landings.sel","discards.sel")%in%slots)){
      sel=harvest(stk)[,-nyr]
      sel=sel%/%apply(sel[ac(range(stk)["minfbar"]:range(stk)["maxfbar"])],2,mean)
      
      if ("landings.sel"%in%slots){
        landings.sel(object)  =propagate(landings.sel(object),dims(ref)$iter)
        landings.sel(object)[]=sel*landings.n(stk)[,-nyr]/catch.n(stk)[,-nyr]}
      
      if ("discards.sel"%in%slots){
        discards.sel(object)  =propagate(discards.sel(object),dims(ref)$iter)
        discards.sel(object)[]=sel*discards.n(stk)[,-nyr]/catch.n(stk)[,-nyr]}
    }
  }
  
  refpts(object)=computeRefpts(object)
  
  res=cbind(as.data.frame(biomass.obs(object)[,-nyr],drop=T)[,c("year","data")],
            as.data.frame(ssb.obs(    object)[,-nyr],drop=T)[,c("data")],
            as.data.frame(catch.obs(  object)[,-nyr],drop=T)[,"data"],
            as.data.frame(rec.obs(    object)[,-nyr],drop=T)[,"data"],
            as.data.frame(biomass.obs(object)[,-nyr]-
                          ssb.obs(    object)[,-nyr],drop=T)[,"data"],
            as.data.frame(prd,drop=T)[,"data"],
            model.frame(FLBRP:::refpts(object)[,"yield"])[,"biomass"])
  
  names(res)=c("year","biomass","ssb","catch","rec","juve","obs","hat")
  
  res[,c("year","biomass","ssb","juve","rec","catch","obs","hat")]}


setGeneric('production',   function(object,stk,...) standardGeneric('production'))

setMethod('production', signature(object='FLBRP',stk='missing'),
          function(object,stk) productionFn(object,stk))

setMethod('production', signature(object='FLBRP',stk='FLStock'),
         function(object,stk,slots=c("landings.sel","discards.sel",
                                     "stock.wt","landings.wt","discards.wt",
                                     "m","mat",
                                     "harvest.spwn","m.spwn"))
         productionFn(object,stk,slots))

## ----albn-om,fig.width=6,fig.height=8------------------------------------
ggplot(ldply(om, function(x) as.data.frame(FLQuants(x,
                              Recruits=rec,
                              SSB     =ssb,
                              Biomass =stock,
                              F       =fbar,
                              Fapex   =fapex,
                              harvest =hrate,
                              Catch   =catch))))+
         geom_line(aes(year,data,col=.id))+
         facet_grid(qname~.,scale="free_y")+
         theme_bw(10)

## ----albn-om-sr,fig.height=6,fig.width=6---------------------------------
srs=FLSRs(llply(om,as.FLSR,model="bevholt"))
srs=FLSRs(llply(srs, fmle,control=list(silent=TRUE)))
names(srs)=names(om)

plot(srs[["Base"]])

## ----albn-om-sr-fits,fig.height=6,fig.width=6----------------------------
sr1=ldply(srs, function(x) 
  model.frame(FLQuants(x,"fitted","ssb","rec","residuals"),drop=TRUE))

ggplot(sr1)+
  geom_line(aes(ssb,fitted),col="red")+
  geom_point(aes(ssb,rec))+
  facet_wrap(~.id,ncol=3)+
  xlab("SSB")+ylab("Recruits")

## ----albn-om-sr-rsdl5----------------------------------------------------
ggplot(sr1,aes(ssb,residuals))+
  geom_hline(aes(yintercept=0),col="red")+
  geom_point()+
  geom_smooth()+
  facet_wrap(~.id,ncol=3)+
  xlab("SSB")+ylab("Residuals")

## ----albn-om-sr-rsdl2----------------------------------------------------
ggplot(sr1,aes(year,residuals))+
  geom_hline(aes(yintercept=mean(residuals)),col="red")+
  geom_point()+
  geom_smooth()+
  facet_wrap(~.id,ncol=3)+
  xlab("Year")+ylab("Residuals")

## ----albn-om-sr-rsdl3----------------------------------------------------
ggplot(sr1,aes(fitted,residuals))+
  geom_hline(aes(yintercept=mean(residuals)),col="red")+
  geom_point()+
  geom_smooth()+
  facet_wrap(~.id,ncol=3)+
  xlab("Fitted")+ylab("Residuals")

## ----albn-om-sr-rsdl4,fig.width=6,fig.height=8---------------------------
qq=ddply(sr1,.(.id),with,as.data.frame(qqnorm(residuals,plot=F)))

ggplot(qq,aes(x,y))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~.id,ncol=3)+
  xlab("Theoretical")+ylab("Sample")

## ----albn-om-sr-ar-------------------------------------------------------
acf=ddply(sr1,.(.id),with,data.frame(lag=rec[-1],rec=rec[-length(rec)]))

ggplot(acf,aes(rec,lag))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  facet_wrap(~.id,ncol=3)+
  xlab(expression(Rec[t]))+ylab(expression(Rec[t+1]))

## ----albn-om-biol,fig.width=6,fig.height=5-------------------------------
names(eql)=names(om)
biol=ldply(eql, function(x) as.data.frame(x[["catch.sel","m","mat","stock.wt"]]))
biol=transform(biol,
               Quantity=factor(qname,labels=c("Selectivity","M","Maturity","Mass")))
ggplot(biol)+
  geom_line(aes(age,data,group=.id,col=.id))+ 
  facet_wrap(~Quantity,scale="free_y")+ 
  theme_bw()+
  scale_colour_discrete(name="Scenario")+
  xlab("Age") +ylab("")

## ----albn-om-eql,fig.width=6,fig.height=5--------------------------------
eql=FLBRPs(mlply(names(srs),function(i) 
  brp(FLBRP(om[[i]],sr=srs[[i]],nyears=dims(om[[i]])$year))))
names(eql)=names(om)
dat=subset(plot(eql,refpts=FALSE)$dat,pnl=="Equilibrium Yield v. SSB")

ggplot(data=dat,aes(x,y,col=.id))+
  geom_line()+
  theme(legend.position="bottom")+
  
  theme_bw(10)

## ----albn-om-prd,fig.height=4--------------------------------------------
bd=FLBRP2biodyn(eql[["Base"]])
#bd=as(eql[["Base"]],"biodyn")
bd=mpb:::fwd(bd,catch=catch(bd)[,-c(1,82)],maxF=0.3)

cpue =stock(om[["Base"]])
catch=catch(om[["Base"]])

fbar(eql[["Base"]])=fbar(eql[["Base"]])/max(fbar(eql[["Base"]]))*.5
eql[["Base"]]=brp(eql[["Base"]])

dat=model.frame(FLQuants(bd,"stock","catch"))

mpb:::plotPrdfn(bd,biomass=FLQuant(seq(0,max(params(bd)['k']),length.out=101)))+
  geom_path( aes(stock,catch),data=dat)+
  geom_path(aes(stock,catch),
            data=model.frame(FLQuants(bd,"stock","catch")))+
  geom_path(aes(stock,catch),
             data=subset(model.frame(FLQuants(eql[["Base"]],"stock","catch")),
                         stock>0),col="blue")+
  theme_bw(10)

## ----albn-om-mp----------------------------------------------------------
dat=
  rbind.fill(cbind(What="MP", as.data.frame(FLQuants(bd,
                      "catch",
                      "harvest"=function(x) catch(x)/stock(x),"stock"))),
             cbind(What="OM", as.data.frame(FLQuants(om[["Base"]],
                        "catch",
                        "harvest"=function(x) catch(x)/stock(x),"stock"))))


ggplot(dat)+
  geom_line(aes(year,data,col=What))+
  facet_grid(qname~.,scale="free")+
  theme_bw()

## ------------------------------------------------------------------------
fapexAge<-function(object){
  tmp=harvest(object)
  tmp[]=fapex(object)
  tmp=FLQuant(ages(tmp)[harvest(object)==tmp],
              dimnames=dimnames(fapex(object)))
  tmp}

bry=eql[["Base"]]
landings.sel(bry)=propagate(iter(landings.sel(bry),1),82)
landings.sel(bry)[]=harvest(om[["Base"]])

tmp=model.frame(computeRefpts(bry)["msy",   c(1:2,4:5)])
tm2=model.frame(computeRefpts(bry)["msy","ssb"]/computeRefpts(bry)["virgin","ssb"])
tm2[,"quantity"]="lro"
tmp=rbind.fill(tmp,tm2,
          transform(as.data.frame(fapexAge(om[["Base"]]),drop=T),
          msy=data,iter=year-1929,quantity="Fapex Age")[,c(3:5)])
tmp=rbind.fill(tmp,
               transform(as.data.frame(computeRefpts(bry)["msy","biomass"]/
                                       computeRefpts(bry)["virgin","biomass"]),
                         msy=data,quantity="shape")[,c(2:3,5)])
tmp=rbind.fill(tmp,
               transform(as.data.frame(computeRefpts(bry)["msy","yield"]/
                                       computeRefpts(bry)["msy","biomass"]),
                         msy=data,quantity="hrate")[,c(2:3,5)])

bdt=maply(1:82,function(i) FLBRP2biodyn(iter(bry,i))@params[c("r","k","p")])

bdt=melt(bdt)
names(bdt)=c("iter","quantity","msy")
bdt=rbind(tmp,bdt)
names(bdt)[1]="data"
#unique(bdt$quantity)

bdt=transform(bdt,quantity=factor(quantity,
            levels=c("yield","harvest","hrate","Fapex Age",
                         "ssb","biomass","lro",
                         "r","k","p","shape"),
            labels=c("MSY","Fmsy","Hmsy","FApex Age",
                     "SSBmsy","Bmsy","Lmsy",
                     "r","k","p","Shape")))

## ------------------------------------------------------------------------
prd=ldply(eql,function(OM) mpb:::production(object=OM))
stR=ddply(prd,.(.id),with,rod(FLQuant(rec)))

ggplot(stR)+
  geom_line(aes(year,rec),data=prd,col="red")+
  geom_polygon(aes(year+1929,data,group=regime),fill="blue",alpha=.4)+
  facet_wrap(~.id,ncol=2)+
  theme_bw()

## ------------------------------------------------------------------------
rgm=ldply(eql,function(OM) {
              res=mpb:::production(object=OM)[,c("year","obs","hat")]
              rod(as.FLQuant(transform(res,data=obs/hat)[,c("year","data")]))})

prd=ldply(eql,function(OM) {
              res=mpb:::production(object=OM)[,c("year","obs","hat")]
              transform(res,data=obs/hat)[,c("year","data")]})
              
ggplot(rgm)+
  geom_line(aes(year,data),data=prd,col="red")+
  geom_polygon(aes(year,data,group=regime),fill="blue",alpha=.4)+
  facet_wrap(~.id,ncol=3)+
  theme_bw()

## ------------------------------------------------------------------------
bench=mdply(names(eql), function(OM){
      
    bry=eql[[OM]]
    landings.sel(bry)=propagate(iter(landings.sel(bry),1),82)
    landings.sel(bry)[]=harvest(om[["Base"]])
    
    
    tmp=model.frame(computeRefpts(bry)["msy",   c(1:2,4:5)])
    tm2=model.frame(computeRefpts(bry)["msy","ssb"]/computeRefpts(bry)["virgin","ssb"])
    tm2[,"quantity"]="lro"
    tmp=rbind.fill(tmp,tm2,
              transform(as.data.frame(fapexAge(om[["Base"]]),drop=T),
              msy=data,iter=year-1929,quantity="Fapex Age")[,c(3:5)])
    tmp=rbind.fill(tmp,
                   transform(as.data.frame(computeRefpts(bry)["msy","biomass"]/
                                           computeRefpts(bry)["virgin","biomass"]),
                             msy=data,quantity="shape")[,c(2:3,5)])
    tmp=rbind.fill(tmp,
                   transform(as.data.frame(computeRefpts(bry)["msy","yield"]/
                                           computeRefpts(bry)["msy","biomass"]),
                             msy=data,quantity="hrate")[,c(2:3,5)])
    
    bdt=maply(1:82,function(i) FLBRP2biodyn(iter(bry,i))@params[c("r","k","p")])
    
    bdt=melt(bdt)
    names(bdt)=c("iter","quantity","msy")
    bdt=rbind(tmp,bdt)
    names(bdt)[1]="data"
    
    bdt=transform(bdt,quantity=factor(quantity,
                levels=c("yield","harvest","hrate","Fapex Age",
                             "ssb","biomass","lro",
                             "r","k","p","shape"),
                labels=c("MSY","Fmsy","Hmsy","FApex Age",
                         "SSBmsy","Bmsy","Lmsy",
                         "r","k","p","Shape")))
    
    bdt})

## ----albn-om-bench3,fig.width=7,fig.height=8-----------------------------
ggplot(bench)+
  geom_line(aes(as.numeric(as.numeric(iter)+1929),data,col=X1))+
  facet_wrap(~quantity,scale="free",ncol=2)+
  xlab("Year")+ylab("")+
  theme_bw(10)

## ----albn-om-bench4,fig.width=7,fig.height=8-----------------------------
ggplot(ddply(subset(bench,as.numeric(iter)>30&quantity!="p"),
        .(quantity,X1), transform, data=(data-mean(data))/mean(data)))+
  geom_line(aes(as.numeric(as.numeric(iter)+1929),data,col=X1))+
  facet_wrap(~quantity,scale="free",ncol=2)+
  xlab("Year")+ylab("")+
  theme_bw(10)+
  scale_y_continuous(labels=percent)+
  geom_hline(aes(yintercept=0),col="grey")

## ----albn-sp,fig.width=6,fig.height=8------------------------------------
prd1=ldply(eql,mpb:::production)
prd2=mdply(names(eql),function(OM) production(object=eql[[OM]],om[[OM]]))
names(prd1)[1]="OM"
prd2=transform(prd2,OM=names(eql)[X1])[,-1]

ggplot(prd1)+
  geom_point(aes(year,obs))+
  geom_line( aes(year,hat),col="blue")+
  geom_line( aes(year,hat),data=prd2,col="red")+
  theme_bw()+facet_wrap(~OM,ncol=3)   

## ----albn-om-prd2,fig.width=6,fig.height=6-------------------------------
ggplot(prd1)+
  geom_histogram(aes(obs/biomass))+
  facet_wrap(~OM)+
  xlab("Surplus Production")+
  theme_bw()

## ----albn-om-prd3,fig.width=6,fig.height=6-------------------------------
ggplot(prd1)+
  geom_histogram(aes((obs-hat)/biomass))+
  geom_vline(aes(xintercept=0),col="red")+
  facet_wrap(~OM)+
  xlab("Surplus Production")+
  theme_bw()

## ----albn-om-ccf,fig.width=8,fig.height=6,eval=FALSE---------------------
## par(mfrow=c(2,2),mar=c(2,2,4,2))
## with(dat,acf(rec))
## with(dat,acf(obs))
## with(dat,ccf(rec,obs))

## ----fig.width=8,fig.height=4,eval=FALSE---------------------------------
## ggplot(spectra(dat$Recruitment))+
##   geom_point(aes(f,mx))+
##   theme_bw()
## 
## ggplot(spectra(dat$Production))+
##   geom_point(aes(f,mx)sa)+
##   theme_bw()
## 
## **Figure** Spectra for production

## ------------------------------------------------------------------------
refs=transform(ldply(eql, function(x) FLBRP:::refpts(x)["msy",1:5,drop=T]))

ts=ldply(om, function(x) model.frame(FLQuants(x,"stock","ssb","catch","fbar","rec"),drop=T))
ts=merge(ts,refs,by=".id")
ts=transform(ts,stock  = stock/biomass,
                ssb    = ssb.x/ssb.y,
                catch  = catch/yield,
                harvest=(catch/stock)/(yield/biomass),
                f      =harvest/harvest,
                rec    =rec.x/rec.y)[,c(".id","year","stock","ssb","catch","harvest","f","rec")]

## ----albn-om-harvest,fig.width=6-----------------------------------------
quads<- rbind(data.frame(x=c(-Inf,-Inf,Inf,Inf), y=c(-Inf,  1,  1,-Inf), 
                         fill=as.factor("green")),
              data.frame(x=c(-Inf,-Inf,Inf,Inf), y=c(   1,Inf,Inf,   1), 
                         fill=as.factor("red")))

ggplot(ts)+
  geom_polygon(data=quads,aes(x,y,fill=fill)) +
  geom_line(aes(year,harvest,group=.id))+
  facet_wrap(~.id,ncol=3)+
  scale_fill_manual(values=c("lightgreen","pink"), guide="none") +  
  ylab(expression(Harvest/H[MSY]))+
  xlab("Year")

## ----albn-om-biomass,fig.width=6-----------------------------------------
ggplot(ts)+
  geom_polygon(data=quads,aes(x,y,fill=fill)) +
  geom_line(aes(year,stock,group=.id))+
  facet_wrap(~.id,ncol=3)+
  scale_fill_manual(values = c("pink","lightgreen"), guide="none") +  
  ylab(expression(Biomass/B[MSY]))+
  xlab("Year")

## ----albn-om-k2pp--------------------------------------------------------
kobePhase(ts)+
  geom_path(aes(stock,harvest,group=.id))+
  geom_point(aes(stock,harvest),col="blue",size=3,data=subset(ts,year==2011))+
  facet_wrap(~.id,ncol=3)+
  ylab(expression(Harvest/H[MSY]))+
  xlab(expression(Biomass/B[MSY]))

## ----eval=FALSE----------------------------------------------------------
## d_ply(prd,.(.id),with,spectrum(juve))

## ----fig101--------------------------------------------------------------
refpts=FLBRP:::refpts

stks=FLStocks(llply(eql,function(eq){
  eq=eql[[1]] 
  fbar(eq) =propagate(FLQuant(rep(1,1001)),3)%*%refpts(eq)["msy","harvest"]
  f        =fbar(eq)[,-1]
  f[,,,,,1]=0.1
  f[,,,,,3]=f[,,,,,2]*2
  
  stk=as(eq,"FLStock")
  stk=FLash:::fwd(stk,f=f,sr=eq,sr.residuals=rlnorm(1,rec(iter(stk,1))*0,.3))[,-(1:200)]
  stk}))

## ----fig111--------------------------------------------------------------
dat=ldply(stks,function(x) mdply(1:3,function(F)
          as.data.frame(spectrum(iter(rec(x),F), log = "dB", ci = 0.8,plot=FALSE)[c("freq","spec")])))

dat=ddply(subset(dat,freq>0.05),.(.id,X1),transform,val=spec/max(spec))
ggplot(dat,aes(freq,val,col=X1))+
  geom_smooth(se=FALSE)+
  facet_wrap(~.id,ncol=3)+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Frequency")+ylab("")

## ----fig121--------------------------------------------------------------
dat=ldply(stks,function(x) mdply(1:3,function(F)
          as.data.frame(spectrum(iter(ssb(x),F), log = "dB", ci = 0.8,plot=FALSE)[c("freq","spec")])))

dat=ddply(subset(dat,freq>0.05),.(.id,X1),transform,val=spec/max(spec))
ggplot(dat,aes(freq,val,col=X1))+
  geom_smooth(se=FALSE)+
  facet_wrap(~.id,ncol=3)+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Frequency")+ylab("")

## ----fig131--------------------------------------------------------------
dat=ldply(stks,function(x) mdply(1:3,function(F)
          as.data.frame(spectrum(iter(catch(x),F), log = "dB", ci = 0.8,plot=FALSE)[c("freq","spec")])))

dat=ddply(subset(dat,freq>0.05),.(.id,X1),transform,val=spec/max(spec))
ggplot(dat,aes(freq,val,col=X1))+
  geom_smooth(se=FALSE)+
  facet_wrap(~.id,ncol=3)+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Frequency")+ylab("")

## ----fig141--------------------------------------------------------------
dat=ldply(stks,function(x) mdply(1:3,function(F)
          as.data.frame(spectrum(iter(stock(x),F), log = "dB", ci = 0.8,plot=FALSE)[c("freq","spec")])))

dat=ddply(subset(dat,freq>0.05),.(.id,X1),transform,val=spec/max(spec))
ggplot(dat,aes(freq,val,col=X1))+
  geom_smooth(se=FALSE)+
  facet_wrap(~.id,ncol=3)+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Frequency")+ylab("")

## ----fig151--------------------------------------------------------------
dat=ldply(stks,function(x) mdply(1:3,function(F)
          as.data.frame(spectrum(iter(stock(x)-ssb(x),F), log = "dB", ci = 0.8,plot=FALSE)[c("freq","spec")])))

dat=ddply(subset(dat,freq>0.05),.(.id,X1),transform,val=spec/max(spec))
ggplot(dat,aes(freq,val,col=X1))+
  geom_smooth(se=FALSE)+
  facet_wrap(~.id,ncol=3)+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Frequency")+ylab("")

## ----fig161--------------------------------------------------------------
source('~/Desktop/flr/git/mp/R/FLBRP-production.R')

dat=ldply(eql,function(x) mdply(1:3,function(F)
          as.data.frame(spectrum(iter(mpb:::production(x),F), log = "dB", ci = 0.8,plot=FALSE)[c("freq","spec")])))

dat=ddply(subset(dat,freq>0.05),.(.id,X1),transform,val=spec/max(spec))
ggplot(dat,aes(freq,val,col=X1))+
  geom_smooth(se=FALSE)+
  facet_wrap(~.id,ncol=3)+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Frequency")+ylab("")

