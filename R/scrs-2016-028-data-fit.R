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

##Please change the dirs to those of your own choice
dirMy   ="/home/laurie/Desktop/scrs-2016/papers/scrs-2016-028"
dirDat  =file.path(dirMy,"data")
dirKobe ="/home/laurie/MEGAsync/mse/albn/inputs/aspic"

scen=expand.grid(stk="albn",method="aspic",  run=paste("run",c(1:7),   sep=""),
                              stringsAsFactors=F)
scen=mdply(scen,function(stk,method,run) 
  data.frame(dir=file.path(dirKobe,run)))

## ----inputs--------------------------------------------------------------
asp=aspics(file.path(scen$dir,"aspic.inp"))
u  =ldply(asp, function(x) index(x))
u  =u[!duplicated(u[,c("name","year")]),c("name","year","index")]

## Prettify names
nms=c("Chinese-Tai\n All","Chinese-Tai\n Early","Chinese-Tai\n Late",
                          "Japanese LL\n Early","Japanese LL\n Late", 
      "Troll Combined")
names(nms)=sort(unique(u$name))
nm2=c("Chinese-Tai Early","Chinese-Tai Late","Chinese-Tai All",
      "Japanese LL Early","Japanese LL Late", 
      "Troll Combined")
names(nm2)=sort(unique(u$name))

u  =transform(u,name=factor(name,levels=sort(unique(u$name)),labels=nm2))

cpue=FLQuants(dlply(u,.(name),with, 
                 as.FLQuant(data.frame(year=year,data=index))))
names(cpue)=nm2

## ----cpue,fig.height=7.5-------------------------------------------------
library(gam)

mpb:::plotIndex(asp[1:7])

## ----cpueCor-------------------------------------------------------------
cr=cor(cast(u,year~name,value="index")[,-1],
       use="pairwise.complete.obs")

dimnames(cr)=list(names(cpue), names(cpue))
cr[is.na(cr)]=0
corrplot(cr,diag=F,order="hclust",addrect=2)  +          
             theme(legend.position="bottom")  

## ----cpueCcf-------------------------------------------------------------
names(cpue)=nms
cc=mdply(expand.grid(a=names(cpue),b=names(cpue)),
         function(a,b){
           #print(paste(a,b))
           res=model.frame(mcf(FLQuants(cpue[c(a,b)])))
           res=subset(res,!is.na(res[,7])&!is.na(res[,8]))

           if (dim(res)[1]>10){
             res=data.frame(lag=-10:10,data=ccf(res[,7],res[,8],plot=F,
                                                lag.max=10)$acf)
             return(res)}else{return(NULL)}}
           )

# cc=transform(subset(cc,a%in%c("us","jll1","jll2")&b%in%c("us","jll1","jll2")),
#                     a=factor(a,levels=rev(names(cpue))))

ggplot(cc)+
  geom_linerange(aes(x=lag,ymin=0,ymax=data))+
  facet_grid(a~b)+
  geom_vline(aes(xintercept=0))+
  theme_bw()    

## ----bdFit---------------------------------------------------------------
bd=biodyn("pellat",
          params=FLPar(r=0.3,k=5.5e5,b0=.95,p=0.001),
          catch =catch(asp[[1]]))

bds=mpb:::biodyns(list("Chinese-Tai"=bd,
                      "Japanese LL"=bd,
                      "Troll"      =bd))

params(bds[[1]])=FLPar(r=.1,k=3.2e6,b0=0.95,p=0.001)
setParams(bds[[1]])=cpue[2:3]
setControl(bds[[1]])=params(bds[[1]])
control(bds[[1]])["p","phase"]=-1
bds[[1]]=fit(bds[[1]],cpue[2:3])

params(bds[[2]])=FLPar(r=0.1,k=0.95e6,b0=0.95,p=0.001)
setParams(bds[[2]])=cpue[4:5]
setControl(bds[[2]])=params(bds[[2]])
control(bds[[2]])["p","phase"]=-1
bds[[2]]=fit(bds[[2]],cpue[4:5])

params(bds[[3]])=FLPar(r=.3,k=0.7e6,b0=0.95,p=0.001)
setParams(bds[[3]])=cpue[6]
setControl(bds[[3]])=params(bds[[3]])
control(bds[[3]])["p","phase"]=-1
bds[[3]]=fit(bds[[3]],cpue[6])

plot(bds)+
  theme_bw()    

## ----bdResidual----------------------------------------------------------
ggplot(ldply(bds, function(x) x@diags))+
    geom_hline(aes(yintercept=0),col="orange")+
    geom_point(aes(year,residual),fill="salmon",col="grey",pch=21)+
    facet_wrap(~.id,ncol=1)+
    theme_bw()    

## ----aspFit--------------------------------------------------------------
names(cpue)=nm2
asp=aspics(list("Chinese-Tai"=mpb:::biodynCpue2aspic(bds[[1]],cpue[2:3]),
                "Japanese LL"=mpb:::biodynCpue2aspic(bds[[2]],cpue[4:5]),
                "Troll"      =mpb:::biodynCpue2aspic(bds[[3]],cpue[6])))

asp[[1]]@model=factor("FOX")
asp[[2]]@model=factor("FOX")
asp[[3]]@model=factor("FOX")

asp[[1]]=fit(asp[[1]])
asp[[2]]=fit(asp[[2]])
asp[[3]]=fit(asp[[3]])

plot(asp)+
  theme_bw()

## ----aspResidual---------------------------------------------------------
ggplot(ldply(asp, function(x) x@diags))+
    geom_hline(aes(yintercept=0),col="orange")+
    geom_point(aes(year,residual),fill="salmon",col="grey",pch=21)+
    facet_wrap(~.id,ncol=1)+
    theme_bw()    

## ------------------------------------------------------------------------
save(bds,cpue,file=file.path(dirDat,"fit.RData"))

