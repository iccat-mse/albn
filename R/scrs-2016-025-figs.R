## ---- echo=FALSE---------------------------------------------------------
library(knitr)

opts_chunk$set(comment=NA, fig.width =6, 
                           fig.height=5,
                           fig.path  ="../tex/",
                           warning=FALSE, 
                           message=FALSE, 
                           error  =FALSE, 
                           echo   =FALSE,
                           cache  =TRUE)

iFig=0

## -----init,cache=FALSE---------------------------------------------------
library(plyr)
library(FLCore)
library(kobe)
library(mpb)
library(ggplotFL)

#source('~/Desktop/flr/git/mp/R/biodyn-hcr.R')

## ------------------------------------------------------------------------
reflevel<-rbind(data.frame(x=c(-Inf,-Inf,Inf,Inf), y=c(-Inf,  1,  1,-Inf), 
                         fill=as.factor("bottom")),
                data.frame(x=c(-Inf,-Inf,Inf,Inf), y=c(   1,Inf,Inf,   1), 
                         fill=as.factor("top")))

## ----prj-----------------------------------------------------------------
bd  =sim()
bd30=window(sim(),end=30)
bds =mpb:::biodyns()

bds["stock"]       =bd30
bds["hcr"]         =fwd(bd30,harvest=hcr(bd30,yr=29,pyr=31:47))
bds["bound F 5%"]  =fwd(bd30,harvest=hcr(bd30,yr=29,pyr=31:47, bndF  =c(0.90,1.1)))
bds["bound TAC 5%"]=fwd(bd30,catch  =hcr(bd30,yr=29,pyr=31:47, tac=T,bndTac=c(0.90,1.1)))

names(bds)=c("Stock","HCR","Bound F","Bound TAC")

plot(bds[1:2])+
  theme_bw()+
  theme(legend.position="bottom")  

## ------------------------------------------------------------------------
dat=ldply(bds,function(x){
  model.frame(mcf(
      FLQuants(stock  =stock(  x)%/%refpts(x)["bmsy"],
               harvest=harvest(x)%/%refpts(x)["fmsy"],
               catch  =catch(  x)%/%refpts(x)[ "msy"])),
  drop=TRUE)
  })

## ----fig.height=4--------------------------------------------------------
ggplot(dat)+
  geom_polygon(data=reflevel,aes(x,y,fill=fill)) +
  scale_fill_manual(values=c("pink","lightgreen"), guide="none") +  
  geom_line(aes(year,stock,col=.id),size=1.5)+
  theme(legend.position="bottom")+
  geom_hline(aes(yintercept=0.80),col="grey")+
  geom_hline(aes(yintercept=0.40),col="red")+
  xlab("Year")+ylab(expression(B/B[MSY]))

## ----fig.height=4--------------------------------------------------------
ggplot(dat)+
  geom_polygon(data=reflevel,aes(x,y,fill=fill)) +
  scale_fill_manual(values=c("lightgreen","pink"), guide="none") +  
  geom_line(aes(year,harvest,col=.id),size=1.5)+
  theme(legend.position="bottom")+
  geom_hline(aes(yintercept=0.70),col="grey")+
  geom_hline(aes(yintercept=0.10),col="red")+
  xlab("Year")+ylab(expression(F/F[MSY]))

## ----fig.height=4--------------------------------------------------------
ggplot(dat)+
  geom_polygon(data=reflevel,aes(x,y,fill=fill)) +
  scale_fill_manual(values=c("pink","lightgreen"), guide="none") +  
  geom_line(aes(year,catch,col=.id),size=1.5)+
  theme(legend.position="bottom")+  
  xlab("Year")+ylab(expression(Yield/MSY))

## ----fig.height=6,fig.width=6--------------------------------------------
hcr=data.frame(stock=c(0,.4,.8,.8,2),harvest=c(.1,.1,.7,.7,.7))
kobePhase(dat)+
  geom_hline(aes(yintercept=0.70),col="grey",size=.25)+
  geom_hline(aes(yintercept=0.10),col="red",size=.25)+
  geom_vline(aes(xintercept=0.80),col="grey",size=.25)+
  geom_vline(aes(xintercept=0.40),col="red",size=.25)+
  geom_line(aes(stock,harvest),data=hcr,col="brown",size=1.2)+
  geom_path(aes(stock,harvest,col=.id),size=1.2)+
  theme(legend.position="bottom")

## ------------------------------------------------------------------------
bd  =sim()
bd30=window(sim(),end=30)
bds =mpb:::biodyns()

bds=mpb:::biodyns(list(stock=bd30,hcr=bd30,"bound F 5%"=bd30,"bound TAC 5%"=bd30))

for (i in seq(29,47,3)){
  bds[["hcr"]]         =fwd(bds[["hcr"]],         
                            harvest=mpb:::hcr(bds[["hcr"]],yr=i,pyr=i+2:4))
  bds[["bound F 5%"]]  =fwd(bds[["bound F 5%"]],  
                            harvest=mpb:::hcr(bds[["hcr"]],yr=i,pyr=i+2:4,       bndF  =c(0.90,1.1)))
  bds[["bound TAC 5%"]]=fwd(bds[["bound TAC 5%"]],
                            catch  =mpb:::hcr(bds[["hcr"]],yr=i,pyr=i+2:4, tac=T,bndTac=c(0.90,1.1)))}

plot(bds)+
  theme_bw()+
  theme(legend.position="bottom")  

## ------------------------------------------------------------------------
dat=ldply(bds,function(x){
  model.frame(mcf(
      FLQuants(stock  =stock(  x)%/%refpts(x)["bmsy"],
               harvest=harvest(x)%/%refpts(x)["fmsy"],
               catch  =catch(  x)%/%refpts(x)[ "msy"])),
  drop=TRUE)
  })

## ----fig.height=4--------------------------------------------------------
ggplot(dat)+
  geom_polygon(data=reflevel,aes(x,y,fill=fill)) +
  scale_fill_manual(values=c("pink","lightgreen"), guide="none") +  
  geom_line(aes(year,stock,col=.id),size=1.5)+
  theme(legend.position="bottom")+
  geom_hline(aes(yintercept=0.80),col="grey")+
  geom_hline(aes(yintercept=0.40),col="red")+
  xlab("Year")+ylab(expression(B/B[MSY]))

## ----fig.height=4--------------------------------------------------------
ggplot(dat)+
  geom_polygon(data=reflevel,aes(x,y,fill=fill)) +
  scale_fill_manual(values=c("lightgreen","pink"), guide="none") +  
  geom_line(aes(year,harvest,col=.id),size=1.5)+
  theme(legend.position="bottom")+
  geom_hline(aes(yintercept=0.70),col="grey")+
  geom_hline(aes(yintercept=0.10),col="red")+
  xlab("Year")+ylab(expression(F/F[MSY]))

## ----fig.height=4--------------------------------------------------------
ggplot(dat)+
  geom_polygon(data=reflevel,aes(x,y,fill=fill)) +
  scale_fill_manual(values=c("pink","lightgreen"), guide="none") +  
  geom_line(aes(year,catch,col=.id),size=1.5)+
  theme(legend.position="bottom")+  
  xlab("Year")+ylab(expression(Yield/MSY))

## ----fig.height=6,fig.width=6--------------------------------------------
hcr=data.frame(stock=c(0,.4,.8,.8,2),harvest=c(.1,.1,.7,.7,.7))
kobePhase(dat)+
  geom_hline(aes(yintercept=0.70),col="grey",size=.25)+
  geom_hline(aes(yintercept=0.10),col="red",size=.25)+
  geom_vline(aes(xintercept=0.80),col="grey",size=.25)+
  geom_vline(aes(xintercept=0.40),col="red",size=.25)+
  geom_line(aes(stock,harvest),data=hcr,col="brown",size=1.2)+
  geom_path(aes(stock,harvest,col=.id),size=1.2)+
  theme(legend.position="bottom")

