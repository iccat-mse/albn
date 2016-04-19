## Albacore MSE"

harvest=function(x) catch(x)/stock(x)[,dimnames(catch)$year]

rFn<-function(x){
  f   =c(FLBRP:::refpts(x)["crash","harvest"])
  if (is.na(f)) 
      f=c(refpts(x)["msy","harvest"])*3
                                  
  r=try(log(lambda(leslie(x,f))))
  if (class(r)=="try-error") 
      r=NA  
  r}

priorFn<-function(x,p=NULL){
  
  # production function based ob S/R
  msy =c(FLBRP:::refpts(x)["msy",   "yield"])
  bmsy=c(FLBRP:::refpts(x)["msy",   "biomass"])
  
  r=rFn(x)
  k=c(FLBRP:::computeRefpts(x)["virgin","biomass"])
  
  if (is.null(p))
  p   =optimise(function(p,bmsy,k) 
                (bmsy-k*(1/(1+p))^(1/p))^2, c(0.001,5), 
                bmsy=bmsy,k=k)$minimum

  kPrime=bmsy/((1/(1+p))^(1/p))
  rPrime=msy/(k*(1/(1+p))^(1/p+1))
  
  cbind(p=p,
        r     =r,     k     =k,
        rPrime=rPrime,kPrime=kPrime,
        msy=msy,bmsy=bmsy,fmsy=msy/bmsy)}

fitBd<-function(om,eql,priors,
                fixQ=FALSE,prr="",
                rMult=1,b0Mult=1,kMult=1,
                rtn=function(x,y){
                 cbind(model.frame(mcf(FLQuants(
                            omF=catch(x)/stock(x),
                            omB=stock(x),
                            mpF=harvest(y),
                            mpB=stock(y))),drop=T),
                            ll =unlist(c(x@ll)))})
  {
  
  bd=FLBRP2biodyn(eql)
  catch(bd)=catch(om)
  
  params(bd)[c("r","k","p","b0")]=unlist(c(priors[,c("rPrime","k","p","b0")]))
  params(bd)["r" ]=params(bd)["r" ]*rMult
  params(bd)["k" ]=params(bd)["k" ]*kMult
  params(bd)["b0"]=params(bd)["b0"]*b0Mult
  
  bd@stock=FLQuant(c(params(bd)["k"]*params(bd)["b0"]),
                   dimnames=list(age="all",year=dimnames(stock(om))$year))
  bd=window(bd,end=dims(stock(om))$maxyear-1)
  
  bd=mp:::fwd(bd,catch=catch(bd))
  
  u =window((stock(om)[,- dim(stock(om))[2]]+
             stock(om)[,-1])/2,start=20)

  setParams(bd)=u
  setControl(bd)=params(bd)
  
  if (fixQ)
    control(bd)["q1",c("phase")]=c(-1)
  else
    control(bd)["q1",c("phase")]=c(1)
  
  if (prr=="q"){
    bd@priors["q1","weight"]=1
    bd@priors["q1","a"]     =1 
  }else if (prr=="bmsy"){
    bd@priors["bmsy","weight"]=1
    bd@priors["bmsy","a"]     =priors[,"bmsy"]
  }else if (prr=="fmsy"){
    bd@priors["fmsy","weight"]=1
    bd@priors["fmsy","a"]     =priors[,"fmsy"]
  }else if (prr=="msy"){
    bd@priors["msy","weight"]=1
    bd@priors["msy","a"]     =priors[,"msy"]
  }else if (prr=="r"){
    bd@priors["r","weight"]=1
    bd@priors["r","a"]     =priors[,"r"]
  }else if (prr=="k"){
    bd@priors["k","weight"]=1
    bd@priors["k","a"]     =priors[,"k"]
  }else if (prr=="rPrime"){
    bd@priors["r","weight"]=1
    bd@priors["r","a"]     =priors[,"rPrime"]
  }else if (prr=="kPrime"){
    bd@priors["k","weight"]=1
    bd@priors["k","a"]     =priors[,"kPrime"]}

  control(bd)[c("p","b0"),"phase"]=-1
  
  bd=fit(bd,u)

  res=rtn(om,bd)

  if (FALSE) 
    ggplot(res)+
      geom_line(aes(year,omB,group=i))+
      geom_line(aes(year,mpB,group=i),col="red")
  
  res}

goodFn<-function(det){
  
  rms=ddply(det,.(i,rMult,fixQ,prr),
            function(x) sum(((x$mpB-x$omB)/x$omB)^2,na.rm=T)/length(x$omB))
  
  rms =rms[do.call(order,rms[,c("i","V1")]),]
  good=rms[!duplicated(rms[,"i"]),]
  
  good[do.call(order,list(V1=good[,"V1"])),]}

pMeasure=function(stk,brp,proxy="msy"){
  
  res=FLQuants(stock  =stock(stk)%/%FLBRP:::refpts(brp)[proxy,"biomass"],
               ssb    =ssb(  stk)%/%FLBRP:::refpts(brp)[proxy,"ssb"],
               rec    =rec(  stk)%/%FLBRP:::refpts(brp)[proxy,"rec"],
               catch  =catch(stk)%/%FLBRP:::refpts(brp)[proxy,"yield"],
               fbar   =fbar( stk)%/%FLBRP:::refpts(brp)[proxy,"harvest"],
               harvest=(catch(stk)/stock(stk))%/%(FLBRP:::refpts(brp)[proxy,"yield"]/FLBRP:::refpts(brp)[proxy,"biomass"]))
  
  model.frame(res,drop=T)}

pm<-function(om,eql){
  dat=pMeas(om,eql)
  
  res=with(dat, cbind(kobeSmry(stock,harvest),
                      dRate(catch,0),
                      dRate(catch,0.05),
                      dRate(catch,0.10),
                      diags:::av(catch),
                      diags:::av(harvest)))
  
  names(res)[10:14]=c("catch","catch5","catch10","aavCatch","aavF")
  res}


## ------------------------------------------------------------------------
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

fitBd<-function(om,eql,priors,
                fixQ =FALSE,
                prr  ="",
                rMult=1,b0Mult=1,kMult=1,
                cpue =c("stock","catch","juvenile","mature"),
                hrate=TRUE,
                rtn=function(x,y,z){
                  cbind(model.frame(mcf(FLQuants(
                    ##Absolute
                    omF=catch(x)/stock(x)[,dimnames(catch(x))$year],
                    omB=stock(x),
                    mpF=catch(y)/stock(y)[,dimnames(catch(y))$year],
                    mpB=stock(y),
                    
                    ##Relative
                    omF.=(catch(x)/stock(x)[,dimnames(catch(x))$year])%/%(FLBRP:::refpts(z)["msy","yield"]/FLBRP:::refpts(z)["msy","biomass"]),
                    omB.=stock(x)%/%FLBRP:::refpts(z)["msy","biomass"],
                    mpF.=(catch(y)/stock(y)[,dimnames(catch(y))$year])%/%refpts(y)["fmsy"],
                    mpB.=stock(y)%/%refpts(y)["bmsy"])),drop=T))})
{
  bd=FLBRP2biodyn(eql)
  catch(bd)=catch(om)

  params(bd)[c("r","k","p","b0")]=unlist(c(priors[,c("rPrime","k","p","b0")]))
  params(bd)["r" ]=params(bd)["r" ]*rMult
  params(bd)["k" ]=params(bd)["k" ]*kMult
  params(bd)["b0"]=params(bd)["b0"]*b0Mult
  
  bd@stock=FLQuant(c(params(bd)["k"]*params(bd)["b0"]),
                   dimnames=list(age="all",year=dimnames(stock(om))$year))
  bd=window(bd,end=dims(stock(om))$maxyear-1)
  
  bd=mp:::fwd(bd,catch=catch(bd))
  
  if (hrate)
    effort=catch(om)%/%stock(om)
  else
    effort=fbar(om)
  
  u=switch(substr(cpue[1],1,1),
           s=(stock(om)[,- dim(stock(om))[2]]+stock(om)[,-1])/2,
           c=catch(om)%/%effort,
           j=FLQuant(apply(catch.n(om)*catch.wt(om)*(1-mat(om))%/%effort,c(2,6), sum)),
           m=FLQuant(apply(catch.n(om)*catch.wt(om)*mat(om)%/%effort,    c(2,6), sum)))
  
  setParams(bd)=u
  setControl(bd)=params(bd)

  control(bd)["q1","val"]=1
  if (fixQ)
    control(bd)["q1","phase"]=-1
  else
    control(bd)["q1","phase"]= 2
    
  if (prr=="q"){
    bd@priors["q1",   "weight"]=0.3
    bd@priors["q1","a"]        =1 
  }else if (prr=="bmsy"){
    bd@priors["bmsy", "weight"]=0.3
    bd@priors["bmsy","a"]      =priors[,"bmsy"]
  }else if (prr=="fmsy"){
    bd@priors["fmsy", "weight"]=0.3
    bd@priors["fmsy","a"]      =priors[,"fmsy"]
  }else if (prr=="msy"){
    bd@priors["msy",  "weight"]=0.3
    bd@priors["msy","a"]       =priors[,"msy"]
  }else if (prr=="r"){
    bd@priors["r",    "weight"]=0.3
    bd@priors["r","a"]         =priors[,"r"]
  }else if (prr=="k"){
    bd@priors["k",    "weight"]=0.3
    bd@priors["k","a"]         =priors[,"k"]
  }else if (prr=="rPrime"){
    bd@priors["r",    "weight"]=0.3
    bd@priors["r","a"]         =priors[,"rPrime"]
  }else if (prr=="kPrime"){
    bd@priors["k",    "weight"]=0.3
    bd@priors["k","a"]         =priors[,"kPrime"]}

  bd=fit(bd,FLQuants(u))
  
  res=rtn(om,bd,eql)
  
  res}

goodFn<-function(det){
  
  rms=ddply(det,.(i,rMult,fixQ,prr),
            function(x) sum(((x$mpB-x$omB)/x$omB)^2,na.rm=T)/length(x$omB))
  
  rms =rms[do.call(order,rms[,c("i","V1")]),]
  good=rms[!duplicated(rms[,"i"]),]
  
  good[do.call(order,list(V1=good[,"V1"])),]}

mpStats=function(x,y){
  
  cbind(
    model.frame(refpts(y)),
    model.frame(FLQuants(stock  =stock(  y)[,ac(range(y)["maxyear"])]/refpts(y)["bmsy"],
                         harvest=harvest(y)[,ac(range(y)["maxyear"])]/refpts(y)["fmsy"],
                         catch  =catch(  y)[,ac(range(y)["maxyear"])]/refpts(y)["msy"]),drop=TRUE))
}


stats=function(x,y){
  model.frame(mcf(FLQuants(
    omF=catch(x)/stock(x),
    omB=stock(x),
    mpF=harvest(y),
    mpB=stock(y))),drop=T)}

