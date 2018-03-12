source('R/prelimStan.r')

if(!exists('allModels')) allModels <- FALSE

#### first main model
if(allModels | exists('runMain')) if(allModels|runMain){
 main <- stan('R/psmod.stan',data=sdat,warmup=1500,chains=10,iter=5000)

## cat('\n\n\n\n',rep('-',40),'\n','MAIN MODEL',rep('-',40),'\n\n\n')
## print(main,c('a0','a1','b0','b1'),c(0.05,0.95))

 save(main,sdat,file='output/mainMod.RData'); rm(main); gc();
}

if(allModels | exists('runInt') ) if(allModels| runInt){
 ### interactions model
 sdatInt <- makeStanDat(dat,advance,xInteract=TRUE)
 xint <- stan('R/psmod.stan',pars=c('a1','b0','b1'),data=sdatInt,warmup=1500,iter=5000,chains=4)
## cat('\n\n\n\n',rep('-',40),'\n','INTERACTIONS',rep('-',40),'\n\n\n')
## print(xint,c('a0','a1','b0','b1'),c(0.05,0.95))
 save(xint,sdatInt,file='output/xInteractions.RData'); rm(xint); gc()
}

if(allModels | exists('runHard') ) if(allModels| runHard){
### only sections with hints
 source('R/hardSections.r')
#cat('\n\n\n\n',rep('-',40),'\n','HARD SECTIONS',rep('-',40),'\n\n\n')
#print(hard,c('a0','a1','b0','b1'),c(0.05,0.95))
 rm(hard);gc()
}
### no teacher

if(allModels | exists('runNoTeach') ) if(allModels| runNoTeach){
 noTeach <- stan('R/psmodNoTeacher.stan',pars=c('a1','b0','b1'),,data=sdat,warmup=1500,iter=3000,chains=4)
 save(noTeach,sdat,file='output/modNoTeacher.RData'); rm(noTeach); gc()
}

### pooled usage data
if(allModels | exists('runPooledU')) if(allModels| runPooledU){
 pooledU <- stan('R/pooledU.stan',pars=c('a1','b0','b1'),data=sdat,warmup=1500,iter=3000,chains=4)
 save(pooledU,sdat,file='output/pooledU.RData'); rm(pooledU); gc()
}

## 2PL
if(allModels | exists('run2pl')) if(allModels|run2pl){
 stanMod2pl <- stan('R/psmod2pl.stan',pars=c('a1','b0','b1','studEff'),data=sdat,warmup=1500,iter=5000,chains=4)
 save(stanMod2pl,sdat,file='output/stanMod2pl.RData'); rm(stanMod2pl); gc()
}

## 3PL
if(allModels | exists('run3pl')) if(allModels|run3pl){
 stanMod3pl <- stan('R/psmod3pl.stan',pars=c('a1','b0','b1','studEff'),data=sdat,warmup=1500,iter=5000,chains=4)
 save(stanMod3pl,sdat,file='output/stanMod3pl.RData'); rm(stanMod3pl); gc()
}

## mbarModel
if(allModels | exists('runMbar')) if(allModels|runMbar){
source('R/prelimMbar.r')
 mbarMod <- stan('R/psmodObs.stan',data=sdatObs); save(mbarMod,sdatObs,file='output/mbarModel.RData')
 rm(mbarMod); gc()
}

## BC model
if(allModels | exists('runBC')) if(allModels|runBC){
 bcMod <- stan('R/psmodBC.stan',data=sdat,iter=4000,warmup=1500,chains=6); save(bcMod,sdat,file='output/bcModel.RData')
 rm(bcMod);gc()
}

## full data
if(allModels | exists('runFull')) if(allModels|runFull){
 totDatFull <- dataPrep(datOrig,advanceOrig,discard=FALSE)
 sdatFull <- makeStanDat(totDatFull$dat,totDatFull$advance)
 full <- stan('R/psmod.stan',pars=c('a1','b0','b1'),data=sdat); save(full,file='output/fullMod.RData')
}

## raw scores
if(allModels | exists('runRaw')) if(allModels|runRaw){
 hs2 <- read.csv('~/Box Sync/CT/data/RANDstudyData/H2_algebra_rcal_20121119_fieldid.csv')
 datRaw <- dat
 datRaw$Y <- hs2$t2score[match(datRaw$field_id,hs2$field_id)]
 sdatRaw <- makeStanDat(datRaw,advance)
 raw <- stan('R/psmod.stan',pars=c('a1','b0','b1'),data=sdatRaw,iter=4000,warmup=1500,chains=6); save(raw,sdatRaw,file='output/rawMod.RData')
 rm(raw);gc()
}

if(allModels | exists('runRawBC')) if(allModels|runRawBC){
 hs2 <- read.csv('~/Box Sync/CT/data/RANDstudyData/H2_algebra_rcal_20121119_fieldid.csv')
 datRaw <- dat
 datRaw$Y <- hs2$t2score[match(datRaw$field_id,hs2$field_id)]
 datRaw$Y <- datRaw$Y/sd(datRaw$Y)
 sdatRaw <- makeStanDat(datRaw,advance)
 rawbc <- stan('R/psmodBC.stan',pars=c('a1','b0','b1','lambda'),data=sdatRaw,iter=4000,warmup=1500,chains=6); save(rawbc,sdatRaw,file='output/rawModBC.RData')
 rm(rawbc);gc()
}

if(allModels | exists('runRawSimp')) if(allModels|runRawSimp){
 hs2 <- read.csv('~/Box Sync/CT/data/RANDstudyData/H2_algebra_rcal_20121119_fieldid.csv')
 datRaw <- dat
 datRaw$Y <- hs2$t2score[match(datRaw$field_id,hs2$field_id)]
 #datRaw <- subset(datRaw,Y>0) ## delete 8 cases, 6 ctl, 2 trt
 sdatRaw <- makeStanDat(datRaw,advance)
 bc <- with(sdatRaw,MASS::boxcox(Y[Y>0]~X[Y>0,]+Z[Y>0]+factor(teacher[Y>0])))
 lambda <- bc$x[which.max(bc$y)] #0.5455
 sdatRaw <- within(sdatRaw,{
     Y <- (Y^lambda-1)/lambda
     sdPooled <- sqrt(((sum(Z)-1)*var(Y[Z==1])+(sum(Z==0)-1)*var(Y[Z==0]))/(nstud-2))
     Y <- Y/sdPooled })
 rawbc2 <- stan('R/psmod.stan',pars=c('a1','b0','b1'),data=sdatRaw,iter=4000,warmup=1500,chains=6); save(rawbc2,sdatRaw,file='output/rawModBC2.RData')

}

if(allModels | exists('runMS')) if(allModels|runMS)
 source('R/ms.r')

## fake
if(allModels | exists('runFake')) if(allModels|runFake){
 source('R/fakeModelsStan.r')
}

if(allModels | exists('runFakeBS')) if(allModels|runFakeBS){
 source('R/fakeModelsStanBS.r')
}

if(allModels | exists('runFakeBSBC')) if(allModels|runFakeBSBC){
 source('R/fakeModelsStanBSBC.r')
}
