source('R/prelimStan.r')

#### first main model
if(exists('runMain')&runMain){
 main <- stan('R/psmod.stan',data=sdat,warmup=1500,chains=10,iter=5000)

## cat('\n\n\n\n',rep('-',40),'\n','MAIN MODEL',rep('-',40),'\n\n\n')
## print(main,c('a0','a1','b0','b1'),c(0.05,0.95))

 save(main,sdat,file='output/stanMod.RData'); rm(main); gc();
}

if(exists('runInt') & runInt){
 ### interactions model
 sdatInt <- makeStanDat(dat,advance,xInteract=TRUE)
 xint <- stan('R/psmod.stan',data=sdatInt,warmup=1500,iter=5000,chains=4)
## cat('\n\n\n\n',rep('-',40),'\n','INTERACTIONS',rep('-',40),'\n\n\n')
## print(xint,c('a0','a1','b0','b1'),c(0.05,0.95))
 save(xint,sdatInt,file='output/xInteractions.RData'); rm(xint); gc()
}

if(exists('runHard') & runHard){
### only sections with hints
 source('R/hardSections.r')
#cat('\n\n\n\n',rep('-',40),'\n','HARD SECTIONS',rep('-',40),'\n\n\n')
#print(hard,c('a0','a1','b0','b1'),c(0.05,0.95))
 rm(hard);gc()
}
### no teacher

if(exists('runNoTeach') & runNoTeach){
 noTeach <- stan('R/psmodNoTeacher.stan',data=sdat,warmup=1500,iter=3000,chains=4)
 save(noTeach,sdat,file='output/modNoTeacher.RData'); rm(noTeach); gc()
}

### pooled usage data
if(exists('runPooledU')& runPooledU){
 stanModPooledU <- stan('R/pooledU.stan',data=sdat,warmup=1500,iter=3000,chains=4)
 save(pooledU,sdat,file='output/pooledU.RData'); rm(pooledU); gc()
}

## 2PL
if(exists('run2pl')&run2pl){
 stanMod2pl <- stan('R/psmod2pl.stan',data=sdat,warmup=1500,iter=5000,chains=4)
 save(stanMod2pl,sdat,file='output/stanMod2pl.RData'); rm(stanMod2pl); gc()
}

## 3PL
if(exists('run3pl')&run3pl){
 stanMod3pl <- stan('R/psmod3pl.stan',data=sdat,warmup=1500,iter=5000,chains=4)
 save(stanMod3pl,sdat,file='output/stanMod3pl.RData'); rm(stanMod3pl); gc()
}

## mbarModel
if(exists('runMbar')&runMbar){
source('R/prelimMbar.r')
 mbarMod <- stan('R/psmodObs.stan',data=sdatObs); save(mbarMod,sdatObs,file='output/mbarModel.RData')
 rm(mbarMod); gc()
}

## BC model
if(exists('runBC')&runBC){
 bcMod <- stan('R/psmodBC.stan',data=sdat); save(bcMod,sdat,file='output/bcModel.RData')
 rm(bcMod);gc()
}

## full data
if(exists('runFull')&runFull){
 totDatFull <- dataPrep(datOrig,advanceOrig,discard=FALSE)
 sdatFull <- makeStanDat(totDatFull$dat,totDatFull$advance)
 full <- stan('R/psmod.stan',data=sdat); save(full,file='output/fullMod.RData')
}

## fake
if(exists('runFake')&runFake){
 source('R/fakeModelsStan.r')
}

if(exists('runFakeBS')&runFakeBS){
 source('R/fakeModelsStanBS.r')
}
