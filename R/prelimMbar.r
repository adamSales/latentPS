### only look at year 2
### and only at students who weren't in  year 1
library(splines)
library(rstan)
library(dplyr)
memory.limit(50000)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#save(list=ls(),file=paste0('prev',Sys.time(),'.RData'))
#rm(list=ls())

#load('~/Box Sync/CT/data/RANDstudyData/HSdata.RData')
#load('~/Box Sync/CT/data/sectionLevelUsageData/advanceData.RData')


## ### estimated eta using only treatment group
## library(jagstools)
## load('~/Google Drive/CTmodels/fullUsage.RData')
## usageMod <- mod
## U <- jagsresults(usageMod,'studEff')[,1]


dataPrepObs <- function(dat,advance,discard=TRUE){
 dat <- dat[!dat$field_id%in%dat$field_id[dat$year==1],]
 print(table(dat$year))

 dat <- droplevels(dat)
### look at promotion
### delete "CP" sections. Logic: these are missing data, since we don't get to see if student
### would have graduated. Student random effects calculated using observed sections. FIML?
 advance <- droplevels(subset(advance, status%in%c('graduated','promoted')))

 advance <- advance[advance$field_id%in%dat$field_id[dat$treatment==1],]

 ### just look at algebra I sections--- likely different advance patterns in other curricula
 ### make sure to keep algebra i units that are also part of other curricula
 algUnit <- unique(advance$unit[advance$curriculum=='algebra i'])
 advance <- subset(advance,unit%in%algUnit)

 advance$grad <- advance$status=='graduated'

 advance <- droplevels(advance)


 ### discard some pairs
 ## discard treatment schools with no usage data
 ## (do a robustness check with everything left in afterwards)
 percUse <- function(scl) mean(unique(dat$field_id[dat$schoolid2==scl])%in%advance$field_id)
 obsUse <- vapply(unique(dat$schoolid2[dat$treatment==1]),percUse,1)

 obsUse <- unique(dat$schoolid2[dat$treatment==1])[obsUse>0.1]
 obsUse <- c(as.character(obsUse),as.character(unique(dat$schoolid2[dat$treatment==0])))

 if(discard)
  dat <- dat[dat$schoolid2%in%obsUse,]

 ## discard pairs with only a treatment or a control school
 pairTrtTab <- with(dat,table(pair,treatment))
 trtVar <- apply(pairTrtTab,1,prod)
 dat <- dat[trtVar[dat$pair]>0,]

 advance <- advance[advance$field_id%in%dat$field_id,]

 aaa <- aggregate(advance$grad,by=list(section=advance$section),FUN=mean)
 aaa$n <- as.vector(table(advance$section))

 if(discard) advance <- subset(advance,section%in%aaa$section[aaa$n>100 & aaa$x<1] & year==2)

 ## delete duplicate sections
 advance <- advance%>%group_by(field_id,section)%>%summarize(grad=min(grad,na.rm=TRUE))

 advance <- droplevels(advance)
 dat <- droplevels(dat)

 mbar <- advance%>%group_by(field_id)%>%summarize(nsec=n(),ngrad=sum(grad,na.rm=TRUE),mbar=mean(grad,na.rm=TRUE))
                                        #aggregate(advance$grad,by=list(field_id=advance$field_id),FUN=mean,na.rm=TRUE)
 dat <- merge(dat,mbar,all.x=TRUE)

 dat
}

makeStanDatObs <- function(dat){
 stanDat <- list()
 dat <- droplevels(dat)

 stanDat$nstudC <- nrow(dat)-sum(dat$treatment)
 stanDat$nstudTO <- sum(dat$treatment==1 & !is.na(dat$mbar))
 stanDat$nstudTM <- sum(dat$treatment==1 & is.na(dat$mbar))
 stanDat$nteacher <- nlevels(dat$teachid2)
 stanDat$nschool <- nlevels(dat$schoolid2)
 stanDat$npair <- nlevels(dat$pair)

 Z <- dat$treatment
 Z[Z==1 & is.na(dat$mbar)] <- 2
 teacher <- as.numeric(dat$teachid2)
 school <- as.numeric(dat$schoolid2)
 pair <- as.numeric(dat$pair)

 stanDat$teacherTO <- teacher[Z==1]
 stanDat$teacherTC <- teacher[Z==0]
 stanDat$teacherTM <- teacher[Z==2]
 stanDat$teacherC <- teacher[Z==0]
 stanDat$pairTO <- pair[Z==1]
 stanDat$pairTM <- pair[Z==2]
 stanDat$pairC <- pair[Z==0]
 stanDat$schoolTO <- school[Z==1]
 stanDat$schoolTM <- school[Z==2]
 stanDat$schoolC <- school[Z==0]


 X <- scale(model.matrix(~poly(xirt,2)+race+sex+spec+state,data=dat)[,-1])
 stanDat$XtO <- X[Z==1,]
 stanDat$XtM <- X[Z==2,]
 stanDat$Xc <- X[Z==0,]

 stanDat$ncov <- ncol(X)

 stanDat$YtO <- dat$Y[Z==1]
 stanDat$Yc <- dat$Y[Z==0]
 stanDat$YtM <- dat$Y[Z==2]

 stanDat$MbarTO <- dat$mbar[Z==1]

 stanDat
}

datObs <- dataPrepObs(dat,advance)

sdatObs <- makeStanDatObs(datObs)

#NOT RUN:
#mod <- stan('~/gitRepos/ctaiAdvance/psmodObs.stan',data=sdat); save(mod,sdat,file='fittedModels/mbarModel.RData')
#printStan(mod)

