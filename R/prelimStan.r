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

load('~/Box Sync/CT/data/RANDstudyData/HSdata.RData')
load('~/Box Sync/CT/data/sectionLevelUsageData/advanceData.RData')


## ### estimated eta using only treatment group
## library(jagstools)
## load('~/Google Drive/CTmodels/fullUsage.RData')
## usageMod <- mod
## U <- jagsresults(usageMod,'studEff')[,1]

dataPrep <- function(dat,advance,discard=TRUE){
 dat <- dat[!dat$field_id%in%dat$field_id[dat$year==1],]
 #print(table(dat$year))

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

### take out bridge to algebra, algebra ii, geometry
 advance <- advance[-grep('ii|geo|bridge',advance$curriculum),]

 advance$grad <- advance$status=='graduated'

 advance <- droplevels(advance)


 ### discard some pairs
 ## discard treatment schools with no usage data
 ## (do a robustness check with everything left in afterwards)
 percUse <- function(scl)
    length(intersect(unique(dat$field_id[dat$schoolid2==scl]),unique(advance$field_id)))/
        length(unique(dat$field_id[dat$schoolid2==scl]))

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

 advance <- advance%>%group_by(field_id,section)%>%summarize(grad=min(grad,na.rm=TRUE))

 advance <- droplevels(advance)
 dat <- droplevels(dat)

 list(dat=dat,advance=advance)
}

makeStanDat <- function(dat,advance,xInteract=FALSE){
 stanDat <- list()
 advance <- droplevels(advance)
 dat <- droplevels(dat)


 stanDat$nsecWorked <- nrow(advance)
 stanDat$nstud <- nrow(dat)
 stanDat$nteacher <- nlevels(dat$teachid2)
 stanDat$nschool <- nlevels(dat$schoolid2)
 stanDat$nsec <- nlevels(advance$section)
 stanDat$npair <- nlevels(dat$pair)

 stanDat$teacher <- as.numeric(as.factor(dat$teachid2))
 stanDat$pair <- as.numeric(dat$pair)
 stanDat$school <- as.numeric(dat$school)
 stanDat$studentM <- seq(stanDat$nstud)[match(advance$field_id,dat$field_id)]
 stanDat$section <- as.numeric(advance$section)

 stanDat$grad <- as.numeric(advance$grad)

 X <- model.matrix(~poly(xirt,2)+race+sex+spec+state,data=dat)[,-1]
 if(xInteract) X <- model.matrix(~poly(xirt,2)*(race+sex+spec)+(race+sex+spec)^2+state,data=dat)[,-1]
 stanDat$X <- scale(X)
 stanDat$ncov <- ncol(X)


 stanDat$Z <- as.numeric(dat$treatment)

 stanDat$Y <- dat$Y

 stanDat
}

totDat <- dataPrep(dat,advance)
datOrig <- dat
advanceOrig <- advance
dat <- totDat$dat
advance <- totDat$advance
rm(totDat);gc()

sdat <- makeStanDat(dat,advance)



