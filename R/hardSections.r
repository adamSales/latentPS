library(dplyr)
load('~/Box Sync/CT/data/RANDstudyData/HSdata.RData')
load('~/Box Sync/CT/data/sectionLevelUsageData/advanceData.RData')
load('~/Box Sync/CT/data/problemLevelUsageData/probLevelData.RData')

## ### estimated eta using only treatment group
## library(jagstools)
## load('~/Google Drive/CTmodels/fullUsage.RData')
## usageMod <- mod
## U <- jagsresults(usageMod,'studEff')[,1]


dataPrep <- function(dat,advance,discard=TRUE){
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

 advance <- droplevels(advance)
 dat <- droplevels(dat)

 list(dat=dat,advance=advance)
}

totDat <- dataPrep(dat,advance)
datOrig <- dat
advanceOrig <- advance
dat <- totDat$dat
advance <- totDat$advance
rm(totDat);gc()


x$unit <- tolower(x$unit)
x$section <- tolower(x$section)
x$version <- tolower(x$version)

advance$unit <- tolower(advance$unit)
advance$section <- tolower(advance$section)
advance$version <- tolower(advance$version)

x <- x[x$field_id%in%advance$field_id,]

x$ts1 <- as.POSIXct(x$ts1,format="%m/%d/%y %H:%M")

Mode <- function(x){
    tab <- table(x)
    names(tab)[which.max(tab)]
}

errHintSum <- x%>%group_by(field_id,version,section,unit)%>%summarize(
                                                                nprob=n(),
                                                                nerr=sum(nerrs1>0,na.rm=TRUE),
                                                                perr=mean(nerrs1>0,na.rm=TRUE),
                                                                nhint=sum(nhints1>0,na.rm=TRUE),
                                                                phint=mean(nhints1>0,na.rm=TRUE),
                                                                ncurr=n_distinct(curriculum),
    curr=Mode(curriculum),
    start= min(ts1,na.rm=TRUE)
    )

adv2 <- merge(advance,errHintSum,all.x=TRUE)

adv2 <- adv2%>%arrange(start)%>%group_by(field_id)%>%mutate(secOrd=rank(start))

##keep1 <- adv2%>%filter(nerr>0|nhint>0)
keep1 <- adv2%>%filter(nhint>0)

sec100 <- keep1%>%group_by(section)%>%summarize(nn=n())
adv <- adv2%>%filter((nhint>0) & section%in%sec100$section[sec100$nn>=100])

save(adv,file='output/advanceHardHint.RData')

source('R/prelimStan.r')

adv$section <- factor(adv$section)
sdatHard <- makeStanDat(dat,adv)

hard <- stan('R/psmod.stan',data=sdatHard,iter=3000)
save(hard,sdatHard,file='output/hardSections.RData')
