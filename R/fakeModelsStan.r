source('R/prelimStan.r')

################################################
## Fake Models
##############################################

### use data, E[studEff] from main model
load('output/mainMod.RData')
dat$U <- colMeans(extract(main,'studEff')[[1]])
### U has a couple outliers:
dat$U[dat$U< -4] <- -4

effs <- extract(main,c('b0','b1'))

rm(main); gc()

##### make fake data:

### first delete control schools
datF <- subset(dat,treatment==1)

### now double the dataset
datF$schoolid2 <- as.character(datF$schoolid2)
datF$teachid2 <- as.character(datF$teachid2)
dat2 <- datF
dat2$schoolid2 <- paste0(dat2$schoolid2,'Fake')
dat2$teachid2 <- paste0(dat2$teachid2,'Fake')
dat2$field_id <- dat2$field_id*100+99
dat2$treatment <- 0

## ### take bootstrap samples within each school
## for(scl in unique(datF$schoolid2)){
##     ind <- which(datF$schoolid2==scl)
##     datF[ind,] <- datF[sample(ind,length(ind),replace=TRUE),]
##     dat2[ind,] <- dat2[sample(ind,length(ind),replace=TRUE),]
## }

datF <- rbind(datF,dat2)

datF$schoolid2 <- factor(datF$schoolid2)
datF$teachid2 <- factor(datF$teachid2)

### now delete usage data for "control" group
advanceF <- advance[advance$field_id%in%datF$field_id[datF$treatment==1],]

#########################################################
### Compile Stan data as before
#####################################################
sdatF <- makeStanDat(datF,advanceF)


noEff <- stan('R/psmod.stan',data=sdatF,iter=3000,chains=6)

##cat('\n\n\n\n',rep('-',40),'\n','TRUTH: NO EFFECT\n',rep('-',40),'\n\n\n')
##print(noEff,c('a0','a1','b0','b1'),c(0.05,0.95))
save(list=ls(),file='output/noEffect.RData'); rm(noEff); gc()




######## constant effect
########
sdatFConst <- within(sdatF,{
 te <- rnorm(sum(Z),0.18,0.1)
 Y[Z==1] <- Y[Z==1]+ te
 }
)

constEff <- stan('R/psmod.stan',data=sdatFConst,iter=3000,chains=6)

##cat('\n\n\n\n',rep('-',40),'\n','TRUTH: CONSTANT EFFECT ATE=0.18\n',rep('-',40),'\n\n\n')
##print(constEff,c('a0','a1','b0','b1'),c(0.05,0.95))
save(constEff,sdatF,file='output/constEff.RData')




##################################
######### linear effect ##########
##################################
datLin <- datF
datLin$te <- mean(effs$b0)+mean(effs$b1)*datF$U
datLin$Y[datLin$treatment==1] <- datLin$Y[datLin$treatment==1]+datLin$te[datLin$treatment==1]

sdatFlin <- makeStanDat(datLin,advanceF)

linEff <- stan('R/psmod.stan',data=sdatFlin,iter=3000,chains=6)
##cat('\n\n\n\n',rep('-',40),'\n','TRUTH: LINEAR EFFECT b1=',round(mean(effs$b1),2),'\n',rep('-',40),'\n\n\n')
##print(linEff,c('a0','a1','b0','b1'),c(0.05,0.95))

save(linEff,sdatFlin,file='output/linEff.RData');rm(linEff);gc();





########### quadratic effects
datQuad <- within(datF,{
    te <- -(U-mean(U))^2
    te <- te/sd(te)*0.1
    te <- te-mean(te)+0.13
    Y[treatment==1] <- Y[treatment==1]+te[treatment==1]
})

sdatFquad <- makeStanDat(datQuad,advanceF)



quadEff <- stan('R/psmod.stan',data=sdatFquad,iter=3000,chains=6)

##cat('\n\n\n\n',rep('-',40),'\n','TRUTH: QUADRATIC EFFECT',rep('-',40),'\n\n\n')
##print(constEff,c('a0','a1','b0','b1'),c(0.05,0.95))



save(quadEff,sdatFquad,file='output/quadEff.RData'); rm(quadEff); gc();





