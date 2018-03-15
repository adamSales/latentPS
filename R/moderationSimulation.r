library(lme4)
library(rstan)

load('fittedModels/stanMod.RData')


U <- extract(main,'studEff')[[1]]

set.seed(613)
U <- U[sample(1:nrow(U),1000),]

multImp <- apply(U,1,
                 function(u) summary(
                                 lmer(Y~treatment*u+poly(xirt,2)+race+sex+spec+state+(1|schoolid2)+(1|teachid2),
                                      data=dat))$coef)

multImp <- t(multImp)
u <- U[nrow(U),]
mod <- lmer(Y~treatment*u+poly(xirt,2)+race+sex+spec+state+(1|schoolid2)+(1|teachid2),data=dat)
ccc <- summary(mod)$coef
colnames(multImp) <- c(rownames(ccc),paste0(rownames(ccc),'SE'),paste0(rownames(ccc),'Tstat'))

print('Est. Effect:')
mean(multImp[,'treatment:u'])

hist(multImp[,'treatment:u'],freq=FALSE)
lines(density(extract(main,'b1')[[1]]))


print('SD of estimated effects')
sd(multImp[,'treatment:u'])

print('multiple imputation SD')
sqrt(var(multImp[,'treatment:u'])+mean(multImp[,'treatment:uSE']^2))

save(multImp,file='fittedModels/multImp.RData')



### try naive way (aka "what does Y buy us?")
load('glmerRasch.RData')

library(rstan)
source('src/prelimStan.r')

sdatRasch <- within(sdat,{
 nstudT <- sum(Z)
 nstudC <- sum(Z==0)
 teacherT <- as.numeric(as.factor(teacher[Z==1]))
 nteacherT <- max(teacherT)
 schoolT <- as.numeric(as.factor(school[Z==1]))
 nschoolT <- max(schoolT)
 teacherC <- as.numeric(as.factor(teacher[Z==0]))
 nteacherC <- max(teacherC)
 schoolC <- as.numeric(as.factor(school[Z==0]))
 nschoolC <- max(schoolC)
 Xtrt <- X[Z==1,]
 Xctl <- X[Z==0,]
 studentM <- as.numeric(as.factor(studentM))
})

sdatRaschY <- within(sdat,{
 nstudT <- sum(Z)
 nstudC <- sum(Z==0)
 teacherT <- as.numeric(as.factor(teacher[Z==1]))
 nteacherT <- max(teacherT)
 schoolT <- as.numeric(as.factor(school[Z==1]))
 nschoolT <- max(schoolT)
 teacherC <- as.numeric(as.factor(teacher[Z==0]))
 nteacherC <- max(teacherC)
 schoolC <- as.numeric(as.factor(school[Z==0]))
 nschoolC <- max(schoolC)
 Xtrt <- cbind(X[Z==1,],Y[Z==1])
 Xctl <- cbind(X[Z==0,],Y[Z==0])
 ncov <- ncov+1
 studentM <- as.numeric(as.factor(studentM))
})

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

rm <- stan('src/raschT.stan',data=sdatRasch,chains=3,iter=2000)
save(rm,file='stanRasch.RData')

rmY <- stan('src/raschT.stan',data=sdatRaschY,chains=3,iter=2000)
save(rmY,file='stanRaschY.RData')

rmDraws <- extract(rm)
U <- cbind(rmDraws$studEffT,rmDraws$studEffC)

outDat <- with(sdatRasch,data.frame(
    Y=c(Y[Z==1],Y[Z==0]),
    teacher=c(teacherT,teacherC),
    school=c(schoolT,schoolC),
    Z=c(rep(1,nstudT),rep(0,nstudC))))

X <- with(sdatRasch,rbind(Xtrt,Xctl))

library(lme4)
multImp <- apply(U[seq(1,nrow(U),length=1000),],1,
                 function(u) summary(
                     lmer(Y~Z*u+X+(1|teacher)+(1|school),
                          data=outDat))$coef)
u <- U[1,]
nm <- summary(lmer(Y~Z*u+X+(1|teacher)+(1|school),
                                      data=outDat))$coef
rownames(multImp) <- outer(rownames(nm),colnames(nm),'paste')
save(multImp,file='multImp2.RData')

mean(multImp['Z:u Estimate',])
sqrt(var(multImp['Z:u Estimate',])+mean(multImp['Z:u Std. Error',]^2))

rmYDraws <- extract(rmY)
UY <- cbind(rmYDraws$studEffT,rmYDraws$studEffC)

multImpY <- apply(UY,1,
                 function(u) summary(
                     lmer(Y~Z*u+X+(1|teacher)+(1|school),
                          data=outDat))$coef)

rownames(multImpY) <- outer(rownames(nm),colnames(nm),'paste')
save(multImpY,file='multImp2Y.RData')

mean(multImpY['Z:u Estimate',])
sqrt(var(multImpY['Z:u Estimate',])+mean(multImpY['Z:u Std. Error',]^2))
