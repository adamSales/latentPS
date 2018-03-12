#### PPMC
### use Q3 (levy et al 2009)

print(load('output/mainModel.RData'))
studEff <- rstan::extract(main,'studEff')[[1]]
secEff <- rstan::extract(main,'secEff')[[1]]




stud2 <- sort(unique(sdat$studentM))[table(sdat$studentM)>1]
studM2 <- sdat$studentM[sdat$studentM%in%stud2]
sec2 <- sdat$section[sdat$studentM%in%stud2]
grad2 <- sdat$grad[sdat$studentM%in%stud2]

studEffObs <- studEff[,stud2]

studEffObs <- studEffObs[seq(1,nrow(studEffObs),length.out=1000),]
secEffThin <- secEff[seq(1,nrow(secEff),length.out=1000),]

nstud <- ncol(studEffObs)
nsec <- ncol(secEffThin)

studObs <- as.numeric(as.factor(studM2))
Xobs <- matrix(NA,nrow=nstud,ncol=nsec)
for(i in 1:length(sec2)) Xobs[studObs[i],sec2[i]] <- grad2[i]


ncomb <- prod(dim(Xobs))


ex <- function(iter){
    linPred <- outer(studEffObs[iter,],secEffThin[iter,],'+')
    linPred[is.na(Xobs)] <- NA
    plogis(linPred)
}

xrep <- function(prob){
    matrix(rbinom(ncomb,1,prob),nrow(prob),ncol(prob))
}

Q3rep <- function(prob){
    Xrep <- xrep(prob)
    e <- Xrep-prob
    cor(e,use='pairwise')
}

Q3obs <- function(prob){
    e <- Xobs-prob
    cor(e,use='pairwise')
}


oneIter <- function(iter){
    prob <- ex(iter)
    diff <- Q3rep(prob)-Q3obs(prob)
    diff <- diff>0
    diff
}

ppps <- function(seed=613){
    set.seed(seed)
    diffs <- vapply(1:1000,oneIter,matrix(TRUE,134,134))
    byPair <- apply(diffs,c(1,2),mean,na.rm=TRUE)
    print(mean(diffs,na.rm=TRUE))
    invisible(byPair)
}


PPP <- ppps()
save(PPP,file='output/PPPdim.RData')
