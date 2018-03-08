## list of tables \& numbers quoted in paper
##


source('R/prelimStan.r')

## covariate means/balance
miss <- NULL
for(i in c('race','sex','spec','xirt')) miss <- rbind(miss,
 c(sum(is.na(covs[[i]])),mean(is.na(covs[[i]])),error[i,'error']))
miss <- as.data.frame(miss)
miss$`Error Type` <- c('PFC','PFC','PFC','SRMSE')
rownames(miss) <- c('Race/Ethnicity','Sex','Special Education','Pretest')
names(miss)[1:3] <- c('# Missing','% Missing','Imputation Error')
miss[,2] <- as.integer(round(miss[,2]*100))
miss[,1] <- as.integer(miss[,1])
miss['Pretest','Imputation Error'] <- sqrt(miss['Pretest','Imputation Error'])/sd(covs$xirt,na.rm=TRUE)

library(RItools) ### this is development version from github

covs2 <- covs[rownames(dat),] ### same as dat but with NAs

balTest <- balanceTest(treatment~race+sex+spec+xirt+cluster(schoolid2)+strata(pair),data=dat,report=c('adj.means','std.diffs','chisquare.test'),include.NA.flags=T)

## use "unstrat" for group means cuz no one knows/cares about "adjusted means"
trtMeans <- round(balTest$results[,'Treatment','Unstrat']*100)
ctlMeans <- round(balTest$results[,'Control','Unstrat']*100)
stdDiff <- sprintf('%.2f',balTest$results[,'std.diff','pair'])
names(stdDiff) <- dimnames(balTest$results)[[1]]

impErr <- sprintf('%.2f',miss[,3])
names(impErr) <- rownames(miss)

results <- list()
results$balP <- sprintf('%.2f',balTest$overall['pair','p.value'])


cat('
\\begin{tabular}{cccrccc}
  \\hline
  \\hline
 &\\% Miss.& Imp. Err.&Levels& Ctl.& Trt.& Std. Diff.\\\\
  \\hline
  \\hline
\\multirow{3}{*}{Ethnicity}&\\multirow{3}{*}{',miss['Race/Ethnicity',2],'\\%}
  &\\multirow{3}{*}{',impErr['Race/Ethnicity'],'}& White/Asian & ',
ctlMeans['raceWhiteAsian'],'\\% &',trtMeans['raceWhiteAsian'],'\\% & ',stdDiff['raceWhiteAsian'],'\\\\
&&&Black/Multi &',
ctlMeans['raceBlackMulti'],'\\% &',trtMeans['raceBlackMulti'],'\\% & ',stdDiff['raceBlackMulti'],'\\\\
 & &&Hispanic/Nat.Am. & ',
ctlMeans['raceHispAIAN'],'\\% &',trtMeans['raceHispAIAN'],'\\% & ',stdDiff['raceHispAIAN'],' \\\\
\\hline
\\multirow{2}{*}{Sex}&\\multirow{2}{*}{',miss['Sex',2],'\\%}&\\multirow{2}{*}{',impErr['Sex'],'}& Female & ',ctlMeans['sexF'],'\\% &',trtMeans['sexF'],'\\% & ',stdDiff['sexF'],'\\\\
  &&&Male &',
ctlMeans['sexM'],'\\% &',trtMeans['sexM'],'\\% & ',stdDiff['sexM'],'\\\\
\\hline
\\multirow{3}{*}{Sp. Ed.}&\\multirow{3}{*}{',miss['Special Education',2],'\\%}&\\multirow{3}{*}{',impErr['Special Education'],'}&Typical &',
ctlMeans['spectypical'],'\\% &',trtMeans['spectypical'],'\\% & ',stdDiff['spectypical'],'\\\\
&&&Spec. Ed & ',
ctlMeans['specspeced'],'\\% &',trtMeans['specspeced'],'\\% & ',stdDiff['specspeced'],' \\\\
&&&Gifted &',
ctlMeans['specgifted'],'\\% &',trtMeans['specgifted'],'\\% & ',stdDiff['specgifted'],'\\\\
\\hline
Pretest&',miss['Pretest',2],'\\%&',impErr['Pretest'],'&   &',
sprintf('%.2f',balTest$results['xirt','Control','Unstrat']),'& ',
sprintf('%.2f',balTest$results['xirt','Treatment','Unstrat']),'& ',stdDiff['xirt'],'\\\\
   \\hline
&&&\\multicolumn{4}{c}{Overall Covariate Balance: p=',results$balP,'}\\\\
\\hline
\\hline
\\end{tabular}\n',sep='',file='output/covariateTable.tex')


### control students w/ usage data
results$numCtlUse <- sum(dat$field_id[dat$treatment==0]%in%advanceOrig$field_id)
results$perCtlUse <- round(results$numCtlUse/sum(dat$treatment==0)*100)

### percent of treatment students with observed mastery data
results$propObs <- round(mean(dat$field_id[dat$treatment==1]%in%advance$field_id)*100)

### sample sizes
results$ntot <- nrow(dat)
results$ntrt <- sum(dat$treatment)
results$nctl <- sum(dat$treatment==0)
results$nteach <- length(unique(dat$teachid2))
results$nschool <- length(unique(dat$schoolid2))
results$nSecWorked <- prettyNum(nrow(advance),big.mark=',')
results$mastPer <- round(mean(advance$grad)*100)

### mbar model results
load('output/mbarModel.RData')
sumMbar <- summary(mbarMod,par=c('b0','b1'))[[1]]
results$mbarB0int <- paste0('[',paste(round(sumMbar['b0',c('2.5%','97.5%')],1),collapse=','),']')
results$mbarB1int <- paste0('[',paste(round(sumMbar['b1',c('2.5%','97.5%')],1),collapse=','),']')
results$perMastHalf <- round(mean(sdatObs$MbarTO>.5)*100)

### comparing mbar to eta
load('output/mainModel.RData')
draws <- extract(main)
secDiff <- colMeans(-draws$secEff)
sss <- secDiff[sdat$sec]
mDiff <- aggregate(sss,list(stud=sdat$studentM),mean)
mbar <- aggregate(sdat$grad,list(sdat$studentM),mean)

mbarDiffDat <- data.frame(mbar=mbar$x,mDiff=mDiff$x)

results$mbarDiffCor <- sprintf('%.2f',with(mbarDiffDat,cor(mbar,mDiff,method='spearman')))

## covariates
perExp1 <- with(draws,apply(betaU,1,function(b) var(sdat$X%*%b)))
perExp <- perExp1/(perExp1 +draws$sigU^2)
results$perExp <- mean(perExp)

## main results
results$perRunNeg <- mean(draws$b1<0)

## compared to ATE
trtEff <- sweep(sweep(draws$studEff,1,draws$b1,'*'),1,draws$b0,'+')
ate <- mean(trt)

sdEta <- apply(draws$studEff,1,sd)
b1Std <- draws$b1/sdEta/ate

### save it all
attach(results)
save(list=names(results),file='output/results.RData')
detach(results)
