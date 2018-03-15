library(ggplot2)
library(tikzDevice) ## allows latex code in figure
library(rstan)

options( tikzLatexPackages = c(
getOption( "tikzLatexPackages" ),
"\\usepackage{amsmath,amsfonts}"
    ))

load('data/HSdata.RData')
load('data/advanceData.RData')

source('R/prelimStan.r')
pooledSD <- with(sdat, sqrt(((sum(Z)-1)*var(Y[Z==1])+(sum(1-Z)-1)*var(Y[Z==0]))/(nstud-2)))
source('R/prelimMbar.r')

##################
### load main model results
################
load('output/mainMod.RData')
load('output/mbarModel.RData')

###############3
### summaries from main model
#################
draws <- extract(main)

### for "multImp" and "trtEff"
set.seed(613)
U <- draws$studEff
Usamp <- U[sample(1:nrow(U),1000),]

### for sampleSizeEta & etaDiff
draw <- 1000
U <- U[,sort(unique(sdat$studentM))]
eta <- U[draw,]
etasd <- apply(U,2,sd)

### for "usageModel"
sdEta <- sqrt(mean(apply(draws$studEff,1,var)))
Eeta <- colMeans(draws$studEff)

#draws$studEff <- Usamp

summMain <- summary(main)[[1]]

save(draws,draw,eta,etasd,sdEta,Eeta,summMain,Usamp,file='output/smallMain.RData')




########################
###  Mbar model figure
########################
samps <- extract(mbarMod)
drawMb <- which.min(abs(samps$b1-mean(samps$b1)))
plotDatObs <- with(sdatObs,data.frame(Y=c(YtO,YtM,Yc)/pooledSD,mbar=c(MbarTO,samps$MbarTM[drawMb,],samps$MbarC[drawMb,]),Z=c(rep(1,nstudTO),rep(1,nstudTM),rep(0,nstudC))))
plotDatObs$treat <- ifelse(plotDatObs$Z==1,'Treatment','Control')
plotDatObs$slope <- (samps$a1[drawMb]+ifelse(plotDatObs$treat=='Control',0,samps$b1[drawMb]))/pooledSD
plotDatObs <- within(plotDatObs,int <- mean(Y[treat=='Control'])-mean(slope[treat=='Control'])*mean(mbar[treat=='Control'])+ifelse(treat=='Control',0,samps$b0[drawMb]))

#plotDatObs <- within(plotDatObs, int <- int-( mean(int+slope*mbar)-mean(Y)))
plotDatObs <- plotDatObs[order(plotDatObs$treat),]
plotDatObs$treat2 <- plotDatObs$treat

tikz(file = "figure/mbarModel.tex",
  standAlone = T,
  width  = 6, height  = 6)
print(ggplot(plotDatObs,aes(mbar,Y,fill=treat,group=treat,color=treat))+geom_point(size=1)+
      geom_abline(aes(intercept=int,slope=slope,color=treat),size=4)+
      geom_abline(aes(intercept=int,slope=slope),size=2)+
    scale_colour_manual(values=c('red','blue'))+
    labs(group=NULL,fill=NULL,alpha=NULL)+xlab('$\\bar{m}_T$')+
    ylab('Posttest Score')+theme(legend.position='top',text=element_text(size=15))+
    guides(color = guide_legend(title=NULL,override.aes=list(alpha=1,size=3),keywidth=3),linetype=guide_legend(title=NULL,keywidth=1,override.aes=list(size=1))))#override.aes=list(size=2)))
dev.off()
setwd('figure'); tools::texi2dvi('mbarModel.tex', pdf = T, clean = T); setwd('..')

###################
### mbar vs sample size
###################

## smart jittering:
datObs$mbarJ <- datObs$mbar
datObs$nsecJ <- datObs$nsec
tab <- with(datObs,table(mbar,nsec))
mult <- which(tab>1,arr.ind=TRUE)
ms <- sort(unique(datObs$mbar))
ns <- sort(unique(datObs$nsec))
for(i in 1:nrow(mult)){
    w <- which(datObs$mbar==ms[mult[i,'mbar']] & datObs$nsec==ns[mult[i,'nsec']])
    s <- length(w)
    if(s>1){
        width=min(s*0.002,0.01)
        height=min(s*0.2,2)
        datObs$nsecJ[w] <- datObs$nsecJ[w]+runif(s,-width,width)
        datObs$mbarJ[w] <- datObs$mbarJ[w]+runif(s,-width,width)
    }
}



tikz(file='figure/mbarSampleSize.tex',
     standAlone=T,
     width=3,height=3)
print(ggplot(datObs,aes(mbarJ,nsecJ))+geom_point()+xlab('$\\bar{m}$')+ylab('$n_{sec}$')+theme(text=element_text(size=15)))
dev.off()
setwd('figure'); tools::texi2dvi('mbarSampleSize.tex', pdf = T, clean = T); setwd('..')

####################3
### eta vs sample size
######################
sdatLat <- sdat
nsec <- as.vector(table(sdatLat$studentM))
## etaDraws <- extract(main,'studEff')[[1]][,sort(unique(sdatLat$studentM))]
## eta <- etaDraws[draw,]
## etasd <- apply(etaDraws,2,sd)

plotDat <- data.frame(nsec=nsec,eta=eta,etasd=etasd)

tikz(file='figure/etaSampleSize.tex',
     standAlone=T,
     width=3,height=3)
print(ggplot(plotDat,aes(eta,nsec,size=1/etasd))+geom_point()+ylab(NULL)+#ylab('$n_{sec}$')+
    labs(size='$1/\\text{SE}(\\eta_T)$')+scale_size(range=c(.5,2))+guides(size=FALSE)+xlab('$\\eta_T$')+
    theme(text=element_text(size=15))+
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()))#+ggtitle('One Posterior Draw')#+xlab('$\\mathbb{E}\\eta$')
dev.off()
setwd('figure'); tools::texi2dvi('etaSampleSize.tex', pdf = T, clean = T); setwd('..')

#######################
#### mbar vs difficulty
####################3
secDiff <- -draws$secEff[draw,]
sss <- secDiff[sdatLat$sec]
mDiff <- aggregate(sss,list(stud=sdatLat$studentM),mean)
mbar <- aggregate(sdatLat$grad,list(sdatLat$studentM),mean)

mbarDiffDat <- data.frame(mbar=mbar$x,mDiff=mDiff$x)



tikz(file='figure/mbarDiff.tex',
     standAlone=T,
     width=3,height=3)
print(ggplot(mbarDiffDat,aes(mbar,mDiff))+geom_point()+xlab('$\\bar{m}$')+ylab('Avg. Sec. Difficulty')+
    theme(text=element_text(size=15)))
dev.off()
setwd('figure'); tools::texi2dvi('mbarDiff.tex', pdf = T, clean = T); setwd('..')

#######################
#### eta vs difficulty
####################3
plotDat$mDiff <- mDiff$x
tikz(file='figure/etaDiff.tex',
     standAlone=T,
     width=3,height=3)
print(ggplot(plotDat,aes(eta,mDiff,size=1/etasd))+geom_point()+#ylab('$Avg. Section Difficulty$')+
    labs(size='$1/\\text{SE}(\\eta_T)$')+scale_size(range=c(.5,2))+guides(size=FALSE)+xlab('$\\eta_T$')+
    theme(text=element_text(size=15))+
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()))#+ggtitle('One Posterior Draw')#+xlab('$\\mathbb{E}\\eta$')
dev.off()
setwd('figure'); tools::texi2dvi('etaDiff.tex', pdf = T, clean = T); setwd('..')


###################
#### usage coef
#################
library(coefplot)

## coefs <- summary(main,'betaU')[[1]]
coefs <- summMain[grep('betaU',rownames(summMain)),]
rownames(coefs) <- colnames(sdatLat$X)

## sqrt of average (over the draws) of the variance of eta
## sdEta <- sqrt(mean(apply(draws$studEff,1,var)))

coefs <- coefs/apply(sdatLat$X,2,sd)/sdEta

coefName <- c('Black/\n Multiracial','Hispanic/\n Native American','Male','Special Ed.','Gifted')
coefName <- factor(coefName,levels=coefName)
cpdf <- data.frame(Value=coefs[3:7,1],Coefficient=coefName,HighInner=coefs[3:7,'75%'],LowInner=coefs[3:7,'25%'],HighOuter=coefs[3:7,'97.5%'],LowOuter=coefs[3:7,'2.5%'],Model='model')

coefplot.data.frame(cpdf,title=NULL,xlab=expression(hat(beta)[std]),ylab=NULL,lwdOuter=0.5,lwdInner=1.5)+theme(text=element_text(size=15))
ggsave('figure/usageCoef.pdf',width=6,height=3)


##################
### pretest vs eta
##################

pdf('figure/pretestEta.pdf',width=6,height=6)
omar <- par()$mar
par(mar=omar+c(0,1,-1,0))
## Eeta <- colMeans(draws$studEff)
plot(dat$xirt,Eeta/sdEta,col=ifelse(dat$treatment==1,'blue','red'),
     xlab='Pretest (std)',ylab=expression(paste('E[',eta[T],'|x]')),cex.lab=1.25)
X <- scale(model.matrix(~poly(xirt,2),data=dat)[,-1])
xpred <- (X[,1]*mean(draws$betaU[,1])+X[,2]*mean(draws$betaU[,2]))/sdEta
set.seed(613)
samp <- sample(1:4000,100)
for(ss in samp){
    xpredS <- (X[,1]*draws$betaU[ss,1]+X[,2]*draws$betaU[ss,2])/sdEta
    lines(sort(dat$xirt),xpredS[order(dat$xirt)],col=adjustcolor('pink',0.5))
}
lines(sort(dat$xirt),xpred[order(dat$xirt)],lwd=2)
legend('bottomright',legend=c('Treated','Control (Imputed)','Model (Avg.)','Model (draws)'),col=c('blue','red','black','pink'),pch=c('o','o','.','.'),lwd=c(0.00,0.00,2,2))
dev.off()


##############################
#### Potential Outcomes Plot
###########################
a0 <- rnorm(length(draws$a1),mean(sdatLat$Y[sdatLat$Z==0]),sd(sdatLat$Y[sdatLat$Z==0])/sqrt(sum(sdatLat$Z==0)))
a1 <- draws$a1
b0 <- draws$b0
b1 <- draws$b1

studEff95 <- quantile(Usamp,c(0.025,0.975))
xx <- seq(studEff95[1],studEff95[2],length=100)
Yc <- outer(a1,xx)
Yc <- sweep(Yc,1,a0,'+')

YcUp <- apply(Yc,2,function(x) quantile(x,0.975))
YcDown <- apply(Yc,2,function(x) quantile(x,0.025))


Yt <- outer(a1+b1,xx)
Yt <- sweep(Yt,1,a0+b0,'+')
YtUp <- apply(Yt,2,function(x) quantile(x,0.975))
YtDown <- apply(Yt,2,function(x) quantile(x,0.025))

pdf('figure/potentialOutcomes.pdf',width=6,height=6)
curve(mean(a0)+mean(a1)*x,from=min(xx), to=max(xx),lwd=2,col='red',xlab=expression(eta[T]),ylab=expression(paste('E[',Y[Z],'|',eta[T],']',sep='')),ylim=range(c(YtDown,YcDown,YtUp,YcUp)),cex.lab=1.25)

curve(mean(a0)+mean(b0)+(mean(b1)+mean(a1))*x,add=TRUE,lwd=2,col='blue')
polygon(c(xx,rev(xx)),c(YcUp,rev(YcDown)),col=adjustcolor('red',0.1))
polygon(c(xx,rev(xx)),c(YtUp,rev(YtDown)),col=adjustcolor('blue',0.1))

legend('topleft',legend=c(expression(Y[C]),expression(Y[T])),col=c('red','blue'),lwd=2)
dev.off()

############################
######## main effect plot
#########################


pdMod <- function(mod,row=1,column=1,func){
    draws <- extract(mod)
    samp <- seq(1,length(draws$b1),length=1000)
    Usamp <- draws$studEff[samp,]
    iqr <- apply(Usamp,1,IQR)
    studEff95 <- quantile(Usamp,c(0.025,0.975))
    Usamp[Usamp<studEff95[1] | Usamp>studEff95[2]] <- NA
    trtEff <- sweep(sweep(Usamp,1,draws$b1[samp],'*'),1,draws$b0[samp],'+')



    if(missing(func)){
        func <- function(x) mean(draws$b0)+mean(draws$b1)*x
        knownTruth <- FALSE
    } else knownTruth <- TRUE
    truth <- curve(func,from=studEff95[1],to=studEff95[2],n=length(samp)/3)
    avg <- curve(mean(draws$b0)+x*mean(draws$b1),
                 from=studEff95[1],to=studEff95[2],n=length(samp)/3)
    postDraw <- curve(mean(draws$b0)+x*mean(draws$b1),
                      from=studEff95[1],to=studEff95[2],n=length(samp)-length(truth$x)-length(avg$x))
    x <- c(postDraw$x,truth$x,avg$x)
    y <- c(postDraw$y,truth$y,avg$y)
    if(knownTruth) truthOrAvg <- c(rep('Posterior\nDraws',length(postDraw$x)),rep('True\nEffect',length(truth$x)),rep('Posterior\nAverage',length(avg$x))) else
     truthOrAvg <- c(rep('Posterior\nDraws',length(postDraw$x)),rep('Posterior\nAverage',length(avg$x)+length(truth$x)))

#    if(knownTruth) title <- paste('True Effe

    pd <- data.frame(b0=draws$b0[samp],b1=draws$b1[samp],id=1:length(samp),row=row,column=column,xmin=studEff95[1],xmax=studEff95[2],ymin=min(trtEff,na.rm=T),ymax=max(trtEff,na.rm=T),x=x,y=y,
                     truthOrAvg=truthOrAvg,
                     iqr=iqr)
    pd
}

pdMain <- pdMod(main)
## pdMain <- within(pdMain,
## {
##     b0 <- b0/pooledSD
##     b1 <- b1/pooledSD*iqr
##     xmin <- xmin/mean(iqr)
##     xmax <- xmax/mean(iqr)
##     ymin <- ymin/pooledSD
##     ymax <- ymax/pooledSD
## }
## )
tikz('figure/mainEffects.tex', standAlone=T,
     width=6,height=5)
print(ggplot(pdMain)+
    geom_abline(aes(intercept=b0,slope=b1,group=id),color='red')+
    coord_cartesian(xlim=c(min(pdMain$xmin),max(pdMain$xmax)),
                    ylim=c(min(pdMain$ymin),max(pdMain$ymax)),expand=FALSE)+
    geom_line(aes(x=x,y=y,group=truthOrAvg,linetype=truthOrAvg,color=truthOrAvg,alpha=truthOrAvg),size=1.5)+
    xlab('$\\eta_T$')+ylab('$\\hat{\\tau}(\\eta_T)$')+
    labs(group=NULL,color=NULL,linetype=NULL)+
    scale_color_manual(values=c('black','red','black'))+scale_linetype_manual(values=c('solid','solid','dotted'))+
    scale_alpha_manual(values=c(1,0,1),guide=FALSE)+theme(legend.position='top')+
    theme(text=element_text(size=15),legend.key.width=unit(.5,'in')))
dev.off()
setwd('figure'); tools::texi2dvi('mainEffects.tex', pdf = T, clean = T); setwd('..')

#################
### mbar vs eta
#################
mbar$Eeta <- colMeans(draws$studEff[,mbar$Group.1])
mbar$nsec <- as.vector(table(sdatLat$studentM)[as.character(mbar$Group.1)])
mbar$nsec2 <- ifelse(mbar$nsec<50,mbar$nsec,50)
mbar$mbar <- mbar$x

tikz('figure/mbarVsEta.tex',standAlone=TRUE,width=6,height=6)
print(ggplot(mbar,aes(Eeta,mbar, size=nsec2))+geom_point()+labs(y='$\\bar{m}$',
           x='$E[\\eta_T]$',
           title=paste0('spearman rho=',round(with(mbar,cor(mbar,Eeta,method='spearman')),2)))+
    scale_size('\\# Sections',breaks=c(1,5,25,50),labels=c(1,5,25,'50+'),range=c(1,3))+theme(text=element_text(size=15)))
dev.off()
setwd('figure'); tools::texi2dvi('mbarVsEta.tex', pdf = T, clean = T); setwd('..')

######################
### results from fake models
#####################
estEff <- list()
load('output/noEffect.RData')
pd <- pdMod(noEff,1,1,function(x) x*0)
estEff$noEff <- summary(noEff,par=c('b0','b1'))[[1]]
rm(noEff);gc()

load('output/constEff.RData')
pd <- rbind(pd,pdMod(constEff,1,2,function(x) x*0+0.18))
estEff$constEff <- summary(constEff,par=c('b0','b1'))[[1]]
rm(constEff); gc()

load('output/linEff.RData')
pd <- rbind(pd,pdMod(linEff,2,1,function(x) mean(effs$b0)+x*mean(effs$b1)))
estEff$linEff <- summary(linEff,par=c('b0','b1'))[[1]]
rm(linEff);gc()

load('output/quadEff.RData')
U <- datF$U
mux <- mean(U)
te <- -(U-mux)^2
sigte <- sd(te)
mute <- mean(te/sigte*0.1)


pd <- rbind(pd,pdMod(quadEff,2,2,function(x) -(x-mux)^2*0.1/sigte-mute+0.13))
estEff$quadEff <- summary(quadEff,par=c('b0','b1'))[[1]]
rm(quadEff);gc()

pd$title <- NA
pd <- within(pd,{
    title[row==1 & column==1] <- paste0('$\\tau=0$\n$\\hat{\\tau}=',sprintf("%.2f",estEff$noEff['b0',1]),
                                        ifelse(estEff$noEff['b1',1]>0,'+',''),
                                        sprintf("%.2f",estEff$noEff['b1',1]),'\\eta_T$')
    title[row==1 & column==2] <- paste0('$\\tau=0.18+\\epsilon$\n$\\hat{\\tau}=',sprintf("%.2f",estEff$constEff['b0',1]),
                                        ifelse(estEff$constEff['b1',1]>0,'+',''),
                                        sprintf("%.2f",estEff$constEff['b1',1]),'\\eta_T$')
    title[row==2 & column==1] <- paste0('$\\tau=0.13-0.02\\eta_T$\n$\\hat{\\tau}=',sprintf("%.2f",estEff$linEff['b0',1]),
                                        ifelse(estEff$linEff['b1',1]>0,'+',''),
                                        (sprintf("%.2f",estEff$linEff['b1',1])),'\\eta_T$')
    title[row==2 & column==2] <- paste0('$\\tau=-0.04\\eta_T^2+0.01\\eta_T+0.02$\n$\\hat{\\tau}=',
                                        sprintf("%.2f",estEff$quadEff['b0',1]),
                                        ifelse(estEff$quadEff['b1',1]>0,'+',' '),
                                        sprintf("%.2f",estEff$quadEff['b1',1]),'\\eta_T$')})

pd <- within(pd, {
    title <- factor(title,levels=c(title[row==1 & column==1][1],
                                   title[row==1 & column==2][1],
                                   title[row==2 & column==1][1],
                                   title[row==2 & column==2][1]))})

tikz('figure/fakePlots.tex',standAlone=TRUE,width=6,height=6)
print(ggplot(pd)+
    geom_abline(aes(intercept=b0,slope=b1,group=id),color='red')+
    coord_cartesian(xlim=c(min(pd$xmin),max(pd$xmax)),ylim=c(min(pd$ymin),max(pd$ymax)),expand=FALSE)+
    geom_line(aes(x=x,y=y,group=truthOrAvg,linetype=truthOrAvg,color=truthOrAvg,alpha=truthOrAvg),size=1.5)+
    facet_wrap(~title,ncol=2)+xlab('$\\eta_T$')+ylab('$\\hat{\\tau}(\\eta_T)$')+
    labs(group=NULL,color=NULL,linetype=NULL)+
    #theme(strip.background = element_blank(),strip.text.x = element_blank(),strip.text.y=element_blank())+
    scale_color_manual(values=c('black','red','black'))+scale_linetype_manual(values=c('solid','solid','dotted'))+scale_alpha_manual(values=c(1,0,1),guide=FALSE)+theme(legend.position='top')+theme(text=element_text(size=15))+theme(legend.key.width=unit(.5,'in')))
    dev.off()
setwd('figure'); tools::texi2dvi('fakePlots.tex', pdf = T, clean = T); setwd('..')


#########################################
### eta_T vs Y
#########################################

draw <- which.min(abs(draws$b1-mean(draws$b1)))
plotDat <- with(sdatLat,data.frame(Y=scale(Y,center=mean(Y[Z==0]),scale=pooledSD),
                                   eta=scale(studEff[draw,],scale=IQR(draws$studEff[draw,])),
                                   Z=Z))

plotDat$treat <- ifelse(plotDat$Z==1,'Treatment','Control')
plotDat$slope <- (draws$a1[draw]+ifelse(plotDat$treat=='Control',0,draws$b1[draw]))*IQR(studEff[draw,])/pooledSD
plotDat$int <- ifelse(plotDat$treat=='Control',0,draws$b0[draw]/pooledSD)

#plotDat <- within(plotDat, int <- int-( mean(int+slope*eta)-mean(plotDat$Y)))
#plotDat <- plotDat[order(plotDat$treat),]
plotDat$treat2 <- plotDat$treat

tikz(file = "figure/etaModel.tex",
  standAlone = T,
  width  = 6, height  = 6)
print(ggplot(plotDat,aes(eta,Y,fill=treat,group=treat,color=treat))+geom_point(size=1)+
      coord_cartesian(xlim=quantile(plotDat$eta,c(0.005,0.995)),ylim=quantile(plotDat$Y,c(0.005,0.995)))+
      geom_abline(aes(intercept=int,slope=slope,color=treat),size=4,alpha=1)+#+scale_alpha_discrete(range=c(0.4,.8))+
      geom_abline(aes(intercept=int,slope=slope),color='black',size=2,alpha=1)+#+scale_alpha_discrete(range=c(0.4,.8))+
    scale_colour_manual(values=c('red','blue'))+
    labs(group=NULL,fill=NULL,alpha=NULL)+xlab('$\\eta_T$')+
    ylab('Posttest Score')+theme(legend.position='top',text=element_text(size=15))+
    guides(color = guide_legend(title=NULL,override.aes=list(alpha=1,size=3),keywidth=3),linetype=guide_legend(title=NULL,keywidth=1,override.aes=list(size=1))))
dev.off()
setwd('figure'); tools::texi2dvi('etaModel.tex', pdf = T, clean = T); setwd('..')

##########################
### posterior predictive dimensionality
###########################
load('output/PPPdim.RData')
pdf('figure/dimensionalityHist.pdf')
hist(PPP[upper.tri(PPP)],main='Posterior Predictve P-values')
dev.off()

### MDS
## distance from p-values: 1-p
d <- 1-PPP
diag(d) <- 0
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
#pdf('figure/pppMDS.pdf')
qplot(x, y, xlab="MDS Coordinate 1", ylab="MDS Coordinate 2",
     main="")#,pch=16)#,type='n')
#text(x,y,
                                        #dev.off()
ggsave('figure/pppMDS.pdf',width=4,height=4)


##### posterior predictive checks: figures
library(bayesplot)
### make Yrep
set.seed(613)
samp <- sample(1:nrow(draws$studEff),size=1000)
draws1k <- lapply(draws, function(x) if(NCOL(x)==1) x[samp] else x[samp,])
Ymean <- with(c(sdatLat,draws1k),
              teacherEffY[,teacher]+schoolEffY[,school]+pairEffect[,pair]+sweep(studEff,1,a1,'*')+
              sweep(sweep(studEff,1,b1,'*'),2,Z,'*')+betaY%*%t(X))

Yrep <- sapply(1:1000,function(i) rnorm(sdatLat$nstud,Ymean[i,],draws1k$sigY[i,sdatLat$Z+1]))

## overall Y
ppc_dens_overlay(sdatLat$Y,t(Yrep))+ggtitle('Pooled Treatment and Control')
ggsave('figure/ppcYoverall.jpg')

ppc_dens_overlay(sdatLat$Y[sdatLat$Z==1],t(Yrep)[,sdatLat$Z==1])+ggtitle('Treatment Group')
ggsave('figure/ppcYtrt.jpg')

ppc_dens_overlay(sdatLat$Y[sdatLat$Z==0],t(Yrep)[,sdatLat$Z==0])+ggtitle('Control Group')
ggsave('figure/ppcYctl.jpg')

### residual plots
fitd <- rowMeans(Yrep)
qplot(fitd,sdatLat$Y-fitd,xlab='Fitted Values',ylab='Residuals',main='Overall')
ggsave('figure/residPlotOverall.jpg')

qplot(fitd[sdatLat$Z==1],sdatLat$Y[sdatLat$Z==1]-fitd[sdatLat$Z==1],xlab='Fitted Values',ylab='Residuals',main='Treatment Group')
ggsave('figure/residPlotTrt.jpg')

qplot(fitd[sdatLat$Z==0],sdatLat$Y[sdatLat$Z==0]-fitd[sdatLat$Z==0],xlab='Fitted Values',ylab='Residuals',main='Control Group')
ggsave('figure/residPlotCtl.jpg')


pdf('figure/qqPlots.pdf', height=3,width=6)
par(mfrow=c(1,2))
qqnorm(sdatLat$Y[sdatLat$Z==1]-fitd[sdatLat$Z==1],main='Treatment Group')
qqline(sdatLat$Y[sdatLat$Z==1]-fitd[sdatLat$Z==1])
qqnorm(sdatLat$Y[sdatLat$Z==0]-fitd[sdatLat$Z==0],main='Control Group')
qqline(sdatLat$Y[sdatLat$Z==0]-fitd[sdatLat$Z==0])
dev.off()


### usage model
set.seed(613)
samp <- sample(1:1000,9)
lp <- with(c(sdatLat,draws1k),studEff[,studentM]+secEff[,section])
prob <- exp(lp)/(1+exp(lp))
p <- ppc_error_binned(sdatLat$grad,prob[samp,])
ggplot2::ggsave('figure/binnedplot.pdf',p)

### ppc m-bar
mbarRep <- apply(prob,1,function(p) aggregate(rbinom(length(sdatLat$grad),1,p),by=list(sdatLat$studentM),FUN=mean)$x)
mbar <- aggregate(sdatLat$grad,by=list(sdatLat$studentM),FUN=mean)$x
ppc_dens_overlay(mbar,t(mbarRep))
ggsave('figure/mbarPPC.jpg')

#####################
### plot results from robustness models
#####################
mainb1 <- summary(main,par='b1')[[1]]['b1',]
robustness <- data.frame(Value=mainb1['mean'],Coefficient='Main Model',HighInner=mainb1['75%'],LowInner=mainb1['25%'],
                                      HighOuter=mainb1['97.5%'],LowOuter=mainb1['2.5%'],Model='model')
ests <- c(mainb1['mean'],mainb1['sd'])
modFiles <- c('xInteractions','modNoTeacher','pooledU','stanMod2pl','stanMod3pl','bcModel','hardSections','fullMod','rawMod','rawModBC2')

for(mmm in modFiles){
    modName <- load(paste0('output/',mmm,'.RData'))
    modName <- modName[1]
    modSumm <- summary(get(modName),par='b1')[[1]]['b1',]
    ests <- rbind(ests,c(modSumm['mean'],modSumm['sd']))
    robustness <- rbind(robustness,
                        data.frame(Value=modSumm['mean'],Coefficient=modName,HighInner=modSumm['75%'],
                                   LowInner=modSumm['25%'],HighOuter=modSumm['97.5%'],LowOuter=modSumm['2.5%'],
                                   Model='model'))
    rm(list=modName); gc()
}

#robustness <- subset(robustness,Coefficient!='rawMod')

modelNames <- c('Main Model','Covariate Interations','No Teacher Effects','Pooled Rasch Model','2PL Mastery', '3PL Mastery','Power Transform Y','Exclude Sections w/o Hints','Include All Data','Raw Scores',
                'Raw Scores (BC Trans)')
robustness$Coefficient <- factor(modelNames,levels=rev(modelNames))

coefplot.data.frame(subset(robustness,modelNames!='Raw Scores'),xlab=expression(hat(b)[1]),ylab=NULL,lwdOuter=0.5,lwdInner=1.5)+theme(text=element_text(size=15))

rob2 <- subset(robustness,modelNames!='Raw Scores')
rob2$estX <- min(rob2$LowOuter)-0.05
rob2$seX <- min(rob2$LowOuter)-0.02
rob2$est <- paste0('  ',round(ests[-which(robustness$Coefficient=='Raw Scores'),'mean'],3),' (',round(ests[-which(robustness$Coefficient=='Raw Scores'),'sd'],3),')')
#rob2$se <- paste0(

ggplot(rob2,aes(x=Value,y=Coefficient))+geom_point(size=2)+
    geom_errorbarh(aes(xmin=LowInner,xmax=HighInner),size=1.5,height=0)+
    geom_errorbarh(aes(xmin=LowOuter,xmax=HighOuter),size=0.5,height=0)+
    geom_text(aes(x=estX,label=est),hjust='left')+
                                        #geom_text(aes(x=seX,label=se))+
    theme_minimal()+
    ylab('')+xlab(expression(b[1]))+scale_x_continuous(breaks=c(-0.1,-0.05,0,0.05))+geom_vline(xintercept=0,linetype='dotted')+theme(text=element_text(size=15))

ggsave('figure/robustness.pdf',height=8,width=6)
