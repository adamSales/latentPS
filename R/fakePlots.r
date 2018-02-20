library(tikzDevice)
library(ggplot2)
library(rstan)

#set.seed(613)
#sample(1:9000,1000)

samp <- seq(1,9000,length=1000)

pdMod <- function(mod,row=1,column=1,func){
    draws <- extract(mod)
    Usamp <- draws$studEff[samp,]
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

    pd <- data.frame(b0=draws$b0[samp],b1=draws$b1[samp],id=1:length(samp),row=row,column=column,xmin=studEff95[1],xmax=studEff95[2],ymin=min(trtEff,na.rm=T),ymax=max(trtEff,na.rm=T),x=x,y=y,
                     truthOrAvg=truthOrAvg)
    pd
}

load('output/noEffect.RData')
pd <- pdMod(noEff,1,1,function(x) x*0)
#rm(noEff);gc()

load('output/constEff.RData')
pd <- rbind(pd,pdMod(constEff,1,2,function(x) x*0+0.18))
#rm(constEff); gc()

load('output/linEff.RData')
pd <- rbind(pd,pdMod(linEff,2,1,function(x) mean(effs$b0)+x*mean(effs$b1)))
#rm(linEff);gc()

load('output/quadEff.RData')
U <- datF$U
mux <- mean(U)
te <- -(U-mux)^2
sigte <- sd(te)
mute <- mean(te/sigte*0.1)


pd <- rbind(pd,pdMod(quadEff,2,2,function(x) -(x-mux)^2*0.1/sigte-mute+0.13))
#rm(quadEff);gc()

tikz('fakePlots.tex',standAlone=TRUE,width=6,height=6)
ggplot(pd)+
    geom_abline(aes(intercept=b0,slope=b1,group=id),color='red')+
    coord_cartesian(xlim=c(min(pd$xmin),max(pd$xmax)),ylim=c(min(pd$ymin),max(pd$ymax)),expand=FALSE)+
    geom_line(aes(x=x,y=y,group=truthOrAvg,linetype=truthOrAvg,color=truthOrAvg,alpha=truthOrAvg),size=1.5)+
    facet_grid(row~column)+xlab('$\\eta_T$')+ylab('$\\hat{\\tau}(\\eta_T)$')+
    labs(group=NULL,color=NULL,linetype=NULL)+
    theme(strip.background = element_blank(),strip.text.x = element_blank(),strip.text.y=element_blank())+
    scale_color_manual(values=c('black','red','black'))+scale_linetype_manual(values=c('solid','solid','dotted'))+scale_alpha_manual(values=c(1,0,1),guide=FALSE)+theme(legend.position='top')+theme(text=element_text(size=20))+theme(legend.key.width=unit(.5,'in'))
dev.off()
tools::texi2dvi('fakePlots.tex', pdf = T, clean = T)

pdMain <- pdMod(main)

tikz('mainEffects.tex', standAlone=T,
     width=6,height=5)
ggplot(pdMain)+
    geom_abline(aes(intercept=b0,slope=b1,group=id),color='red')+
    coord_cartesian(xlim=c(min(pd$xmin),max(pd$xmax)),ylim=c(min(pd$ymin),max(pd$ymax)),expand=FALSE)+
    geom_line(aes(x=x,y=y,group=truthOrAvg,linetype=truthOrAvg,color=truthOrAvg,alpha=truthOrAvg),size=1.5)+
    xlab('$\\eta_T$')+ylab('$\\hat{\\tau}(\\eta_T)$')+
    labs(group=NULL,color=NULL,linetype=NULL)+
    scale_color_manual(values=c('black','red','black'))+scale_linetype_manual(values=c('solid','solid','dotted'))+scale_alpha_manual(values=c(1,0,1),guide=FALSE)+theme(legend.position='top')+theme(text=element_text(size=20),legend.key.width=unit(.5,'in'))
dev.off()
tools::texi2dvi('mainEffects.tex', pdf = T, clean = T)
