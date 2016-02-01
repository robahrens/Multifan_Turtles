rm(list = ls())
setwd("/Users/robertahrens/ADMBprojects/BrazilTurtles/ADMB")
source("read.admb.R")
report=read.admb("turtlelf")

par(mfcol=c(2,2),mai=c(0.35,0.35,0.2,0.2),omi=c(0.3,0.3,0.4,0.1))
for(i in 1:length(report$lfdata[,1]))
{
	plot(report$lbinmids,report$lfdata[i,],type="n",xlab="",ylab="",yaxs="i")
	rect(report$lbinmids-0.5,0,report$lbinmids+0.5,report$lfdata[i,],col="grey",border="dark grey")
	abline(v=report$meanl[i,],lty=2)
	x=seq(range(report$lbinmids)[1],range(report$lbinmids)[2],length=100)
	points(report$lbinmids,report$pnl[i,],type="l")
	nobs=sum(report$lfdata[i,])
	for(j in 1:length(report$meanl[1,]))
	{
		mn=report$meanl[i,j]
		sd=report$sdl[i,j]
		y=dnorm(x,mn,sd)
		y=y/sum(y)*report$npage[i,j]*nobs
		points(x,y,type="l",col="white",lwd=2)
	
	}
	yy=20
	xx=qnorm(c(0.025,0.975),report$pmeanlone,report$psdmeanlone);
	if(i==1)points(report$pmeanlone,20,pch=19)
	if(i==1)points(xx,c(yy,yy),type="l")

	box()
	if(i==1)
	{
		mtext("Caripace Length (cm)",side=1,line=0.5,outer=TRUE)
		mtext("Frequency",side=2,line=0.5,outer=TRUE)
		mtext("Quarter 1",side=3,line=0)
		text(75,35,substitute(fage==a,list(a=round(report$fage,0))),adj=0)
		text(75,33,substitute(nage==a,list(a=round(report$nage,0))),adj=0)
		text(75,31,substitute(k==a,list(a=round(report$vbk,3))),adj=0)
		text(75,29,substitute(L[infinity]==a,list(a=round(report$linf,2))),adj=0)
		text(75,27,substitute(t[o]==a,list(a=round(report$tzero,2))),adj=0)
		text(75,25,substitute(AIC==a,list(a=round(report$AIC,2))),adj=0)
		text(75,23,substitute(AICc==a,list(a=round(report$AICc,2))),adj=0)
		text(75,21,substitute(BIC==a,list(a=round(report$BIC,2))),adj=0)
	}
	if(i==2)mtext("Quarter 2",side=3,line=0)
	if(i==3)mtext("Quarter 3",side=3,line=0)
	if(i==4)mtext("Quarter 4",side=3,line=0)
}
