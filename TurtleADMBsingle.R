rm(list = ls())
setwd("/Users/robertahrens/ADMBprojects/BrazilTurtles/ADMB")
source("read.admb.R")
report=read.admb("turtlelf")

	plot(report$lbinmids,report$lfdata,type="n",xlab="",ylab="",yaxs="i")
	rect(report$lbinmids-0.5,0,report$lbinmids+0.5,report$lfdata,col="grey",border="dark grey")
	abline(v=report$meanl,lty=2)
	x=seq(range(report$lbinmids)[1],range(report$lbinmids)[2],length=100)
	points(report$lbinmids,report$pnl,type="l")
	nobs=sum(report$lfdata)
	for(j in 1:length(report$meanl))
	{
		mn=report$meanl[j]
		sd=report$sdl[j]
		y=dnorm(x,mn,sd)
		y=y/sum(y)*report$npage[j]*nobs
		points(x,y,type="l",col="white",lwd=2)
	
	}
	box()
		mtext("Caripace Length (cm)",side=1,line=0.5,outer=TRUE)
		mtext("Frequency",side=2,line=0.5,outer=TRUE)
		mtext("Quarter 1",side=3,line=0)
		text(75,150,substitute(fage==a,list(a=round(report$fage,0))),adj=0)
		text(75,145,substitute(nage==a,list(a=round(report$nage,0))),adj=0)
		text(75,140,substitute(k==a,list(a=round(report$vbk,3))),adj=0)
		text(75,135,substitute(L[infinity]==a,list(a=round(report$linf,2))),adj=0)
		text(75,130,substitute(t[o]==a,list(a=round(report$tzero,2))),adj=0)
		text(75,125,substitute(AIC==a,list(a=round(report$AIC,2))),adj=0)
		text(75,120,substitute(AICc==a,list(a=round(report$AICc,2))),adj=0)
		text(75,115,substitute(BIC==a,list(a=round(report$BIC,2))),adj=0)
