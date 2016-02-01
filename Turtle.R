vbk=c(0.045,0.0955,0.0749,0.0543,0.044,0.051,0.072,0.074,0.076,0.059,0.037,0.12,0.077,0.066,0.051)
exp(mean(log(vbk)))
sd(log(vbk))
hist(log(rlnorm(1000,mean(log(vbk)), sd(log(vbk)))))

linf=c(92.0,102,115,118.5,113.0,105.5,116.6,111.9,96.1,112.5,94.6,95.63,99,102.7)
exp(mean(log(linf)))
sd(log(linf))
hist(log(rlnorm(1000,mean(log(linf)), sd(log(linf)))))

lone=c(17.7,15)
exp(mean(log(lone)))
sd(log(lone))
hist(log(rlnorm(1000,mean(log(lone)), sd(log(lone)))))


lzero

setwd("/Users/robertahrens/Documents/Rob/A_Projects/BrazilTurtles")
data=read.table(file="lfdata.csv",header=TRUE,sep=",")
bin=data[,4]-data[,4]%%5+2.5
data=cbind(data,bin)
ii=which(data$YEAR<2012)
data=data[ii,]
xtabs(~bin+MONTH+YEAR,data)
ii=which(colSums(xx),arr.ind=TRUE)

