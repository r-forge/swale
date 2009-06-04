#############################################
# PLOT FUNCTIONS							#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

plot.iterates <- function(avgdat,yhat,fhat,dfhat,abhat,grad,dif,maxiter) 
{
	
	mn=min(c(avgdat,yhat))
	mx=max(c(avgdat,yhat))
	
	par(las=1)
	plot(NA,NA,xlab='time',ylab='mV',xlim=c(1,length(avgdat)),ylim=c(mn,mx),bty='n',main='fit')
	lines(avgdat,lwd=3,lty=1,col='gray')
	lines(yhat,lwd=3,lty=1,col='black')
	
	mn=min(c(fhat,dfhat))
	mx=max(c(fhat,dfhat))
	
	plot(NA,NA,xlab='time',ylab='mV',xlim=c(1,length(avgdat)),ylim=c(mn,mx),bty='n',main='fhat - dfhat')
	for(i in 1:ncol(fhat)) {
		lines(dfhat[,i],lwd=3,lty=1,col='gray')
		lines(fhat[,i],lwd=3,lty=1,col=1)
		
	}
	
	
	mn=min(grad)
	mx=max(grad[-1])
	
	plot(NA,NA,xlab='iterate',ylab='gradient/1e3',xlim=c(1,maxiter),ylim=c(mn,mx),bty='n',main='gradient',axes=F)
	axis(1)
	axis(2,at=axTicks(2),labels=axTicks(2)/1000)
	lines(grad,lwd=2,col=1)
	
	
	mn=min(dif)
	mx=max(dif)
	plot(NA,NA,xlab='iterate',ylab='fit/1e3',xlim=c(1,maxiter),ylim=c(mn,mx),bty='n',main='objective',axes=F)
	axis(1)
	axis(2,at=axTicks(2),labels=axTicks(2)/1000)
	lines(dif,lwd=2,col=1)
	
	
	mn=min(abhat[,1:ncol(fhat)])
	mx=max(abhat[,1:ncol(fhat)])
	plot(NA,NA,xlab='trials',ylab='amp',xlim=c(1,length(abhat[,1])),ylim=c(mn,mx),bty='n',main='amplitude',axes=F)
	axis(1)
	axis(2)
	for(i in 1:ncol(fhat)) {
		points(abhat[,i],pch=19,col=i)
	}
	
	mn=min(abhat[,(ncol(fhat)+1):ncol(abhat)])
	mx=max(abhat[,(ncol(fhat)+1):ncol(abhat)])
	plot(NA,NA,xlab='trials',ylab='amp',xlim=c(1,length(abhat[,1])),ylim=c(mn,mx),bty='n',main='latency',axes=F)
	axis(1)
	axis(2)
	p=1
	for(i in (ncol(fhat)+1):ncol(abhat)) {
		points(abhat[,i],pch=19,col=p)
		p=p+1
	}
	
}

plot.trials <- function(dat,sol) 
{
	x11(width=8,height=8)
	layout(1)
	
	avgdata=apply(dat$data,2,mean)
	
	mn=min(c(dat$data,sol$internal$yhat))
	mx=max(c(dat$data,sol$internal$yhat))
	
	for(tr in 1:nrow(dat$data)) {
		
		par(las=1)
		plot(NA,NA,xlab='time',ylab='mV',xlim=c(1,ncol(dat$data)),ylim=c(mn,mx),bty='n',main=paste('trial',tr))
		lines(dat$data[tr,],lwd=2,lty=1,col='gray')
		lines(apply(sol$internal$yhat,2,mean),lwd=2,lty=2,col='red')
		lines(sol$internal$yhat[tr,],lwd=3,lty=1,col='black')
		
		text(ncol(dat$data)/2,max(apply(sol$internal$yhat,2,mean))+30,paste('AMP=',round(sol$internal$abhat[tr,1],4),'@ LAT=',round(sol$solution$b[tr],0),sep=''))
		arrows(which.max(apply(sol$internal$yhat,2,mean)),mean(apply(sol$internal$yhat,2,mean)),which.max(apply(sol$internal$yhat,2,mean))+round(sol$solution$b[tr],0),mean(apply(sol$internal$yhat,2,mean)),lwd=2,length=.08,col='red')
		
		
		mean(apply(sol$internal$yhat,2,mean))
		cat(paste(paste('AMP=',round(sol$internal$abhat[tr,1],4),'@ LAT=',round(sol$solution$b[tr],0),sep='')))
		
		browser()
	}
	
}


plot.solution <- function(dat,sol)
{
	x11(width=7,height=7)
	layout(matrix(1:4,2,byrow=T))
	
	avgdata=apply(dat$data,2,mean)
	
	mn=min(c(avgdata,sol$solution$avgy))
	mx=max(c(avgdata,sol$solution$avgy))
	
	par(las=1)
	plot(NA,NA,xlab='time',ylab='mV',xlim=c(1,ncol(dat$data)),ylim=c(mn,mx),bty='n',main='avgerage fit')
	lines(avgdata,lwd=3,lty=1,col='gray')
	lines(sol$solution$avgy,lwd=3,lty=1,col='black')
	
	mn=min(c(avgdata,sol$solution$fhat))
	mx=max(c(avgdata,sol$solution$fhat))
	plot(NA,NA,xlab='time',ylab='mV',xlim=c(1,ncol(dat$data)),ylim=c(mn,mx),bty='n',main='normalized waveform')
	lines(sol$solution$fhat,lwd=3,lty=1,col=1)
	
	
	mn=min(c(sol$solution$fhat,sol$solution$dfhat))
	mx=max(c(sol$solution$fhat,sol$solution$dfhat))
	
	plot(NA,NA,xlab='time',ylab='mV',xlim=c(1,ncol(dat$data)),ylim=c(mn,mx),bty='n',main='estimated functions')
	lines(sol$solution$dfhat,lwd=2,lty=1,col='gray')
	lines(sol$solution$fhat,lwd=3,lty=1,col=1)
	
	mn=min(sol$solution$dfhat)
	mx=max(sol$solution$dfhat)
	plot(NA,NA,xlab='time',ylab='mV',xlim=c(1,ncol(dat$data)),ylim=c(mn,mx),bty='n',main='normalized derivative')
	lines(sol$solution$dfhat,lwd=3,lty=1,col=1)
	
}

plot.sims <- function(dat,sol)
{
	
	x11(width=7,height=7,title='sims')
	
	layout(matrix(1:4,2))
	
	par(las=1)
	plot(dat$amps,sol$internal$abhat[,1],xlab='real amps',ylab='est amps',bty='n',main='amplitude')
	plot(dat$lats,sol$internal$abhat[,2]/sol$internal$abhat[,1],xlab='real lats',ylab='est lats',bty='n',main='latency')
	hist(sol$internal$abhat[,1],main='amplitude')
	hist(sol$internal$abhat[,2]/sol$internal$abhat[,1]*-1,main='latency')
	
	cat('amplitude\n')
	print(cor.test(dat$amps,sol$internal$abhat[,1]))
	cat('latency\n')
	print(cor.test(dat$lats,sol$internal$abhat[,2]/sol$internal$abhat[,1]))

	
}

plot.freq <- function(vec,Hz=512)
{
	
	Nyq <- Hz/2
	blocklen <- length(vec)*(1/Hz)
	
	vecfft <- fft(vec)
	magn <- Mod(vecfft)
	
	magn <- magn[1:(length(magn)/2)]
	
	xax <- c(1:length(magn)/blocklen)
	
	x11(width=3,height=6,title='Low-pass filter settings')
	layout(1:2)
	par(las=1)
	plot(vec,type='l',bty='n',lwd=2,xlab='time samples',ylab='mV',main='Filter (spatial)')
	plot(xax,magn,type='l',lwd=2,bty='n',xlab='frequencies Hz',ylab='Magnitude',main='Frequency plot')
	
	
}


plot.data <- function(dat,each=FALSE)
{
	x11(width=7,height=7,title='data')
	
	mn=min(dat$data)
	mx=max(dat$data)
	
	par(las=1)
	plot(NA,NA,xlim=c(1,ncol(dat$data)),ylim=c(mn,mx),xlab='time (samples)',ylab='mv',bty='n')
	
	for(i in 1:nrow(dat$data)) lines(dat$data[i,],lwd=1,col=i)
	
	lines(apply(dat$data,2,mean),lwd=4,col=1)
		
	
}

