#############################################
# PLOT FUNCTIONS							#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################


#CONTAINS
#plotSingleTrials
#plotSolution


plotSingleTrials <-
function(swalesol,what=c('all','sum'),which=numeric(0)) 
#plot singletrial data
{
	#layout(1)
	
	what = match.arg(what,c('all','sum'))
	
	#data
	eegdat = .eeg.data.data(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	sampRate = .eeg.data.sampRate(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	
	#calculate absolute maxima and minima for plotwindow
	rng = range(c(.swale.solution.model(swalesol),eegdat))
	
	#set discard
	if(length(.swale.solution.discard(swalesol))==0) .swale.solution.discard(swalesol)=rep(0,nrow(eegdat))
	
	#waves amps and lats
	waveform = .swale.internal.waves(.swale.solution.internal(swalesol))
	TP = .basis.matrix(.swale.internal.basis(.swale.solution.internal(swalesol)))
	dTP = .basis.deriv(.swale.internal.basis(.swale.solution.internal(swalesol)))
	amps = .swale.internal.amps(.swale.solution.internal(swalesol))
	lats = .swale.internal.lats(.swale.solution.internal(swalesol))
	
	swalemod = .swale.solution.model(swalesol)
		
	#create plot
	for(trial in 1:nrow(eegdat)) {
		if(!.swale.solution.discard(swalesol)[trial]) {
			par(ask=T,las=1)
			plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=rng,xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main=paste('trial [',trial,']',sep=''))
			axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
			axis(2)
			
			#plot avg data
			lines(eegdat[trial,],lty=1,lwd=3,col=gray(.5))
			
			legs = character(0)
			cols = numeric(0)
			sumwaves = numeric(ncol(eegdat))
			
			#plot waves
			for(wave in 1:ncol(.swale.solution.waveform(swalesol))) {
				#st.wave = (TP%*%waveform[,wave]*amps[trial,wave]+(dTP%*%waveform[,wave])*lats[trial,wave]) #check on function
				st.wave = (.swale.solution.waveform(swalesol)[,wave]+.swale.solution.derivwave(swalesol)[,wave]*.swale.solution.latency(swalesol)[trial,wave])*.swale.solution.amplitude(swalesol)[trial,wave]
				#sumwaves = sumwaves + st.wave #check on function
				
				if(what=='all') {
					lines(st.wave,lty=2,lwd=2,col=wave+1)
					legs=c(legs,paste('Wave(',wave,'): ',round(.swale.solution.latency(swalesol)[trial,wave]*(1/sampRate)*1000),' @ ',round(.swale.solution.amplitude(swalesol)[trial,wave],2),sep=''))
					cols=c(cols,wave+1)
				}
			}
			
			#plot cumulativemodel 
			lines(swalemod[trial,],lwd=5,col=1)
			#lines(sumwaves,col=gray(1),lty=1,lwd=1) #check on function
			
			#plot legend
			legs = c('Data','Model',legs)
			cols = c(gray(.5),1,cols)
			legend(1,rng[2],legend=legs,col=cols,lwd=3,lty=1,bty='n')
			
		} #discardloop
	}
	#layout(1)
	#par(ask=F)
}

plotSolution <-
function(swalesol) 
#plot singletrial data
{
	
	#make layout
	lmat = matrix(1,3,3)
	lmat[,3]=c(2,3,6)
	lmat[3,]=c(4,5,6)
	layout(lmat)
	
	#data
	eegdat = .eeg.data.data(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	sampRate = .eeg.data.sampRate(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	
	#averages
	avgdat = apply(eegdat,2,mean)
	avgmod = apply(.swale.solution.model(swalesol),2,mean)
	
	#plot averge and averagemodel
	rng = range(c(avgmod,avgdat))
	par(ask=F,las=1)
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=rng,xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main=paste(.eeg.data.channel(.swale.internal.eeg.data(.swale.solution.internal(swalesol))),':',.eeg.data.condition(.swale.internal.eeg.data(.swale.solution.internal(swalesol))),' - ','Model fit (AIC): ',round(.swale.solution.aic(swalesol)),sep=''))
	axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	axis(2)
	lines(avgdat,lty=1,lwd=3,col=gray(.5))
	lines(avgmod,lty=1,lwd=3,col=gray(0))
	
	
	#plot waves
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=c(-2,2),xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main='Waveform')
	axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	axis(2)
	for(wave in 1:ncol(.swale.solution.waveform(swalesol))) {
		dmwave = (.swale.solution.waveform(swalesol)[,wave]-mean(.swale.solution.waveform(swalesol)[,wave]))
		dmwave = dmwave / max(abs(dmwave))
		lines(dmwave,lwd=2,lty=1,col=wave+1)
	}
	
	#plot derivs
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=c(-2,2),xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main='Derivative')
	axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	axis(2)
	for(wave in 1:ncol(.swale.solution.derivwave(swalesol))) {
		dmwave = (.swale.solution.derivwave(swalesol)[,wave]-mean(.swale.solution.derivwave(swalesol)[,wave]))
		dmwave = dmwave / max(abs(dmwave))
		lines(dmwave,lwd=2,lty=1,col=wave+1)
	}
	
	
	#plot amp densities
	disc = which(.swale.solution.discard(swalesol)!=0)
	if(length(disc)>0) {
		amps = matrix(.swale.solution.amplitude(swalesol)[-which(.swale.solution.discard(swalesol)!=0),],,ncol(.swale.solution.amplitude(swalesol)))
	} else amps = .swale.solution.amplitude(swalesol)
		
	allden = try(density(amps),silen=T)
	if(class(allden)!='try-error') {
		plot(NA,NA,xlim=range(allden$x),ylim=range(allden$y)*1.5,xlab='',ylab='',main='Amplitude distributions',bty='n',axes=F)
		axis(1)
		axis(2,at=axTicks(2))
		for(wave in 1:ncol(.swale.solution.amplitude(swalesol))) {
			lines(allden$x,density(amps[,wave])$y,col=wave+1,lwd=3)
		}
	} else {
		plot(NA,NA,xlim=c(0,1),ylim=c(0,1),bty='n',axes=F,xlab='',ylab='')
	}

	#plot latency densities
	disc = which(.swale.solution.discard(swalesol)!=0)
	if(length(disc)>0) {
		lats = matrix(.swale.solution.latency(swalesol)[-which(.swale.solution.discard(swalesol)!=0),],,ncol(.swale.solution.latency(swalesol)))
	} else lats = .swale.solution.latency(swalesol)
	
	allden = try(density(lats),silen=T)	
	if(class(allden)!='try-error') {
		plot(NA,NA,xlim=range(allden$x),ylim=range(allden$y)*1.5,xlab='',ylab='',main='Latency distributions',bty='n',axes=F)
		axis(1)
		axis(2,at=axTicks(2))
		for(wave in 1:ncol(.swale.solution.latency(swalesol))) {
			lines(allden$x,density(lats[,wave])$y,col=wave+1,lwd=3)
		}
	} else {
		plot(NA,NA,xlim=c(0,1),ylim=c(0,1),bty='n',axes=F,xlab='',ylab='')
	}
	
	#plot discarded data
	rng = range(c(.swale.solution.model(swalesol),eegdat))
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=rng,xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main=paste('Discarded trials (',length(which(.swale.solution.discard(swalesol)!=0)),')',sep=''))
	axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	axis(2)
	if(length(.swale.solution.discard(swalesol))>0) {
		for(i in which(.swale.solution.discard(swalesol)!=0)) {
			lines(eegdat[i,],col=i+1,lwd=1)
		}
		
		if(length(which(.swale.solution.discard(swalesol)!=0))>1) lines(apply(eegdat[which(.swale.solution.discard(swalesol)!=0),],2,mean),lty=1,lwd=3,col=gray(0)) 
	}  
	
	#layout(1)
}



eeg.plot <- 
function(dat) 
{
	require(colorRamps)
	
		#if(min(dat)<0) {
		layout(cbind(rbind(matrix(1,8,10),matrix(3,2,10)),rep(2,10)))
		par(las=1,mar=c(1,4,4,1)+0.1,mgp=c(3,1,0))
		image(t(dat),ylab='trial',axes=F,col=matlab.like(64))
		axis(2,at=seq(0,1,.25),labels=round(seq(0,dim(dat)[1],dim(dat)[1]/4),0))
		#axis(1,at=seq(0,1,.1),labels=rep('',11))
		
		legbar=seq(round(min(dat),0),round(max(dat),0))
		legbar=matrix(legbar,1,length(legbar))
		par(mar=c(5.1, 0, 4.1, 4.1),mgp=c(2,1,0))
		image(legbar,col=matlab.like(64),axes=F,xlim=c(0,.01))
		#axis(4,at=seq(0,1,.1),labels=c(round(seq(min(legbar),0,abs(min(legbar))/5),0),round(seq(0,max(legbar),max(legbar)/5),0)[-1]))
		mtext(expression(paste(mu,'V',sep='')),4,cex=0.7,outer=F,adj=-1.5)
		
		avgdat=apply(dat,2,mean)
		par(mar=c(5.1, 4.1, .1, 1.1),mgp=c(3,1,0))
		#image(matrix(avgdat,length(avgdat),1),col=matlab.like(64),axes=F)
		plot(1:length(avgdat),avgdat,bty='n',axes=F,xlab='',ylab='',type='l',lwd=2)
		axis(1,pos=0,at=axTicks(1),label=rep('',length(axTicks(1))))
		axis(2)
		#axis(1)
		mtext('time in units',1,cex=0.7,padj=2.5)
		
		#layout(1)
		
	#} else { cat('No negatives in dataset -nyi-\n') }
	
}



plotTraces <-
function(swalesol) 
{
	#data
	eegdat = .eeg.data.data(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	sampRate = .eeg.data.sampRate(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	
	#averages
	avgdat = apply(eegdat,2,mean)
	
	yrng = range(eegdat)
	par(las=1,ask=F)
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=yrng,bty='n',xlab='time(ms)',ylab=expression(paste(mu,'V',sep='')),axes=F,main='Single Trials')
	axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	axis(2)
	
	for(i in 1:nrow(eegdat)) lines(eegdat[i,],lwd=1,col=i)
		
	lines(avgdat,lwd=4,col=1)
	#layout(1)

}

