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
	if(length(which)==0) trialvec = 1:nrow(eegdat) else trialvec = which
	
	colvec = c(4,3,5:ncol(waveform))
	
	for(trial in trialvec) {
		if(!.swale.solution.discard(swalesol)[trial]) {
			par(ask=T,las=1)
			plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=rng,xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main=paste('trial [',trial,']',sep=''))
			#axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
			prettyTime(1,c(1,ncol(eegdat)),Hz=sampRate)
			axis(2)
			
			#plot avg data
			lines(eegdat[trial,],lty=1,lwd=3,col=gray(.5))
			
			legs = character(0)
			cols = numeric(0)
			sumwaves = numeric(ncol(eegdat))

			#plot waves
			for(wave in 1:ncol(.swale.solution.waveform(swalesol))) {
				
				#st.wave = (TP%*%waveform[,wave]*amps[trial,wave]+(dTP%*%waveform[,wave])*lats[trial,wave]) #check on function
				st.wave = (.swale.solution.waveform(swalesol)[,wave]+.swale.solution.derivwave(swalesol)[,wave]*(.swale.internal.lats(.swale.solution.internal(swalesol))[trial,wave]/.swale.internal.amps(.swale.solution.internal(swalesol))[trial,wave]))*.swale.internal.amps(.swale.solution.internal(swalesol))[trial,wave]
				#sumwaves = sumwaves + st.wave #check on function
				
				if(what=='all') {
					lines(st.wave,lty=2,lwd=2,col=colvec[wave])
					if(!any(c(is.na(.swale.solution.amplitude(swalesol)[trial]),is.na(is.na(.swale.solution.latency(swalesol)[trial]))))) {
						points(.swale.solution.latency(swalesol)[trial,wave],.swale.solution.amplitude(swalesol)[trial,wave],col=colvec[wave],pch=19,cex=1.3)
					} #incorrect amp/lat estimation
					#legs=c(legs,paste('Wave(',wave,'): ',round(.swale.solution.latency(swalesol)[trial,wave]*(1/sampRate)*1000),' @ ',round(.swale.solution.amplitude(swalesol)[trial,wave],2),sep=''))
					#cols=c(cols,wave+1)
				}
			
			}
			
			#plot cumulativemodel 
			lines(swalemod[trial,],lwd=2,col=1)
			#lines(sumwaves,col=gray(1),lty=1,lwd=1) #check on function
			
			#plot legend
			#legs = c('Data','Model',legs)
			#cols = c(gray(.5),1,cols)
			#legend(1,rng[2],legend=legs,col=cols,lwd=3,lty=1,bty='n')
			
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
	#axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	prettyTime(1,c(1,ncol(eegdat)),Hz=sampRate)
	axis(2)
	lines(avgdat,lty=1,lwd=3,col=gray(.5))
	lines(avgmod,lty=1,lwd=3,col=gray(0))
	
	
	#plot waves
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=c(-2,2),xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main='Waveform')
	#axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	prettyTime(1,c(1,ncol(eegdat)),Hz=sampRate)
	axis(2)
	for(wave in 1:ncol(.swale.solution.waveform(swalesol))) {
		dmwave = (.swale.solution.waveform(swalesol)[,wave]-mean(.swale.solution.waveform(swalesol)[,wave]))
		dmwave = dmwave / max(abs(dmwave))
		lines(dmwave,lwd=2,lty=1,col=wave+1)
	}
	
	#plot derivs
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=c(-2,2),xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main='Derivative')
	#axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	prettyTime(1,c(1,ncol(eegdat)),Hz=sampRate)	
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
	
	hs = hist(amps,plot=F)
	yl = max(hs$counts*1.2)
	xl = range(amps)
	
	hist(amps[,1],xlim=xl,ylim=c(0,yl),xlab='Amplitude',ylab='Frequency',col=2,main='Amplitude distribution')
	if(ncol(.swale.solution.latency(swalesol))>1) {
		for(wave in 2:ncol(.swale.solution.latency(swalesol))) {
			hist(amps[,wave],add=T,col=wave+1)
		}
	}
	
	#plot latency densities
	disc = which(.swale.solution.discard(swalesol)!=0)
	if(length(disc)>0) {
		lats = matrix(.swale.solution.latency(swalesol)[-which(.swale.solution.discard(swalesol)!=0),],,ncol(.swale.solution.latency(swalesol)))
	} else lats = .swale.solution.latency(swalesol)
	
	hs = hist(lats,plot=F)
	yl = max(hs$counts*1.2)
	xl = range(lats)
	
	hist(lats[,1],xlim=xl,ylim=c(0,yl),xlab='Latency',ylab='Frequency',col=2,main='Latency distribution')
	if(ncol(.swale.solution.latency(swalesol))>1) {
		for(wave in 2:ncol(.swale.solution.latency(swalesol))) {
			hist(lats[,wave],add=T,col=wave+1)
		}
	}
	
	#plot discarded data
	rng = range(c(.swale.solution.model(swalesol),eegdat))
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=rng,xlab='time (ms)',ylab=expression(paste(mu,'V',sep='')),bty='n',axes=F,main=paste('Discarded trials (',length(which(.swale.solution.discard(swalesol)!=0)),')',sep=''))
	#axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	prettyTime(1,c(1,ncol(eegdat)),Hz=sampRate)
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
function(dat,sampRate=512,main='',order=NULL,latency=NULL) 
{
	require(colorRamps)
	
	if(!is.null(order)) {
		dat = dat[order,]
		latency = matrix(latency[order,],,ncol(latency))
	}
	
	layout(cbind(rbind(matrix(1,8,10),matrix(3,2,10)),rep(2,10)))
	par(las=1,mar=c(1,4,4,1)+0.1,mgp=c(3,1,0))
	image(t(dat),ylab='trial',axes=F,col=matlab.like(64),main=main)
	axis(2,at=seq(0,1,.25),labels=round(seq(0,dim(dat)[1],dim(dat)[1]/4),0))
	axis(1,at=seq(0,1,.1),labels=rep('',11))
	
	if(!is.null(latency)) {
		over = matrix(NA,dim(dat)[1],dim(dat)[2])
		for(i in 1:nrow(over)) {
			for(j in 1:ncol(latency)) {
				over[i,latency[i,j]] = 1
			}
		}
		image(t(over),add=T,col=1)
	}
	
	legbar=seq(round(min(dat),0),round(max(dat),0))
	legbar=matrix(legbar,1,length(legbar))
	par(mar=c(5.1, 0, 4.1, 4.1),mgp=c(2,1,0))
	image(legbar,col=matlab.like(64),axes=F,xlim=c(0,.01))
	axis(4,at=seq(0,1,.1),labels=c(round(seq(min(legbar),0,abs(min(legbar))/5),0),round(seq(0,max(legbar),max(legbar)/5),0)[-1]))
	mtext(expression(paste(mu,'V',sep='')),4,cex=0.7,outer=F,adj=-1.5)
		
	avgdat=apply(dat,2,mean)
	par(mar=c(5.1, 4.1, .1, 1.1),mgp=c(3,1,0))
	#image(matrix(avgdat,length(avgdat),1),col=matlab.like(64),axes=F)
	plot(1:length(avgdat),avgdat,bty='n',axes=F,xlab='',ylab=expression(paste(mu,'V',sep='')),type='l',lwd=2)
	#axis(1,pos=0,at=axTicks(1),label=round(axTicks(1)*(1/sampRate)*1000))
	prettyTime(1,c(1,length(avgdat)),Hz=sampRate)	
	axis(2)
	mtext('time (ms)',1,cex=0.7,padj=2.5)
	
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
	#axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	prettyTime(1,c(1,ncol(eegdat)),Hz=sampRate)
	axis(2)
	
	for(i in 1:nrow(eegdat)) lines(eegdat[i,],lwd=1,col=i)
		
	lines(avgdat,lwd=4,col=1)
	#layout(1)

}

plotTraceAverage <-
function(simdata,swalesol,what='signal',main='',yrng=NULL) 
{
	#data
	eegdat = .eeg.data.data(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	sampRate = .eeg.data.sampRate(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	
	#averages
	avgdat = apply(eegdat,2,mean)
	
	if(is.null(yrng)) yrng = range(apply(simdata$data,2,mean)*1.1)
	par(las=1,ask=F)
	plot(NA,NA,xlim=c(1,ncol(eegdat)),ylim=yrng,bty='n',xlab='time(ms)',ylab=expression(paste(mu,'V',sep='')),axes=F,main=main)
	#axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/sampRate)*1000))
	#axis(1,at=c(77,154,230,307),labels=c(150,300,450,600))
	prettyTime(1,c(1,ncol(eegdat)),Hz=sampRate)
	axis(2)
	#if(what=='signal') lines(apply(simdata$data,2,mean),lwd=2,col=1)
	if(what=='all') {
		lines(apply(simdata$data,2,mean),lwd=2,col=gray(.5))
		lines(swalesol@waveform,lwd=3,col=1)
	}
	
	#for(i in 1:nrow(eegdat)) lines(eegdat[i,],lwd=1,col=i)
	
	#lines(avgdat,lwd=4,col=1)
	#layout(1)
	
}

eeg.plot.raster <- 
function(dat,sampRate=512,main='',order=NULL,labs=T) 
{
	require(colorRamps)
	
	par(las=1)
	if(!is.null(order)) {
		dat = dat[order,]
	}
	image(1:ncol(dat),1:nrow(dat),t(dat),axes=F,col=matlab.like(64),main=main,xlab='time(ms)',ylab='')
	if(labs) axis(2,c(1,nrow(dat)/2,nrow(dat)),label=c(1,nrow(dat)/2,nrow(dat)))
	if(labs) prettyTime(1,c(1,ncol(dat)),Hz=sampRate) #axis(1,at=c(77,154,230,307),labels=c(150,300,450,600))
	
}


plotPeakMethods <-
function(solution) 
#plot the peak method function + ranges
{
	quartz(width=9,height=3)
	layout(matrix(1:3,1,3,byrow=T))
	avg = apply(.swale.solution.model(solution),2,mean)
	ml = setMaxLat(solution)
	
	control = .swale.solution.control(solution)
	
	plot(avg,type='l',lwd=2,lty=1,col=1,ylim=range(avg)*1.2,main='Peak selection (average)')
	st = .05
	for(i in 1:ncol(.swale.solution.waveform(solution))) {
		ranges = mean(.control.peak.windows(control)[[i]])+ml[[i]]
		pm = peakModel(matrix(avg,1,),.control.peak.windows(control)[[i]])
		arrows(.control.peak.windows(control)[[i]][1],max(avg)*(1+st),.control.peak.windows(control)[[i]][2],max(avg)*(1+st),angle=90,code=3,length=.05,lwd=2,col=i+1,lty=1)
		for(j in 1:length(pm[[1]]$lat)) {
			points(pm[[1]]$lat[j],pm[[1]]$amps[j],col=i+1,pch=19,cex=1.2)
		}
	
	}
	
	#plot(avg,type='l',lwd=2,lty=1,col=1,ylim=range(avg)*1.2)
	st = .05
	for(i in 1:ncol(.swale.solution.waveform(solution))) {
		tr = .swale.solution.waveform(solution)[,i]
		plot(tr,type='l',lwd=2,lty=1,col=1,ylim=range(tr)*1.2,main=paste('Min/max within range wave',i))
		ranges = mean(.control.peak.windows(control)[[i]])+ml[[i]]
		pm = peakModel(matrix(tr,1,),ranges)
		arrows(ranges[1],max(tr)*(1+st),ranges[2],max(tr)*(1+st),angle=90,code=3,length=.05,lwd=2,col=i+1,lty=2)
		#arrows(.control.peak.windows(control)[[i]][1],max(tr)*(1+st),.control.peak.windows(control)[[i]][2],max(tr)*(1+st),angle=90,code=3,length=.05,lwd=2,col=i+1,lty=1)
		
		for(j in 1:length(pm[[1]]$lat)) {
			points(pm[[1]]$lat[j],pm[[1]]$amps[j],col=i+1,pch=19,cex=1.2)
		}
		
		points(.swale.solution.latency(solution)[,i],rep(min(tr)*1.05,length(.swale.solution.latency(solution)[,i])),col=i+1)
			
	}
}

prettyTime <- 
function(axisnum=1,samprange=NULL,Hz=512,ticks=length(axTicks(axisnum)),every=NULL,scale=c('ms','s'),plot=T) 
{
	#make pretty division
	tsrange = range(samprange*(1/Hz)*1000)
		
	secs = tsrange[2]%/%1000
	dsec = tsrange[2]%/%100
	msec = tsrange[2]%/%10
	mmsec =tsrange[2]%/%1
		
	if(secs>0) {
		if(round(msec/10) > secs) secs = secs + 1
		if(is.null(every)) pticks=seq(tsrange[1]%/%1000,secs,length.out=ticks) else seq(tsrange[1]%/%1000,secs,every)
		ats = pticks/(1/Hz*1)
		labs = pticks*1000
		#cat('secs')
	} else if(dsec>0) {
		if(round(msec/10) > dsec) dsec = dsec + 1
		if(is.null(every)) pticks=seq(tsrange[1]%/%100,dsec,length.out=ticks) else pticks=seq(tsrange[1]%/%100,dsec,every)
		ats = pticks/(1/Hz*10)
		labs = pticks*100
		#cat('dsec')
	} else if(msec>0) {
		if(round(mmsec/10) > msec) msec =msec + 1
		if(is.null(every)) pticks=seq(tsrange[1]%/%10,msec,length.out=ticks) else pticks=seq(tsrange[1]%/%100,msec,every)
		ats = pticks/(1/Hz*100)
		labs = pticks*10
		#cat('msec')
	} else if(msec>0) {
		if(round((tsrange[2]%/%.1)/10) > msec) msec =msec + 1
		if(is.null(every)) pticks=seq(tsrange[1]%/%1,mmsec,length.out=ticks) else pticks=seq(tsrange[1]%/%100,mmsec,every)
		ats = pticks/(1/Hz*1000)
		labs = pticks*1
		#cat('mmsec')
	}
	
	if(plot) axis(axisnum,at=ats,label=labs)	
	
	return(invisible(list(ticks=pticks,at=ats,label=labs)))
		
}
