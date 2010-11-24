#############################################
# SIMTOOLS FUNCTIONS		 				#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#CONTAINS
#simulateEEGsignal
#makeNoise
#calculateSNR
#estimateSNR
#add discards

simulateEEGsignal <-
function(signal=c(150,-200),peakwidth=15,mVscale=500,trials=40,samples=350,sample.freq=512,amp.dist='norm',amp.sd=.5,amp.mean=1,amp.range=c(-2,2),lat.dist='norm',lat.sd=40,lat.mean=0,lat.range=c(-100,100),snr=10,noisemethod=c('white'),linkmethod='none',plot=TRUE)
# simulate EEG datasets
{
	funcall = list(signal=signal,peakwidth=peakwidth,mVscale=mVscale,trials=trials,samples=samples,sample.freq=sample.freq,amp.dist=amp.dist,amp.sd=amp.sd,amp.mean=amp.mean,amp.range=amp.range,lat.dist=lat.dist,lat.sd=lat.sd,lat.mean=lat.mean,lat.range=lat.range,snr=snr,noisemethod=noisemethod,linkmethod=linkmethod,plot=plot)
	
	#in HZ
	ms = samples*(1/sample.freq)*1000
	
	#show simulated information
	if(plot) {
		cat('simulate EEG data ',date(),'\n')
		cat(' num peaks : ',length(signal),' w(',peakwidth,')','\n',sep='')
		cat(' located at:',abs(signal),'\n')
		cat(' trials    : ',trials,'\n',sep='')
		cat(' samples   : ',samples,' (',round(ms),' ms)','\n',sep='')
		cat(' amplitude : ',amp.dist,'(',amp.mean,',',amp.sd,')',' r[',amp.range[1],',',amp.range[2],']\n',sep='')
		cat(' latency   : ',lat.dist,'(',lat.mean,',',lat.sd,')',' r[',lat.range[1],',',lat.range[2],']\n',sep='')
		cat(' latencySD : ',lat.sd*(1/sample.freq)*1000,' ms\n',sep='')
		cat(' input SNR : ',snr,'\n',sep='')
		cat(' noise meth:',noisemethod,'\n')
	}
	
	#create amplitude matrix (fill and flip)
	nwave = length(signal)
	amps = matrix(NA,trials,nwave)
	for(wave in 1:nwave) {
		if(amp.dist=='norm') amps[,wave] = rnorm(trials,amp.mean,amp.sd)
		amps[,wave][amps[,wave]>amp.range[2]]=amp.range[2]
		amps[,wave][amps[,wave]<amp.range[1]]=amp.range[1]
	}
	
	#create latency matrix
	lats = matrix(NA,trials,nwave)
	for(wave in 1:nwave) {
		if(lat.dist=='norm') lats[,wave] = rnorm(trials,lat.mean,lat.sd)
		lats[,wave][lats[,wave]>lat.range[2]]=lat.range[2]
		lats[,wave][lats[,wave]<lat.range[1]]=lat.range[1]
	}
	
	#link amps and lats (reference is first wave)
	if(linkmethod=='all') {
		if(ncol(amps)>1) {
			for(wave in 2:ncol(amps)) {
				amps[,wave]=amps[,1]
				lats[,wave]=lats[,1]
			}
		}
	}
	
	#create and fill datamatrix
	data = matrix(NA,trials,samples)
	for(i in 1:nrow(data)) {
		f = numeric(samples)
		for(j in 1:nwave) {
			if(signal[j]<0) neg=-1 else neg=1
			f = f + dnorm(1:samples,(abs(signal[j])-lats[i,j]),peakwidth)*amps[i,j]*mVscale*neg
		}
		data[i,] = f
	}
	
	#make noise
	noise = makeNoise(data,snr,noisemethod)
	
	#add noise
	dataplusnoise = data + noise
	
	#check SNRS
	snrinfo = list(input.snr=snr,calc.h=calculateSNR(data,noise,'H'),est.h=estimateSNR(dataplusnoise,'H'),calc.ft=calculateSNR(data,noise,'FT'),est.ft=estimateSNR(dataplusnoise,'FT'))
	if(plot) {
		cat(' *** SNR info ****\n')
		cat(' Calculated SNR (max avg/sdnoise) = ',snrinfo$calc.h,'\n')
		cat(' Estimated SNR (Handy method)     = ',snrinfo$est.h,'\n')
		
		cat(' Calculated SNR (ss/sn)           = ',snrinfo$calc.ft,'\n')
		cat(' Estimated SNR (Fein/Turetsky)    = ',snrinfo$est.ft,'\n')
		cat(' ***\n')
	}
	
	#return
	return(list(data=dataplusnoise,signal=data,noise=noise,amps=amps,lats=lats,snr=snrinfo,call=funcall))
	
}


makeNoise <-
function(data,snr,noisemethod)
{
	
	#create and fill noise matrix
	noise = matrix(NA,nrow(data),ncol(data))
	#trial_sdnoise=apply(data,1,function(x) {max((x))/(snr)})
	trial_sdnoise=rep(max(apply(data,2,mean))/snr,nrow(data))
	
	if(snr==Inf) trial_sdnoise = rep(0,nrow(data))
	
	for(i in 1:nrow(noise)) 	noise[i,] = rnorm(ncol(data),0,trial_sdnoise[i])
	
	#make correlated noise
	if(noisemethod[1]=='pink') {
		noise = lowpass(noise,filsd=as.numeric(noisemethod[2]),fillen=100,filmean=50)
	}
	
	return(noise)
	
}



calculateSNR <- function(data,noise,meth='FT') 
#calculate SNR based on data and noise matrices
{
	meth = match.arg(meth,c('FT','H'))
	
	#Fein & Turetsky
	if(meth=='FT') {
		#noise power
		sn=0
		for(j in 1:nrow(data)) sn = sn + sum(noise[j,])^2
		sn = sn * (1/(ncol(data)*(nrow(data)-1)))
		
		#signal power
		avg <- apply(data,2,mean)
		ss = (1/(ncol(data)))*sum(avg^2)
		
		#signal-to-noise ratio
		snr = ss/sn
		
		return(snr)
	}
	
	#Handy
	if(meth=='H') {
		avg <- apply(data,2,mean)
		maxavg <- max(abs(avg))
		
		#noise SD (signal-avg)
		sdnoise = sd(as.vector(noise))
		
		#signal to noise ratio
		snr = maxavg/sdnoise 
		
		return(snr)
	}
	
	return(NULL)
}


plotAmpLat <- function(dat,sol,plot=T)
#plot Amplitude and Latency data and estimates
{
	
	disc = sol@discard
	
	quartz(title='Amplitude')
	layout(1);par(las=1)
	plot(NA,NA,xlim=range(dat$amps),ylim=range(sol@amplitude,na.rm=T),xlab='Data',ylab='Estimate',axes=T,bty='n')
	abline(0,1)
	if(length(which(disc>0))>0) {
		points(dat$amps[which(disc==0)],sol@amplitude[which(disc==0)],pch=19,col=1)
		points(dat$amps[which(disc==1)],sol@amplitude[which(disc==1)],pch='X',col=gray(.5))	
		ampcor = cor.test(dat$amps[which(disc==0)],sol@amplitude[which(disc==0)])
	} else {
		points(dat$amps,sol@amplitude,pch=19,col=1)
		ampcor = cor.test(dat$amps,sol@amplitude)
	}
	
	quartz(title='Latency')
	layout(1);par(las=1)
	plot(NA,NA,xlim=range(dat$lats),ylim=range(sol@latency,na.rm=T),xlab='Data',ylab='Estimate',axes=T,bty='n')
	abline(0,1)
	if(length(which(disc>0))>0) {
		points(dat$lats[which(disc==0)],sol@latency[which(disc==0)],pch=19,col=1)
		points(dat$lats[which(disc==1)],sol@latency[which(disc==1)],pch='X',col=gray(.5))
		latcor = cor.test(dat$lats[which(disc==0)],sol@latency[which(disc==0)])
	} else {
		points(dat$lats,sol@latency,pch=19,col=1)
		latcor = cor.test(dat$lats,sol@latency)
	}
	
	return(invisible(list(ampcor=ampcor,latcor=latcor)))
	
}


