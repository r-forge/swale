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
function(signal=c(125,-225),sigmeth='shift',peakwidth=15,mVscale=500,trials=100,samples=350,sample.freq=512,amp.dist='norm',amp.sd=0.5,amp.mean=1,amp.range=c(-2,2),lat.dist='norm',lat.sd=40,lat.mean=0,lat.range=c(-100,100),snr=10,noisemethod=c('white'),linkmethod='none',plot=TRUE,poly=20)
# simulate EEG datasets
{
	funcall = list(signal=signal,sigmeth=sigmeth,peakwidth=peakwidth,mVscale=mVscale,trials=trials,samples=samples,sample.freq=sample.freq,amp.dist=amp.dist,amp.sd=amp.sd,amp.mean=amp.mean,amp.range=amp.range,lat.dist=lat.dist,lat.sd=lat.sd,lat.mean=lat.mean,lat.range=lat.range,snr=snr,noisemethod=noisemethod,linkmethod=linkmethod,plot=plot,poly=poly)
	
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
		if(amp.dist=='lnorm') amps[,wave] = rlnorm(trials,amp.mean-1,amp.sd/5) 
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
	
	#rescale amps to mean 1
	for(wave in 1:nwave) {
		amps[,wave] = amps[,wave]/mean(amps[,wave])
	}

	#create and fill datamatrix
	data = matrix(NA,trials,samples)
	for(i in 1:nrow(data)) {
		f = numeric(samples)
		
		for(j in 1:nwave) {
			if(signal[j]<0) neg=-1 else neg=1
			if(sigmeth=='shift') f = f + dnorm(1:samples,(abs(signal[j])-lats[i,j]),peakwidth)*amps[i,j]*mVscale*neg
			if(sigmeth=='model') {
				if(linkmethod=='none') {
					
					#make basis object
					basis = makePoly(poly,samples,'FD')
					s = matrix(dnorm(1:samples,(abs(signal[j])),peakwidth),,1)*mVscale*neg
					
					fb = (t(basis@matrix)%*%s)
					fs = basis@matrix%*%fb
					ds = basis@deriv%*%fb
					f = f + amps[i,j]*(fs+ds*lats[i,j])
					
				}
			}
		}
		
		if(sigmeth=='model') {
			if(linkmethod=='all') {
				
				#make basis object
				basis = makePoly(poly,samples,'FD')
				s=numeric(length(f))
				for(j in 1:nwave) {
					if(signal[j]<0) neg=-1 else neg=1
					s = s + matrix(dnorm(1:samples,(abs(signal[j])),peakwidth),,1)*mVscale*neg
				}
				
				fb = (t(basis@matrix)%*%s)
				fs = basis@matrix%*%fb
				ds = basis@deriv%*%fb
				f = f + amps[i,1]*(fs+ds*lats[i,1])
			}
		}
		
		data[i,] = f
		
	}
	
	#make noise
	noise = makeNoise(data,snr,noisemethod)
	
	#add noise
	if(noisemethod[1]=='ar') {
		dataplusnoise = data
		mbet = max(apply(data,2,mean))
		mnoise = mean(apply(noise,1,sd))
		for(i in 1:nrow(data)) {
			data[i,]=data[i,]/mbet
			data[i,]=data[i,]*mnoise*snr	
		}
	}	
	dataplusnoise = data + noise
	
	#check SNRS
	nsd = dsd = numeric(nrow(data))
	for(i in 1:nrow(data)) {
		rng = which(abs(data[i,])>1e-3)
		dsd[i] = sd(data[i,rng])
		nsd[i] = sd(noise[i,rng])
	}	
	snrQQ = (dsd)/(nsd)
	snrQQdB = (20*log(dsd/nsd,base=10))

	snrinfo = list(input.snr=snr,calc.h=calculateSNR(data,noise,'H'),est.h=estimateSNR(dataplusnoise,'H'),calc.ft=calculateSNR(data,noise,'FT'),est.ft=estimateSNR(dataplusnoise,'FT'),snrQQ=snrQQ,snrQQdB=snrQQdB)
	
	
	if(plot) {
		cat(' *** SNR info ****\n')
		cat(' Calculated SNR (max avg/sdnoise) = ',snrinfo$calc.h,'\n')
		#cat(' Estimated SNR (Handy method)     = ',snrinfo$est.h,'\n')
		
		cat(' Calculated SNR (FT: ss/sn)       = ',snrinfo$calc.ft,'\n')
		#cat(' Estimated SNR (Fein/Turetsky)    = ',snrinfo$est.ft,'\n')
		
		#cat(' Calculated SNR (QQ)              = ',mean(snrinfo$snrQQ),'\n')
		#cat(' Calculated SNR 20log10(QQ)       = ',mean(snrinfo$snrQQdB),'\n')
		
		snrHeaven(data,noise)
		
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
	trial_sdnoise=rep(max(apply(data,2,mean))/snr,nrow(data))
	
	if(snr==Inf) trial_sdnoise = rep(0,nrow(data))
	
	if(noisemethod[1]=='white') {
		for(i in 1:nrow(noise)) 	noise[i,] = rnorm(ncol(data),0,trial_sdnoise[i])
	}
	
	#make correlated noise
	if(noisemethod[1]=='lp') {
		noise = lowpass(noise,filsd=as.numeric(noisemethod[2]),fillen=100,filmean=50)
	}
	
	#make correlated noise
	if(noisemethod[1]=='ar') {
		order = attr(noisemethod,'order')
		ar = attr(noisemethod,'ar')
		for(i in 1:nrow(noise)) noise[i,] = arima.sim(list(order=c(order,0,0),ar=ar),n=ncol(data))
		
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
	
	plot(NA,NA,xlim=range(dat$amps),ylim=range(sol@amplitude,na.rm=T),xlab='Data',ylab='Estimate',axes=T,bty='n',main='Amplitude')
	abline(0,1)
	if(length(which(disc>0))>0) {
		points(dat$amps[which(disc==0),1],sol@amplitude[which(disc==0)],pch=19,col=1)
		points(dat$amps[which(disc==1),1],sol@amplitude[which(disc==1)],pch='X',col=gray(.5))	
		ampcor = cor.test(dat$amps[which(disc==0),1],sol@amplitude[which(disc==0)])
	} else {
		points(dat$amps[,1],sol@amplitude,pch=19,col=1)
		ampcor = cor.test(dat$amps[,1],sol@amplitude)
	}
	
	plot(NA,NA,xlim=range(dat$lats),ylim=range(sol@latency,na.rm=T),xlab='Data',ylab='Estimate',axes=T,bty='n',main='Latency')
	abline(0,1)
	if(length(which(disc>0))>0) {
		points(dat$lats[which(disc==0),1],sol@latency[which(disc==0)],pch=19,col=1)
		points(dat$lats[which(disc==1),1],sol@latency[which(disc==1)],pch='X',col=gray(.5))
		latcor = cor.test(dat$lats[which(disc==0),1],sol@latency[which(disc==0)])
	} else {
		points(dat$lats[,1],sol@latency,pch=19,col=1)
		latcor = cor.test(dat$lats[,1],sol@latency)
	}
	
	return(invisible(list(ampcor=ampcor,latcor=latcor)))
	
}


snrHeaven <- function(data,noise=NULL)
{
		
	if(is.null(noise)) {
		cat('Estimated noise\n')
		noise = data
		avg = apply(data,2,mean)
		for(i in 1:nrow(data)) noise[i,] = data[i,]-avg
	} else cat('Real noise\n') 
	
	
	rng = apply(data,1,function(x) which(abs(x)>1e-3))
	
	vardata = varnoise = numeric(nrow(data))
	for(i in 1:nrow(data)) {
		vardata[i] = var(data[i,rng[[i]]])
		varnoise[i] = var(noise[i,rng[[i]]]) 
	}
	
	meanvardata = mean(vardata)
	
	iyer = mean(meanvardata/varnoise)
	ivan = 10*log(meanvardata/mean(varnoise),base=10)	
	truc = mean( 20*log(sqrt(vardata)/sqrt(varnoise),base=10) )	
	quia = mean( sqrt(vardata)/sqrt(varnoise) )
			
	cat('iyer [.2    1] =', iyer, '\n')
	cat('ivan [-27 -18] =', ivan, 'dB\n')
	cat('truc [-34  27] =', truc, 'dB\n')
	cat('quia [.5    3] =', quia, '\n')
	
	return(invisible(NULL))
}

swaleEEG <-
		function(DATA,which=1,channel='Fz',npoly=16,deriv='FD',control=new('control')) 
{
	
	data = new('eeg.data')
	
	if(is.numeric(which)) {
		
		if(class(DATA[[which]])!='bvadata') stop('Input must be of class BVADATA (see EEG package)')
		
		.eeg.data.condition(data) = names(DATA)[which]
		wc = which(dimnames(DATA[[which]]$x())[[2]]==channel)
		.eeg.data.channel(data) = dimnames(DATA[[which]]$x())[[2]][wc]
		.eeg.data.data(data) = t(DATA[[which]]$x()[,wc,])
		.eeg.data.sampRate(data) = DATA[[which]]$header$`Common Infos`$SamplingInterval ## NOT WORKING!!!!!!!!!!!!!!
		
	} else {
		if(class(DATA[[match(which,names(DATA))]])!='bvadata') stop('Input must be of class BVADATA (see EEG package)')
		
		.eeg.data.condition(data) = which
		wc = which(dimnames(DATA[[match(which,names(DATA))]]$x())[[2]]==channel)
		.eeg.data.channel(data) = dimnames(DATA[[match(which,names(DATA))]]$x())[[2]][wc]
		.eeg.data.data(data) = t(DATA[[match(which,names(DATA))]]$x()[,wc,])
		.eeg.data.sampRate(data) = DATA[[match(which,names(DATA))]]$header$`Common Infos`$SamplingInterval ## NOT WORKING!!!!!!!!!!!!!!
	}
	
	.eeg.data.trials(data) = nrow(.eeg.data.data(data))
	.eeg.data.samples(data) = ncol(.eeg.data.data(data))
	
	cat(.eeg.data.condition(data),'\n')
	
	.eeg.data.data(data) = detrend(.eeg.data.data(data))
	
	basis = makePoly(npoly,.eeg.data.samples(data),'FD')
	swaledat = new('swale.internal',eeg.data=data,basis=basis)
	
	
	.control.start.value(control) = makeStart(swaledat)
	solution = iterate(swaledat,control)
	solution = calcSolution(solution)
	
	return(solution)
	
}

plotDatEst <- function(estlat,datlat,model,data)
{
	par(las=1,ask=T)
	for(i in 1:nrow(data)) {
		plot(NA,NA,ylim=range(c(as.vector(model),as.vector(data))),xlim=c(1,ncol(data)),bty='n',xlab='time',ylab='mV',main=paste('trial [',i,']'))	
		lines(data[i,],col='darkgray',lwd=2,lty=1)
		lines(model[i,],col='black',lwd=2,lty=1)
		points(datlat[i],model[i,datlat[i]],col=2,pch=21,bg=2,cex=1)
		points(estlat[i],model[i,estlat[i]],col=4,pch=21,lwd=2)
	}
	
}




