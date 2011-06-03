#############################################
# eegdata FUNCTIONS     				    #
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#CONTAINS
#cleanTrials
#setMaxLat
#modelPeakPick

cleanTrials <-
function(swalesol,lat.thres=c(-100,100),amp.thres=c(1e-3,1e2)) 
#clean op the single-trials based on latency and/or amplitude estimates
{

	if(length(grep('control',slotNames(new('swale.solution'))))!=0) swalesol = new('swale.solution',swalesol,discard=integer(0),control=new('control')) else swalesol = new('swale.solution',swalesol,discard=integer(0))
	
	amps = .swale.solution.amplitude(swalesol)
	lats = .swale.solution.latency(swalesol)
	
	rm.lat = rm.amp = numeric(0)
	
	for(i in 1:ncol(amps)) {
		rm.lat = c(rm.lat,which(lats[,i]>lat.thres[2] | lats[,i]<lat.thres[1]))
		rm.amp = c(rm.amp,which(amps[,i]>amp.thres[2] | amps[,i]<amp.thres[1]))		
	}
	
	discard = sort(unique(c(rm.lat,rm.amp)))
	discardvec = numeric(nrow(amps))
	discardvec[discard]=1
	
	.swale.solution.discard(swalesol) = discardvec
	
	return(swalesol)

}

setMaxLat <-
function(swalesol,window=list(NULL),prec.fac=list(1.5),method=list('abs'),discarded=list(10),plot=F) 
#determine the maximum latency a waveform + derivative can model
{
	
	if(length(grep('control',slotNames(new('swale.solution'))))!=0) swalesol = new('swale.solution',swalesol,discard=integer(0),control=new('control')) else swalesol = new('swale.solution',swalesol,discard=integer(0))
	
	ranges = vector('list',ncol(.swale.solution.waveform(swalesol)))
	
	for(wave in 1:ncol(.swale.solution.waveform(swalesol)))  {
	
		method[[wave]] = match.arg(method[[wave]],c('abs','max','min','close','first','last'))
		nsamp = .eeg.data.samples(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
		
		if(is.null(window[[wave]])) window[[wave]] = c(1,nsamp) 
		
		if(method[[wave]]=='abs') absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='max') absmax = which.max((.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='min') absmax = which.min((.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='close') absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='first') absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='last') absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1

		vseq = seq(-nsamp,nsamp,1)
		maxes = numeric(length(vseq))
		
		p=1
		for(s in vseq) {
			trial = .swale.solution.waveform(swalesol)[,wave]+.swale.solution.derivwave(swalesol)[,wave]*s
			trial[1:discarded[[wave]]] = 0
			trial[(length(trial)-discarded[[wave]]):(length(trial))] = 0
			if(.swale.solution.waveform(swalesol)[absmax,wave]<0) maxes[p] = which.min(trial) else maxes[p] = which.max(trial)	
			p=p+1
		}
	
		latrange = ( range(maxes) - absmax ) * prec.fac[[wave]]
	
		if(plot) {
			quartz('Latency Ranges',width=8,height=4)
			layout(matrix(1:2,,2))
			plot(.swale.solution.waveform(swalesol)[,wave]+.swale.solution.derivwave(swalesol)[,wave]*latrange[1],type='l',lwd=2,bty='n',xlab='time',ylab='mV',main=latrange[1])
			plot(.swale.solution.waveform(swalesol)[,wave]+.swale.solution.derivwave(swalesol)[,wave]*latrange[2],type='l',lwd=2,bty='n',xlab='time',ylab='mV',main=latrange[2])
			
		}
	
		ranges[[wave]] = latrange
	}
	
	return(ranges)
}


summarizeModel <-
function(swale.solution,exp.peak=list(NULL),method=list('max'),maxlatfac=list(1.5),discarded=list(10),plot=F,output=T) 
#show and plot model summary (maxlat,latency+amplitude solution)
{
	numwave = ncol(.swale.solution.waveform(swale.solution))
	outpeak = vector('list',numwave)
	
	intern = .swale.solution.internal(swale.solution)
	
	if(numwave!=length(exp.peak)) stop('exp.peak must be a list of length numwaves.')
	
	ml = setMaxLat(swale.solution,window=exp.peak,prec.fac=maxlatfac,method=method,discarded=discarded,plot=F)
	
	for(wave in 1:numwave) {
		
		method[[wave]] = match.arg(method[[wave]],c('abs','max','min','close','first','last'))
		
		wavemod = intern
		.swale.internal.waves(wavemod) = matrix(.swale.internal.waves(intern)[,wave],,1)
		.swale.internal.amps(wavemod) = matrix(.swale.internal.amps(intern)[,wave],,1)
		.swale.internal.lats(wavemod) = matrix(.swale.internal.lats(intern)[,wave],,1)
		
		wavemodel = model(wavemod)
		peaks = peakModel(.swale.internal.model(wavemodel),window=mean(exp.peak[[wave]])+ml[[wave]],plot=plot)
		
		if(output) {
			cat('--------------------------------------------\n')
			cat('        model summary wave [',wave,']\n')
			cat('--------------------------------------------\n')
			cat('  Expected peak window  :',exp.peak[[wave]],'\n')
			cat('                latency :',round(mean(exp.peak[[wave]])),'\n')
			cat('                range   :',round(mean(exp.peak[[wave]])+ml[[wave]]),'\n')
			cat('--------------------------------------------\n')
			cat('  Selection method      :',method[[wave]],'\n')
			cat('--------------------------------------------\n')
		}
				
		selpeak = list(amps=numeric(length(peaks)),lats=numeric(length(peaks)))
		
		for(trials in 1:length(peaks)) {
			
			if(method[[wave]]=='abs') peak = which.max(abs(peaks[[trials]]$amps))
			if(method[[wave]]=='max') {
				peak = which.max(peaks[[trials]]$amps)
				if(length(peak)!=0) if(peaks[[trials]]$amps[peak]<0) peak = numeric(0) 	
			}
			if(method[[wave]]=='min') {
				peak = which.min(peaks[[trials]]$amps)
				if(length(peak)!=0) if(peaks[[trials]]$amps[peak]>0) peak = numeric(0) 	
			}
			if(method[[wave]]=='close') peak = which.min(abs(peaks[[trials]]$lat-mean(exp.peak[[wave]])))
			if(method[[wave]]=='first') peak = 1
			if(method[[wave]]=='last') peak = length(peaks[[trials]]$amps)

			if(length(peak)!=0) {
				selpeak$amps[trials] = peaks[[trials]]$amps[peak]
				selpeak$lats[trials] = peaks[[trials]]$lat[peak]
			} else {
				selpeak$amps[trials] = NA
				selpeak$lats[trials] = NA
			}
			
			
			if(output) {
				cat(sprintf(' [%3.0d] >> %6.2f @ %4.0f \n',trials,as.numeric(selpeak$amps[trials]),as.numeric(selpeak$lats[trials])))
			}
			
		} #trial loop
	
		outpeak[[wave]] = selpeak
		if(output) cat('--------------------------------------------\n')
	} #waveloop
	
	return(invisible(outpeak))
}

