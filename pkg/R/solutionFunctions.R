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
function(swalesol,window=NULL,prec.fac=1.5,plot=F) 
#determine the maximum latency a waveform + derivative can model
{
	
	if(length(grep('control',slotNames(new('swale.solution'))))!=0) swalesol = new('swale.solution',swalesol,discard=integer(0),control=new('control')) else swalesol = new('swale.solution',swalesol,discard=integer(0))
	
	
	nsamp = .eeg.data.samples(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
	
	if(is.null(window)) window = c(1,nsamp) 
	absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[1]:window[2]]))+window[1]-1
	vseq = seq(-nsamp,nsamp,1)
	maxes = numeric(length(vseq))
	
	p=1
	for(s in vseq) {
		trial = .swale.solution.waveform(swalesol)+.swale.solution.derivwave(swalesol)*s
		if(.swale.solution.waveform(swalesol)[absmax]<0) maxes[p] = which.min(trial) else maxes[p] = which.max(trial)	
		p=p+1
	}

	latrange = ( range(maxes) - absmax ) * prec.fac

	if(plot) {
		quartz('Latency Ranges',width=8,height=4)
		layout(matrix(1:2,,2))
		plot(.swale.solution.waveform(swalesol)+.swale.solution.derivwave(swalesol)*latrange[1],type='l',lwd=2,bty='n',xlab='time',ylab='mV',main=latrange[1])
		plot(.swale.solution.waveform(swalesol)+.swale.solution.derivwave(swalesol)*latrange[2],type='l',lwd=2,bty='n',xlab='time',ylab='mV',main=latrange[2])
		
	}
	
	
	return(latrange)
}


summarizeModel <-
function(solution1,solution2) 
{
	cat('AIC 1: ',.swale.solution.aic(solution1),'\n')
	cat('AIC 2: ',.swale.solution.aic(solution2),'\n')
	if(.swale.solution.aic(solution2)<.swale.solution.aic(solution1)) {
		cat('AIC indicates multiple waveforms\n')
		
		cors = swalecor(solution2)
		
		if(cors$amplitude$p.value[1,2]<.05) cat('Amplitude fixed.\n') else cat('Amplitude free\n')
		if(cors$latency$p.value[1,2]<.05) cat('Latency fixed.\n') else cat('Latency free\n')
		
		if(cors$amplitude$p.value[1,2]<.05 | cors$latency$p.value[1,2]<.05) {
			cat('Treat as one waveform model!\n')	
			return(1)
		} else { return(2) }
	
	
	} else {
		cat('AIC indicates one waveform\n')
		return(1)
	}



}

