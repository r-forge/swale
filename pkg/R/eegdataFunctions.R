#############################################
# eegdata FUNCTIONS     				    #
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#CONTAINS
#cleanTrials
#peakPick

cleanTrials <-
function(swalesol,lat.thres=100,amp.thres=c(1e-3,1e2)) 
#clean op the single-trials based on latency and/or amplitude estimates
{

	if(length(grep('control',slotNames(new('swale.solution'))))!=0) swalesol = new('swale.solution',swalesol,discard=integer(0),control=new('control')) else swalesol = new('swale.solution',swalesol,discard=integer(0))
	
	amps = .swale.solution.amplitude(swalesol)
	lats = .swale.solution.latency(swalesol)
	
	rm.lat = rm.amp = numeric(0)
	
	for(i in 1:ncol(amps)) {
		rm.lat = c(rm.lat,which(abs(lats[,i])>lat.thres))
		rm.amp = c(rm.amp,which(abs(amps[,i])>amp.thres[2] | abs(amps[,i])<amp.thres[1]))		
	}
		
	discard = sort(unique(c(rm.lat,rm.amp)))
	discardvec = numeric(nrow(amps))
	discardvec[discard]=1
	
	.swale.solution.discard(swalesol) = discardvec
	
	return(swalesol)

}


peakPick <- 
function(object,window) 
#peak pick data or swale model
{

	if(class(object)=='swale.solution') {
		eegdat = .swale.solution.model(object)
		window = .swale.internal..swale.solution.internal(object)
		
	}
	if(class(object)=='eeg.data') eegdat = .eeg.data.data(object)
	
	peakmat = matrix(NA,nrow(eegdat),nrow(window))
	
	#latency PP
	for(trial in 1:nrow(eegdat)) {
		for(peak in 1:nrow(window)) {
			peakmat[trial,peak] = window[peak,1] + which.max(abs(eegdat[trial,window[peak,1]:window[peak,2]]))
		}
	}

	slopemat = matrix(NA,nrow(eegdat),nrow(window)-1)
	#caluclate slopes
	for(trial in 1:nrow(eegdat)) {
		for(slope in 2:nrow(window)) {
			slopemat[trial,slope-1] = (eegdat[trial,peakmat[trial,slope]]-eegdat[trial,peakmat[trial,slope-1]])/(peakmat[trial,slope]-peakmat[trial,slope-1])
		}
	}
	
	return(list(peaks=peakmat,slopes=slopemat))

}



