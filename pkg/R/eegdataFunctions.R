#############################################
# eegdata FUNCTIONS     				    #
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#CONTAINS
#peakPick

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



