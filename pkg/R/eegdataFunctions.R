#############################################
# eegdata FUNCTIONS     				    #
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#CONTAINS
#peakPick

peakPick <- 
function(eegdat,window) 
#peak pick data or swale model
{

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

peakModel <-
function(data,window=NULL,plot=F) 
#peakPick using derivative of single-trial models
{
	
	trials = nrow(data)
	peakvec = vector('list',trials)
	
	if(plot) par(ask=T)
	
	for(i in 1:trials) {
		
		deriv = fdDeriv(matrix(data[i,],,1))
		peaks = zerocross(deriv)
		minmax = attr(peaks,'minmax')
		
		if(!is.null(window)) {
			rm = c(which(peaks<window[1]),which(peaks>window[2]))
			if(length(rm)>0) {
				peaks = peaks[-rm]
				minmax = minmax[-rm]
				attr(peaks,'minmax') = minmax
			}
		}
		if(plot) {
			plot(data[i,],lwd=2,col=1,type='l')
			points(peaks,data[i,peaks],col=2,lwd=2)
			mm = which(attr(peaks,'minmax')==1)
			mm = mm[which.max(abs(data[i,peaks[mm]]))]			
			if(length(mm)>0) text(peaks[mm],data[i,peaks[mm]],'MAX')
			
			mm = which(attr(peaks,'minmax')==-1)
			mm = mm[which.max(abs(data[i,peaks[mm]]))]
			if(length(mm)>0) text(peaks[mm],data[i,peaks[mm]],'MIN')
		}
		
		peakvec[[i]] = list(lat=peaks,amps=data[i,peaks])
		
	}
	
	return(peakvec)	
}

zerocross <- 
function(vec) 
#define zero-crossings
{
	len = length(vec)
	sig = sign(vec)
	zc = mm = numeric(0)
	
	for(i in 1:(len-1)) {
		
		if(sig[i]==1 & sig[i+1]==-1) {
			zc = c(zc,(i-1+which.min(abs(c(sig[i],sig[i+1])))))
			mm = c(mm,1)
		}
		if(sig[i]==-1 & sig[i+1]==1) {
			zc = c(zc,(i-1+which.min(abs(c(sig[i],sig[i+1])))))
			mm = c(mm,-1)	
		}
	
	}
	
	attr(zc,'minmax') = mm
	#browser()
	return(zc)
	
}


