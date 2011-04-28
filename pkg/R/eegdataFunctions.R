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
function(data,window=NULL,plot=F,simdat=NULL,simeeg=NULL) 
#peakPick using derivative
{
	
	trials = nrow(data)
	
	peakvec = vector('list',trials)
	
	for(i in 1:trials) {
		
		deriv = fdDeriv(matrix(data[i,],,1))
		peaks = zerocross(deriv)
		
		if(!is.null(window)) {
			rm = c(which(peaks<window[1]),which(peaks>window[2]))
			if(length(rm)>0) peaks = peaks[-rm]
		}
		
		if(plot) {
			par(ask=T)
			
			if(!is.null(simeeg)) {
				plot(simeeg[i,],col='gray',type='l',lty=1,bty='n',main='peaks')
				lines(data[i,],lwd=2,col=1)
				rsq = (sum(simeeg[i,]^2)-sum((simeeg[i,]-data[i,])^2))/sum(simeeg[i,]^2)
				text(200,5,round(rsq,6))
				
			} else plot(data[i,],type='l',lwd=2,col=1,bty='n',main='peaks')
			
			points(peaks,data[i,peaks],pch=19,col=2)
			points(window,data[i,window],pch=19,col=4)
			
			if(!is.null(simdat)) {
				points(150-simdat$lats[i,1],data[i,round(150-simdat$lats[i,1])],col=3,pch=21,cex=1.5)
				
			}
		}
	
		peakvec[[i]] = list(lat=peaks,amps=data[i,peaks])
		
	}
	
	return(peakvec)	
}


zerocross <- function(vec) 
{
	len = length(vec)
	sig = sign(vec)
	zc = numeric(0)
	
	for(i in 1:(len-1)) {
		
		if(sig[i]==1 & sig[i+1]==-1) zc = c(zc,(i-1+which.min(abs(c(sig[i],sig[i+1])))))
		if(sig[i]==-1 & sig[i+1]==1) zc = c(zc,(i-1+which.min(abs(c(sig[i],sig[i+1])))))
	
	}
	
	return(zc)
	
}


