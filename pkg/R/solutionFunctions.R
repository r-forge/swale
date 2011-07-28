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
function(swalesol,plot=F) 
#determine the maximum latency a waveform + derivative can model ** NOW FIXED RANGE
{
	
	#if(length(grep('control',slotNames(new('swale.solution'))))!=0) swalesol = new('swale.solution',swalesol,discard=integer(0),control=new('control')) else swalesol = new('swale.solution',swalesol,discard=integer(0))
	
	control = .swale.solution.control(swalesol)
	ranges = vector('list',ncol(.swale.solution.waveform(swalesol)))
	window = .control.peak.windows(control)

	if(length(.control.max.lat(control)[[1]])>1) {
		ranges = .control.max.lat(control)
	} else {
		for(i in 1:length(window)) {
			rng = (window[[i]][1]*.control.max.lat(control)[[i]]-window[[i]][1])
			abs = c((rng)*-1,(rng))
			ranges[[i]] = abs
		}
	}
		
	return(ranges)
}


summarizeModel <- 
function(swale.solution,plot=F,output=T) 
#show and plot model summary (maxlat,latency+amplitude solution)
{
	numwave = ncol(.swale.solution.waveform(swale.solution))
	intern = .swale.solution.internal(swale.solution)
	control = .swale.solution.control(swale.solution)
	
	exp.peak = .control.peak.windows(control)
	method = .control.peak.method(control)
	maxlatfac = .control.max.lat(control)
	
	outpeak = vector('list',numwave)
	
	if(numwave!=length(exp.peak)) stop('exp.peak must be a list of length numwaves.')

	ml = setMaxLat(swale.solution,plot=F)
	
	for(wave in 1:numwave) {
		
		method[[wave]] = match.arg(method[[wave]],c('abs','max','min','close','firstmax','firstmin','lastmax','lastmin'))
		
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
				minmax = which(attr(peaks[[trials]]$lat,'minmax')==1)
				peak = minmax[which.max(abs(peaks[[trials]]$amps[minmax]))]
			}
			if(method[[wave]]=='min') {
				minmax = which(attr(peaks[[trials]]$lat,'minmax')==-1)
				peak = minmax[which.max(abs(peaks[[trials]]$amps[minmax]))] 	
			}
			if(method[[wave]]=='close') peak = which.min(abs(peaks[[trials]]$lat-mean(exp.peak[[wave]])))
			if(method[[wave]]=='firstmax') {
				minmax = which(attr(peaks[[trials]]$lat,'minmax')==1)
				peak = minmax[1]
			}			
			if(method[[wave]]=='firstmin') {
				minmax = which(attr(peaks[[trials]]$lat,'minmax')==-1)
				peak = minmax[1]
			}
			if(method[[wave]]=='lastmax') {
				minmax = which(attr(peaks[[trials]]$lat,'minmax')==1)
				peak = minmax[length(minmax)]
			}
			if(method[[wave]]=='lastmin') {
				minmax = which(attr(peaks[[trials]]$lat,'minmax')==-1)
				peak = minmax[length(minmax)]
			}
			
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

