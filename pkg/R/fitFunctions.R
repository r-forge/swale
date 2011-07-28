#############################################
# FIT function FUNCTIONS				    #
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#CONTAINS
#iterate
#makeStart
#calcSolution
#iterateDiscard
#iterateSplit
#autoSplit

iterate <-
function(swaledat,control=new('control'),posGradStop=F) 
#iterate a f/ab estimation (with or without split)
{
	
	
	#check if startingvalues are sane
	if(nrow(.control.start.value(control))!=ncol(.basis.matrix(.swale.internal.basis(swaledat)))) stop('Length of startingvalues do not match with transform matrix.')

	#set startingvalues
	.swale.internal.waves(swaledat) = .control.start.value(control)
	.swale.internal.model(swaledat) = matrix(0,.eeg.data.trials(.swale.internal.eeg.data(swaledat)),.eeg.data.samples(.swale.internal.eeg.data(swaledat)))
	
	#do iterations
	exitIterate = F
	iterNum = 0
	gradient = NA
	objective = numeric(0)	
	
	#set old to swaledat
	swaledat_old = swaledat
	
	#plot iteration info
	cat('Intial rss [',round(.swale.internal.rss(rss(swaledat))),'], iterations started at: ',date(),' Running...',sep='')
		
	while(!exitIterate) {
		
		iterNum = iterNum + 1
		
		#order if one waveform is estimated
		if(.control.split.type(control)=='none') {
			swaledat_new = estimate.ab(swaledat_old)
			swaledat_new = estimate.f(swaledat_new)
			swaledat_new = model(swaledat_new)
			swaledat_new = rss(swaledat_new)
		}
		
		
		#calculate objective and gradient
		objective = c(objective,.swale.internal.rss(swaledat_new))
		if(iterNum>1) {
			gradient = c(gradient,.swale.internal.rss(swaledat_new)-objective[iterNum-1])
		} else {
			if(.control.output(control)) cat('\n')
		}

		#show gradient and objective information
		if(.control.output(control)) cat(sprintf('%3d: %10.0f ~ (%16.9f)',iterNum,objective[iterNum],gradient[iterNum]),'\n')
		
		#check if gradient is smaller than convergence or maxIter is reached
		if(iterNum>2) if(abs(gradient[iterNum])<.control.iter.convergence(control)) exitIterate=TRUE
		if(posGradStop) if(iterNum>1) if(gradient[iterNum]>0) exitIterate=TRUE
		if(iterNum==.control.iter.limit(control)) exitIterate=TRUE
		
		if(exitIterate==FALSE) swaledat_old = swaledat_new
		
	}
	
	#set swale internal objects 
	.swale.internal.rss(swaledat_old) = objective[iterNum-1]
	.swale.internal.gradient(swaledat_old) = gradient[iterNum-1]
	
	solution = new('swale.solution',internal=swaledat_old,control=control)
	cat('done\n')	
	cat(' final rss [',round(.swale.internal.rss(swaledat_old)),']\n',sep='')
	cat('       aic [',round(aic(swaledat_old)),']\n',sep='')
	
	return(solution)

}

makeStart <-
function(swaledat,method='mean') 
#make startingvalues
{
	if(class(swaledat)!='swale.internal') stop('Input must be of class \'swale.internal\'')
	method = match.arg(method,'mean')
	
	if(method=='mean') start = t(.basis.matrix(.swale.internal.basis(swaledat)))%*%apply(.eeg.data.data(.swale.internal.eeg.data(swaledat)),2,mean)
		
	return(start)
}


calcSolution <-
function(swaledat,summarize = TRUE) 
#calculate solution of a swale.solution object
{
	if(class(swaledat)!='swale.solution') stop('Input must be of class \'swale.solution\'')
	
	solution = swaledat
	swaledat = .swale.solution.internal(swaledat)
	
	.swale.solution.amplitude(solution) = .swale.internal.amps(swaledat)
	.swale.solution.latency(solution) = .swale.internal.lats(swaledat)/.swale.internal.amps(swaledat)
	.swale.solution.model(solution) = .swale.internal.model(swaledat)
	.swale.solution.waveform(solution) = .basis.matrix(.swale.internal.basis(swaledat))%*%.swale.internal.waves(swaledat)
	.swale.solution.derivwave(solution) = .basis.deriv(.swale.internal.basis(swaledat))%*%.swale.internal.waves(swaledat)
	.swale.solution.aic(solution) = aic(swaledat)
	
	if(summarize) {
		#calc amplitude and latency
		control = .swale.solution.control(solution)
		peaks = .control.peak.windows(control)
		nsplit = length(peaks)
		
		sm = summarizeModel(solution,plot=F,output=F)
		
		for(i in 1:nsplit) {
			.swale.solution.amplitude(solution)[,i] = sm[[i]]$amps
			.swale.solution.latency(solution)[,i] = sm[[i]]$lats
		}
	}
	
	.swale.solution.discard(solution) = numeric(nrow(.swale.solution.amplitude(solution))) 
	.swale.solution.discard(solution)[unique(which(is.na(cbind(.swale.solution.amplitude(solution),.swale.solution.latency(solution))),arr=T)[,1])] = 1
	
	return(solution)
}




multiSplit <-
function(swalesol,control=new('control'),latRange=NULL,ampRange=c(1e-06,Inf),posGradStop=F,stepsize=1,fixstep=NULL)
{
	
	swaledat = .swale.solution.internal(swalesol)
	window = .control.peak.windows(control)
	
	#set estimation objects
	TP = .basis.matrix(.swale.internal.basis(swaledat))
	f = .swale.internal.waves(swaledat)
	ft = TP%*%f
	
	#search for maxima within window
	nc = length(window)
	if(nc>2) {warning('multiSplit not available for more than two splits');return(NULL)}
	
	mx = numeric(nc)
	for(i in 1:nc) 	mx[i] = window[[i]][1] + which.max(abs(ft[window[[i]][1]:window[[i]][2]]))
	
	steps = seq(mx[1],mx[2],stepsize)
	
	if(!is.null(fixstep)) steps = fixstep
	
	aicvec = numeric(length(steps))
	
	cat('[multisplit] fitting',length(steps),'splits...')
	
	if(.control.output(control)) cat('\n')
	
	for(i in 1:length(steps)) 
	{
		
		if(.control.output(control)) cat( 'fitting split',i,'@',steps[i],'\n')
		
		#split and estimate
		swalesol_new = split.f.fix(.swale.solution.internal(swalesol),steps[i])
		swalesol_new = estimate.ab(swalesol_new)
		
		#recalculate model/rss/aic
		swalesol_new = model(swalesol_new)
		swalesol_new = rss(swalesol_new)
		
		#recalculate solution
		swalesol_new = new('swale.solution',internal=swalesol_new,control=control)	
		#swalesol_new = calcSolution(swalesol_new)
		
		aicvec[i] = aic(swalesol_new) 
		
	}
	
	if(.control.output(control)) cat('\n')

	cat('done\n')
	#get best fit
	names(aicvec) = steps
	bestaic = which.min(aicvec)
	
	#split and estimate
	swalesol_new = split.f.fix(.swale.solution.internal(swalesol),steps[bestaic])
	swalesol_new = estimate.ab(swalesol_new)
	
	#recalculate model/rss/aic
	swalesol_new = model(swalesol_new)
	swalesol_new = rss(swalesol_new)
	
	#recalculate solution
	swalesol_new = new('swale.solution',internal=swalesol_new,control=control)	
	swalesol_new = calcSolution(swalesol_new)
	
	.swale.solution.aic(swalesol_new) = aic(swalesol_new)
	
	cat(' best AIC (',aicvec[bestaic],') at split ',names(aicvec[bestaic]),'\n',sep='')
	
	
	return(list(solution=swalesol_new,aicvec=aicvec))
	
}

