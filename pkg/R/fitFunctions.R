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
#swaleEEG

iterate <-
function(swaledat,control=new('control')) 
#iterate a f/ab estimation (with or without split)
{
	
	#plot iteration info
	cat('SWALE\n')
	cat('Iterations started:',date(),'\n')
	
	#check if startingvalues are sane
	if(nrow(.control.start.value(control))!=ncol(.basis.matrix(.swale.internal.basis(swaledat)))) stop('Length of startingvalues do not match with transform matrix.')

	#set startingvalues
	.swale.internal.waves(swaledat) = .control.start.value(control)
	
	#do iterations
	exitIterate = F
	iterNum = 0
	gradient = NA
	objective = numeric(0)	
	
	swaledat_old = swaledat
	
	while(!exitIterate) {
		
		iterNum = iterNum + 1
		
		#order if one waveform is estimated
		if(.control.split.type(control)=='none') {
			swaledat_new = estimate.ab(swaledat_old)
			swaledat_new = estimate.f(swaledat_new)
			swaledat_new = model(swaledat_new)
			swaledat_new = rss(swaledat_new)
		}
		
		#order if multiple waveforms are estimated
		if(.control.split.type(control)=='window') {
			
			if(iterNum==1) {
				cuts = split.f.one(swaledat_old,.control.split.data(control))
				swaledat_new = split.f.fix(swaledat_old,cuts)
			}
			
			if(iterNum>1) {
				swaledat_new = average.f(swaledat_old)
				swaledat_new = split.f.fix(swaledat_new,cuts)
			}
			
			swaledat_new = estimate.ab(swaledat_new)
			swaledat_new = model(swaledat_new)
			swaledat_new = rss(swaledat_new)
		}
		
		#calculate objective and gradient
		objective = c(objective,.swale.internal.rss(swaledat_new))
		if(iterNum>1) gradient = c(gradient,.swale.internal.rss(swaledat_new)-objective[iterNum-1])
		
		#show gradient and objective information
		if(.control.output(control)) cat(sprintf('%3d: %10.0f ~ (%16.9f)',iterNum,objective[iterNum],gradient[iterNum]),'\n')
		
		#check if gradient is smaller than convergence or maxIter is reached
		if(iterNum>2) if(abs(gradient[iterNum])<.control.iter.convergence(control)) exitIterate=TRUE
		if(iterNum>1) if(gradient[iterNum]>0) exitIterate=TRUE
		if(iterNum==.control.iter.limit(control)) exitIterate=TRUE
		
		if(exitIterate==FALSE) swaledat_old = swaledat_new
		
	}
	
	#set swale internal objects 
	.swale.internal.rss(swaledat_old) = objective[iterNum-1]
	.swale.internal.gradient(swaledat_old) = gradient[iterNum-1]
	
	solution = new('swale.solution',internal=swaledat_old,control=control)
	
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
function(swaledat) 
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
	
	#cleanTrials
	#solution = cleanTrials(solution)
	
	return(solution)
	

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
