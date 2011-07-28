#############################################
# OLD FUNCTIONS								#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

iterateOld <-
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
		
		#order if multiple waveforms are estimated
		if(.control.split.type(control)=='window') {
			
			if(iterNum==1) {
				cuts = split.f.one(swaledat_old,.control.peak.windows(control))
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


calcSolutionOld <-
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
	
	return(solution)
}

iterateDiscard <-
		function(swaledat,control=new('control'),latRange=NULL,ampRange=c(1e-06,Inf),posGradStop=F) 
#iterate and discard trials until good fit is reached
{
	stopDisc = FALSE 
	removals = numeric(0)
	
	while(!stopDisc) {
		#iterate
		solution = iterate(swaledat,control,posGradStop)
		solution = calcSolution(solution)
		
		#check solution
		if(is.null(latRange)) latRange = setMaxLat(solution,prec.fac=1.5,plot=F)
		solution = cleanTrials(solution,lat.thres=latRange,amp.thres=ampRange)
		
		#determine discards
		if(sum(.swale.solution.discard(solution))==0) {
			stopDisc = TRUE
			break()
		} else {
			
			if(length(which(.swale.solution.discard(solution)!=0))>=(nrow(.eeg.data.data(.swale.internal.eeg.data(.swale.solution.internal(solution))))-1)) { 
				warning('Stopping iterations no valid model!')
				stopDisc=TRUE
				break()
			}
			
			#set data 
			oldat = .swale.internal.eeg.data(.swale.solution.internal(solution))
			
			#make new data with discards removed
			newdat = .eeg.data.data(.swale.internal.eeg.data(.swale.solution.internal(solution)))
			
			remrows = which(.swale.solution.discard(solution)!=0)
			newdat = newdat[-remrows,]
			newdat = detrend(newdat,meantrend=F)
			
			#make new data element
			data = new('eeg.data')
			.eeg.data.data(data) = as.matrix(newdat$data)
			.eeg.data.trend(data) = as.matrix(newdat$trend)
			.eeg.data.trials(data) = nrow(newdat$data)
			.eeg.data.samples(data) = .eeg.data.samples(oldat)
			.eeg.data.sampRate(data) = .eeg.data.sampRate(oldat)
			.eeg.data.channel(data) = .eeg.data.channel(oldat)
			.eeg.data.condition(data) = .eeg.data.condition(oldat)
			
			swaledat = new('swale.internal',eeg.data=data,basis=.swale.internal.basis(.swale.solution.internal(solution)))
			.control.start.value(control) = makeStart(swaledat)
			
		} #discardsloop
		removals = c(removals,list(remrows))
		cat('[discard] Removed',length(which(.swale.solution.discard(solution)!=0)),'trials. Refitting.\n')
		
	} #main iteration loop
	
	cat('All trials valid.\n')
	return(list(solution=solution,removals=removals))
	
}

iterateSplit <-
		function(swalesol,control=new('control'),latRange=NULL,ampRange=c(1e-06,Inf),posGradStop=F) 
#iterate and split based on removal of estimated waves
{
	#split up the data	
	cuts = split.f.one(.swale.solution.internal(swalesol),.control.peak.windows(control))
	swalesol_new = split.f.fix(.swale.solution.internal(swalesol),cuts)
	swalesol_new = estimate.ab(swalesol_new)
	
	#set estimation objects
	TP = .basis.matrix(.swale.internal.basis(swalesol_new))
	dTP = .basis.deriv(.swale.internal.basis(swalesol_new))
	a = .swale.internal.amps(swalesol_new)
	b = .swale.internal.lats(swalesol_new)
	f = .swale.internal.waves(swalesol_new)
	origdat = .eeg.data.data(.swale.internal.eeg.data(swalesol_new))
	trials = nrow(origdat)
	
	
	solutionlist = vector('list',ncol(f))
	
	#reverse estimate
	for(waves in 1:ncol(f)) {
		
		#set correct matrices
		newdat = origdat
		cTP = TP
		cdTP = dTP
		ca = as.matrix(a[,-waves])
		cb = as.matrix(b[,-waves])
		cf = as.matrix(f[,-waves])
		
		#remove wavetrend from data 
		fhat = matrix(0,trials,nrow(cTP))
		for(i in 1:ncol(ca)) {
			for(trial in 1:trials) {
				newdat[trial,] = newdat[trial,] - ((cTP%*%cf[,i])*(ca[trial,i])+(cdTP%*%cf[,i])*(cb[trial,i]))
			}
		}		
		
		#set swaleobject newdata
		.eeg.data.data(.swale.internal.eeg.data(swalesol_new)) = newdat
		control_new = control
		.control.split.type(control_new) = 'none'
		.control.peak.windows(control_new) = NULL
		
		cat('[split] Fitting wave [',waves,']\n',sep='')
		
		#iterate onewave
		solution = iterate(swalesol_new,control_new,posGradStop=posGradStop)	
		solution = calcSolution(solution)
		solutionlist[[waves]] = list(solution=solution,data=newdat)
		
	}
	
	#rebuild solution
	soldouble = new('swale.internal')
	.swale.internal.eeg.data(soldouble) = .swale.internal.eeg.data(.swale.solution.internal(swalesol))
	.swale.internal.basis(soldouble) = .swale.internal.basis(.swale.solution.internal(swalesol))
	
	#concatenate waves and amps and lats
	waves = amps = lats = numeric(0)
	for(i in 1:length(solutionlist)) {
		waves = cbind(waves,.swale.internal.waves(.swale.solution.internal(solutionlist[[i]]$solution)))
		amps = cbind(amps,.swale.internal.amps(.swale.solution.internal(solutionlist[[i]]$solution)))
		lats = cbind(lats,.swale.internal.lats(.swale.solution.internal(solutionlist[[i]]$solution)))
	}
	.swale.internal.waves(soldouble) = waves
	.swale.internal.amps(soldouble) = amps
	.swale.internal.lats(soldouble) = lats
	
	#recalculate model/rss/aic
	soldouble = model(soldouble)
	soldouble = rss(soldouble)
	
	#recalculate solution
	soldouble = new('swale.solution',internal=soldouble,control=control)	
	soldouble = calcSolution(soldouble)
	
	#return
	return(soldouble)
}


autoSplit <-
		function(swalesol,control=new('control'),latRange=NULL,ampRange=c(1e-06,Inf),posGradStop=F,gradlim=.1)
#split reiterated splitted solution again
{
	stopSplit = FALSE
	
	newsol = swalesol
	i=1
	
	obj = gradient = numeric(0)
	while(!stopSplit) {
		if(i>1) .swale.internal.waves(.swale.solution.internal(newsol)) = matrix(apply(.swale.internal.waves(.swale.solution.internal(soldouble)),1,sum),,1)	
		soldouble = iterateSplit(newsol,control=control)
		obj = c(obj,soldouble@internal@rss)
		
		if(i>1) {
			cat('[gradient]',obj[i]-obj[i-1],'\n')
			gradient = c(gradient,obj[i]-obj[i-1])
			if(abs(gradient[i-1])<=gradlim) stopSplit=T
		}
		
		i=i+1
	}
	
	return(soldouble)
	
}

testSplits <-
		function(swaledat,control=new('control')) 
#calculate and test which (amplitude or latency) contribute to the differences
{
	if(class(swaledat)=='swale.solution') swaledat = .swale.solution.internal(swaledat) else if(class(swaledat)!='swale.internal') stop('input must be of class internal or solution')
	
	numwave = ncol(.swale.internal.waves(swaledat))
	numtrial = .eeg.data.trials(.swale.internal.eeg.data(swaledat))
	
	if(numwave<=1) warning('Only one wave in model, new models might give same results!')
	
	#allfix
	allfix = swaledat
	.swale.internal.amps(allfix) = matrix(apply(.swale.internal.amps(swaledat),1,mean),numtrial,numwave)
	.swale.internal.lats(allfix) = matrix(apply(.swale.internal.lats(swaledat),1,mean),numtrial,numwave)
	allfix = model(allfix)
	allfix = rss(allfix)
	allfix = new('swale.solution',internal=allfix,control=control)	
	allfix = calcSolution(allfix)
	.swale.solution.aic(allfix) = aic(allfix,adjustAmp=numwave,adjustLat=numwave)
	
	#ampfix
	ampfix = swaledat
	.swale.internal.amps(ampfix) = matrix(apply(.swale.internal.amps(swaledat),1,mean),numtrial,numwave)
	ampfix = model(ampfix)
	ampfix = rss(ampfix)
	ampfix = new('swale.solution',internal=ampfix,control=control)	
	ampfix = calcSolution(ampfix)
	.swale.solution.aic(ampfix) = aic(ampfix,adjustAmp=numwave)
	
	#latfix
	latfix = swaledat
	.swale.internal.lats(latfix) = matrix(apply(.swale.internal.lats(swaledat),1,mean),numtrial,numwave)
	latfix = model(latfix)
	latfix = rss(latfix)
	latfix = new('swale.solution',internal=latfix,control=control)	
	latfix = calcSolution(latfix)
	.swale.solution.aic(latfix) = aic(latfix,adjustLat=numwave)
	
	return(list(allfix=allfix,ampfix=ampfix,latfix=latfix))
	
}

setMaxLatOld <-
		function(swalesol,window=list(NULL),prec.fac=list(1.5),method=list('abs'),discarded=list(10),plot=F,manual=NULL) 
#determine the maximum latency a waveform + derivative can model
{
	
	if(length(grep('control',slotNames(new('swale.solution'))))!=0) swalesol = new('swale.solution',swalesol,discard=integer(0),control=new('control')) else swalesol = new('swale.solution',swalesol,discard=integer(0))
	
	ranges = vector('list',ncol(.swale.solution.waveform(swalesol)))
	
	for(wave in 1:ncol(.swale.solution.waveform(swalesol)))  {
		
		method[[wave]] = match.arg(method[[wave]],c('abs','max','min','close','firstmax','firstmin','lastmax','lastmin'))
		
		nsamp = .eeg.data.samples(.swale.internal.eeg.data(.swale.solution.internal(swalesol)))
		
		if(is.null(window[[wave]])) window[[wave]] = c(1,nsamp) 
		
		if(method[[wave]]=='abs') absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='max') absmax = which.max((.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='min') absmax = which.min((.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='close') absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='firstmax' | method[[wave]]=='firstmin') absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		if(method[[wave]]=='lastmax' | method[[wave]]=='lastmin') absmax = which.max(abs(.swale.solution.waveform(swalesol)[window[[wave]][1]:window[[wave]][2],wave]))+window[[wave]][1]-1
		
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
	
	#manual override
	if(!is.null(manual)) ranges = manual
	
	return(ranges)
}
