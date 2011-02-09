#############################################
# Basis function FUNCTIONS				    #
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################


#CONTAINS
#estimate.f
#average.f
#split.f.window
#split.f.fix
#estimate.ab
#model
#rss
#aic

estimate.f <- 
function(swaledat) 
#estimate waveform (polynomial coeficients)
{
	if(class(swaledat)!='swale.internal') stop('Input must be of class \'swale.internal\'')
	
	#set estimation objects
	Y = .eeg.data.data(.swale.internal.eeg.data(swaledat))
	TP = .basis.matrix(.swale.internal.basis(swaledat))
	dTP = .basis.deriv(.swale.internal.basis(swaledat))
	a = .swale.internal.amps(swaledat)
	b = .swale.internal.lats(swaledat)
	
	#estimate f
	X = TP%x%a + dTP%x%b
	f = solve(t(X)%*%X)%*%t(X)%*%as.vector(Y)
	.swale.internal.waves(swaledat) = matrix(f,,ncol(a),byrow=T)
	
	return(swaledat)
	
}

average.f <- 
function(swaledat) 
#calculate summed waveform based on mean amps and lats
{
	if(class(swaledat)!='swale.internal') stop('Input must be of class \'swale.internal\'')
	
	#set estimation objects
	TP = .basis.matrix(.swale.internal.basis(swaledat))
	dTP = .basis.deriv(.swale.internal.basis(swaledat))
	a = .swale.internal.amps(swaledat)
	b = .swale.internal.lats(swaledat)
	f = .swale.internal.waves(swaledat)
	trials = nrow(.eeg.data.data(.swale.internal.eeg.data(swaledat)))
	
	#calculate average f
	fhat = matrix(0,trials,nrow(TP))
	for(i in 1:ncol(a)) {
		for(trial in 1:trials) {
		 	fhat[trial,] = fhat[trial,] + ((TP%*%f[,i])*(a[trial,i])+(dTP%*%f[,i])*(b[trial,i]))
		}
	}
	fhatmean = as.vector(apply(fhat,2,mean))

	.swale.internal.waves(swaledat) = t(TP)%*%fhatmean
	
	return(swaledat)
	
}

split.f.window <- 
function(swaledat,window,plot=FALSE) 
#cut up data into different waveforms
#checks for maxima within a window and cuts in the midpoint between maxima
{
	if(class(swaledat)!='swale.internal') stop('Input must be of class \'swale.internal\'')
	
	#set estimation objects
	TP = .basis.matrix(.swale.internal.basis(swaledat))
	f = .swale.internal.waves(swaledat)
	ft = TP%*%f
	
	#search for maxima within window
	nc = nrow(window)
	mx = numeric(nc)
	for(i in 1:nc) 	mx[i] = window[i,1] + which.max(abs(ft[window[i,1]:window[i,2]]))
	
	#estimate cutpoints
	mx = sort(mx)
	cutpoints = numeric(nc-1)
	for(i in 1:(nc-1)) cutpoints[i] = round((mx[i]+mx[i+1])/2)
	
	#cut-up-nicely
	fs = matrix(ft,nrow(ft),nc)
	end = length(ft)
	
	#first and last
	fs[cutpoints[1]:end,1] = fs[cutpoints[1],1]
	fs[1:cutpoints[length(cutpoints)],nc] = fs[cutpoints[length(cutpoints)],nc]

	#cut in between
	if(length(cutpoints)>1) {
		for(i in 2:(nc-1)) {
			fs[1:cutpoints[i-1],i] = fs[cutpoints[i-1],i]
			fs[cutpoints[i]:end,i] = fs[cutpoints[i],i]
		}
	}
	
	#set waves matrix
	waves = matrix(NA,ncol(TP),nc)
	for(i in 1:nc) 	waves[,i] = t(TP)%*%fs[,i]	
	.swale.internal.waves(swaledat) = waves
	
	
	if(plot) {
		sumwave = apply(fs,1,sum)
		yl = c(min(cbind(fs,sumwave)),max(cbind(fs,sumwave)))
		
		par(las=1)
		
		plot(sumwave,type='l',lty=2,xlab='samples',ylab='mV',ylim=yl,bty='n',col=gray(.5),lwd=3)
		lines(apply(.eeg.data.data(.swale.internal.eeg.data(swaledat)),2,mean),col=gray(.2),lwd=5)
		for(i in 1:ncol(fs)) lines(fs[,i],lwd=3,lty=1,col=i)
		
	}
	
	return(swaledat)	

}

split.f.one <- 
function(swaledat,window) 
#cut up data into different waveforms
#checks for maxima within a window and cuts in the midpoint between maxima
{
	if(class(swaledat)!='swale.internal') stop('Input must be of class \'swale.internal\'')
	
	#set estimation objects
	TP = .basis.matrix(.swale.internal.basis(swaledat))
	f = .swale.internal.waves(swaledat)
	ft = TP%*%f
	
	#search for maxima within window
	nc = nrow(window)
	mx = numeric(nc)
	for(i in 1:nc) 	mx[i] = window[i,1] + which.max(abs(ft[window[i,1]:window[i,2]]))
	
	#estimate cutpoints
	mx = sort(mx)
	cutpoints = numeric(nc-1)
	for(i in 1:(nc-1)) cutpoints[i] = round((mx[i]+mx[i+1])/2)
	
	return(cutpoints)	
	
}

split.f.fix <- 
function(swaledat,cutpoints) 
#cut up data into different waveforms
#checks for maxima within a window and cuts in the midpoint between maxima
{
	if(class(swaledat)!='swale.internal') stop('Input must be of class \'swale.internal\'')
	
	#set estimation objects
	TP = .basis.matrix(.swale.internal.basis(swaledat))
	f = .swale.internal.waves(swaledat)
	ft = TP%*%f
	nc = length(cutpoints)+1
	
	#cut-up-nicely
	fs = matrix(ft,nrow(ft),nc)
	end = length(ft)
	
	#first and last
	fs[cutpoints[1]:end,1] = fs[cutpoints[1],1]
	fs[1:cutpoints[length(cutpoints)],nc] = fs[cutpoints[length(cutpoints)],nc]
	
	#cut in between
	if(length(cutpoints)>1) {
		for(i in 2:(nc-1)) {
			fs[1:cutpoints[i-1],i] = fs[cutpoints[i-1],i]
			fs[cutpoints[i]:end,i] = fs[cutpoints[i],i]
		}
	}
	
	#set waves matrix
	waves = matrix(NA,ncol(TP),nc)
	for(i in 1:nc) 	waves[,i] = t(TP)%*%fs[,i]	
	.swale.internal.waves(swaledat) = waves
	
	return(swaledat)	
	
}



estimate.ab <- 
function(swaledat) 
#estimate amplitude and latency
{
	if(class(swaledat)!='swale.internal') stop('Input must be of class \'swale.internal\'')
	
	#set estimation objects
	Y = .eeg.data.data(.swale.internal.eeg.data(swaledat))
	TP = .basis.matrix(.swale.internal.basis(swaledat))
	dTP = .basis.deriv(.swale.internal.basis(swaledat))
	f = .swale.internal.waves(swaledat)
	I = diag(nrow(Y))

	#estimate ab
	X = cbind( I%x%(TP%*%f) , I%x%(dTP%*%f) )
	ab = as.vector(t(Y))%*%X%*%solve(t(X)%*%X)
	
	#rearrange vectors to matrices
	nf = ncol(f)
	ab1 = matrix(ab,,nf,byrow=T)
	.swale.internal.amps(swaledat) = matrix(ab1[1:nrow(Y),],,nf)
	.swale.internal.lats(swaledat) = matrix(ab1[(nrow(Y)+1):nrow(ab1),],,nf)
	
	#scale amplitude to mean 1
	ma = apply(.swale.internal.amps(swaledat),2,mean)
	for(i in 1:ncol(.swale.internal.amps(swaledat))) .swale.internal.amps(swaledat)[,i] = .swale.internal.amps(swaledat)[,i] / ma[i] 
	
	return(swaledat)
	
}

model <- 
function(swaledat)
#calculate model estimates
{
	if(class(swaledat)!='swale.internal') stop('Input must be of class \'swale.internal\'')
	
	#set estimation objects
	TP = .basis.matrix(.swale.internal.basis(swaledat))
	dTP = .basis.deriv(.swale.internal.basis(swaledat))
	a = .swale.internal.amps(swaledat)
	b = .swale.internal.lats(swaledat)
	f = .swale.internal.waves(swaledat)
	
	#calculate model estimates
	Y = (TP%x%a + dTP%x%b)%*%as.vector(t(f))
	.swale.internal.model(swaledat) = matrix(Y,nrow(a),nrow(TP),byrow=F)
	
	return(swaledat)
	
}

rss <- 
function(swaledat) 
#calculate rss
{
	if(class(swaledat)=='swale.solution') swaledat = .swale.solution.internal(swaledat) else if(class(swaledat)!='swale.internal') stop('input must be of class internal or solution')
	
	#set estimation objects
	Y = .eeg.data.data(.swale.internal.eeg.data(swaledat))
	yhat = .swale.internal.model(swaledat)
	
	#calculate rss (innerproduct of residuals)
	.swale.internal.rss(swaledat) = as.vector(t(as.vector(Y)-as.vector(yhat))%*%((as.vector(Y)-as.vector(yhat))))
	
	return(swaledat)
}


aic <-
function(swaledat,correct=F,adjustWave=1,adjustAmp=1,adjustLat=1,adjustN=0) 
#calculate modelfit using AIC
{
	if(class(swaledat)=='swale.solution') swaledat = .swale.solution.internal(swaledat) else if(class(swaledat)!='swale.internal') stop('input must be of class internal or solution')
	
	k = (.basis.num.funcs(.swale.internal.basis(swaledat))*ncol(.swale.internal.waves(swaledat))*adjustWave)+length(as.vector(.swale.internal.amps(swaledat)))/adjustAmp+length(as.vector(.swale.internal.lats(swaledat)))/adjustLat
	n = length(as.vector(.eeg.data.data(.swale.internal.eeg.data(swaledat))))-adjustN
	rss = .swale.internal.rss(swaledat)
	
	aicvalue = try( 2*k + n*(log(rss/n)) )
	
	
	if(correct) aicvalue = aicvalue + ( (2*k*(k+1)) / (n-k-1) )
		
	return(aicvalue)
}


bic <-
function(swaledat) 
#calculate modelfit using BIC
{
	if(class(swaledat)=='swale.solution') swaledat = .swale.solution.internal(swaledat) else if(class(swaledat)!='swale.internal') stop('input must be of class internal or solution')
	
	k = .basis.num.funcs(.swale.internal.basis(swaledat))*ncol(.swale.internal.waves(swaledat))+length(as.vector(.swale.internal.amps(swaledat)))+length(as.vector(.swale.internal.lats(swaledat)))
	n = length(as.vector(.eeg.data.data(.swale.internal.eeg.data(swaledat))))
	resids = as.vector(.eeg.data.data(.swale.internal.eeg.data(swaledat)) - .swale.internal.model(swaledat))
	rss = .swale.internal.rss(swaledat)
	
	#bicvalue = try( n*(log(var(resids)/n)) + k*log(n) )
	bicvalue = try( n*(log(rss/n)) - k*log(n) )
	
	return(bicvalue)
}


Rsq <-
function(swaledat) 
#calculate Rsquared value for a model
{
	if(class(swaledat)=='swale.solution') swaledat = .swale.solution.internal(swaledat) else if(class(swaledat)!='swale.internal') stop('input must be of class internal or solution')
	
	#set estimation objects
	Y = .eeg.data.data(.swale.internal.eeg.data(swaledat))
	yhat = .swale.internal.model(swaledat)
	
	#calculate rss (innerproduct of residuals)
	SSres = as.vector(t(as.vector(Y)-as.vector(yhat))%*%((as.vector(Y)-as.vector(yhat)))) 
	SStot = as.vector(t(as.vector(Y))%*%((as.vector(Y)))) 
	
	Rsq = 1 - (SSres/SStot)
	
	return(Rsq)
}

swalecor <-
function(swalesol) 
#calculate correlations between amplitude and latency parameters in a solution
{
	ncors = ncol(.swale.solution.waveform(swalesol))
	
	ampcor = ampcor.p = latcor = latcor.p = matrix(NA,ncors,ncors)
	
	valid = which(.swale.solution.discard(swalesol)!=1)
	
	if(length(valid)<2) {warning('Not enough valid trials!')} else {
		
		for(row in 1:ncors) {
			for(col in 1:ncors) {
				ampcor[row,col] = cor.test(.swale.solution.amplitude(swalesol)[valid,row],.swale.solution.amplitude(swalesol)[valid,col])$estimate
				ampcor.p[row,col] = cor.test(.swale.solution.amplitude(swalesol)[valid,row],.swale.solution.amplitude(swalesol)[valid,col])$p.value
				latcor[row,col] = cor.test(.swale.solution.latency(swalesol)[valid,row],.swale.solution.latency(swalesol)[valid,col])$estimate
				latcor.p[row,col] = cor.test(.swale.solution.latency(swalesol)[valid,row],.swale.solution.latency(swalesol)[valid,col])$p.value
			}
		}
		
		
	}
	
	return(list(amplitude=list(estimate=ampcor,p.value=ampcor.p),latency=list(estimate=latcor,p.value=latcor.p)))
}


