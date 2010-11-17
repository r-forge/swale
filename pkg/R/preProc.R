#############################################
# PREPROCESSING FUNCTIONS					#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#CONTAINS
#lowpass
#detrend
#dmean
#s2t
#estimateSNR

lowpass <- 
function(data,filsd=10,fillen=100,filmean=50) 
#apply lowpass filter to data
{
	
	data_sm <- matrix(0,nrow(data),ncol(data))	
	
	len=ncol(data)
	fl=round(fillen/2)
	filt=dnorm(1:100,mean=50,sd=filsd)
	
	for(i in 1:nrow(data)) 	{
		x <- convolve(data[i,],filt,type='open')
		data_sm[i,] <- x[-c(1:(fl-1),(length(x)-(fl-1)):length(x))]
	}
		
	return(data_sm)	
	
}

detrend <- 
function(data) 
#remove linear trend from data
{
	
	trenddata = data
	
	for(i in 1:nrow(data)) {
		
		y <- data[i,]
		x <- seq(1:length(y))
		
		trend <- lm(y ~ x)
		
		tr=coef(trend)[2]*x
		tr=tr+coef(trend)[1]
		
		data[i,] <- data[i,] - tr 
		trenddata[i,] = tr
	}
	
	return(list(data=data,trend=trenddata))
}

demean <- 
function(data)
#demean data
{
	for(i in 1:nrow(data)) data[i,]=data[i,]-mean(data[i,])
	return(data)
}


s2t <- 
function(s,Hz=512,unit='ms',quiet=T) 
#convert samples to time (ms or s)
{
	
	unit = match.arg(unit,c('ms','s'))
	
	if(unit=='ms') t = s*((1/Hz)*1000)
	if(unit=='s') t = s*(1/Hz)
	
	if(!quiet) cat(t,unit,'\n')
	return(invisible(t))
}


estimateSNR <- function(data,meth='FT')
#estimate snr (using fein & turetsky method and Handy method)
{
	
	meth = match.arg(meth,c('FT','H'))
	
	#Fein & Turetsky
	if(meth=='FT') {
		
		#average signal
		avg <- apply(data,2,mean)
		
		#noise power
		sn=0
		for(j in 1:nrow(data)) sn = sn + sum((data[j,]-avg)^2)
		sn = sn * (1/(ncol(data)*(nrow(data)-1)))
		
		#signal power
		ss = (1/(ncol(data)))*sum(avg^2)-((1/nrow(data))*sn)
		
		#signal-to-noise ratio
		snr = ss/sn
		
		return(snr)
	}
	
	#Handy
	if(meth=='H') {
		
		avg <- apply(data,2,mean)
		maxavg <- max(abs(avg))
		
		#noise SD (signal-avg)
		for(i in 1:nrow(data)) data[i,] = data[i,]-avg
		sdnoise = sd(as.vector(data))
		
		#signal to noise ratio
		snr = maxavg/sdnoise 
		
		return(snr)
	}
	
	return(NULL)
}


interval2Hz <- function(sampInt) 
{
	return(1/sampInt*1e+6)	
	
}

hz2interval <- function(Hz) 
{
	return(1/Hz*1e+6)	
	
}

