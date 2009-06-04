#############################################
# PREPROCESSING FUNCTIONS					#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

lowpass <- function(dat,filsd=10,fillen=100,filmean=50,replace=TRUE) 
{
	
	cat('lowpass filtering...')
	dat$data_sm <- matrix(0,nrow(dat$data),ncol(dat$data))	
	
	len=ncol(dat$data)
	fl=round(fillen/2)
	filt=dnorm(1:100,mean=50,sd=filsd)
	
	plot.freq(filt)
	
	
	for(i in 1:nrow(dat$data)) 	{
		x <- convolve(dat$data[i,],filt,type='open')
		dat$data_sm[i,] <- x[-c(1:(fl-1),(length(x)-(fl-1)):length(x))]
		
	}
	
	if(replace) dat$data <- dat$data_sm
	cat('ok\n')
	return(dat)	
	
}

detrend <- function(dat) 
{
	
	cat('detrending...')
	
	for(i in 1:nrow(dat$data)) {
		
		y <- dat$data[i,]
		x <- seq(1:length(y))
		
		trend <- lm(y ~ x)
		
		tr=coef(trend)[2]*x
		tr=tr+coef(trend)[1]
		
		dat$data[i,] <- dat$data[i,] - tr 
		
	}
	cat('ok\n')
	
	return(dat)
	
}

demean <- function(dat)
{
	
	cat('demeaning...')
	
	for(i in 1:nrow(dat$data)) dat$data[i,]=dat$data[i,]-mean(dat$data[i,])
		
	cat('ok\n')
	
	return(dat)
	
	
}