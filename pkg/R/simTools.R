#############################################
#  SIMTOOLS FUNCTIONS		 				#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################


makeSims <- function(datinp,asd=.5,lsd=10,tn=40,sql=350,snr=10,snrmethod=2) {
	
	
	y <- datinp
	x <- seq(1:length(y))
	
	trend <- lm(y ~ x)
	
	tr=coef(trend)[2]*x
	tr=tr+coef(trend)[1]
	
	datinp <- datinp - tr 
	
	
	#amp creation
	arng=c(1e-10,1e+10)
	amps=abs(rnorm(tn,mean=1,sd=asd))
	amps[amps<arng[1]]=arng[1]
	amps[amps>arng[2]]=arng[2]
	#amps[1]=1
	
	#lat creation
	lrng=c(-3000,3000)
	mid=which.max(datinp)
	lats=round(rnorm(tn,mean=0,sd=lsd))
	lats[lats<lrng[1]]=lrng[1]
	lats[lats>lrng[2]]=lrng[2]
	#lats[1]=0
	
	#shift latencies and do amplitude
	dat=matrix(0,tn,sql)

	for(i in 1:dim(dat)[1]) {

		start=mid-(sql/2)
		start=start-lats[i]
		stop=start+(sql)
		
		if(start<1) {start=1;stop=start+sql}
		
		dat[i,]=datinp[seq(start,stop-1)]
		
		n=which(is.na(dat[i,]))
		if(length(n)>0) {
			#cat(i,'\n');
			#cat(dat[i,],'\n');
			dat[i,n]=rep(dat[i,n[1]-1],length(n));
			#cat(dat[i,],'\n')
		}
		
		dat[i,]=dat[i,]*amps[i] 
		
	}

	
	ml=round(mean(lats))	
	
	
	#signal to noise ratio 
	avg=apply(dat,2,mean)
	
	if(snrmethod==1) {
		if(snr==0) sdnoise=0 else sdnoise=max(avg)/(snr/sqrt(tn))
		cat('Set SNR =',snr,'\n')
		cat('Per trial SNR =',(snr/sqrt(tn)),'\n')
		sdnoise=rep(sdnoise,tn)
	}
	
	if(snrmethod==2) {
		sdnoise=numeric(tn)
		for(i in 1:tn) sdnoise[i]=max(avg)/(snr)
		cat('Set per trial SNR =',snr,'\n')
	}
	

	#cat('sdnoise',sdnoise,'\n')
	#add noise
	for(i in 1:tn) 	dat[i,]=dat[i,]+rnorm(length(dat[i,]),0,sdnoise[i])
	cat(' est SNR (mean) =',max(avg)/sd(apply(dat,2,mean)-avg),'\n')
	cat('\n')
	
	
	avgnoise <- apply(dat,2,mean)
	
	snrs=numeric(0)
	
	for(i in 1:nrow(dat)) 	{
		
		x <- dat[i,]-avgnoise
		snrs=c(snrs,max(avg)/sd(x))
		
	}
	
	#return data
	return(list(data=dat,amps=amps,lats=lats,inp=datinp,avginp=avg,stsnr=snrs))
}

makeSims2 <- function(asd=.5,lsd=10,tn=40,sql=350,snr=10,snrmethod=2) {
	
	#amp creation
	#arng=c(1e-10,1e+10)
	arng=c(1e-1,1e+1)
	
	amps1=abs(rnorm(tn,mean=1,sd=asd))
	amps2=abs(rnorm(tn,mean=1,sd=asd))
	
	amps1[amps1<arng[1]]=arng[1]
	amps1[amps1>arng[2]]=arng[2]
	
	amps2[amps2<arng[1]]=arng[1]
	amps2[amps2>arng[2]]=arng[2]
	
	amps1=amps2 ## SET TO SAME AMP
	
	#lat creation
	lrng=c(-75,75)
	mid=round(sql/2)
	
	lats1=round(rnorm(tn,mean=0,sd=lsd))
	lats2=round(rnorm(tn,mean=0,sd=lsd))
	
	lats1[lats1<lrng[1]]=lrng[1]
	lats1[lats1>lrng[2]]=lrng[2]
	
	lats2[lats2<lrng[1]]=lrng[1]
	lats2[lats2>lrng[2]]=lrng[2]
	
	
	cn=500
	
	#shift latencies and do amplitude
	dat=matrix(0,tn,sql)
	
	
	for(i in 1:nrow(dat)) {
		
		f1=dgamma((1:sql),(110-lats2[i]))*cn*(amps1[i]/4) ## SET TO SAME AMP AND SAME LAT
		#f2=dnorm(1:sql,(200-lats2[i]),sd=20)*cn*amps2[i]
		
		#f1a=dgamma((1:sql),(100-lats1[i]))*cn*amps1[i]*1
		#f1b=dgamma((1:sql),(100-lats1[i]))*cn*amps1[i]*1
		#f1=(f1a+f1b)/2
		#l1=which.max(f1)
		#f2=dgamma((1:sql),(l1-25-lats2[i]))*cn*amps2[i]*-1
		#f2=dgamma((1:sql),(220-lats2[i]))*cn*amps2[i]*-1
		
		#dat[i,]=(f1+f2)/2
		dat[i,]=f1
		#browser()
	}
	
	datinp=(((dgamma((1:sql),100)*cn+dgamma((1:sql),100)*cn)/2)+dgamma((1:sql),220)*-cn)/2
	
	
	#signal to noise ratio 
	avg=apply(dat,2,mean)
	
	if(snrmethod==1) {
		if(snr==0) sdnoise=0 else sdnoise=max(avg)/(snr/sqrt(tn))
		cat('Set SNR =',snr,'\n')
		cat('Per trial SNR =',(snr/sqrt(tn)),'\n')
		sdnoise=rep(sdnoise,tn)
	}
	
	if(snrmethod==2) {
		sdnoise=numeric(tn)
		for(i in 1:tn) sdnoise[i]=max(avg)/(snr)
		cat('Set per trial SNR =',snr,'\n')
	}
	
	
	#add noise
	for(i in 1:tn) 	dat[i,]=dat[i,]+rnorm(length(dat[i,]),0,sdnoise[i])
	cat('Est SNR =',max(avg)/sd(apply(dat,2,mean)-avg),'\n')
	cat('\n')
	
	avgnoise <- apply(dat,2,mean)
	
	snrs=numeric(0)
	
	for(i in 1:nrow(dat)) 	{
		
		x <- dat[i,]-avgnoise
		snrs=c(snrs,max(avg)/sd(x))
		
	}
	
	#return data
	return(list(data=dat,amps=amps2,lats=lats2,inp=datinp,avginp=avg,stsnr=snrs))
}



makeSims3 <- function(datinp,T,asd=.5,lsd=.1,tn=40,sql=350,snr=10,snrmethod=2)
{
	
	#abhat
	avec=rnorm(tn,1,asd)
	bvec=rnorm(tn,0,lsd)
	#avec=rep(1,tn)
	#bvec=rep(0,tn)
	
	abvec=bvec*avec
	
	
	abhat = cbind(avec,abvec)
	
	#fhat
	mid=which.max(datinp)
	start=mid-(sql/2)
	start=start-0
	stop=start+(sql)
	datinp=datinp[seq(start,stop-1)]
	fhat = t(T)%*%datinp
	
	TD=dP(T,an=TRUE)
	
	#data with correct model
	dat = model(T,TD,fhat,abhat)
		
	
	#signal to noise ratio 
	avg=apply(dat,2,mean)
	
	if(snrmethod==1) {
		if(snr==0) sdnoise=0 else sdnoise=max(avg)/(snr/sqrt(tn))
		cat('Set SNR =',snr,'\n')
		cat('Per trial SNR =',(snr/sqrt(tn)),'\n')
		sdnoise=rep(sdnoise,tn)
	}
	
	if(snrmethod==2) {
		sdnoise=numeric(tn)
		for(i in 1:tn) sdnoise[i]=max(avg)/(snr)
		cat('Set per trial SNR =',snr,'\n')
	}
	
	
	#cat('sdnoise',sdnoise,'\n')
	#add noise
	for(i in 1:tn) 	dat[i,]=dat[i,]+rnorm(length(dat[i,]),0,sdnoise[i])
	cat(' est SNR (mean) =',max(avg)/sd(apply(dat,2,mean)-avg),'\n')
	cat('\n')
	
	
	avgnoise <- apply(dat,2,mean)
	
	snrs=numeric(0)
	
	for(i in 1:nrow(dat)) 	{
		
		x <- dat[i,]-avgnoise
		snrs=c(snrs,max(avg)/sd(x))
		
	}
	
	#return data
	return(list(data=dat,amps=avec,lats=abvec,inp=datinp,avginp=avg,stsnr=snrs))
		
}





est.snr <- function(dat)
{
	
	avg <- apply(dat$data,2,mean)
	
	snrs=numeric(0)
	
	for(i in 1:nrow(dat$data)) 	{
		
		x <- dat$data[i,]-avg
		snrs=c(snrs,max(avg)/sd(x))
				
	}
	
	snrmean <- max(dat$avginp)/sd(avg-dat$avginp)
	
		
	cat('estimated SNR:',mean(snrs),'\n')
	cat('estimated SNR mean:',snrmean,'\n')
	
	
	return(snrs)
}


Delta <- function(a) 
{
	
	x=diag(length(a))
	diag(x)=1
	diag(x[-1,])=-1
	
	return(x%*%a)	
}







