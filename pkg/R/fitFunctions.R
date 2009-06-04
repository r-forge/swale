#############################################
# FIT FUNCTIONS							    #
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#LEGEND#
# Y = data matrix (trial x samples)
# T = transformation matrix (samples x npoly)
# D = delta matrix (npoly x npoly)
# note: T%*%D may be replaced with analytical derivatives of f. TD = T%*%D (samples x npoly)
# I = identity matix (trial x trial)
# ab = amplitude (trial x 1) , latency (tria x 1) -> (trial x 2)
# f = waveform (npoly x 1)

#estimate waveform (polynomial coeficients)
estimate.f <- function(Y,T,TD,ab) 
{

	X = T%x%ab[,1] + TD%x%ab[,2]
	
	f = solve(t(X)%*%X)%*%t(X)%*%as.vector(Y)

	return(f)
		
}

#estimate amplitude and latency
estimate.ab <- function(Y,T,TD,I,f) 
{

	X = cbind( I%x%(T%*%f) , I%x%(TD%*%f) )
	
	ab = as.vector(t(Y))%*%X%*%solve(t(X)%*%X)

	ab = matrix(ab,,2)
	
	ab[,1]=ab[,1]/mean(ab[,1])
	
	return(ab)
	
}

#estimate model
model <- function(T,TD,f,ab)
{
	
	Y = (T%x%ab[,1] + TD%x%ab[,2])%*%f
	Y = matrix(Y,nrow(ab))
	
	return(Y)
		
}

#make orthogonal polynomials
P <- function(len=1,degree=1) 
{
	return(invisible(poly(seq(0,1,len=len),degree=degree)))
}

#make derivatives of polynomials
dP <- function(P,analytical=T)
{
	
	if(analytical) {
		dP=cbind(rep(1,nrow(P)),P)
		dP=dP[,-ncol(dP)]
		dP=dP*matrix(seq(1,ncol(P)),nrow(P),ncol(P),byrow=T)
	} else {
		x=diag(nrow(P))
		diag(x)=1
		diag(x[-1,])=-1
		dP=x%*%P	
	}
	
	return(dP)
	
}

#calculate rss
rss <- function(Y,yhat) 
{
	return(t(as.vector(Y)-as.vector(yhat))%*%((as.vector(Y)-as.vector(yhat))))
}

iterate <- function(Y,npoly=12,gradtol=c('fix',1e-06),burnin=1,maxiter=100,plot=TRUE,analytical.deriv=TRUE) 
{
	
	#init
	T = P(ncol(Y),npoly) 								#polynomials
	TD = dP(T,analytical=analytical.deriv)				#deriv poly
	I = diag(nrow(Y))		

	#start
	fstart = t(T)%*%apply(Y,2,mean)
	
	obj=numeric(0)
	
	if(plot){x11(width=7,height=5);layout(matrix(1:6,2,byrow=FALSE))}

	#burn-in	
	cat('burn-in (',burnin,')...',sep='')
	
	abhat = estimate.ab(Y,T,TD,I,fstart) 
	
	for(bi in 1:burnin) {
		fhat = estimate.f(Y,T,TD,abhat)
		abhat = estimate.ab(Y,T,TD,I,fhat) 
		fhat = estimate.f(Y,T,TD,abhat)
	}
	
	cat('ok\n')
	
	#calc modelfit after burn-in
	yhat = model(T,TD,fhat,abhat)
	obj = rss(Y,yhat)
	
	#set gradlim
	if(gradtol[1]=='fix') gradlim=as.numeric(gradtol[2]) else gradlim=obj*(as.numeric(gradtol[2])/100)
	cat('gradient limit:',gradlim,'\n')
	grad=obj
	
	#start iterates
	st_t=Sys.time()
	cat('starting iterations (',as.character(st_t),')\n')
	cat(sprintf('%3d: %10.0f',1,obj),'\n')
	
	#iterate init
	i=2
	dif=obj
	
	#iterate
	while(abs(grad[i-1])>gradlim) {	
		
		#estimate
		fhat = estimate.f(Y,T,TD,abhat)
		abhat = estimate.ab(Y,T,TD,I,fhat) 
		fhat = estimate.f(Y,T,TD,abhat)
		
		#modelfit
		yhat = model(T,TD,fhat,abhat)
		obj = rss(Y,yhat)

		#add vectors
		dif <- c(dif,obj)
		grad <- c(grad,dif[i]-dif[i-1])	
		
		#print add and continue
		cat(sprintf('%3d: %10.0f ~ (%16.6f)',i,obj,grad[i]),'\n')
		
		if(plot) plot.iterates(apply(Y,2,mean),apply(yhat,2,mean),T%*%fhat,TD%*%fhat,abhat,grad,dif,maxiter)
				
		i=i+1
		
		if(i>=maxiter) break()
		
		
	}
	
	
	yhat = model(T,TD,fhat,abhat)
	obj = rss(Y,yhat)
	
	dif <- c(dif,obj)
	
	cat('\n')
	cat('Final solution',dif[length(dif)],'\n')
	
	en_t=Sys.time()
	tt_t=difftime(en_t,st_t,units='min')
	
	cat('Iterated',i,'times. Process took',round(tt_t,1),'minutes.\n')

	
	internal <- list(yhat=yhat,T=T,TD=TD,fstart=fstart,pfhat=fhat,abhat=abhat,rss=dif[length(dif)],grad=grad[length(grad)])
	sol <- list(internal=internal,solution=list())
	
	return(sol)
	
}

solution <- function(sol)
{
	
	cat('calculating parameters...')
	
	#if(length(which(sol$internal$abhat[,1]==0))>0) sol$internal$abhat=sol$internal$abhat[-which(sol$internal$abhat[,1]==0),]
	
	
	sol$solution$a <- (sol$internal$abhat[,1])
	sol$solution$b <- (sol$internal$abhat[,2]/sol$internal$abhat[,1])*-1
	
	sol$solution$fhat <- (sol$internal$T%*%sol$internal$pfhat)*mean(sol$internal$abhat[,1])
	sol$solution$dfhat <- (sol$internal$TD%*%sol$internal$pfhat)*mean(sol$internal$abhat[,2])
	
	sol$solution$avgy <- sol$solution$fhat+sol$solution$dfhat
				
	sol$npar <- length(as.vector(sol$pfhat))+length(as.vector(sol$abhat))
	n <- length(as.vector(sol$yhat))

	sol$aic <- aic(n,sol$npar,sol$rss)
	
	cat('ok\n')
	return(sol)
}


aic <- function(n,k,rss)
{
	
	return( try( 2*k+(n*(log(2*pi*rss/n)+1)) ) )
	
}


peak.pick <- function(sol)
{
	
	zero <- which.max(sol$solution$avgy)
	
	trials <- nrow(sol$internal$yhat)
	
	lats <- numeric(trials)
	for(tr in 1:nrow(sol$internal$yhat)) {
		
		lats[tr] <- which.max(sol$solution$a[tr]*sol$solution$fhat+sol$solution$b[tr]*sol$solution$dfhat)-zero
		
	}
	
	return(lats)
}


peak.pick.data <- function(dat)
{
	
	avg = apply(dat$data,2,mean)
	zero = which.max(avg)
	trials = nrow(dat$data)
	
	lats <- numeric(trials)
	
	for(tr in 1:trials) {
		
		lats[tr] <- which.max(dat$data[tr,])-zero
		
	}
	
	return(lats)
	
	
}
