#############################################
# METHODS									#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################


setGeneric('plot',package='graphics',where=.GlobalEnv)

# eeg.data.methods
setMethod('plot',signature(x='eeg.data',y='missing'),
		function(x) {
			
			mn=min(x@data)
			mx=max(x@data)
			
			par(las=1)
			plot(NA,NA,xlim=c(1,ncol(x@data)),ylim=c(mn,mx),xlab='time (ms)',ylab='mV',bty='n',main=paste(x@channel,': ',x@condition,sep=''),axes=F)
			axis(2)
			axis(1,at=axTicks(1),labels=round(axTicks(1)*(1/x@sampRate)*1000))
			
			for(i in 1:nrow(x@data)) lines(x@data[i,],lwd=1,col=i)
			
			lines(apply(x@data,2,mean),lwd=5,col=1,lty=1)
			
		}
)

setMethod('show',signature(object='eeg.data'),
		function(object) {
			cat('eeg.data [',object@condition,']\n')
			cat('trials      : ',object@trials,'\n',sep='')
			cat('samples     : ',object@samples,' @',object@sampRate,'Hz\n',sep='')
			cat('channel     : ',object@channel,sep='')
			cat('\n')
		}
)


#swalesolution methods
setMethod('plot',signature(x='swale.solution',y='missing'),function(x) plotSolution(x))
