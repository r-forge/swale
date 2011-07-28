#############################################
# METHODS									#
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################


setGeneric('plot',package='graphics')
setGeneric('summary',package='base')

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
			cat('[ eeg.data ]\n')
			cat(' condition   :',object@condition,'\n')
			cat(' trials      : ',object@trials,'\n',sep='')
			cat(' samples     : ',object@samples,' @',object@sampRate,'Hz\n',sep='')
			cat(' channel     : ',object@channel,sep='')
			cat('\n')
		}
)


setMethod('show',signature(object='swale.solution'),
		function(object) {
			cat('[ swale solution ]\n')
			cat(' Number of waves: ',ncol(object@internal@waves),'\n')
			cat(' Basisfunctions : ',object@internal@basis@num.funcs,'\n')
			cat(' Minimum        : ',object@internal@rss,'\n')
			cat(' Rsquared       : ',Rsq(object),'\n')
			cat(' AIC            : ',object@aic,'\n')
			cat(' BIC            : ',bic(object),'\n')
			cat(' ','\n')
			show(object@internal@eeg.data)
			cat('\n')
		}
)

setMethod('summary','swale.solution',
		function(object) {
			cat('[ swale solution ]\n')
			cat(' Minimum        : ',object@internal@rss,'\n')
			cat(' Rsquared       : ',Rsq(object),'\n')
			cat(' AIC            : ',object@aic,'\n')
			cat(' BIC            : ',bic(object),'\n')
			cat(' ','\n')
			
			ml = setMaxLat(object,F)
			
			nwaves = ncol(object@internal@waves)
			cat(' [summary] ',nwaves,'@',object@internal@basis@num.funcs,'\n')
			cat('  wave [1] (',object@control@peak.windows[[1]][1],',',object@control@peak.windows[[1]][2],')\n',sep='')
			cat('    method          :',object@control@peak.method[[1]],'\n')
			cat('    latency ranges  :',ml[[1]][1],',',ml[[1]][2],'\n')
			cat('    amplitude (IQR) :',round(quantile(object@amplitude[,1],c(.25,.5,.75)),4),'\n')
			cat('    latency   (IQR) :',round(quantile(object@latency[,1],c(.25,.5,.75)),4),'\n')
			if(nwaves>1) {
				for(i in 2:nwaves) {
					cat('  wave [',i,'] (',object@control@peak.windows[[i]][1],',',object@control@peak.windows[[i]][2],')\n',sep='')
					cat('    method          :',object@control@peak.method[[i]],'\n')
					cat('    latency ranges  :',ml[[i]][1],',',ml[[i]][2],'\n')
					cat('    amplitude (IQR) :',round(quantile(object@amplitude[,i],c(.25,.5,.75)),4),'\n')
					cat('    latency   (IQR) :',round(quantile(object@latency[,i],c(.25,.5,.75)),4),'\n')
				}
			}
			cat('\n')
		}
)


#swalesolution methods
setMethod('plot',signature(x='swale.solution',y='missing'),function(x) plotSolution(x))
