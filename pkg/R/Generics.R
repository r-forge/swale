classname <-'version'
funcname <-'.version.version'
standGen <- function(object) standardGeneric('.version.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.version.version<-'
standGen <- function(x, value) standardGeneric('.version.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'version'
funcname <-'.version.build'
standGen <- function(object) standardGeneric('.version.build')
standMethod <- function(object) object@build
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.version.build<-'
standGen <- function(x, value) standardGeneric('.version.build<-')
standMethod <- function(x, value) {x@build<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'version'
funcname <-'.version.update'
standGen <- function(object) standardGeneric('.version.update')
standMethod <- function(object) object@update
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.version.update<-'
standGen <- function(x, value) standardGeneric('.version.update<-')
standMethod <- function(x, value) {x@update<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'version'
funcname <-'.version.svnrev'
standGen <- function(object) standardGeneric('.version.svnrev')
standMethod <- function(object) object@svnrev
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.version.svnrev<-'
standGen <- function(x, value) standardGeneric('.version.svnrev<-')
standMethod <- function(x, value) {x@svnrev<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'basis'
funcname <-'.basis.type'
standGen <- function(object) standardGeneric('.basis.type')
standMethod <- function(object) object@type
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.basis.type<-'
standGen <- function(x, value) standardGeneric('.basis.type<-')
standMethod <- function(x, value) {x@type<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'basis'
funcname <-'.basis.deriv.type'
standGen <- function(object) standardGeneric('.basis.deriv.type')
standMethod <- function(object) object@deriv.type
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.basis.deriv.type<-'
standGen <- function(x, value) standardGeneric('.basis.deriv.type<-')
standMethod <- function(x, value) {x@deriv.type<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'basis'
funcname <-'.basis.num.funcs'
standGen <- function(object) standardGeneric('.basis.num.funcs')
standMethod <- function(object) object@num.funcs
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.basis.num.funcs<-'
standGen <- function(x, value) standardGeneric('.basis.num.funcs<-')
standMethod <- function(x, value) {x@num.funcs<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'basis'
funcname <-'.basis.num.points'
standGen <- function(object) standardGeneric('.basis.num.points')
standMethod <- function(object) object@num.points
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.basis.num.points<-'
standGen <- function(x, value) standardGeneric('.basis.num.points<-')
standMethod <- function(x, value) {x@num.points<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'basis'
funcname <-'.basis.matrix'
standGen <- function(object) standardGeneric('.basis.matrix')
standMethod <- function(object) object@matrix
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.basis.matrix<-'
standGen <- function(x, value) standardGeneric('.basis.matrix<-')
standMethod <- function(x, value) {x@matrix<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'basis'
funcname <-'.basis.deriv'
standGen <- function(object) standardGeneric('.basis.deriv')
standMethod <- function(object) object@deriv
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.basis.deriv<-'
standGen <- function(x, value) standardGeneric('.basis.deriv<-')
standMethod <- function(x, value) {x@deriv<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'basis'
funcname <-'.basis.version'
standGen <- function(object) standardGeneric('.basis.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.basis.version<-'
standGen <- function(x, value) standardGeneric('.basis.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.data'
standGen <- function(object) standardGeneric('.eeg.data.data')
standMethod <- function(object) object@data
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.data<-'
standGen <- function(x, value) standardGeneric('.eeg.data.data<-')
standMethod <- function(x, value) {x@data<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.trend'
standGen <- function(object) standardGeneric('.eeg.data.trend')
standMethod <- function(object) object@trend
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.trend<-'
standGen <- function(x, value) standardGeneric('.eeg.data.trend<-')
standMethod <- function(x, value) {x@trend<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.trials'
standGen <- function(object) standardGeneric('.eeg.data.trials')
standMethod <- function(object) object@trials
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.trials<-'
standGen <- function(x, value) standardGeneric('.eeg.data.trials<-')
standMethod <- function(x, value) {x@trials<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.samples'
standGen <- function(object) standardGeneric('.eeg.data.samples')
standMethod <- function(object) object@samples
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.samples<-'
standGen <- function(x, value) standardGeneric('.eeg.data.samples<-')
standMethod <- function(x, value) {x@samples<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.sampRate'
standGen <- function(object) standardGeneric('.eeg.data.sampRate')
standMethod <- function(object) object@sampRate
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.sampRate<-'
standGen <- function(x, value) standardGeneric('.eeg.data.sampRate<-')
standMethod <- function(x, value) {x@sampRate<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.channel'
standGen <- function(object) standardGeneric('.eeg.data.channel')
standMethod <- function(object) object@channel
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.channel<-'
standGen <- function(x, value) standardGeneric('.eeg.data.channel<-')
standMethod <- function(x, value) {x@channel<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.condition'
standGen <- function(object) standardGeneric('.eeg.data.condition')
standMethod <- function(object) object@condition
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.condition<-'
standGen <- function(x, value) standardGeneric('.eeg.data.condition<-')
standMethod <- function(x, value) {x@condition<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.within'
standGen <- function(object) standardGeneric('.eeg.data.within')
standMethod <- function(object) object@within
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.within<-'
standGen <- function(x, value) standardGeneric('.eeg.data.within<-')
standMethod <- function(x, value) {x@within<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'eeg.data'
funcname <-'.eeg.data.version'
standGen <- function(object) standardGeneric('.eeg.data.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.eeg.data.version<-'
standGen <- function(x, value) standardGeneric('.eeg.data.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.waves'
standGen <- function(object) standardGeneric('.swale.internal.waves')
standMethod <- function(object) object@waves
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.waves<-'
standGen <- function(x, value) standardGeneric('.swale.internal.waves<-')
standMethod <- function(x, value) {x@waves<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.amps'
standGen <- function(object) standardGeneric('.swale.internal.amps')
standMethod <- function(object) object@amps
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.amps<-'
standGen <- function(x, value) standardGeneric('.swale.internal.amps<-')
standMethod <- function(x, value) {x@amps<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.lats'
standGen <- function(object) standardGeneric('.swale.internal.lats')
standMethod <- function(object) object@lats
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.lats<-'
standGen <- function(x, value) standardGeneric('.swale.internal.lats<-')
standMethod <- function(x, value) {x@lats<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.model'
standGen <- function(object) standardGeneric('.swale.internal.model')
standMethod <- function(object) object@model
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.model<-'
standGen <- function(x, value) standardGeneric('.swale.internal.model<-')
standMethod <- function(x, value) {x@model<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.rss'
standGen <- function(object) standardGeneric('.swale.internal.rss')
standMethod <- function(object) object@rss
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.rss<-'
standGen <- function(x, value) standardGeneric('.swale.internal.rss<-')
standMethod <- function(x, value) {x@rss<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.gradient'
standGen <- function(object) standardGeneric('.swale.internal.gradient')
standMethod <- function(object) object@gradient
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.gradient<-'
standGen <- function(x, value) standardGeneric('.swale.internal.gradient<-')
standMethod <- function(x, value) {x@gradient<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.version'
standGen <- function(object) standardGeneric('.swale.internal.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.version<-'
standGen <- function(x, value) standardGeneric('.swale.internal.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.basis'
standGen <- function(object) standardGeneric('.swale.internal.basis')
standMethod <- function(object) object@basis
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.basis<-'
standGen <- function(x, value) standardGeneric('.swale.internal.basis<-')
standMethod <- function(x, value) {x@basis<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.internal'
funcname <-'.swale.internal.eeg.data'
standGen <- function(object) standardGeneric('.swale.internal.eeg.data')
standMethod <- function(object) object@eeg.data
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.internal.eeg.data<-'
standGen <- function(x, value) standardGeneric('.swale.internal.eeg.data<-')
standMethod <- function(x, value) {x@eeg.data<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.waveform'
standGen <- function(object) standardGeneric('.swale.solution.waveform')
standMethod <- function(object) object@waveform
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.waveform<-'
standGen <- function(x, value) standardGeneric('.swale.solution.waveform<-')
standMethod <- function(x, value) {x@waveform<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.derivwave'
standGen <- function(object) standardGeneric('.swale.solution.derivwave')
standMethod <- function(object) object@derivwave
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.derivwave<-'
standGen <- function(x, value) standardGeneric('.swale.solution.derivwave<-')
standMethod <- function(x, value) {x@derivwave<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.amplitude'
standGen <- function(object) standardGeneric('.swale.solution.amplitude')
standMethod <- function(object) object@amplitude
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.amplitude<-'
standGen <- function(x, value) standardGeneric('.swale.solution.amplitude<-')
standMethod <- function(x, value) {x@amplitude<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.latency'
standGen <- function(object) standardGeneric('.swale.solution.latency')
standMethod <- function(object) object@latency
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.latency<-'
standGen <- function(x, value) standardGeneric('.swale.solution.latency<-')
standMethod <- function(x, value) {x@latency<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.model'
standGen <- function(object) standardGeneric('.swale.solution.model')
standMethod <- function(object) object@model
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.model<-'
standGen <- function(x, value) standardGeneric('.swale.solution.model<-')
standMethod <- function(x, value) {x@model<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.latencyRange'
standGen <- function(object) standardGeneric('.swale.solution.latencyRange')
standMethod <- function(object) object@latencyRange
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.latencyRange<-'
standGen <- function(x, value) standardGeneric('.swale.solution.latencyRange<-')
standMethod <- function(x, value) {x@latencyRange<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.pp.latency'
standGen <- function(object) standardGeneric('.swale.solution.pp.latency')
standMethod <- function(object) object@pp.latency
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.pp.latency<-'
standGen <- function(x, value) standardGeneric('.swale.solution.pp.latency<-')
standMethod <- function(x, value) {x@pp.latency<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.pp.amplitude'
standGen <- function(object) standardGeneric('.swale.solution.pp.amplitude')
standMethod <- function(object) object@pp.amplitude
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.pp.amplitude<-'
standGen <- function(x, value) standardGeneric('.swale.solution.pp.amplitude<-')
standMethod <- function(x, value) {x@pp.amplitude<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.aic'
standGen <- function(object) standardGeneric('.swale.solution.aic')
standMethod <- function(object) object@aic
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.aic<-'
standGen <- function(x, value) standardGeneric('.swale.solution.aic<-')
standMethod <- function(x, value) {x@aic<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.discard'
standGen <- function(object) standardGeneric('.swale.solution.discard')
standMethod <- function(object) object@discard
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.discard<-'
standGen <- function(x, value) standardGeneric('.swale.solution.discard<-')
standMethod <- function(x, value) {x@discard<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.internal'
standGen <- function(object) standardGeneric('.swale.solution.internal')
standMethod <- function(object) object@internal
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.internal<-'
standGen <- function(x, value) standardGeneric('.swale.solution.internal<-')
standMethod <- function(x, value) {x@internal<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.control'
standGen <- function(object) standardGeneric('.swale.solution.control')
standMethod <- function(object) object@control
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.control<-'
standGen <- function(x, value) standardGeneric('.swale.solution.control<-')
standMethod <- function(x, value) {x@control<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'swale.solution'
funcname <-'.swale.solution.version'
standGen <- function(object) standardGeneric('.swale.solution.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.swale.solution.version<-'
standGen <- function(x, value) standardGeneric('.swale.solution.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.iter.limit'
standGen <- function(object) standardGeneric('.control.iter.limit')
standMethod <- function(object) object@iter.limit
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.iter.limit<-'
standGen <- function(x, value) standardGeneric('.control.iter.limit<-')
standMethod <- function(x, value) {x@iter.limit<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.iter.convergence'
standGen <- function(object) standardGeneric('.control.iter.convergence')
standMethod <- function(object) object@iter.convergence
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.iter.convergence<-'
standGen <- function(x, value) standardGeneric('.control.iter.convergence<-')
standMethod <- function(x, value) {x@iter.convergence<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.split.type'
standGen <- function(object) standardGeneric('.control.split.type')
standMethod <- function(object) object@split.type
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.split.type<-'
standGen <- function(x, value) standardGeneric('.control.split.type<-')
standMethod <- function(x, value) {x@split.type<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.peak.windows'
standGen <- function(object) standardGeneric('.control.peak.windows')
standMethod <- function(object) object@peak.windows
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.peak.windows<-'
standGen <- function(x, value) standardGeneric('.control.peak.windows<-')
standMethod <- function(x, value) {x@peak.windows<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.peak.method'
standGen <- function(object) standardGeneric('.control.peak.method')
standMethod <- function(object) object@peak.method
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.peak.method<-'
standGen <- function(x, value) standardGeneric('.control.peak.method<-')
standMethod <- function(x, value) {x@peak.method<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.disc.edge'
standGen <- function(object) standardGeneric('.control.disc.edge')
standMethod <- function(object) object@disc.edge
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.disc.edge<-'
standGen <- function(x, value) standardGeneric('.control.disc.edge<-')
standMethod <- function(x, value) {x@disc.edge<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.max.lat'
standGen <- function(object) standardGeneric('.control.max.lat')
standMethod <- function(object) object@max.lat
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.max.lat<-'
standGen <- function(x, value) standardGeneric('.control.max.lat<-')
standMethod <- function(x, value) {x@max.lat<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.start.value'
standGen <- function(object) standardGeneric('.control.start.value')
standMethod <- function(object) object@start.value
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.start.value<-'
standGen <- function(x, value) standardGeneric('.control.start.value<-')
standMethod <- function(x, value) {x@start.value<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.output'
standGen <- function(object) standardGeneric('.control.output')
standMethod <- function(object) object@output
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.output<-'
standGen <- function(x, value) standardGeneric('.control.output<-')
standMethod <- function(x, value) {x@output<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
classname <-'control'
funcname <-'.control.version'
standGen <- function(object) standardGeneric('.control.version')
standMethod <- function(object) object@version
setGeneric(funcname,standGen,package='swale')
setMethod(funcname,classname,standMethod)
slotreplace <-'.control.version<-'
standGen <- function(x, value) standardGeneric('.control.version<-')
standMethod <- function(x, value) {x@version<- value;x}
setGeneric(slotreplace,standGen)
setReplaceMethod(funcname,classname,standMethod)
