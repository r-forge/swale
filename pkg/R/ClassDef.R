#############################################
# SWALE S4 CLASS DEFINITIONS				#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#version
#solution


## Swale version class (version is set here)
setClass(
	Class='version',
	representation=representation(
		version='numeric',
		build='numeric',
		update='numeric',
		svnrev='numeric'
	),
	prototype=prototype(
		version=1,
		build=3,
		update=5,
		svnrev=37
	)#,
	#package='swale'
)

## SWALE basisfunctions
setClass(
		Class='basis',
		representation=representation(
				type='character',
				deriv.type='character',
				num.funcs='numeric',
				num.points='numeric',
				matrix='matrix',
				deriv='matrix',
				version='ANY'
		),
		prototype=prototype(
				version=new('version')
		)#,
		#package='swale'
)

## EEG/MEG data object
setClass(
		Class='eeg.data',
		representation=representation(
				data='matrix',
				trend='matrix',
				trials='numeric',
				samples='numeric',
				sampRate='numeric',
				channel='character',
				condition='character',
				within='character',
				version='ANY'				
		),
		prototype=prototype(
				version=new('version')
		)#,
		#package='swale'
)

## SWALE internal object
setClass(
		Class='swale.internal',
		representation=representation(
				waves='matrix',
				amps='matrix',
				lats='matrix',
				model='matrix',
				rss='numeric',
				gradient='numeric',
				version='ANY',
				basis='ANY',
				eeg.data='ANY'

		),
		prototype=prototype(
				version=new('version')
		)#,
		#package='swale'
)

## SWALE solution object
setClass(
		Class='swale.solution',
		representation=representation(
				waveform='matrix',
				derivwave='matrix',
				amplitude='matrix',
				latency='matrix',
				model='matrix',
				latencyRange='numeric',
				pp.latency='matrix',
				pp.amplitude='matrix',
				aic='numeric',
				discard='numeric',
				internal='ANY',
				control='ANY',
				version='ANY'
		),
		prototype=prototype(
				version=new('version')
		)#,
		#package='swale'
)

setClass(
		Class='control',
		representation=representation(
				iter.limit='numeric',
				iter.convergence='numeric',
				split.type='character', #how to split
				peak.windows='ANY',       #where to split ROWS = SPLITS (col = start,end)
				peak.method = 'ANY',
				disc.edge = 'ANY',
				max.lat = 'ANY',
				start.value = 'matrix',
				output='logical',
				version='ANY'
		),
		prototype=prototype(
				iter.limit = 500,
				iter.convergence = 1e-6,
				split.type = 'none',
				peak.windows = NULL,
				peak.method =  list('max'),
				disc.edge =  list(15),
				max.lat = list(1.5),
				output = TRUE,
				version=new('version')
		)#,
#package='swale'
)

