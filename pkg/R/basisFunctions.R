#############################################
# basis function FUNCTIONS				    #
# simultaneous waveform and amplitude		#
# Copyright(c) 2009 Wouter D. Weeda			#
# University of Amsterdam					#
#############################################

#CONTAINS
#makePoly
#fdDeriv
#anDerivPoly

makePoly <-
function(degree,points,deriv='FD')
{
	deriv=match.arg(deriv,c('FD','AN'))
	
	basis = new('basis')
	
	.basis.type(basis) = 'poly'
	.basis.num.funcs(basis) = degree
	.basis.num.points(basis) = points
	
	polymat = poly(seq(0,1,len=points),degree=degree,raw=F)
	.basis.matrix(basis) = matrix(polymat,points,degree)
	
	if(deriv=='FD') {.basis.deriv(basis) = fdDeriv(.basis.matrix(basis));.basis.deriv.type(basis)='finitedifference'}
	if(deriv=='AN') {.basis.deriv(basis) = anDerivPoly(.basis.matrix(basis));.basis.deriv.type(basis)='analytical'}
	
	return(basis)
}

#make finite-difference derivatives of a matrix
fdDeriv <- function(mat)
{
	x=diag(nrow(mat))
	diag(x)=1
	diag(x[-1,])=-1
	x[1,1]=-1
	x[1,2]=1
	dmat=x%*%mat

	return(dmat)
}

#make analytical derivatives of a matrix
anDerivPoly <- function(mat,norms)
{
	dmat=mat

	dmat = cbind(rep(0,nrow(mat)),mat)
	dmat = dmat[,-ncol(dmat)]
	
	for(i in 1:ncol(dmat)) {
		dmat[,i]=i*dmat[,i]
	}
	
	return(dmat)
}




