toDG = function(coef_mat, count_vec, baseline = FALSE){
    #' Convert MCMC Coefficients to Z-Scores Following Draper-Gittoes
	#'
	#' \code{toDG} reparamatrizes a matrix of MCMC coefficients into Draper-Gittoes
    #' Z-Scores with a linear transformation. 
	#'
	#' @param coef_mat <IXp> coefficient matrix, with I MCMCM iterations.
	#'    The first column should be the overall intercept.
	#' @param count_vec <px1> vector of integers, total number of observations
	#'    per institution or group.
	#' @param baseline Logical; If TRUE, the first coeffient is constrained to be
	#'     0, and the D-G c	onstraint is not satisfied. 

	if (baseline){ #Baseline encoding, first coefficient constrained to be zero.
		coef_mat_transform = as.matrix(coef_mat[ ,-1])
	} else { #No baseline, effect sizes follow D-G constraint
		coef_mat_transform = dgConstraint(coef_mat, count_vec)
	}
	p = dim(coef_mat_transform)[2]
	
	#Convert to z-scores
	z_mat = coef_mat_transform %*% 
		diag(sqrt(1 / diag( cov(coef_mat_transform))), ncol = p)
	if (baseline){z_mat = cbind( 0, z_mat)}	
	return(z_mat)
}

dgConstraint = function(coef_mat, count_vec){
	#' Transform MCMC coefficients to satisfy Draper-Gittoes linear constraint
	#'
	#' Transforms coefficients to satisfy the Draper-Gittoes constraint (dot product of the Z-scores and
	#' the count_vec is constrained to be 0). Written by Daniel Kirsner
	#'
	#' @param coef_mat <IXp> coefficient matrix, with I MCMCM iterations.
	#'    The first column should be the overall intercept.
	#' @param count_vec <px1> vector of integers, total number of observations
	#'    per institution or group.
	
	p = dim(coef_mat)[2]
	n = dim(coef_mat)[1]
	
	if (min(count_vec) <= 0) stop('Must have positive number of observations in each institution')
	
	coef_mat[ ,1] = 0 #Reparametrize
	xform_mat=diag(1,p)-replicate(p,count_vec)/sum(count_vec)
	return(coef_mat %*% xform_mat)
}