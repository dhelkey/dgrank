toDG = function(coef_mat, count_vec, baseline = FALSE){
    #' Convert MCMC Coefficients to Z-Scores Following Draper-Gittoes
	#'
	#' \code{toDG} reparamatrizes a matrix of MCMC coefficients into Draper-Gittoes
    #' Z-Scores with a linear transformation. The dot product of the Z-scores and
	#' the count_vec is constrained to be 0.
	#'
	#' @param coef_mat <IXp> coefficient matrix, with I MCMCM iterations.
	#'    The first column should be the overall intercept.
	#' @param count_vec <px1> vector of integers, total number of observations
	#'    per institution or group.
	#' @param baseline Logical; If TRUE, the first coeffient is constrained to be
	#'     0, and the D-G constraint is not satisfied. 

    p = dim(coef_mat)[2]
	n = dim(coef_mat)[1]
	
    #Construct x-form mat
	if (!baseline){
		coef_mat[ ,1] = 0 #Reparametrize
	    xform_mat=diag(1,p)-replicate(p,count_vec)/sum(count_vec)
		coef_mat_transform = coef_mat %*% xform_mat
	} else {
		coef_mat_transform = as.matrix(coef_mat[ ,-1])
	}
	p_temp = dim(coef_mat_transform)[2]
    #Return z-scores 
    z_mat = coef_mat_transform %*% 
		diag(sqrt(1 / diag( cov(coef_mat_transform))), ncol = p_temp)
	if (baseline){z_mat = cbind( rep(0, n), z_mat)}	
    return(z_mat)
}