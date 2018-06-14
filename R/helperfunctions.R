   pointQuantile = function(x_mat, alpha){
	#Helper function
        #Helper function to obtain estimates and quantiles
        #Input x_mat - (iter X p) matrix
        return(list(point = apply(x_mat, 2, median),
        range = apply(x_mat,2 , quantile, probs = c(alpha, 1-alpha)))
        )
    }
	
posteriorCredible = function(mcmc_mat, alpha = 0.05, est_fun = median){
	#Compute 1-alpha% credible interval and point estimate
	#Input:
	#      mcmc_mat: IxP matrix with iterations stored in rows.
	#Output:
	#      data.frame P rows and three columns (post_est, lower, upper)
	out_frame = data.frame(post_est = apply(mcmc_mat, 2, est_fun))
	out_frame[ ,c('lower', 'upper')] = t(apply(mcmc_mat, 2, quantile,
											probs = c(alpha, 1-alpha)) )                    
	return(out_frame)
}

linCholSolver = function(R, y){
		#' Solve System of Equations w/ Cholesky
		#'
		#' \code{linCholSolver} solves Ax = y for x.
		#' It is first required to take the Cholesky decomposition,
		#' obtaining A = t(R) %*% R. This must be performed outside of this function for speed.
		#' Results should be checked against a known solver, as this function optimizes for speed.
		#'
		#' @param R - Upper triangular matrix of dimension <n x p>; Cholesky decomposition of A.
		#' @param y - Vector of length n.
		q = forwardsolve(t(R), y)
		return(backsolve(R, q))
}