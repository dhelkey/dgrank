probitFitk = function(y, X,  prior_var_vec, iters = 100, sigmasq=1){
	#' MCMC Samples for Probit Regression
	#'
	#' \code{probitFit} generates MCMCM samples of Probit regression, following
	#' the data-augmentation Gibbs sampler of Albert and Chib (1993).
	#' Written by Daniel Kirsner
	#' Requires truncnorm package to sample from truncated normal.
	#' Assumes prior centered at 0, and diagonal prior covariance structure.
	#' Assumes response y is coded as 0-1.
	#'
	#' @param y <nX1> observed data vector
	#' @param X <nXp> design matrix (takes objects of class matrix or sparce Matrix)
	#' @param prior_var_vec <px1> prior variance vector
	#' @param iters Integer; number of desired output iterations
	#' @param sigmasq TODO
    n = length(y)
    p = dim(X)[2]  

    prior_mean = rep(0, p)
    stopifnot(length(prior_var_vec) == p) #Ensure variance vec is the right length
    #create limits for the gibbs sampler truncated normals
    lower=rep(-Inf,n)*(-(y-1))
    lower[is.nan(lower)]<-0
    upper=rep(Inf,n)*y
    upper[is.nan(upper)]<-0

    #TODO assert correct dimension for mean and var vecs...

    #precompute to make regression faster
    XpX=Matrix::crossprod(X)
    IR=backsolve(
    	chol(XpX/1+diag(1/prior_var_vec, ncol = p)),
	diag(p))
    mcmc_betas = matrix(0, nrow = iters, ncol = p)
    mcmc_latent=truncnorm::rtruncnorm(n,lower,upper, as.numeric(X%*%mcmc_betas[1,])) 
    for(i in 2:iters)
      {
        #first perform regression on latent normals
        Xpy=Matrix::crossprod(X,mcmc_latent)
		b1 = Matrix::crossprod(t(IR))
		b2 = Xpy/sigmasq+diag(1/prior_var_vec, ncol = p)%*%prior_mean
        btilde= b1 %*% b2
        mcmc_betas[i,] = t(replicate(1,as.vector(btilde))) + rnorm(p)%*%t(IR)
        #Conditional on the regression coefficients, generate truncated normals
        mcmc_latent=truncnorm::rtruncnorm(n,lower,upper,as.numeric(X%*%mcmc_betas[i,]))
		
		
      }
    return(mcmc_betas)
}
