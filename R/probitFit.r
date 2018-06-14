probitFit = function(y, X,  prior_var_vec, iters = 100, diagnostic = TRUE){
	#Generate MCMCM samples of Probit regression, following
	#Albert and Chib (1993)
	#Rquires truncnorm package to sample from truncated normal	
	#Assumes prior centered at 0, and diagonal prior covariance structure
	#Assumes response y is coded as 0-1

	#Number of coefficients
	n = dim(X)[1]
	p = dim(X)[2]
	B_inv = diag(1/prior_var_vec, ncol = p) #inverse prior variance mat
	
	#Truncation limits, 1 row per observation, colums are lower and upper limits.
	trunc_limits = matrix(0, nrow = n, ncol = 2) 
	trunc_limits[y==0 ,1] = -Inf
	trunc_limits[y==1, 2] = Inf
		

	
	#Precompute Cholesky decomposition (to efficiently sample from multivariate normal)
	XtX = Matrix::crossprod(X)
	Xty = Matrix::crossprod(X,y)
	A = XtX + B_inv
	R = chol(A) #Cholesky decomposition: A = t(R) %*% R

	#Start MCMC algorithm at Least-Squares estimate (guarenteed to exist w/ nonzero prior variance)
	beta_ls = linCholSolver(R, Xty)

	#Generate the normal rancom variates outside of loop
	rnorm_mat = matrix( rnorm(p * iters, mean = 0, sd = 1), 
				ncol = p)
				
	#Helper function to sample Z|... from truncated normal
	sampleZ = function(beta_vec){
		#Samples z given beta from truncated normal distribtion
		truncnorm::rtruncnorm(n, trunc_limits[ ,1], trunc_limits[ ,2], 
							mean =  as.numeric(X %*% beta_vec))							
	}
	
	#MCMC storage, each row is an iteration
	beta_mat = matrix(0, nrow = iters, ncol = p )
	beta_mat[1, ] = as.numeric(beta_ls)
	z_latent = sampleZ(beta_mat[1, ])
	for (i in 2:iters){
		#Gibbs sample beta given z
		Xtz = Matrix::crossprod(X, z_latent)
		beta_z = linCholSolver(R,  Xtz)
		beta_mat[i, ] =  as.numeric(beta_z +
			backsolve(R, rnorm_mat[i, ]))

		#Gibbs sample z
		z_latent = sampleZ(beta_mat[i, ])
	}
	return(beta_mat)
}