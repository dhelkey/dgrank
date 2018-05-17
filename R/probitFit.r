probitFit = function(y, X,  prior_var_vec, iters = 100){
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
	
	sampleZ = function(beta_vec){
		#Samples z given beta 
		rtruncnorm(n, trunc_limits[ ,1], trunc_limits[ ,2], 
							mean =  as.numeric(X %*% beta_vec))
							
	}
	
	#Start MCMC algorithm at Least-Squares estimate (guarenteed to exist w/ nonzero prior variance)
	XtX = crossprod(X)
	Xty = crossprod(X,y)
	beta_ls = solve(B_inv + XtX, Xty)

	#Precompute Cholesky decomposition to sample from multivariate normal
	L_transpose =chol( B_inv + XtX)

	#MCMC storage, each row is an iteration
	beta_mat = matrix(0, nrow = iters, ncol = p )
	beta_mat[1, ] = as.numeric(beta_ls)
	z_latent = sampleZ(beta_mat[1, ])
	for (i in 2:iters){
		#Gibbs sample beta given z
		Xtz = crossprod(X, z_latent)
		beta_z = solve(B_inv + XtX, Xtz)
		std_norm_vec = rnorm(p, mean = 0, sd = 1)
		beta_mat[i, ] =  as.numeric(beta_z + backsolve(L_transpose, std_norm_vec))

		#Gibbs sample z
		z_latent = sampleZ(beta_mat[i, ])
	}
	return(beta_mat)
}