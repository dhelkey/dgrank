bayesianRegression = function(y, X,  prior_var_vec, iters = 100, sigmasq=1){
	#Requires truncnorm package
    #Assumes prior mean centered at 0
    #TODO prior mean and proior cov need input checks
    ##Written by Daniel Kirsner
    n = length(y)
	p = dim(X)[2]  	
	stopifnot(length(prior_var_vec) == p)
	
	#Data augmentation
	X_aug = as.matrix(rbind(X, diag(p)))
	y_aug = c(y, rep(0, p))
	y_aug = as.numeric(y_aug)
	
	bayesian_model=lm(y_aug~0+X_aug) #, weights = 1/prior_var_vec) 
	mcmc_mat=mvrnorm(iters,coef(bayesian_model), vcov(bayesian_model))
    return(mcmc_mat)
}