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