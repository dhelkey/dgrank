toDG = function(coef_mat, count_vec){
    #Reparameteriz coefficients to satisfy D-G constraint and convert to Z-scores
    #Coefficients should be baseline encoded, i.e. first coefficient is intercept
	#
	#Inputs:
	#	coef_mat - institution coefficient matrix with MCMC iterations as the rows
	#	count_vec - Number of observations per institution
	
    coef_mat[ ,1] = 0 #Reparametrize
    p = dim(coef_mat)[2]

    #Construct x-form mat
    xform_mat=diag(1,p)-replicate(p,count_vec)/sum(count_vec)
    coef_mat_transform = coef_mat %*% xform_mat
    
    #Return z-scores 
    z_mat = coef_mat_transform %*% diag(sqrt(1 / diag( cov(coef_mat_transform))))
    return(z_mat)
}


   pointQuantile = function(x_mat, alpha){
	#Helper function
        #Helper function to obtain estimates and quantiles
        #Input x_mat - (iter X p) matrix
        return(list(point = apply(x_mat, 2, median),
        range = apply(x_mat,2 , quantile, probs = c(alpha, 1-alpha)))
        )
    }