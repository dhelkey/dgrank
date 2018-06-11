toDG = function(coef_mat, count_vec, baseline = FALSE){
    #Reparameteriz coefficients to satisfy D-G constraint and convert to Z-scores
    #Coefficients should be baseline encoded, i.e. first coefficient is intercept
	#
	#Inputs:
	#	coef_mat - institution coefficient matrix with MCMC iterations as the rows
	#	count_vec - Number of observations per institution
	#   baseline - If TRUE, the first coefficient is constrained to be zero
	#				and the Draper-Gittoes constraint is not satisfied
    
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