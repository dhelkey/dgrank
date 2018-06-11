modelMatrix = function(x, interactions = FALSE, sparse = TRUE, intercept = FALSE){
    #Creates a design matrix from a data frame, using [0-1] indicator variables 
	#Note, names of categorical variables will change slightly 
	#	(n-1 unique values appended to name)
	
    #Factors use treatment encoding ("one-hot" encoding)
	#	 i.e. indicator variables for 2nd through nth level
	if(is.null(x))return(NULL) #Dealing w/ empty data

	if (sum(!sapply(x, is.finite)) != 0){
		stop('Creating model matrix requires all input elements be finite')}
	
	#Set R to use 1-hot encoding
	options(contrasts = rep("contr.treatment", 2))
    model_formula = as.formula( ~.) #All first order Effects
    if (interactions){model_formula = as.formula( ~(.)^2)} #1st order and interactions	 
   mat = model.matrix(model_formula, data = data.frame(x))
	
    if (!intercept){
		#Remove intercept column, making sure names are preserved
		n_vec = colnames(mat)
		mat = data.frame(mat[ ,-1])
		colnames(mat) = n_vec[-1]
	}
    #Drop any empty levels
	n_vec = colnames(mat)
	nonempty_indices = apply(abs(mat),2, sum) != 0
	mat = data.frame(mat[ ,nonempty_indices])
	colnames(mat) = n_vec[nonempty_indices]

    return( as.matrix(mat) )
}
