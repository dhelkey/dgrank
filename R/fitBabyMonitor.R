fitBabyMonitor = function(minimal_data, num_cat, num_cont, 
						var_intercept,var_inst, var_cat, var_cat_2way, var_cont,
						var_subset = 1,	var_inst_subset_2way = 1,
						subset = FALSE, sparse = TRUE, visualize = FALSE,
                         fit_method = 'probit', burn_in = 100,
                          iters = 1000,  alpha = 0.1,  alpha_star = 0.01,
						verbose = TRUE){					  
	#Inputs:
	#		minimal_data: 
	#		num_cat: Scalar number of categorical variables 
	#		num_cont: Number of continious variables
	#		var_intercept, var_inst, var_cat, var_cat_2way, var_cont, var_subset, var_inst_subset_2way: Prior variances for each type of coefficient
	#		subset: Is analysis to be conducted with a subset variable? (TODO, doesn't work yet)
	#		sparse: Store w/ sparce matrices, requires Matrix package
	#		visualize: (Boolean) If TRUE, display plot
	#		fit_method: Which fit method to use, 'probit' or 'bayesianregression'
	#		burn_in: How many inital MCMC iterations to discard
	# 		iters: Number of desired MCMC iterations to be output
	#		verbose: If TRUE, display diagosnitic messages
	#		alpha: Constructs (1-alpha)% posterior intervals
	#		alpha_star: Flags instituions with (1-alpha_star)% that exclude zero. (TODO, doesn't do anything yet)
	#	
	#Outputs:
	#		inst_mat: Matrix (one row per institution) computing summary statistics, D-G instituiton ranking, and intervals
	#		bayesian_z_ests: D-G instituion ranking (computed in inst_mat, returning here as well for convienence)
	# 		group_labels: A vector of the names of each institution
	#		mcmc_fit: Matrix of MCMC iterations for each coefficient
	#		dg_z: Matrix of commputed z score for each instituion at each MCMC iterations+
	#		coefs: Names of each coefficient (1st is intercept, then institution, then everything else)
	#		prior_var_vec: Vector of prior variances for each coefficient
	#		model_matrix: Design matrix 
	#		dat: List that stores a few variables used by this function, e.g. number of instituions and number of individuals
	
	##Process inputs
	N_full =  dim(minimal_data)[1]
	minimal_data = minimal_data[complete.cases(minimal_data), ]
	
	#Extract outcomes and institutions
	outcome_vec = as.numeric(as.character(minimal_data[ ,1])) #I don't like this
	inst_vec = as.factor(minimal_data[ ,2])

	#If subsetting by a variable (e.g. race)
	var_start_index = 3
	subset_vec = NULL
	if (subset == TRUE){
        subset_vec = as.factor(minimal_data[ ,3])
        var_start_index = 4
    }
    #Categorical risk adjusters
	cat_var_mat = NULL
	if (num_cat > 0){
	    cat_var_mat = as.data.frame(minimal_data[  ,(var_start_index):(			var_start_index + num_cat - 1)])
      #Explicitly turn into factor      
  for (i in 1:num_cat){
       cat_var_mat[ ,i] = as.factor(cat_var_mat[ ,i] )
				}
	}
#Continious risk adjusters
    cont_var_mat = NULL
    if (num_cont > 0){
        cont_var_mat = as.matrix(minimal_data[  ,(var_start_index + num_cat):(var_start_index + num_cat + num_cont - 1)])
     
    }
	
	#Compute institution level summary statistics
    unique_inst_vec = sort(unique(inst_vec))
    inst_mat = data.frame(id = unique_inst_vec,
                         n = sapply(unique_inst_vec, function(id) sum(inst_vec == id)),
                             o_rate = sapply(unique_inst_vec,
                                             function(id) sum( outcome_vec[inst_vec == id]) /  sum(inst_vec == id)))
	
	#Compute a few variables
	N = dim(minimal_data)[1]
	p = length(unique_inst_vec)
	indicator_name = names(minimal_data)[1]
	o_overall = mean(outcome_vec)
	 
	 
	#Send messages
	if (verbose){
			message('Scoring: ', indicator_name ,' (overall rate: ', signif(o_overall,digits = 2), ')')
		}	
	if (N != N_full){message('Warning: ', N_full - N, ' individual records removed due to missing data' ) }
	if (verbose){
		message('Processed data for ', N, ' individuals within ', p, ' institutions.')
		message('Categorical Risk-adjusters: ', paste(colnames(cat_var_mat), collapse = ','))
		message('Continious Risk-adjusters: ', paste(colnames(cont_var_mat), collapse = ','))
	}
	##Construct design matrix additivly and save certain kinds of coefficients
	#Main variables of interest (group, subgroup)
	if (subset){
			#intercept, inst coefficients, subset coefficients, interaction coefficients
			main_mat = modelMatrix( data.frame(inst_vec, subset_vec), 
								   interactions = TRUE, intercept = TRUE, sparse = sparse)
			temp_coefs = colnames(main_mat)
			subset_coefs = temp_coefs[(p + 1):(p + length(levels(subset_vec)) - 1)]
			inst_subset_2way_coefs = temp_coefs[grep(":", temp_coefs)]   
	}   else if (!subset){
			main_mat =  modelMatrix( inst_vec,
									intercept = TRUE, sparse = sparse)
			subset_coefs = inst_subset_2way_coefs = NULL
	}
	all_inst_coefs = colnames(main_mat)[1:p] #1st column is intercept
	inst_coefs = all_inst_coefs[-1]  #Institution 2 to n

	#Categorical risk adjusters

	cat_mat = NULL
	if (num_cat > 0){
		cat_mat =  modelMatrix(cat_var_mat, 
						   interactions = TRUE, sparse = sparse) #Meow
	}
	cat_vars_coefs = colnames(cat_mat)
	cat_vars_2way_coefs = cat_vars_coefs[grep(":", cat_vars_coefs)]
	cat_vars_1way_coefs = setdiff(cat_vars_coefs, cat_vars_2way_coefs)

	#Continious risk adjusters
	cont_mat = NULL
	if (num_cont > 0){
		cont_mat =  modelMatrix(cont_var_mat, sparse = sparse)
	}
	cont_vars_coefs = colnames(cont_mat)

	
	#Construct design matrix
	model_matrix = Matrix(cbind(main_mat, cat_mat, cont_mat), sparse = sparse)  
	coefs = colnames(model_matrix)
	num_coefs = dim(model_matrix)[2]
	##Construct prior variance vector
	prior_var_vec = rep(0, num_coefs)
	prior_var_vec[1] = var_intercept
	prior_var_vec[coefs %in% inst_coefs] = var_inst
	prior_var_vec[coefs %in% subset_coefs] = var_subset
	prior_var_vec[coefs %in% inst_subset_2way_coefs] = var_inst_subset_2way
	prior_var_vec[coefs %in% cat_vars_coefs] = var_cat
	prior_var_vec[coefs %in% cat_vars_2way_coefs] = var_cat_2way
	prior_var_vec[coefs %in% cont_vars_coefs] = var_cont

	
	##Return matrices for null simulations
	no_inst_mat = cbind( rep(1,N),
							 cat_mat, cont_mat)
	no_inst_var_vec = prior_var_vec[!(coefs %in% inst_coefs)]
	
	if (verbose){
	message('Prior variances: intercept: ',var_intercept, ', institution: ',var_inst, 
	', categorical: ',var_cat,'\n, categorical_2way: ', var_cat_2way, ', continious: ', var_cont,
	', subset: ', var_subset, ', institution_subset_2way: ', var_inst_subset_2way)
	}
	#if (!subset){message('Variance values for subset will be ignored')}
	fit_Fun = switch(fit_method, 
					 'probit' = probitFit,
	                'bayesianregression' = bayesianRegression)
	mcmc_fit = fit_Fun(outcome_vec, model_matrix, prior_var_vec, iters = burn_in + iters)
	mcmc_fit = mcmc_fit[ -(1:burn_in),  ]
	colnames(mcmc_fit) = coefs
	#Compute institution level posterior point estimate and posterior credible intervals
	inst_fit = mcmc_fit[ , all_inst_coefs] 
	#Apply D-G transformation
	inst_fit_dg_z = toDG(inst_fit, count_vec = inst_mat$n)
	inst_mat$post_est = apply(inst_fit_dg_z, 2, median)
	inst_mat[ , c('lower', 'upper')] = t(apply(inst_fit_dg_z, 2, quantile, probs = c(alpha, 1- alpha)))
	
	if (visualize){
		inst_order = order(inst_mat$post_est, decreasing = TRUE)
		plot(1:p, inst_mat$post_est[inst_order])
		segments(1:p, inst_mat$lower[inst_order], y1 = inst_mat$upper[inst_order])
	}
	
	#Compute DG nonbaseline subset coding (e.g. for ethnicity)
	#TODO: In progress
	
	#COmpute DG baseline subset coding (e.g. for time)
	#TODO: In progress
		
	returner = list(
				coefs = coefs,
				prior_var_vec = prior_var_vec,
				model_matrix = model_matrix,
				mcmc_fit = mcmc_fit, 
				post_mean_vec = apply(mcmc_fit, 2, mean),
				inst_mat = inst_mat,
				dg_z = inst_fit_dg_z,
				iters = iters,
				p = p,
				bayesian_z_ests = inst_mat$post_est,
				group_labels = as.vector(inst_mat$id),
				no_inst_mat = no_inst_mat, ##For null simulations
				no_inst_var_vec = no_inst_var_vec,
				outcome_vec = outcome_vec,
				all_inst_coefs = all_inst_coefs)
	if (verbose) {message('Fit Complete')}
	return(returner)
}

	
	
	
