fitBabyMonitor = function(minimal_data, num_cat, num_cont, 
						var_intercept,var_inst, var_cat, var_cat_2way, var_cont,
						var_subset = 1,	var_inst_subset_2way = 1,
						subset = FALSE, 
                        fit_method = 'probit', 
						sparse = TRUE, burn_in = 100,
                        iters = 1000,  alpha = 0.05, 
						verbose = TRUE){
	#' Fit Baby-MONITOR for CPQCC/VON
	#'
	#' \code{fitBabyMonitor} comprehensively applies the Baby-MONITOR score
	# to a \code{minimal_data} file. Designed for CPQCC/VON applications.
	#' Returns a large list. Use inst_mat for institution rankings, 
	#' full_subset_mat_baseline and full_subset_mat_nobaseline for subset rankings,
	#' and subset_baseline_mat and subset_nobaseline_mat for subset rankings within institution.
	#'
	#' @param minimal_data Data_frame with a particular format:
	#'
	#' 1st column: Outcome vector (0-1 encoding)
	#'
	#' 2nd column: Institution ID
	#'
	#' 3rd column (if subset == TRUE): 
	#'
	#' Next: num_cat columns of categorical variables (num_cat can equal 0)
	#'
	#' Next: num_cont columns of continuous variables (num_cont can equal 0)
	#'
	#' @param num_cat Scalar number of categorical variables. 
	#' @param num_cont Number of continuous variables.
	#' @param var_intercept Prior variance for intercept parameter.
	#' @param var_inst Prior variance for institution parameters.
	#' @param var_cat Prior variance for categorical parameters.
	#' @param var_cat_2way Prior variance for categorical interaction parameters.
	#' @param var_cont Prior variance for continuous parameters.
	#' @param var_subset Prior variance for subset parameters.
	#' @param var_inst_subset_2way Prior variance for subset interaction parameters.
	#' @param iters Number of MCMC iterations to use.
	#' @param burn_in Number of initial iterations to discard for burn in.
	#' @param sparse Should design_matrix be stored as Sparse matrix? Requires Matrix package.
	#' @param alpha We look at posterior (1-alpha)\% posterior intervals. 
	#' @param verbose If TRUE, display status messages while fitting
	#' @param subset If TRUE, perform subset fitting tasks
	#' @param fit_method Analysis method. Options: 'probit' for Probit Regression, 'bayesianregression' for 
	#' Bayesian linear regression, and 'probitk' for a Probit Regression implementation written by Daniel Kirsner.
	#'
	#' @return
	#'	Returns a large list with the following components:
	#'
	#'		inst_mat: Matrix (one row per institution) computing summary statistics, D-G institution ranking, and intervals
	#'
	#' 	     full_subset_mat_baseline, full_subset_mat_nobaseline: Matrix with rankings and intervals for the various subset categories.
	#'
	#'
	#'	    subset_baseline_mat, subset_nobaseline_mat: Matrix with rankings and intervals for subset categories within institution. 
	#'
	#'
	#'		group_labels: A vector of the names of each institution
	#'
	#'	    mcmc_fit: Matrix of MCMC iterations for each coefficient
	#'		dg_z: Matrix of computed z score for each institution at each MCMC iterations.
	#'
	#'		coefs: Names of each coefficient (1st is intercept, then institution, then everything else)
	#'
	#'		prior_var_vec: Vector of prior variances for each coefficient
	#'
	#'		model_matrix: Design matrix 
	
	##Sort in order of increasing inst_id (2nd column of minimal_data)
	minimal_data = minimal_data[order(minimal_data[ ,2]), ]

	##Process inputs
	N_full =  dim(minimal_data)[1]
	minimal_data = minimal_data[complete.cases(minimal_data), ]
	
	#Extract outcomes and institutions
	outcome_vec = as.numeric(as.character(minimal_data[ ,1]))
	inst_vec = as.factor(minimal_data[ ,2])

	#If subsetting by a variable (e.g. race)
	var_start_index = 3
	subset_vec = NULL
	if (subset == TRUE){
		subset_vec = as.numeric(as.character(minimal_data[ ,3]))
		subset_vec[is.nan(subset_vec) | is.na(subset_vec)]= 99 #Missing subset category code
        subset_vec = as.factor(subset_vec)
        var_start_index = 4
    }
	unique_subset_vec = NULL
	if (subset){unique_subset_vec = sort(unique(subset_vec)) }
	
    #Extract categorical risk adjusters
	cat_var_mat = NULL
	if (num_cat > 0){
		cat_var_locat = var_start_index:(var_start_index + num_cat - 1)
	    cat_var_mat = as.data.frame(minimal_data[  ,cat_var_locat])
		colnames(cat_var_mat) = names(minimal_data)[cat_var_locat]
		#Explicitly turn into factor      
		for (i in 1:num_cat){
		   cat_var_mat[ ,i] = as.factor(cat_var_mat[ ,i] )
		}
	}
	
	#Extract continuous risk adjusters
    cont_var_mat = NULL
    if (num_cont > 0){
		cont_var_locat = (var_start_index + num_cat):(var_start_index + num_cat + num_cont - 1)
        cont_var_mat = as.matrix(minimal_data[  ,cont_var_locat])
		colnames(cont_var_mat) = names(minimal_data)[cont_var_locat]
    }
	
	#Compute institution level summary statistics
    unique_inst_vec = sort(unique(inst_vec))
    inst_mat = data.frame(id = unique_inst_vec,
                         n = sapply(unique_inst_vec, function(id) sum(inst_vec == id)),
                             o_rate = sapply(unique_inst_vec,
                                             function(id) sum( outcome_vec[inst_vec == id]) /  sum(inst_vec == id)))
	
	#Compute summary statistics
	N = dim(minimal_data)[1]
	p = length(unique_inst_vec)
	indicator_name = names(minimal_data)[1]
	o_overall = mean(outcome_vec)
	 
	#Display messages
	if (verbose){
			message('Scoring: ', indicator_name ,' (overall rate: ', signif(o_overall,digits = 2), ')')
		}	
	if (N != N_full){message('Warning: ', N_full - N, ' individual records removed due to missing data' ) }
	if (verbose){
		message('Processed data for ', N, ' individuals within ', p, ' institutions.')
		message('Categorical Risk-adjusters: ', paste(colnames(cat_var_mat), collapse = ','))
		message('Continuous Risk-adjusters: ', paste(colnames(cont_var_mat), collapse = ','))
	}

	##Construct design matrix additivelyand save certain kinds of coefficients
	#Main variables of interest (group, subgroup)
	if (subset){
			#intercept, inst coefficients, subset coefficients, interaction coefficients
			main_mat = modelMatrix( data.frame(inst_vec, subset_vec), 
								   interactions = TRUE, intercept = TRUE, sparse = sparse)
			temp_coefs = colnames(main_mat)
			subset_coefs = temp_coefs[(p + 1):(p + length(levels(subset_vec)) - 1)] #Indexing based on how main_mat is constructed above
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

	#Continuous risk adjusters
	cont_mat = NULL
	if (num_cont > 0){
		cont_mat =  modelMatrix(cont_var_mat, sparse = sparse)
	}
	cont_vars_coefs = colnames(cont_mat)
	
	#Build model additively
	model_matrix = Matrix(cbind(main_mat, cat_mat, cont_mat), sparse = sparse)  
	
	#Save location of coefficients
	coefs = colnames(model_matrix)
	num_coefs = dim(model_matrix)[2]
	
	#Integer location indices of commonly used coefficients
	locat_all_inst = which(coefs %in% all_inst_coefs)
	locat_all_subset = c(1, which(coefs %in% subset_coefs))
	
	##Construct prior variance vector
	prior_var_vec = rep(0, num_coefs)
	prior_var_vec[1] = var_intercept
	prior_var_vec[coefs %in% inst_coefs] = var_inst
	prior_var_vec[coefs %in% subset_coefs] = var_subset
	prior_var_vec[coefs %in% inst_subset_2way_coefs] = var_inst_subset_2way
	prior_var_vec[coefs %in% cat_vars_coefs] = var_cat
	prior_var_vec[coefs %in% cat_vars_2way_coefs] = var_cat_2way
	prior_var_vec[coefs %in% cont_vars_coefs] = var_cont

	##Return design matrices w/o inst. indicators for null simulations
	no_inst_mat = cbind( rep(1,N),
				cat_mat, cont_mat)
	no_inst_var_vec = prior_var_vec[!(coefs %in% inst_coefs)]
	
	if (verbose){
	message('Prior variances: intercept: ',var_intercept, ', institution: ',var_inst, 
	', categorical: ',var_cat,'\n, categorical_2way: ', var_cat_2way, ', continious: ', var_cont,
	', subset: ', var_subset, ', institution_subset_2way: ', var_inst_subset_2way)
	if (!subset){message('Variance values for subset will be ignored')}
	}

	#Fit model with MCMC inference and link function of choice
	fit_Fun = switch(fit_method, 
			'probit' = probitFit,
			'probitk' = probitFitk,
	        'bayesianregression' = bayesianRegression)
	mcmc_fit = fit_Fun(outcome_vec, model_matrix, prior_var_vec, 
					iters = burn_in + iters)
					
	mcmc_fit = mcmc_fit[ -(1:burn_in),  ] #Remove burn-in
	colnames(mcmc_fit) = coefs
	
	
	#Compute institution level posterior point estimate and posterior credible intervals
	inst_fit = mcmc_fit[ ,all_inst_coefs] 
	#Apply D-G transformation
	inst_fit_dg_z = toDG(inst_fit, count_vec = inst_mat$n)
	#Add posterior credible intervals
	inst_mat = cbind(inst_mat,
			posteriorCredible(inst_fit_dg_z, alpha = alpha))
	
	#Compute baseline and non baseline ranks for each subset category
	full_subset_mat_baseline = full_subset_mat_nobaseline = data.frame(
											id = unique_subset_vec,
											n = sapply(unique_subset_vec,
												function(s) sum(subset_vec == s)),
											o_rate = sapply(unique_subset_vec,
												function(s) mean(outcome_vec[subset_vec == s])))

	subset_fit_baseline = subset_fit_nobaseline = NULL
	if (subset){
		##Separate fitting for , just look at subset coefficients w/ risk adjusters
		model_matrix_subset = modelMatrix(subset_vec, intercept = TRUE)
		locat_just_subset = 1:dim(model_matrix_subset)[2] 
		model_matrix_subset = cbind(model_matrix_subset, cat_mat, cont_mat)

		#Build variance vector
		coefs_subset = colnames(model_matrix_subset)
		prior_var_vec_subset = rep(var_subset, length(coefs_subset)) #Subset varince
		prior_var_vec_subset[coefs_subset %in% cat_vars_coefs] = var_cat #Categorical variance
		prior_var_vec_subset[coefs_subset %in% cat_vars_2way_coefs] = var_cat_2way
		prior_var_vec_subset[coefs_subset %in% cont_vars_coefs] = var_cont
		prior_var_vec_subset[1] = var_intercept #intercept

		#Fit the model w/o institutions
		mcmc_fit_subset = fit_Fun(outcome_vec, model_matrix_subset, prior_var_vec_subset,
				iters = burn_in + iters)
		mcmc_fit_subset = mcmc_fit_subset[-(1:burn_in), ]
		colnames(mcmc_fit_subset) = coefs_subset
		
		subset_fit_baseline = toDG(mcmc_fit_subset[ ,locat_just_subset], 
								full_subset_mat_baseline$n, baseline = TRUE)
		subset_fit_nobaseline = toDG(mcmc_fit_subset[ ,locat_just_subset], 
									full_subset_mat_nobaseline$n)
		full_subset_mat_baseline = cbind( full_subset_mat_baseline,
							posteriorCredible(subset_fit_baseline))
		full_subset_mat_nobaseline = cbind( full_subset_mat_nobaseline,
							posteriorCredible(subset_fit_nobaseline))
	}

	#Compute DG nonbaseline subset coding (e.g. for ethnicity)
	subset_list = vector("list", length = p)
	subset_nobaseline_mat = subset_baseline_mat = data.frame(
							inst = NULL, 
							subset_cat = NULL,
							post_est = NULL, 
							lower = NULL, upper = NULL)
	inst_subset_nobaseline_mat = inst_subset_baseline_mat = NULL
	
	
	if (subset){
		# #Count how many of each subset at each instituton
		##NOTE: Any given institution may have not have all subset categories
		subset_table = table( data.frame(inst_vec, subset_vec))
		subset_categories = colnames(subset_table)
	
	##Extract subset variable coefficients for each Inst
		#First institution uses base coefficients due to one-hot encoding.
		mcmc_subset_mat_base = mcmc_fit[ ,locat_all_subset]
		inst_subset_indices = subset_table[1, ] > 0	
		subset_list[[1]]$mcmc_coefs = mcmc_subset_mat_base[ ,inst_subset_indices]
		for (i in 2:p){
			#Locate the correct interaction coefficients
			coef_locat = grep( 
				paste0('inst_vec', inst_mat$id[i],':'),
				coefs) 
			inst_subset_indices = subset_table[i, ] > 0
			
			mcmc_subset_mat_inst =  mcmc_fit[ ,coef_locat]
			
			#Only add intercept if institution has 1st category of subset variable
			if (inst_subset_indices[1]){
				mcmc_subset_mat_inst = cbind(mcmc_fit[ ,1],
												mcmc_subset_mat_inst)}
			#Because of one-hot encoding, 2nd-nth institutions get their coefficients added to base coefficients 
			subset_list[[i]]$mcmc_coefs = mcmc_subset_mat_base[ ,inst_subset_indices]
									+ mcmc_subset_mat_inst
		}
		
		#Baseline and non-baseline encoding for each institution
		for (i in 1:p){
			inst_subset_indices = subset_table[i, ] > 0
			
			#A special case is when there is only one subset category in an institution
			if (sum(inst_subset_indices) == 1){
				subset_list[[i]]$dg_z = subset_list[[i]]$baseline_z = 
								matrix(subset_list[[i]]$mcmc_coefs,
													ncol = 1)
			} else {
				subset_list[[i]]$dg_z = toDG(subset_list[[i]]$mcmc_coefs,
						subset_table[i, inst_subset_indices])
				subset_list[[i]]$baseline_z = toDG(subset_list[[i]]$mcmc_coefs,
						subset_table[i, inst_subset_indices], baseline = TRUE)
			}
			#Now, compute posterior estimate and credible intervals
			subset_mat_temp = data.frame(
				inst = rep(inst_mat$id[i], sum(inst_subset_indices)),
				subset_cat = subset_categories[inst_subset_indices])
			
			nobaseline_temp = cbind(subset_mat_temp,
				posteriorCredible(subset_list[[i]]$dg_z,alpha = alpha))
				
			baseline_temp = cbind(subset_mat_temp, 
					posteriorCredible(subset_list[[i]]$baseline_z, alpha = alpha))
						
			subset_list[[i]]$nonbaseline_mat = nobaseline_temp
			subset_list[[i]]$baseline_mat = baseline_temp
			
			subset_baseline_mat = rbind(subset_baseline_mat, baseline_temp)
			subset_nobaseline_mat = rbind(subset_nobaseline_mat, nobaseline_temp)
			
			#TODO stop storing the intermediate....get rid of subset_list 
			inst_subset_baseline_mat = cbind(inst_subset_baseline_mat, subset_list[[i]]$baseline_z)
			inst_subset_nobaseline_mat = cbind(inst_subset_nobaseline_mat, subset_list[[i]]$dg_z)
		}
	}
	
	##Add Draper Gittoes scoring
	if (!is.null(dim(cat_mat))){
		pcf_cat_vec = apply(cat_mat,1, paste, collapse = ',') #Unique categorical PCF 'categories'
	} else {
		pcf_cat_vec = rep(1, N_full)
	}

	design_based_inst = designBased(0.5, outcome_vec, pcf_cat_vec, inst_vec)
	if (subset){design_based_subset = designBased(0.5, outcome_vec, pcf_cat_vec, subset_vec)}
	
	#Add to existing matrices
	inst_mat$dg_score = design_based_inst$Z
	if (subset){
		full_subset_mat_baseline$dg_score = full_subset_mat_nobaseline$dg_score = design_based_subset$Z
	}
	returner = list(
				outcome_vec = outcome_vec,
				subset = subset,
				inst_vec = inst_vec,
				coefs = coefs,
				prior_var_vec = prior_var_vec,
				model_matrix = model_matrix,
				mcmc_fit = mcmc_fit, 
				post_mean_vec = apply(mcmc_fit, 2, mean),
				inst_mat = inst_mat,
				subset_nobaseline_mat = subset_nobaseline_mat,
				subset_baseline_mat = subset_baseline_mat,
				dg_z = inst_fit_dg_z,
				iters = iters,
				p = p,
				bayesian_z_ests = inst_mat$post_est,
				group_labels = as.vector(inst_mat$id),
				no_inst_mat = no_inst_mat, ##For null simulations
				no_inst_var_vec = no_inst_var_vec,
				outcome_vec = outcome_vec,
				all_inst_coefs = all_inst_coefs,
				subset_list = subset_list,
				full_subset_mat_baseline = full_subset_mat_baseline,
				full_subset_mat_nobaseline = full_subset_mat_nobaseline,
				subset_fit_baseline = subset_fit_baseline,
				subset_fit_nobaseline = subset_fit_nobaseline,
				unique_subset_vec = unique_subset_vec,
				inst_subset_baseline_mat = inst_subset_baseline_mat,
				inst_subset_nobaseline_mat = inst_subset_nobaseline_mat)
	if (verbose) {message('Fit Complete')}
	return(returner)
}

	
	
	


