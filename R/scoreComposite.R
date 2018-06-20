scoreComposite = function(returner_list, alpha = 0.05, type = 'inst'){
#' Composite Scores with Intervals
#'
#'Constructs composite scores as well
#' as (1-alpha)% confidence intervals from scores on multiple performance indicators.
#' Composite scores constructed by taking the average of each score
#'
#' @param returner_list List of lists, each element should be the return value (a list) from fitBabyMonitor
#' @param alpha Constructs (1-alpha)"\%" posterior intervals for the composite score
#' @param type Type of fitting to be done. 'inst' for institution composite scores
#'				'subset_baseline' or 'subset_nobaseline' for subset rankings
#'				'inst_subset_nobaseline' or 'inst_subset_baseline' for subset within institution rankings

##The process is slightly different for institution specific subset ranking
inst_subset = c('inst_subset_nobaseline', 'inst_subset_baseline')
if (type %in% inst_subset){
	
	codeFun = function(x){paste(x, collapse = '-')}
	##First loop through and process
	save_list = list()
	for (i in 1:length(returner_list)){
		if (type == 'inst_subset_nobaseline'){
				unique_id_mat = returner_list[[i]]$subset_nobaseline_mat[ ,c('inst', 'subset_cat')]
				mcmc_mat = returner_list[[i]]$inst_subset_nobaseline_mat
		} else if (type == 'inst_subset_baseline'){
			unique_id_mat = returner_list[[i]]$subset_baseline_mat[ ,c('inst', 'subset_cat')]
			mcmc_mat = returner_list[[i]]$inst_subset_baseline_mat
		}
		unique_id = apply(unique_id_mat, 1, codeFun)
		save_list[[i]] = list( unique_id_mat = unique_id_mat, 
						mcmc_mat = mcmc_mat,
						unique_id = unique_id)
	}

	##Get the overall list of unique ids
	iters = returner_list[[1]]$iters
	
	code_mat = do.call(rbind, lapply(save_list, function(x) x$unique_id_mat))
	code_mat = code_mat[ order(code_mat[ ,1], code_mat[,2]), ]
	id_list = apply(code_mat, 1, codeFun)
	unique_id_list = unique(id_list)
	
	#Set up storage
	composite_mat = matrix(0, nrow = iters, ncol =  length(unique_id_list))
	count_vec = rep(0, length(unique_id_list))
	
	for (i in 1:length(returner_list)){
		indices = unique_id_list %in% save_list[[i]]$unique_id 
		composite_mat[ ,indices] = composite_mat[ ,indices] + save_list[[i]]$mcmc_mat
		count_vec[indices] = count_vec[indices] + 1
	}
	
	if (min(count_vec) != max(count_vec)){
	message('Min number of scores in composite is ', min(count_vec) )
	}
    #Take avg
    composite_mat = t( t(composite_mat) / count_vec)

    #Take median and quaniles to get a range for scores
    point_and_range = pointQuantile(composite_mat, alpha)
	
	#Compile Data
	inst_mat = data.frame(unique_id = unique_id_list,
						point_estimate = point_and_range$point,
						lower = point_and_range$range[1, ],
						upper = point_and_range$range[2, ], 
						n_composite = count_vec)
	
	#COnvert the codes back into inst and subset
	id_indices = match(unique_id_list, id_list)
	inst_mat$inst = code_mat[id_indices, 1]
	inst_mat$subset = code_mat[id_indices, 2]

	return(inst_mat[ ,c('inst', 'subset', 'point_estimate', 'lower', 'upper', 'n_composite')])
}

	#Find the identifying variables and mcmc mat, math is the same in each case
	iters = 'iters'	
	if (type == 'inst'){
		g_lab = 'group_labels'
		mcmc_mat = 'dg_z'
	} else if (type == 'subset_nobaseline'){
		g_lab = 'unique_subset_vec'
		mcmc_mat = 'subset_fit_nobaseline'
	} else if (type == 'subset_baseline'){
		g_lab = 'unique_subset_vec'
		mcmc_mat = 'subset_fit_baseline'
	} else {
		stop("Invalid type, try 'inst', 'subset_nobaseline', or 'subset_baseline'.")
	}

for (i in 1:length(returner_list)){
	#Process elements of list
	temp = returner_list[[i]]
	
	#Ensure group_labels are numeric
	temp[[g_lab]] = as.numeric(as.character(unlist(temp[[g_lab]])))

	#Ensure iters is numeric
	temp$iters = as.numeric(temp[[iters]])

	returner_list[[i]] = temp
}
    
    iter_vec = sapply(returner_list, function(x) x[[iters]])
    stopifnot( length( unique(iter_vec) ) == 1)
    iters = iter_vec[1]    
	
    ##Idenitfy unique set of MCMC iterations
    group_list = unique(unlist(lapply(returner_list, function(x) x[[g_lab]])))
    group_list = group_list[order(group_list)] #Put in numeric order
    p_total = length(group_list) #Number of total institutions
	
    ##Set up storage
    composite_mat = matrix(0, nrow = iters, ncol = p_total)
    count_vec = rep(0, p_total) 

    #Go through indicators and add to matrix
     for (returner in returner_list){
        indices = group_list %in% returner[[g_lab]]
        #Add to dg_z to composite mat and add 1 to count_vec
        composite_mat[ ,indices] = composite_mat[ ,indices] + returner[[mcmc_mat]]
        count_vec[indices] = count_vec[indices] + 1
    }

	if (min(count_vec) != max(count_vec)){
		message('Min number of scores in composite is ', min(count_vec) )
	}
    #Take avg
    composite_mat = t( t(composite_mat) / count_vec)

    #Take median and quaniles to get a range for scores
    point_and_range = pointQuantile(composite_mat, alpha)
    inst_mat = data.frame(id = group_list, point_estimate = point_and_range$point,
                  lower = point_and_range$range[1, ],
                  upper = point_and_range$range[2, ], 
				  n_composite = count_vec)    
    return(inst_mat)
}
