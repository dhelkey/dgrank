#Need to have same number of MCMC iterations in each fit
scoreComposite = function(returner_list, alpha = 0.05, type = 'inst'){
    #Constructs composite scores as well as (1-alpha)% confidence intervals from scores on multiple performance indicators
    #Composite scores constructed by taking the average of each score
    #
    #Inputs:
    #       returner_list: List of lists, each element should be the return value (a list) from fitBabyMonitor
    #       alpha: Constructs (1-alpha)% posterior intervals for the composite score
	#		mode: 'inst' or 'subset_baseline' or 'subset_nobaseline'
    #
    #Outputs:
    #      inst_mat

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
		#print(dim(composite_mat[ ,indices]))
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
                  upper = point_and_range$range[2, ])    
    return(inst_mat)
}
