#Need to have same number of MCMC iterations in each fit
scoreComposite = function(returner_list, alpha = 0.05){
    #Constructs composite scores as well as (1-alpha)% confidence intervals from scores on multiple performance indicators
    #Composite scores constructed by taking the average of each score
    #
    #Inputs:
    #       returner_list: List of lists, each element should be the return value (a list) from fitBabyMonitor
    #       alpha: Constructs (1-alpha)% posterior intervals for the composite score
    #
    #Outputs:
    #      inst_mat


   pointQuantile = function(x_mat, alpha){
	#Helper function
        #Helper function to obtain estimates and quantiles
        #Input x_mat - (iter X p) matrix
        return(list(point = apply(x_mat, 2, median),
        range = apply(x_mat,2 , quantile, probs = c(alpha, 1-alpha)))
        )
    }

##Preprocess Returner List to be safe
for (i in 1:length(returner_list)){
	temp = returner_list[[i]]
	
	#Ensure group_labels are numeric
	temp$group_labels = as.numeric(as.character(unlist(temp$group_labels)))

	#Ensure iters is numeric
	temp$iters = as.numeric(temp$iters)

	returner_list[[i]] = temp
}
    
    iter_vec = sapply(returner_list, function(x) x$iters)
    stopifnot( length( unique(iter_vec) ) == 1)
    iters = iter_vec[1]    

    ##Idenitfy unique set of MCMC iterations
    inst_list = unique(unlist(lapply(returner_list, function(x) x$group_labels)))
    inst_list = inst_list[order(inst_list)] #Put in numeric order
    p_total = length(inst_list) #Number of total institutions

    ##Set up storage
    composite_mat = matrix(0, nrow = iters, ncol = p_total)
    inst_count_vec = rep(0, p_total) 

    #Go through indicators and add to matrix
     for (returner in returner_list){
        indices = inst_list %in% returner$group_labels
		#print(dim(composite_mat[ ,indices]))
        #Add to dg_z to composite mat and add 1 to inst_count_vec
        composite_mat[ ,indices] = composite_mat[ ,indices] + returner$dg_z
        inst_count_vec[indices] = inst_count_vec[indices] + 1
    }
	if (min(inst_count_vec) != max(inst_count_vec)){
		message('Min number of scores in composite is ', min(inst_count_vec) )
	}
    #Take avg
    composite_mat = t( t(composite_mat) / inst_count_vec)

    #Take median and quaniles to get a range for scores
    point_and_range = pointQuantile(composite_mat, alpha)
    inst_mat = data.frame(id = inst_list, point_estimate = point_and_range$point,
                  lower = point_and_range$range[1, ],
                  upper = point_and_range$range[2, ])    
    return(inst_mat)
}
