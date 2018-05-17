processMinimalData = function(minimal_data, num_cat, num_cont, subset = FALSE){
    #Parses data from SAS and computes institution level statistics
	#Removes incomplete data rows
    #
    #Inputs:
    # Same as fitBabyMonitor()
    #
    #Returns a list with elements:
    # N - total number of individuals [scalar]
    # p - number of institutions [scalar]
    # outcome_vec, instid_vec, subset_vec [N]
    # inst_mat - matrix with rows (id, n, o_rate) [p rows, 1 per institution]

    #Remove incomplete data rows
	N_full =  dim(minimal_data)[1]
	minimal_data = minimal_data[ complete.cases(minimal_data), ]
	
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
    cat_var_mat = as.matrix(minimal_data[  ,(var_start_index):(var_start_index + num_cat - 1)])
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
    return( list(outcome_vec = outcome_vec, inst_vec = inst_vec, subset_vec = subset_vec, 
                cat_var_mat = cat_var_mat, cont_var_mat = cont_var_mat,
                inst_mat = inst_mat,
                N = dim(minimal_data)[1], 
				N_full = N_full,
                p = length(unique_inst_vec),
                indicator_name = names(minimal_data)[1],
                o_overall = mean(outcome_vec)))
}