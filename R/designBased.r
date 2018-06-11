designBased = function(gamma, outcome_vec, pcf_cat_vec, instid_vec){
	#Implement design based approach and compute O, E, and z vectors

	#Some of the inputs need to be factors
	pcf_cat_vec = as.factor(pcf_cat_vec)
	instid_vec = as.factor(instid_vec)
	
	#Summarize
	instid_unique = levels(instid_vec)
	n_students = length(outcome_vec)
	n_u = length(levels(instid_vec))
	n_cat = length(levels(pcf_cat_vec))
	p_global = mean(outcome_vec)
	
	#Compute Table 2 in Draper-Gittoes
	y_mat = n_mat = matrix(0, nrow =n_u, ncol = n_cat)

	for (i in 1:n_students){
    inst = as.integer(instid_vec[i])
    categ = as.integer(pcf_cat_vec[i])
    n_mat[inst,categ] = n_mat[inst,categ] + 1
    y_mat[inst, categ] = y_mat[inst, categ] + outcome_vec[i]
	}   
	p_mat = y_mat / n_mat
	p_mat[n_mat == 0] = 0

	#Compute marginals
	pcf_p_vec = apply(y_mat, 2, sum) / apply(n_mat, 2, sum)
	n_u_vec = apply(n_mat, 1, sum)
	n_cat_vec = apply(n_mat, 2, sum)

	
	#Equation (2) in Draper-Gittoes (2004)
	O = apply(y_mat, 1, sum) / apply(n_mat, 1, sum)
	#Equation (3) in D-G
	E = apply(t(t(n_mat) * pcf_p_vec), 1, sum) / n_u_vec
	#Equation (4) in D-G
	D = O - E
	
	##Now compute variance
	gammaFun = function(i,k,j){
		#Equation (7) in D-G
		if (i == k){
			return(n_mat[i,j]/n_u_vec[i] * (1 - (n_mat[i,j]/n_cat_vec[j])))
		} else{
			return( -1 * (n_mat[i,j] * n_mat[k,j]) / (n_u_vec[i] * n_cat_vec[j]) )
		}
	}
	vFun = function(k,j, gamma = 0.5){
		#Equation (11) in D-G
		n = n_mat[k,j]
		pcomb = function(p_local){gamma * p_global + (1-gamma) * p_local}
		if (n > 0){
			p_local = p_mat[k, j]
			p_star = pcomb(p_local)
			return(p_star * (1-p_star)/n)
		} else {
			return(0)
		}
	}
	computeVi = function(i, gamma){
		#Equation (8) in D-G
		v = 0
		for (k in 1:n_u){
			for (j in 1:n_cat){
				v = v + gammaFun(i,k,j)^2 * vFun(k, j, gamma)
			}
		}
		return(v)
	}
	#Now compute w/ lapply
	V = sapply(1:n_u, computeVi, gamma)
	SE = sqrt(V)
	return( data.frame(
			u = instid_unique,
			O = O,
			E = E,
			D = D,
			SE = SE,
			Z = (O - E) / SE,
			n = n_u_vec	
	))
}
