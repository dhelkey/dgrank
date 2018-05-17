nullSimulate = function(minimal_data, num_cat, num_cont,
                var_intercept,var_inst, var_cat, var_cat_2way, var_cont,
                var_subset = 1,	var_inst_subset_2way = 1,
                       visualize = TRUE,
                       mc_iters = 10,
                       iters = 10,
                       burn_in = 10){

    #Leverage the function we already have
    dat = fitBabyMonitor(minimal_data, num_cat,
                            num_cont,
                        1.1,1.2,1.3,1.4,1.5, verbose = FALSE, iters = 2) 
    
    ##Fit with probit IGNORING institution 
    null_fit = probitFit(dat$outcome_vec, dat$no_inst_mat,
                        dat$no_inst_var_vec, iter = mc_iters + burn_in)
    null_fit = null_fit[-(1:burn_in), ] #Rows are MCMC iter, cols are beta coefs

    #COnvert to individual level probabilities (assuming Probit regression)
    prob_mat = pnorm(dat$no_inst_mat %*%  t(null_fit)) #Each row represents a person, each col is an iteration

    #Simulate all at once (might have storeage problems w/ too large values for mc_iters and iters)
    outcome_mat = matrix(rbinom(n = length(prob_mat), size = 1 , prob = prob_mat), ncol = mc_iters)

    #Storage for simulating
    z_mat = matrix(0, nrow = mc_iters, ncol = dat$p)

    for (i in 1:mc_iters){
        outcome_vec = outcome_mat[ ,i]
        mcmc_fit = probitFit(outcome_vec, dat$model_matrix, 
                             dat$prior_var_vec, iters = burn_in + iters )
        mcmc_fit = mcmc_fit[-(1:burn_in), ]
        colnames(mcmc_fit) = dat$coefs
        inst_fit = mcmc_fit[ ,dat$all_inst_coefs]
        #D-G transformation
        inst_fit_dg_z = toDG(inst_fit, count_vec = dat$inst_mat$n)
        z_mat[i, ] = apply(inst_fit_dg_z, 2, median)
    }
    return(z_mat)
}