dataListFromSAS = function(sas7bdat_path, out_data_path){
	#Takes a while (mainly reading from the SAS file)
	library(sas7bdat)

	#This defines how the data is set up
	# use 'NULL' to indicate a lack of risk adjusters (cont or cat)
	inst_identifier = 'hospno'
	indicator_list = list(
		list(
			'outcome' = 'survival',
			'cat_vars' = c('sga_rev1_cat', 'male_cat', 'outborn_cat', 'mult_cat', 'pcare_cat'),
			'cont_vars' = c('z_gaweeks', 'z_gaweeks_sq', 'z_ap5')
			),
		list(
			'outcome' = 'breastmilk_rev',
			'cat_vars' = c('mult_cat', 'outborn_cat', 'pcare_cat', 'csect_cat'),
			'cont_vars' = c('z_gaweeks', 'z_gaweeks_sq', 'z_ap5')
			),
		list(
			'outcome' = 'notcold',
			'cat_vars' = c('sga_rev1_cat', 'outborn_cat', 'mult_cat', 'pcare_cat', 'csect_cat'),
			'cont_vars' = c('z_gaweeks', 'z_gaweeks_sq', 'z_ap5')
			),
		list(
			'outcome' = 'nocld',
			'cat_vars' = c('sga_rev1_cat', 'male_cat', 'outborn_cat', 'mult_cat',
						   'csect_cat', 'pcare_cat'),
			'cont_vars' = c('z_gaweeks', 'z_gaweeks_sq', 'z_ap5')
			),
		list(
			'outcome' = 'nonci_rev',
			'cat_vars' = c('sga_rev1_cat', 'outborn_cat', 'mult_cat'),
			'cont_vars' = c('z_gaweeks', 'z_gaweeks_sq', 'z_ap5')
			),
		list(
			'outcome' = 'nopntx',
			'cat_vars' = c('mult_cat', 'csect_cat', 'sga_rev1_cat', 'male_cat', 'outborn_cat', 'pcare_cat'),
			'cont_vars' =  c('z_gaweeks', 'z_gaweeks_sq', 'z_ap5')
			),
		list(
			'outcome' = 'eyex',
			'cat_vars' = NULL,
			'cont_vars' =  NULL
			),
		list(
			'outcome' = 'ante',
			'cat_vars' = NULL,
			'cont_vars' =  NULL
			)
		)

	#Read data (takes a while)
	dat = read.sas7bdat(sas7bdat_path)
	#Pull data for each indicator
	full_data_list = list()
	for (i in 1:length(indicator_list)){
		temp = indicator_list[[i]]
		#Count number of categorical and continous risk adjusters
		temp$num_cat = length(temp$cat_vars)
		temp$num_cont = length(temp$cont_vars)
		#Build the minimal_data matrix
		outcome_vec = dat[ ,temp$outcome]
		outcome_vec[is.nan(outcome_vec)] = 0
		temp$minimal_data = cbind(outcome_vec,
		   dat[ ,c(inst_identifier, temp$cat_vars, temp$cont_vars)])
		names(temp$minimal_data)[1] = temp$outcome

		#Store in list with the outcome name as the key value
		full_data_list[[temp$outcome]] = temp
	}
	#Save
	save(full_data_list, file = out_data_path)
}

