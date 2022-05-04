*** Predictors of RSBB (B3911) - G0 mother analysis script
*** Created 16/11/2021 by Dan Smith
*** Stata v17.0

*** This script reads in the cleaned G0 mother data, explores associations between exposures, and then conducts an analysis exploring how all the demographic and SES variables are associated with various facets of RSBB.


**********************************************************************************
**** Set working directory, start a log file, and read in cleaned dataset

cd "X:\Groups\ARC\DanS\Descriptive_PredictorsOfRSBB_B3911"

capture log close
log using ".\G0Mother_Results\Desc_RSBB_B3911_G0MotherAnalysis_log", replace text

use "G0Mother_PredictorsOfRSBB_B3911.dta", clear


** As want to use the number of distinct values in a variable later in the script, need to install the user-written 'distinct' package (type 'search distinct' and install package 'dm0042_2')

** Also want to install the 'heatplot' package for displaying correlation plots (plus some dependencies)
*ssc install heatplot, replace
*ssc install palettes, replace
*ssc install colrspace, replace

** Install 'grc1leg' to merge plots together with a single legend
*ssc install grc1leg, replace


** For this analysis, as only have RSBB data for core participants, will drop all non-core pregnancies
tab in_core, m
drop if in_core == 2
drop in_core


**********************************************************************************
*** Descriptive statistics

** Put the RSBB variables at the start of the dataset
order aln d810 d813 d813_grp d816 Y3000 Y3040 Y3040_grp Y3080 Y3080_OccNever Y3080_OccYr Y3153 Y3153_cat Y3160 Y3170 Y3155 Y3155_cat

* Keep just pregnancy RSBB variables (belief in God, religious affiliation and frequency of church attendance)
drop Y3000 Y3040 Y3040_grp Y3080 Y3080_OccNever Y3080_OccYr Y3153 Y3153_cat Y3160 Y3170 Y3155 Y3155_cat

* Also drop the cognitive/psychological variables (as will focus on these in another paper), and also the adverse childhood experiences variables
drop b916-b921 d842 h151b logic_mem-fom_cog_factor1
drop c432 c433

** Using the 'distinct' command (see above), for each variable inspect the number of unque values; if < 10 then display table, while if >= 10 then displays means/SDs/IQRs

* Outcomes
foreach var of varlist d810-d816 {
	quietly distinct `var'
	local unique = r(ndistinct)
	
	display ""
	display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	display "Variable " "`var'" " has " `unique' " values."
	
	if `unique' < 10 {
		tab `var'
		tab `var', m
	}
	else {
		sum `var'
		sum `var', d
	}
}

* Exposures
foreach var of varlist mz028b-partner_ab {
	quietly distinct `var'
	local unique = r(ndistinct)
	
	display ""
	display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	display "Variable " "`var'" " has " `unique' " values."
	
	if `unique' < 10 {
		tab `var'
		tab `var', m
	}
	else {
		sum `var'
		sum `var', d
		quietly misstable sum `var'
		local missing = r(N_eq_dot)
		local percent = (`missing' / _N) * 100
		display "Variable " "`var'" " has " `missing' " missing values ( " `percent' " %). "
	}
}


** Amount of missing data in outcomes and exposures

* Outcomes
misstable sum d810-d816, all

* Exposures
misstable sum mz028b-partner_ab, all



**********************************************************************************
*** Correlations between exposures

** Explore correlations between exposures to see how inter-related these factors are
desc mz028b-partner_ab

* Most variables are either continuous, ordered categories, or binary variables, so will just use standard pearson correlations for these. Only unordered categories are home ownership (owned/mortaged vs renting vs counci/housing association vs other) and marital status (never married vs married vs widowed/divorced/separated). So will exclude these two variables from the correlation matrix and calculate their associations using these variables as outcomes in a multinomial regression and square-rooting the pseudo-R2 value (similar to how the 'rexposome' package in R for ExWAS analyses does)

* First, order variables into categories of 'demographic' and 'socioeconomic/material insecurity' (Y9992 is age at @28 RSBB questionnaire, so will pop at end as not needed here)
order aln-d816 ///
mz028b c800_grp a525_grp a005_grp jan1993ur01ind_grp b032_grp ///
c645a c686a c706a c755_grp c_sc_mgm_grp c_sc_mgf_grp logavinceq jan1993imd2010q5_M jan1993Townsendq5_M a006_grp b594_grp c525 c429a a053 a551 a636 partner_ab ///
Y9992

* Next, rename all these exposures so they are more intuitive and can be read easier on the correlation heatmaps below
rename mz028b ageAtBirth
rename c800_grp nonWhiteEthnic
rename a525_grp maritalStatus
rename a005_grp mobility 
rename jan1993ur01ind_grp rural
rename b032_grp parity
rename c645a education
rename c686a maternalEdu
rename c706a paternalEdu
rename c755_grp highSocClass
rename c_sc_mgm_grp highSocClass_mat
rename c_sc_mgf_grp highSocClass_pat
rename logavinceq income
rename jan1993imd2010q5_M IMD
rename jan1993Townsendq5_M townsendDep
rename c429a poorerChildhood
rename a053 accessToCar
rename a006_grp housing 
rename b594_grp financeDiffs
rename c525 financeDiffsScore
rename a551 crowding
rename a636 neighPercept
rename partner_ab partnerAbsence
rename Y9992 ageAt28

* Associations between demographic factors (excluding marital status; a525_grp) - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
corr ageAtBirth nonWhiteEthnic mobility rural parity

matrix cor_demo = r(C)
matrix list cor_demo

heatplot cor_demo, values(format(%9.3f)) color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) legend(off) xlabel(, angle(45))

* Save heatmap
graph export ".\G0Mother_Results\corr_heatplot_demoOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Mother_Results\corrMatrix_demoOnly.xlsx", replace
putexcel A1=matrix(cor_demo), names nformat(number_d2)

* And now for associations of demographic variables with unordered categorical variable marital status (and save to a CSV file)
capture postclose marital_corrs_demoOnly
postfile marital_corrs_demoOnly str30 variable corr ///
	using ".\G0Mother_Results\marital_corrs_demoOnly.dta", replace

foreach var of varlist ageAtBirth nonWhiteEthnic mobility rural parity {
	quietly mlogit maritalStatus `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on marital status is: " round((sqrt(`r2')), .001)
	
	post marital_corrs_demoOnly ("`var'") (sqrt(`r2'))
}

postclose marital_corrs_demoOnly

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G0Mother_Results\marital_corrs_demoOnly.dta", clear
list, clean
outsheet using ".\G0Mother_Results\marital_corrs_demoOnly.csv", comma replace
restore


** Now repeat for socioeconomic/material insecurity variables (to exclude housing status [housing], as is an unordered categorical variable)
corr education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept partnerAbsence

matrix cor_socio = r(C)

* As lots of entries, the correlation coefficients are hard to read, so will drop these values and include the legend
heatplot cor_socio, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))
	
* Save heatmap
graph export ".\G0Mother_Results\corr_heatplot_socioOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Mother_Results\corrMatrix_socioOnly.xlsx", replace
putexcel A1=matrix(cor_socio), names nformat(number_d2)


* And now for associations of socioeconomic variables with unordered categorical variable home ownership status (and save to a CSV file)
capture postclose housing_corrs_socioOnly
postfile housing_corrs_socioOnly str30 variable corr ///
	using ".\G0Mother_Results\housing_corrs_socioOnly.dta", replace
	
foreach var of varlist education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept partnerAbsence {
	quietly mlogit housing `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on home ownership status is: " round((sqrt(`r2')), .001)

	post housing_corrs_socioOnly ("`var'") (sqrt(`r2'))
}

postclose housing_corrs_socioOnly

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G0Mother_Results\housing_corrs_socioOnly.dta", clear
list, clean
outsheet using ".\G0Mother_Results\housing_corrs_socioOnly.csv", comma replace
restore

	
*** Finally, repeat this on all of the demographic and SES exposures together (excluding unordered cateogorical variables housing status and marital status)
corr ageAtBirth nonWhiteEthnic mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept partnerAbsence

matrix cor_all = r(C)

heatplot cor_all, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle(""))

* Save heatmap
graph export ".\G0Mother_Results\corr_heatplot_all.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Mother_Results\corrMatrix_all.xlsx", replace
putexcel A1=matrix(cor_all), names nformat(number_d2)

	
* And now for associations of all other exposures variables with unordered categorical variables marital status and home ownership status

* Marital status
capture postclose marital_corrs_all
postfile marital_corrs_all str30 variable corr ///
	using ".\G0Mother_Results\marital_corrs_all.dta", replace
	
foreach var of varlist ageAtBirth nonWhiteEthnic mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep housing financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept partnerAbsence {
	quietly mlogit maritalStatus `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on marital status is: " round((sqrt(`r2')), .001)

	post marital_corrs_all ("`var'") (sqrt(`r2'))
}

postclose marital_corrs_all

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G0Mother_Results\marital_corrs_all.dta", clear
list, clean
outsheet using ".\G0Mother_Results\marital_corrs_all.csv", comma replace
restore


* Housing status
capture postclose housing_corrs_all
postfile housing_corrs_all str30 variable corr ///
	using ".\G0Mother_Results\housing_corrs_all.dta", replace
	
foreach var of varlist ageAtBirth nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep financeDiffs poorerChildhood financeDiffsScore accessToCar crowding neighPercept partnerAbsence {
	quietly mlogit housing `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on home ownership status is: " round((sqrt(`r2')), .001)


	post housing_corrs_all ("`var'") (sqrt(`r2'))
}

postclose housing_corrs_all

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G0Mother_Results\housing_corrs_all.dta", clear
list, clean
outsheet using ".\G0Mother_Results\housing_corrs_all.csv", comma replace
restore



************************************************************************************
*** Next, we want to run the actual analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab d810, m

** This will be quite complicated, as want to post results to file, but as exposures differ extracting the results will be variable specific. To adjust for multiple corrections will use conservative bonferroni adjustment when constructing confidence intervals and interpreting p-values - As 23 exposures, will a Bonferroni p-value threshold of 0.05/23 = 0.0022.
display 0.05/23

** We also want to store both estimates adjusting for age (other than for the age-only model), and also the interaction between age and the exposure, to see whether it's moderated by age. Again, this makes the set-up a bit more complicated.

** Create a postfile to post results to, then start the loop - Will create three postfiles; one for coefficients and CIs, another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model), and then a third for the (pseduo-)R2 values to assess model fit - NOTE: Have to store pvalues as 'double' format, else really tiny p-values get coded as '0' (as default if float format, which has minimum value of -3.40282346639e+38 [https://blog.stata.com/2012/04/02/the-penultimate-guide-to-precision/]).
capture postclose mother_belief
postfile mother_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Mother_Results\mother_belief_results.dta", replace

capture postclose mother_belief_lr
postfile mother_belief_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Mother_Results\mother_belief_results_lr.dta", replace
	
capture postclose mother_belief_r2
postfile mother_belief_r2 str30 exposure r2_main r2_int ///
	using ".\G0Mother_Results\mother_belief_results_r2.dta", replace

foreach var of varlist ageAtBirth nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep housing financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept partnerAbsence {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAtBirth' separately first as will be adjusted for in all other models
	if "`var'" == "ageAtBirth" {
		mlogit d810 `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/yes)
		local outcome_level = "Yes (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,1]
		local lci = res[5,1]
		local uci = res[6,1]
		local p = res[4,1]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local age_main = .
		local exp_main = .
		
		post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/not sure)
		local outcome_level = "Not sure (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,3]
		local lci = res[5,3]
		local uci = res[6,3]
		local p = res[4,3]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local age_main = .
		local exp_main = .
		
		post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit d810 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit d810 `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d810 `var', baseoutcome(3) rrr
		local r2_main = e(r2_p)
		
		// As no interaction model for age, will just fill with missing values
		local lr_p_int = .
		local r2_int = .
		
		post mother_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_belief_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "highSocClass_mat" | "`var'" == "highSocClass_pat" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "financeDiffsScore" | "`var'" == "poorerChildhood" | "`var'" == "accessToCar" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" {
		
		mlogit d810 ageAtBirth `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/yes)
		local outcome_level = "Yes (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// Now for interaction model
		mlogit d810 c.ageAtBirth##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/not sure)
		mlogit d810 ageAtBirth `var', baseoutcome(3) rrr
		
		local outcome_level = "Not sure (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
				
		// Now for interaction model
		mlogit d810 c.ageAtBirth##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit d810 ageAtBirth if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit d810 ageAtBirth `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d810 ageAtBirth if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit d810 ageAtBirth `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit d810 c.ageAtBirth##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit d810 c.ageAtBirth##c.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post mother_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_belief_r2 ("`exp'") (`r2_main') (`r2_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Married (ref = Never married)"
			}
			else if "`var'" == "parity" {
				local exp_level = "1 (ref = 0)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,3]
			local lci = res[5,3]
			local uci = res[6,3]
			local p = res[4,3]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Move to the next category of the exposure (category 3)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity" {
				local exp_level = "2 or more (ref = 0)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,3]
			local lci = res[5,3]
			local uci = res[6,3]
			local p = res[4,3]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,3]
			local lci = res[5,3]
			local uci = res[6,3]
			local p = res[4,3]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,17]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,3]
			local lci = res[5,3]
			local uci = res[6,3]
			local p = res[4,3]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures and R2 values
		mlogit d810 ageAtBirth if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d810 ageAtBirth if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post mother_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_belief_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose mother_belief
postclose mother_belief_lr
postclose mother_belief_r2
	

************************************************************************************
*** Now to the next RSBB outcome: Religious affiliation

*** As this is an unordered categorical variable, will again use multinomial logistic model (with 'none' as reference)
tab d813_grp

* Do need to re-order the categories so can simply copy and paste the script from above without having to faff around with editing the cells to take the statistics from
recode d813_grp (1 = 3) (2 = 1) (3 = 2)
label define relig_lb 1 "Christian" 2 "Other" 3 "None", modify
numlabel relig_lb, add
tab d813_grp


*** Now run the loop to save all the results
capture postclose mother_relig
postfile mother_relig str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Mother_Results\mother_relig_results.dta", replace

capture postclose mother_relig_lr
postfile mother_relig_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Mother_Results\mother_relig_results_lr.dta", replace
	
capture postclose mother_relig_r2
postfile mother_relig_r2 str30 exposure r2_main r2_int ///
	using ".\G0Mother_Results\mother_relig_results_r2.dta", replace

foreach var of varlist ageAtBirth nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep housing financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept partnerAbsence {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAtBirth' separately first as will be adjusted for in all other models
	if "`var'" == "ageAtBirth" {
		mlogit d813_grp `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/Chrstian)
		local outcome_level = "Christian (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,1]
		local lci = res[5,1]
		local uci = res[6,1]
		local p = res[4,1]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local age_main = .
		local exp_main = .
		
		post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/Other)
		local outcome_level = "Other (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,3]
		local lci = res[5,3]
		local uci = res[6,3]
		local p = res[4,3]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local age_main = .
		local exp_main = .
		
		post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit d813_grp if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit d813_grp `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d813_grp `var', baseoutcome(3) rrr
		local r2_main = e(r2_p)
		
		// As no interaction model for age, will just fill with missing values
		local lr_p_int = .
		local r2_int = .
		
		post mother_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_relig_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "highSocClass_mat" | "`var'" == "highSocClass_pat" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "financeDiffsScore" | "`var'" == "poorerChildhood" | "`var'" == "accessToCar" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" {
		
		mlogit d813_grp ageAtBirth `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/Chrstian)
		local outcome_level = "Christian (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// Now for interaction model
		mlogit d813_grp c.ageAtBirth##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/Other)
		mlogit d813_grp ageAtBirth `var', baseoutcome(3) rrr
		
		local outcome_level = "Other (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
				
		// Now for interaction model
		mlogit d813_grp c.ageAtBirth##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit d813_grp ageAtBirth if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit d813_grp ageAtBirth `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d813_grp ageAtBirth if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit d813_grp ageAtBirth `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit d813_grp c.ageAtBirth##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit d813_grp c.ageAtBirth##c.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post mother_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_relig_r2 ("`exp'") (`r2_main') (`r2_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Married (ref = Never married)"
			}
			else if "`var'" == "parity" {
				local exp_level = "1 (ref = 0)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,3]
			local lci = res[5,3]
			local uci = res[6,3]
			local p = res[4,3]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Move to the next category of the exposure (category 3)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity" {
				local exp_level = "2 or more (ref = 0)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,3]
			local lci = res[5,3]
			local uci = res[6,3]
			local p = res[4,3]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,3]
			local lci = res[5,3]
			local uci = res[6,3]
			local p = res[4,3]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,17]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,3]
			local lci = res[5,3]
			local uci = res[6,3]
			local p = res[4,3]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures and R2 values
		mlogit d813_grp ageAtBirth if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d813_grp ageAtBirth if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post mother_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_relig_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose mother_relig
postclose mother_relig_lr
postclose mother_relig_r2



************************************************************************************
*** Now to the next RSBB outcome: Attendance at church/place of worship

*** As this is an ordered categorical variable, will use ordinal regression model. For consistency with previous models, will recode so that higher values indicate greater RSBB
tab d816

recode d816 (1 = 3) (2 = 2) (3 = 1) (4 = 0), gen(d816_rev)
label define attend_rev_lb 0 "Not at all" 1 "Min once a year" 2 "Min once a month" 3 "Min once a week"
numlabel attend_rev_lb, add
label values d816_rev attend_rev_lb
tab d816_rev


** Quick test of whether proportional odds assumption been violated in most basic model (with just age at birth). Ah, it has been violated, and gets worse if add in additional predictors. 
ologit d816_rev ageAtBirth, or
brant, detail

ologit d816_rev ageAtBirth i.IMD, or
brant, detail


** So instead will just run multinomial models with 'not at all' as the baseline/reference category (this also means all outcomes are using the same model, making them easier to compare).


*** Now run the loop to save all the results
capture postclose mother_attend
postfile mother_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Mother_Results\mother_attend_results.dta", replace

capture postclose mother_attend_lr
postfile mother_attend_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Mother_Results\mother_attend_results_lr.dta", replace
	
capture postclose mother_attend_r2
postfile mother_attend_r2 str30 exposure r2_main r2_int ///
	using ".\G0Mother_Results\mother_attend_results_r2.dta", replace

foreach var of varlist ageAtBirth nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep housing financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept partnerAbsence {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAtBirth' separately first as will be adjusted for in all other models
	if "`var'" == "ageAtBirth" {
		mlogit d816_rev `var', baseoutcome(0) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/once year)
		local outcome_level = "Min once year (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,3]
		local lci = res[5,3]
		local uci = res[6,3]
		local p = res[4,3]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local age_main = .
		local exp_main = .
		
		post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/once month)
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local age_main = .
		local exp_main = .
		
		post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// Now onto the next reference category (3/once week)
		local outcome_level = "Min once week (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local age_main = .
		local exp_main = .
		
		post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit d816_rev if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit d816_rev `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d816_rev `var', baseoutcome(0) rrr
		local r2_main = e(r2_p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		local r2_int = .
		
		post mother_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_attend_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "highSocClass_mat" | "`var'" == "highSocClass_pat" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "financeDiffsScore" | "`var'" == "poorerChildhood" | "`var'" == "accessToCar" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" {
		
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/once year)
		local outcome_level = "Min once year (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Now for interaction model
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/once month)
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr
		
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/once week)
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr
		
		local outcome_level = "Min once week (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit d816_rev ageAtBirth if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d816_rev ageAtBirth if `var' != ., baseoutcome(0) rrr
		local r2_base = e(r2_p)
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post mother_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_attend_r2 ("`exp'") (`r2_main') (`r2_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Married (ref = Never married)"
			}
			else if "`var'" == "parity" {
				local exp_level = "1 (ref = 0)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity" {
				local exp_level = "2 or more (ref = 0)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (2/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "paternalEdu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures and R2 values
		mlogit d816_rev ageAtBirth if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit d816_rev ageAtBirth if `var' != ., baseoutcome(0) rrr
		local r2_base = e(r2_p)
		mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post mother_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post mother_attend_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose mother_attend
postclose mother_attend_lr
postclose mother_attend_r2



* And close the log file
log close

* Save this file if we want to use it later
save ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", replace



***********************************************************************************
***********************************************************************************
***** Next step is to make some nice plots of these results

**** Start with outcome of belief in God/divine power 
use ".\G0Mother_Results\mother_belief_results_lr.dta", clear

** Convert string exposure var to numeric
count
local n = r(N)

capture drop exp_num
gen exp_num = 0
label define exp_lb 0 "NA", replace
label values exp_num exp_lb
tab exp_num

forvalues i = 1(1)`n' {
	
	* Save the variable name
	local var = exposure in `i'
	display "Variable " `i' " is " "`var'"
	display ""
	
	* Code variable as numeric and add value label
	replace exp_num = `i' if exposure == "`var'"
	label define exp_lb `i' "`var'", modify
}

tab exp_num


*** First display the p-value plots to show overall associations between exposures and outcome

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

** Start with likelihood ratio results comparing null model to model including exposure

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 23 exposures, will do 0.05/23)
local bon_thresh = -log10(0.05/23)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Main effect") ///
	name(belief_main, replace)
	
graph export ".\G0Mother_Results\belief_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Age interaction") ///
	name(belief_int, replace)
	
graph export ".\G0Mother_Results\belief_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/23)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(belief_both, replace)

graph export ".\G0Mother_Results\belief_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\G0Mother_Results\belief_pvalues.dta", replace


**** Now read in the next outcome - religious affiliation
use ".\G0Mother_Results\mother_relig_results_lr.dta", clear

** Convert string exposure var to numeric
count
local n = r(N)

capture drop exp_num
gen exp_num = 0
label define exp_lb 0 "NA", replace
label values exp_num exp_lb
tab exp_num

forvalues i = 1(1)`n' {
	
	* Save the variable name
	local var = exposure in `i'
	display "Variable " `i' " is " "`var'"
	display ""
	
	* Code variable as numeric and add value label
	replace exp_num = `i' if exposure == "`var'"
	label define exp_lb `i' "`var'", modify
}

tab exp_num


*** First display the p-value plots to show overall associations between exposures and outcome

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

** Start with likelihood ratio results comparing null model to model including exposure

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 23 exposures, will do 0.05/23)
local bon_thresh = -log10(0.05/23)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Main effect") ///
	name(relig_main, replace)
	
graph export ".\G0Mother_Results\relig_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Age interaction") ///
	name(relig_int, replace)
	
graph export ".\G0Mother_Results\relig_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/23)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(relig_both, replace)

graph export ".\G0Mother_Results\relig_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\G0Mother_Results\relig_pvalues.dta", replace


**** Now read in the next outcome - church attendance
use ".\G0Mother_Results\mother_attend_results_lr.dta", clear

** Convert string exposure var to numeric
count
local n = r(N)

capture drop exp_num
gen exp_num = 0
label define exp_lb 0 "NA", replace
label values exp_num exp_lb
tab exp_num

forvalues i = 1(1)`n' {
	
	* Save the variable name
	local var = exposure in `i'
	display "Variable " `i' " is " "`var'"
	display ""
	
	* Code variable as numeric and add value label
	replace exp_num = `i' if exposure == "`var'"
	label define exp_lb `i' "`var'", modify
}

tab exp_num


*** First display the p-value plots to show overall associations between exposures and outcome

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

** Start with likelihood ratio results comparing null model to model including exposure

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 23 exposures, will do 0.05/23)
local bon_thresh = -log10(0.05/23)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Main effect") ///
	name(attend_main, replace)
	
graph export ".\G0Mother_Results\attend_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Age interaction") ///
	name(attend_int, replace)
	
graph export ".\G0Mother_Results\attend_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/23)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(attend_both, replace)

graph export ".\G0Mother_Results\attend_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'church attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\G0Mother_Results\attend_pvalues.dta", replace


*** Combine all these datasets together
use ".\G0Mother_Results\belief_pvalues.dta", clear
append using ".\G0Mother_Results\relig_pvalues.dta"
append using ".\G0Mother_Results\attend_pvalues.dta"


** Now look at combined results

* Pregnancy vars main effects
local bon_thresh = -log10(0.05/23)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_main if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num logp_main if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(preg_main, replace)

graph export ".\G0Mother_Results\beliefReligAttend_mainEffects_pvalues.pdf", replace

* Pregnancy vars interaction effects
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if outcome == "Belief" & exp_num != 1, ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Religious affil." & exp_num != 1, ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Church attendance" & exp_num != 1, ///
		col(blue) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)23, valuelabel labsize(small) angle(0)) ///
	title("Age interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(preg_int, replace)

graph export ".\G0Mother_Results\beliefReligAttend_ageInt_pvalues.pdf", replace


** Combine all these graphs together
graph combine preg_main preg_int, ysize(3) xsize(6)

graph export ".\G0Mother_Results\allData_pvalues.pdf", replace

graph close _all


** Save these p-values as CSV files
list outcome exposure lr_p_main lr_p_int in 1/10

keep outcome exp_num lr_p_main lr_p_int
order outcome exp_num lr_p_main lr_p_int

* Convert to wide format (need to edit outcomes strings so have no spaces)
replace outcome = "Attend" if outcome == "Church attendance"
replace outcome = "Religion" if outcome == "Religious affil."
tab outcome

reshape wide lr_p_main lr_p_int, i(exp_num) j (outcome) string

order exp_num lr_p_mainBelief lr_p_intBelief lr_p_mainReligion lr_p_intReligion lr_p_mainAttend lr_p_intAttend

format %9.4f lr_p_mainBelief-lr_p_intAttend

outsheet exp_num-lr_p_intAttend using ".\G0Mother_Results\pvalue_results.csv", comma replace


** And how many exposures were associated with the outcome at both Bonferroni and standard alpha levels?

* Belief in god - main effect
count if lr_p_mainBelief < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainBelief < 0.05
display (r(N) / _N) * 100

* Belief in god - interaction
count if lr_p_intBelief < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_intBelief < 0.05
display (r(N) / (_N - 1)) * 100

* Religious affiliation - main effect
count if lr_p_mainReligion < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainReligion < 0.05
display (r(N) / _N) * 100

* Religious affiliation - interaction
count if lr_p_intReligion < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_intReligion < 0.05
display (r(N) / (_N - 1)) * 100

* Church attendance - main effect
count if lr_p_mainAttend < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainAttend < 0.05
display (r(N) / _N) * 100

* Church attendance - interaction
count if lr_p_intAttend < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_intAttend < 0.05
display (r(N) / (_N - 1)) * 100



********************************************************************************
*** Repeat the plots above, but this time for pseudo-R2 value results
**** Start with outcome of belief in God/divine power 
use ".\G0Mother_Results\mother_belief_results_r2.dta", clear

** Convert string exposure var to numeric
count
local n = r(N)

capture drop exp_num
gen exp_num = 0
label define exp_lb 0 "NA", replace
label values exp_num exp_lb
tab exp_num

forvalues i = 1(1)`n' {
	
	* Save the variable name
	local var = exposure in `i'
	display "Variable " `i' " is " "`var'"
	display ""
	
	* Code variable as numeric and add value label
	replace exp_num = `i' if exposure == "`var'"
	label define exp_lb `i' "`var'", modify
}

tab exp_num


** Start with pseduo-R2 values comparing null model to model including exposure
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Main effect") ///
	name(belief_main, replace)
	
graph export ".\G0Mother_Results\belief_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Age interaction") ///
	name(belief_int, replace)
	
graph export ".\G0Mother_Results\belief_ageInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(belief_both, replace)

graph export ".\G0Mother_Results\belief_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\G0Mother_Results\belief_r2.dta", replace


**** Now read in the next outcome - religious affiliation
use ".\G0Mother_Results\mother_relig_results_r2.dta", clear

** Convert string exposure var to numeric
count
local n = r(N)

capture drop exp_num
gen exp_num = 0
label define exp_lb 0 "NA", replace
label values exp_num exp_lb
tab exp_num

forvalues i = 1(1)`n' {
	
	* Save the variable name
	local var = exposure in `i'
	display "Variable " `i' " is " "`var'"
	display ""
	
	* Code variable as numeric and add value label
	replace exp_num = `i' if exposure == "`var'"
	label define exp_lb `i' "`var'", modify
}

tab exp_num


** Start with pseduo-R2 values comparing null model to model including exposure
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Main effect") ///
	name(relig_main, replace)
	
graph export ".\G0Mother_Results\relig_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Age interaction") ///
	name(relig_int, replace)
	
graph export ".\G0Mother_Results\relig_ageInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(relig_both, replace)

graph export ".\G0Mother_Results\relig_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\G0Mother_Results\relig_r2.dta", replace


**** Now read in the next outcome - church attendance
use ".\G0Mother_Results\mother_attend_results_r2.dta", clear

** Convert string exposure var to numeric
count
local n = r(N)

capture drop exp_num
gen exp_num = 0
label define exp_lb 0 "NA", replace
label values exp_num exp_lb
tab exp_num

forvalues i = 1(1)`n' {
	
	* Save the variable name
	local var = exposure in `i'
	display "Variable " `i' " is " "`var'"
	display ""
	
	* Code variable as numeric and add value label
	replace exp_num = `i' if exposure == "`var'"
	label define exp_lb `i' "`var'", modify
}

tab exp_num


** Start with pseduo-R2 values comparing null model to model including exposure
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Main effect") ///
	name(attend_main, replace)
	
graph export ".\G0Mother_Results\attend_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Age interaction") ///
	name(attend_int, replace)
	
graph export ".\G0Mother_Results\attend_ageInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(attend_both, replace)

graph export ".\G0Mother_Results\attend_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'church attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\G0Mother_Results\attend_r2.dta", replace


*** Combine all these datasets together
use ".\G0Mother_Results\belief_r2.dta", clear
append using ".\G0Mother_Results\relig_r2.dta"
append using ".\G0Mother_Results\attend_r2.dta"


** Now look at combined results

* Pregnancy vars main effects
twoway (scatter exp_num r2_main if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)23, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(preg_main, replace)

graph export ".\G0Mother_Results\beliefReligAttend_mainEffects_r2.pdf", replace

* Pregnancy vars interaction effects
twoway (scatter exp_num r2_int if outcome == "Belief" & exp_num != 1, ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Religious affil." & exp_num != 1, ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Church attendance" & exp_num != 1, ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)23, valuelabel labsize(small) angle(0)) ///
	title("Age interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(preg_int, replace)

graph export ".\G0Mother_Results\beliefReligAttend_ageInt_r2.pdf", replace


** Combine all these graphs together
graph combine preg_main preg_int, ysize(3) xsize(6)

graph export ".\G0Mother_Results\allData_r2.pdf", replace

graph close _all


** Save these pseudo_R2 as CSV files
list outcome exposure r2_main r2_int in 1/10

keep outcome exp_num r2_main r2_int
order outcome exp_num r2_main r2_int

* Convert to wide format (need to edit outcomes strings so have no spaces)
replace outcome = "Attend" if outcome == "Church attendance"
replace outcome = "Religion" if outcome == "Religious affil."
tab outcome

reshape wide r2_main r2_int, i(exp_num) j (outcome) string

order exp_num r2_mainBelief r2_intBelief r2_mainReligion r2_intReligion r2_mainAttend r2_intAttend

format %9.4f r2_mainBelief-r2_intAttend

outsheet exp_num-r2_intAttend using ".\G0Mother_Results\r2_results.csv", comma replace


********************************************************************************
*** Next, want to plot some of the actual coefficient results

** Will read the datasets in then combine together into one single dataset
use ".\G0Mother_Results\mother_belief_results.dta", clear
gen outcome = "Belief"

append using ".\G0Mother_Results\mother_relig_results.dta"
replace outcome = "Relig" if outcome == ""
tab outcome, m

append using ".\G0Mother_Results\mother_attend_results.dta"
replace outcome = "Attend" if outcome == ""
tab outcome, m


** Save these results as CSV files to add to the SI

* Save each result in turn
format coef lci uci coef_int lci_int uci_int %9.3f
format p p_int %9.4f

outsheet exposure-p_int using ".\G0Mother_Results\belief_coefs.csv" if outcome == "Belief", comma replace

outsheet exposure-p_int using ".\G0Mother_Results\relig_coefs.csv" if outcome == "Relig", comma replace

outsheet exposure-p_int using ".\G0Mother_Results\attend_coefs.csv" if outcome == "Attend", comma replace

* Convert format back to default format (so axis on plots display correctly)
format coef lci uci coef_int lci_int uci_int %9.0g
format p p_int %10.0g


**** Now make the actual plots

** First, make a plot for the age results - As some outcomes are on different scales, will just use the results from multinomial regression for all outcomes - Having all variables on the same plot on the same scale makes things easier to visualise.
capture drop level_num
gen level_num = 0
replace level_num = 1 if outcome_level == "Not sure (ref = No)"
replace level_num = 3 if outcome_level == "Christian (ref = None)"
replace level_num = 4 if outcome_level == "Other (ref = None)"
replace level_num = 6 if outcome_level == "Min once year (ref = Not at al"
replace level_num = 7 if outcome_level == "Min once month (ref = Not at a"
replace level_num = 8 if outcome_level == "Min once week (ref = Not at al"

label define level_lb 0 "Belief - Yes (ref = No)" 1 "Belief - Not sure (ref = No)" 3 "Affiliation - Christian (ref = None)" 4 "Affiliation - Other (ref = None)" 6 "Attendance - Min 1/year (ref = Not at all)" 7 "Attendance - Min 1/month (ref = Not at all)" 8 "Attendance - Min 1/week (ref = Not at all)", replace
label values level_num level_lb
tab level_num

twoway (scatter level_num coef if outcome == "Belief" & exposure == "ageAtBirth", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & ///
			exposure == "ageAtBirth", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "ageAtBirth", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & ///
			exposure == "ageAtBirth", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "ageAtBirth", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & ///
			exposure == "ageAtBirth", horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Age and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(age_cat, replace)
		
graph export ".\G0Mother_Results\ageResults.pdf", replace


** Create plot for ethnicity (ref = white)
sum lci uci if exposure == "nonWhiteEthnic" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & ///
		exposure == "nonWhiteEthnic", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & ///
			exposure == "nonWhiteEthnic", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & ///
			exposure == "nonWhiteEthnic", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & ///
			exposure == "nonWhiteEthnic", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & ///
			exposure == "nonWhiteEthnic", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & ///
			exposure == "nonWhiteEthnic", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = White)") ///
		title("Other than White ethnicity and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.5 1 2 5 10 15, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(ethnic_cat, replace)
		
graph export ".\G0Mother_Results\ethnicityResults.pdf", replace


** Create plot for marital status (ref = never married)

* As two exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.2 if exp_level == "Married (ref = Never married)"
replace level_split = level_num + 0.2 if exp_level == "Wid/Div/Sep (ref = Never married)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci uci if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal col(red)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = Never married)") ///
		title("Marital status and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.5 1 2 3 5 8, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Married" 7 "Widowed/Divorced/Separated")) ///
		name(marital_cat, replace)
		
graph export ".\G0Mother_Results\maritalStatusResults.pdf", replace


** Create plot for education (ref = CSE/None)

* As four exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "education" & exp_level == "Vocational (ref = CSE/None)"
replace level_split = level_num - 0.1 if exposure == "education" & exp_level == "O-level (ref = CSE/None)"
replace level_split = level_num + 0.1 if exposure == "education" & exp_level == "A-level (ref = CSE/None)"
replace level_split = level_num + 0.3 if exposure == "education" & exp_level == "Degree (ref = CSE/None)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci uci if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = CSE/None)") ///
		title("Education and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.5 1 2 3 5 8, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Vocational" 7 "O-levels" 13 "A-levels" ///
			19 "Degree") rows(1)) ///
		name(edu_cat, replace)
		
graph export ".\G0Mother_Results\eduResults.pdf", replace
	

** Create plot for occupational social class (ref = lower)
sum lci uci if exposure == "highSocClass" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"highSocClass", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"highSocClass", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"highSocClass", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = lower [III manual/IV/V])") ///
		title("Occupational Social Class and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(1 1.5 2 3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(off) name(socClass_cat, replace)
		
graph export ".\G0Mother_Results\socClassResults.pdf", replace


** Create plot for income
sum lci uci if exposure == "income" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"income", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"income", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"income", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"income", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"income", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"income", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in log income)") ///
		title("Household income and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(1 1.5 2 2.5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(off) name(income_cat, replace)
		
graph export ".\G0Mother_Results\incomeResults.pdf", replace


** Create plot for IMD (ref = 1/least deprived)

* As four exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "IMD" & exp_level == "2 (ref = 1/Least dep.)"
replace level_split = level_num - 0.1 if exposure == "IMD" & exp_level == "3 (ref = 1/Least dep.)"
replace level_split = level_num + 0.1 if exposure == "IMD" & exp_level == "4 (ref = 1/Least dep.)"
replace level_split = level_num + 0.3 if exposure == "IMD" & exp_level == "5/Most dep. (ref = 1/Least dep.)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci uci if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"2 (ref = 1/Least dep.)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"2 (ref = 1/Least dep.)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"2 (ref = 1/Least dep.)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"2 (ref = 1/Least dep.)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"2 (ref = 1/Least dep.)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"2 (ref = 1/Least dep.)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"3 (ref = 1/Least dep.)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"3 (ref = 1/Least dep.)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"3 (ref = 1/Least dep.)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"3 (ref = 1/Least dep.)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"3 (ref = 1/Least dep.)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"3 (ref = 1/Least dep.)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"4 (ref = 1/Least dep.)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"4 (ref = 1/Least dep.)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"4 (ref = 1/Least dep.)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"4 (ref = 1/Least dep.)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"4 (ref = 1/Least dep.)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"4 (ref = 1/Least dep.)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", horizontal col(green)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = 1/Least Deprived)") ///
		title("IMD and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.3 0.5 0.7 1 1.5 2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(order(1 "2" 7 "3" 13 "4" 19 "5/Most dep.") rows(1)) ///
		name(imd_cat, replace)
		
graph export ".\G0Mother_Results\imdResults.pdf", replace


** Create plot for housing (ref = owned/mortgaged)

* As three exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.25 if exposure == "housing" & exp_level == "Rent (ref = Own/Mortgage)"
replace level_split = level_num - 0 if exposure == "housing" & exp_level == "Council/HA (ref = Own/Mortgage)"
replace level_split = level_num + 0.25 if exposure == "housing" & exp_level == "Other (ref = Own/Mortgage)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci uci if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = Owned/Mortgaged)") ///
		title("Housing status and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.3 0.5 1.5 2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Rented" 7 "Council/HA" 13 "Other") rows(1)) ///
		name(housing_cat, replace)
		
graph export ".\G0Mother_Results\housingResults.pdf", replace


** Create plot for mobility (ref = 1/least deprived)

* As five exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "mobility" & exp_level == "1 move (ref = 0 moves)"
replace level_split = level_num - 0.15 if exposure == "mobility" & exp_level == "2 moves (ref = 0 moves)"
replace level_split = level_num if exposure == "mobility" & exp_level == "3 moves (ref = 0 moves)"
replace level_split = level_num + 0.15 if exposure == "mobility" & exp_level == "4 moves (ref = 0 moves)"
replace level_split = level_num + 0.3 if exposure == "mobility" & exp_level == "5 + moves (ref = 0 moves)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci uci if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"1 move (ref = 0 moves)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"1 move (ref = 0 moves)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"1 move (ref = 0 moves)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"1 move (ref = 0 moves)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"1 move (ref = 0 moves)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"1 move (ref = 0 moves)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"2 moves (ref = 0 moves)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"2 moves (ref = 0 moves)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"2 moves (ref = 0 moves)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"2 moves (ref = 0 moves)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"2 moves (ref = 0 moves)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"2 moves (ref = 0 moves)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"3 moves (ref = 0 moves)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"3 moves (ref = 0 moves)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"3 moves (ref = 0 moves)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"3 moves (ref = 0 moves)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"3 moves (ref = 0 moves)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"3 moves (ref = 0 moves)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"4 moves (ref = 0 moves)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"4 moves (ref = 0 moves)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"4 moves (ref = 0 moves)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"4 moves (ref = 0 moves)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"4 moves (ref = 0 moves)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"4 moves (ref = 0 moves)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"5 + moves (ref = 0 moves)", col(orange) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"5 + moves (ref = 0 moves)", horizontal col(orange)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"5 + moves (ref = 0 moves)", col(orange) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"5 + moves (ref = 0 moves)", horizontal col(orange)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"5 + moves (ref = 0 moves)", col(orange) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"5 + moves (ref = 0 moves)", horizontal col(orange)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = 0 moves)") ///
		title("Residential mobility and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.5 1 1.5 2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "1" 7 "2" 13 "3" 19 "4" 25 "5+") rows(1)) ///
		name(mobility_cat, replace)

graph export ".\G0Mother_Results\mobilityResults.pdf", replace


graph close _all


****** And now for some interaction plots

** Create plot for education by age interaction (ref = CSE/None)

* As four exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "education" & exp_level == "Vocational (ref = CSE/None)"
replace level_split = level_num - 0.1 if exposure == "education" & exp_level == "O-level (ref = CSE/None)"
replace level_split = level_num + 0.1 if exposure == "education" & exp_level == "A-level (ref = CSE/None)"
replace level_split = level_num + 0.3 if exposure == "education" & exp_level == "Degree (ref = CSE/None)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci_int uci_int if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = CSE/None)") ///
		title("Education*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.85 0.9 0.95 1 1.05 1.1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Vocational" 7 "O-levels" 13 "A-levels" ///
			19 "Degree") rows(1)) ///
		name(edu_int, replace)
		
graph export ".\G0Mother_Results\eduResults_int.pdf", replace
	

** Create plot for occupational social class by age interaction (ref = lower)
sum lci_int uci_int if exposure == "highSocClass" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == ///
			"highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == ///
			"highSocClass", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == ///
			"highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == ///
			"highSocClass", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == ///
			"highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == ///
			"highSocClass", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = lower [III manual/IV/V])") ///
		title("Occupational Class*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.9 0.95 1 1.05, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(socClass_int, replace)
		
graph export ".\G0Mother_Results\socClassResults_int.pdf", replace


** Create plot for income by age interaction
sum lci_int uci_int if exposure == "income" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == ///
			"income", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == ///
			"income", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == ///
			"income", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == ///
			"income", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == ///
			"income", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == ///
			"income", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in log income)") ///
		title("Household income*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.9 0.92 0.94 0.96 0.98 1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 1, valuelabel labsize(small) angle(0)) ///
		legend(off) name(income_int, replace)
		
graph export ".\G0Mother_Results\incomeResults_int.pdf", replace


** Create plot for housing (ref = owned/mortgaged)

* As three exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.25 if exposure == "housing" & exp_level == "Rent (ref = Own/Mortgage)"
replace level_split = level_num - 0 if exposure == "housing" & exp_level == "Council/HA (ref = Own/Mortgage)"
replace level_split = level_num + 0.25 if exposure == "housing" & exp_level == "Other (ref = Own/Mortgage)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci_int uci_int if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = Owned/Mortgaged)") ///
		title("Housing status*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.9 1 1.1 1.2 1.3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Rented" 15 "Council/HA" 29 "Other") rows(1)) ///
		name(housing_int, replace)
		
graph export ".\G0Mother_Results\housingResults_int.pdf", replace


graph close _all



*********************************************************************************
** For the multinomial regression results, as interpretation not intuitive, could convert to predicted probabilities using the 'margins' command? (see: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/). Some - pretty rough - predicted probability plots are below.

* Age only model
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d810 ageAtBirth, rrr baseoutcome(3)

clear
set obs 30
egen ageAtBirth = fill(15 16)
gen d810 = 1
predict p1, outcome(1)
sum p1

replace d810 = 2
predict p2, outcome(2)
sum p1 p2

replace d810 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 ageAtBirth) ///
	(line p2 ageAtBirth) ///
	(line p3 ageAtBirth), ///
	legend(order(1 "Believer" 2 "Not sure believes" 3 "Non-believer")) ///
	name(age, replace)
	
	
* As above, but trying to add CIs to the predicted probability plots
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d810 ageAtBirth, rrr baseoutcome(3)
margins, at(ageAtBirth = (15(1)44))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAtBirth = fill(15 16)
gen prob_yes = .
gen lci_yes = .
gen uci_yes = .
gen prob_notSure = .
gen lci_notSure = .
gen uci_notSure = .
gen prob_no = .
gen lci_no = .
gen uci_no = .

forvalues i = 1(1)`n' {
	replace prob_yes = res[1,`i'] if _n == `i'
	replace lci_yes = res[5,`i'] if _n == `i'
	replace uci_yes = res[6,`i'] if _n == `i'
	replace prob_notSure = res[1,`n' + `i'] if _n == `i'
	replace lci_notSure = res[5,`n' + `i'] if _n == `i'
	replace uci_notSure = res[6,`n' + `i'] if _n == `i'
	replace prob_no = res[1,`n' + `n' + `i'] if _n == `i'
	replace lci_no = res[5,`n' + `n' + `i'] if _n == `i'
	replace uci_no = res[6,`n' + `n' + `i'] if _n == `i'
}

list, clean

twoway (line prob_yes ageAtBirth, col(black)) ///
	(rarea lci_yes uci_yes ageAtBirth, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure ageAtBirth, col(red)) ///
	(rarea lci_notSure uci_notSure ageAtBirth, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no ageAtBirth, col(blue)) ///
	(rarea lci_no uci_no ageAtBirth, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Belief in God/divine power", size(medium)) ///
	legend(order(1 "Believer" 3 "Not sure" 5 "Non-believer") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(age_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religion
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d813_grp ageAtBirth, rrr baseoutcome(3)
margins, at(ageAtBirth = (15(1)44))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAtBirth = fill(15 16)
gen prob_xian = .
gen lci_xian = .
gen uci_xian = .
gen prob_other = .
gen lci_other = .
gen uci_other = .
gen prob_none = .
gen lci_none = .
gen uci_none = .

forvalues i = 1(1)`n' {
	replace prob_xian = res[1,`i'] if _n == `i'
	replace lci_xian = res[5,`i'] if _n == `i'
	replace uci_xian = res[6,`i'] if _n == `i'
	replace prob_other = res[1,`n' + `i'] if _n == `i'
	replace lci_other = res[5,`n' + `i'] if _n == `i'
	replace uci_other = res[6,`n' + `i'] if _n == `i'
	replace prob_none = res[1,`n' + `n' + `i'] if _n == `i'
	replace lci_none = res[5,`n' + `n' + `i'] if _n == `i'
	replace uci_none = res[6,`n' + `n' + `i'] if _n == `i'
}

list, clean

twoway (line prob_xian ageAtBirth, col(black)) ///
	(rarea lci_xian uci_xian ageAtBirth, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other ageAtBirth, col(red)) ///
	(rarea lci_other uci_other ageAtBirth, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none ageAtBirth, col(blue)) ///
	(rarea lci_none uci_none ageAtBirth, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(medium)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(age_relig, replace)
	
	
* Attend church
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d816_rev ageAtBirth, rrr baseoutcome(0)
margins, at(ageAtBirth = (15(1)44))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen ageAtBirth = fill(15 16)
gen prob_no = .
gen lci_no = .
gen uci_no = .
gen prob_yr = .
gen lci_yr = .
gen uci_yr = .
gen prob_mth = .
gen lci_mth = .
gen uci_mth = .
gen prob_wk = .
gen lci_wk = .
gen uci_wk = .

forvalues i = 1(1)`n' {
	replace prob_no = res[1,`i'] if _n == `i'
	replace lci_no = res[5,`i'] if _n == `i'
	replace uci_no = res[6,`i'] if _n == `i'
	replace prob_yr = res[1,`n' + `i'] if _n == `i'
	replace lci_yr = res[5,`n' + `i'] if _n == `i'
	replace uci_yr = res[6,`n' + `i'] if _n == `i'
	replace prob_mth = res[1,`n' + `n' + `i'] if _n == `i'
	replace lci_mth = res[5,`n' + `n' + `i'] if _n == `i'
	replace uci_mth = res[6,`n' + `n' + `i'] if _n == `i'
	replace prob_wk = res[1,`n' + `n' + `n' + `i'] if _n == `i'
	replace lci_wk = res[5,`n' + `n' + `n' + `i'] if _n == `i'
	replace uci_wk = res[6,`n' + `n' + `n' + `i'] if _n == `i'
}

list, clean

twoway (line prob_no ageAtBirth, col(black)) ///
	(rarea lci_no uci_no ageAtBirth, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr ageAtBirth, col(red)) ///
	(rarea lci_yr uci_yr ageAtBirth, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth ageAtBirth, col(blue)) ///
	(rarea lci_mth uci_mth ageAtBirth, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk ageAtBirth, col(green)) ///
	(rarea lci_wk uci_wk ageAtBirth, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Attendance at place of worship", size(medium)) ///
	legend(order(1 "Not at all" 3 "1/yr" 5 "1/mth" 7 "1/wk") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(age_attend, replace)
	
	
* Combine these plots together
graph combine age_bel age_relig age_attend, iscale(0.5) rows(2)

graph export ".\G0Mother_Results\agePredProbs_combined.pdf", replace

graph close _all


* Ethnicity predicted probs (at mean age)
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d810 ageAtBirth nonWhiteEthnic, rrr baseoutcome(3)

margins, atmeans at(nonWhiteEthnic = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)

clear
set obs `n'
egen v1 = fill(1 2)
gen prob = .
gen lci = .
gen uci = .

forvalues i = 1(1)`n' {
	replace prob = res[1,`i'] if v1 == `i'
	replace lci = res[5,`i'] if v1 == `i'
	replace uci = res[6,`i'] if v1 == `i'
}

list, clean

gen v1_grp = v1 + 0.3 if inlist(v1, 1, 3, 5)
replace v1_grp = v1 - 0.3 if inlist(v1, 2, 4, 6)
	
twoway (scatter prob v1_grp if inlist(v1, 1, 3, 5), col(black) msize(vsmall) msym(D)) ///
	(rspike lci uci v1_grp if inlist(v1, 1, 3, 5), col(black)) ///
	(scatter prob v1_grp if inlist(v1, 2, 4, 6), col(red) msize(vsmall) msym(D)) ///
	(rspike lci uci v1_grp if inlist(v1, 2, 4, 6), col(red)), ///
	title("Belief in God/Divine Power - Predicted Probabilities") ///
	xlabel(1.5 "Believe" 3.5 "Not sure believe" 5.5 "Not believe") ///
	xtitle("") ytitle("Predicted probability") xscale(range(1 6)) ///
	legend(order(1 "White" 3 "Other than White")) ///
	name(ethnicity, replace)
	

* Ethnicity by age interaction predicted probs
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d810 c.ageAtBirth##c.nonWhiteEthnic, rrr baseoutcome(3)

clear
set obs 60
egen ageAtBirth = fill(15 15 16 16)
egen nonWhiteEthnic = fill(0 1 0 1)
gen d810 = 1
predict p1, outcome(1)
sum p1

replace d810 = 2
predict p2, outcome(2)
sum p1 p2

replace d810 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 ageAtBirth if nonWhiteEthnic == 0, col(black)) ///
	(line p1 ageAtBirth if nonWhiteEthnic == 1, col(red)) ///
	(line p2 ageAtBirth if nonWhiteEthnic == 0, col(black) lpattern(longdash)) ///
	(line p2 ageAtBirth if nonWhiteEthnic == 1, col(red) lpattern(longdash)) ///
	(line p3 ageAtBirth if nonWhiteEthnic == 0, col(black) lpattern(shortdash)) ///
	(line p3 ageAtBirth if nonWhiteEthnic == 1, col(red) lpattern(shortdash)), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	legend(order(1 "White believer" 2 "Other than White believer" ///
	3 "White not sure" 4 "Other than White not sure" 5 "White non-believer" ///
	6 "Other than White non-believer") cols(2) size(small)) ///
	name(eth_int, replace)
	
graph export ".\G0Mother_Results\athnicAgePredProbs.pdf", replace
	
	
* As above, but trying to add CIs to the predicted probability plots
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d810 c.ageAtBirth##c.nonWhiteEthnic, rrr baseoutcome(3)
margins, at(nonWhiteEthnic = (0 1) ageAtBirth = (15(1)44))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAtBirth = fill(15 15 16 16)
egen nonWhiteEthnic = fill(0 1 0 1)
gen prob_yes = .
gen lci_yes = .
gen uci_yes = .
gen prob_notSure = .
gen lci_notSure = .
gen uci_notSure = .
gen prob_no = .
gen lci_no = .
gen uci_no = .

forvalues i = 1(1)`n' {
	replace prob_yes = res[1,`i'] if _n == `i'
	replace lci_yes = res[5,`i'] if _n == `i'
	replace uci_yes = res[6,`i'] if _n == `i'
	replace prob_notSure = res[1,`n' + `i'] if _n == `i'
	replace lci_notSure = res[5,`n' + `i'] if _n == `i'
	replace uci_notSure = res[6,`n' + `i'] if _n == `i'
	replace prob_no = res[1,`n' + `n' + `i'] if _n == `i'
	replace lci_no = res[5,`n' + `n' + `i'] if _n == `i'
	replace uci_no = res[6,`n' + `n' + `i'] if _n == `i'
}

list, clean

* This is too messy and is really rather horrible to look at...
twoway (line prob_yes ageAtBirth if nonWhiteEthnic == 0, col(black)) ///
	(rarea lci_yes uci_yes ageAtBirth if nonWhiteEthnic == 0, ///
	lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yes ageAtBirth if nonWhiteEthnic == 1, col(red)) ///
	(rarea lci_yes uci_yes ageAtBirth if nonWhiteEthnic == 1, ///
	lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_notSure ageAtBirth if nonWhiteEthnic == 0, col(black) lpattern(longdash)) ///
	(rarea lci_notSure uci_notSure ageAtBirth if nonWhiteEthnic == 0, ///
	lcol(black) lwidth(vthin) fcol(black%20) lpattern(longdash)) ///
	(line prob_notSure ageAtBirth if nonWhiteEthnic == 1, col(red) lpattern(longdash)) ///
	(rarea lci_notSure uci_notSure ageAtBirth if nonWhiteEthnic == 1, ///
	lcol(black) lwidth(vthin) fcol(red%20) lpattern(longdash)) ///
	(line prob_no ageAtBirth if nonWhiteEthnic == 0, col(black) lpattern(longdash)) ///
	(rarea lci_no uci_no ageAtBirth if nonWhiteEthnic == 0, ///
	lcol(black) lwidth(vthin) fcol(black%20) lpattern(longdash)) ///
	(line prob_no ageAtBirth if nonWhiteEthnic == 1, col(red) lpattern(shortdash)) ///
	(rarea lci_no uci_no ageAtBirth if nonWhiteEthnic == 1, ///
	lcol(black) lwidth(vthin) fcol(red%20) lpattern(shortdash)), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	legend(order(1 "White believer" 3 "Other than White believer" ///
	5 "White not sure" 7 "Other than White not sure" 9 "White non-believer" ///
	11 "Other than White non-believer") cols(2) size(small)) ///
	name(eth_int2, replace)


* Also try and make some nice predicted probability plots for continuous exposure income
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum income, d
local inc_25 = r(p25)
local inc_50 = r(p50)
local inc_75 = r(p75)

mlogit d810 ageAtBirth income, rrr baseoutcome(3)
margins, atmeans at(income = (`inc_25' `inc_50' `inc_75'))

matrix res = r(table)
matrix list res

local n = colsof(res)

clear
set obs `n'
egen v1 = fill(1 1 1 2 2 2)
gen prob = .
gen lci = .
gen uci = .

forvalues i = 1(1)`n' {
	replace prob = res[1,`i'] if _n == `i'
	replace lci = res[5,`i'] if _n == `i'
	replace uci = res[6,`i'] if _n == `i'
}

list, clean

gen v1_grp = v1
replace v1_grp = v1 - 0.2 if inlist(_n, 1, 4, 7)
replace v1_grp = v1 + 0.2 if inlist(_n, 3, 6, 9)
	
twoway (scatter prob v1_grp if inlist(_n, 1, 4, 7), col(black) msize(vsmall) msym(D)) ///
	(rspike lci uci v1_grp if inlist(_n, 1, 4, 7), col(black)) ///
	(scatter prob v1_grp if inlist(_n, 2, 5, 8), col(red) msize(vsmall) msym(D)) ///
	(rspike lci uci v1_grp if inlist(_n, 2, 5, 8), col(red)) ///
	(scatter prob v1_grp if inlist(_n, 3, 6, 9), col(blue) msize(vsmall) msym(D)) ///
	(rspike lci uci v1_grp if inlist(_n, 3, 6, 9), col(blue)), ///
	title("Belief in God/Divine Power - Predicted Probabilities") ///
	xlabel(1 "Believe" 2 "Not sure believe" 3 "Not believe") ///
	xtitle("") ytitle("Predicted probability") xscale(range(0.5 3.5)) ///
	legend(order(1 "Lower IQR income" 3 "Median income" 5 "Upper IQR income") ///
	rows(1) size(small)) ///
	name(income, replace)


* And now interaction between age and income
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum income, d
local inc_25 = r(p25)
local inc_50 = r(p50)
local inc_75 = r(p75)

mlogit d810 ageAtBirth income, rrr baseoutcome(3)
margins, at(income = (`inc_25' `inc_50' `inc_75') ageAtBirth = (15(1)44))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAtBirth = fill(15 15 15 16 16 16)
egen income = fill(`inc_25' `inc_50' `inc_75' `inc_25' `inc_50' `inc_75')
gen d810 = 1
predict p1, outcome(1)
sum p1

replace d810 = 2
predict p2, outcome(2)
sum p1 p2

replace d810 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 ageAtBirth if income == `inc_25', col(black)) ///
	(line p1 ageAtBirth if income == `inc_50', col(red)) ///
	(line p1 ageAtBirth if income == `inc_75', col(blue)) ///
	(line p2 ageAtBirth if income == `inc_25', col(black) lpattern(longdash)) ///
	(line p2 ageAtBirth if income == `inc_50', col(red) lpattern(longdash)) ///
	(line p2 ageAtBirth if income == `inc_75', col(blue) lpattern(longdash)) ///
	(line p3 ageAtBirth if income == `inc_25', col(black) lpattern(shortdash)) ///
	(line p3 ageAtBirth if income == `inc_50', col(red) lpattern(shortdash)) ///
	(line p3 ageAtBirth if income == `inc_75', col(blue) lpattern(shortdash)), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	legend(order(1 "Low income - yes" 2 "Median income - yes" ///
	3 "High income - yes" 4 "Low income - not sure" 5 "Median income - not sure" ///
	6 "High income - not sure" 7 "Low income - no" 8 "Median income - no" ///
	9 "High income - no") cols(3) size(vsmall)) ///
	name(inc_int, replace)


** And now interaction between age and education
use ".\G0Mother_Results\G0Mother_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

* Belief in God
mlogit d810 c.ageAtBirth##i.education, rrr baseoutcome(3)
margins i.education, at(ageAtBirth = (20 30 40))
capture drop p1 p2 p3
predict p1 p2 p3
sort ageAtBirth
twoway (line p1 ageAtBirth if education == 1) ///
	(line p1 ageAtBirth if education == 2) ///
	(line p1 ageAtBirth if education == 3) ///
	(line p1 ageAtBirth if education == 4) ///
	(line p1 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Belief in God - Yes", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(yes, replace)
	
twoway (line p2 ageAtBirth if education == 1) ///
	(line p2 ageAtBirth if education == 2) ///
	(line p2 ageAtBirth if education == 3) ///
	(line p2 ageAtBirth if education == 4) ///
	(line p2 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Belief in God - Not sure", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(notSure, replace)
	
twoway (line p3 ageAtBirth if education == 1) ///
	(line p3 ageAtBirth if education == 2) ///
	(line p3 ageAtBirth if education == 3) ///
	(line p3 ageAtBirth if education == 4) ///
	(line p3 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Belief in God - No", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(no, replace)
	
grc1leg yes notSure no, ycommon

graph export ".\G0Mother_Results\eduAgeIntPredProbs_belief.pdf", replace

* Religious affiliation
mlogit d813_grp c.ageAtBirth##i.education, rrr baseoutcome(3)
margins i.education, at(ageAtBirth = (20 30 40))
capture drop p1 p2 p3
predict p1 p2 p3
sort ageAtBirth
twoway (line p1 ageAtBirth if education == 1) ///
	(line p1 ageAtBirth if education == 2) ///
	(line p1 ageAtBirth if education == 3) ///
	(line p1 ageAtBirth if education == 4) ///
	(line p1 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Religious affiliation - Christian", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(Xian, replace)
	
twoway (line p2 ageAtBirth if education == 1) ///
	(line p2 ageAtBirth if education == 2) ///
	(line p2 ageAtBirth if education == 3) ///
	(line p2 ageAtBirth if education == 4) ///
	(line p2 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Religious affiliation - Other", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(other, replace)
	
twoway (line p3 ageAtBirth if education == 1) ///
	(line p3 ageAtBirth if education == 2) ///
	(line p3 ageAtBirth if education == 3) ///
	(line p3 ageAtBirth if education == 4) ///
	(line p3 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Religious affiliation - None", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(none, replace)
	
grc1leg Xian other none, ycommon

graph export ".\G0Mother_Results\eduAgeIntPredProbs_relig.pdf", replace

* Chuch attendance
mlogit d816_rev c.ageAtBirth##i.education, rrr baseoutcome(0)
margins i.education, at(ageAtBirth = (20 30 40))
capture drop p1 p2 p3
predict p1 p2 p3 p4
sort ageAtBirth
twoway (line p1 ageAtBirth if education == 1) ///
	(line p1 ageAtBirth if education == 2) ///
	(line p1 ageAtBirth if education == 3) ///
	(line p1 ageAtBirth if education == 4) ///
	(line p1 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Church attendance - Not at all", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(never, replace)
	
twoway (line p2 ageAtBirth if education == 1) ///
	(line p2 ageAtBirth if education == 2) ///
	(line p2 ageAtBirth if education == 3) ///
	(line p2 ageAtBirth if education == 4) ///
	(line p2 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Church attendance - Min 1/year", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(year, replace)
	
twoway (line p3 ageAtBirth if education == 1) ///
	(line p3 ageAtBirth if education == 2) ///
	(line p3 ageAtBirth if education == 3) ///
	(line p3 ageAtBirth if education == 4) ///
	(line p3 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Church attendance - Min 1/month", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(month, replace)
	
twoway (line p4 ageAtBirth if education == 1) ///
	(line p4 ageAtBirth if education == 2) ///
	(line p4 ageAtBirth if education == 3) ///
	(line p4 ageAtBirth if education == 4) ///
	(line p4 ageAtBirth if education == 5), ///
	xscale(range(13 47)) xlabel(15(5)45, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age at birth") ytitle("Predicted probability") ///
	title("Church attendance - Min 1/week", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(week, replace)
	
grc1leg never year month week, ycommon

graph export ".\G0Mother_Results\eduAgeIntPredProbs_attend.pdf", replace

graph close _all
	

********************************************************************************	
*** Example code demonstrating collider bias (figure S28)
clear
set obs 64

* Create 'age' variable (0 = under 50;  1 = over 50)
egen age = fill(0 1 0 1)
tab age

* Create 'religious belief' variable (0 = no belief; 1 = belief; 75% of people aged over 50 believe in God, but only 25% of people aged under 50 believe in God)
egen belief = fill(0 0 0 1 0 1 1 1 0 0 0 1 0 1 1 1)
tab belief
tab belief age, col

* Create 'married' variable (0 - not married; 1 = married; 25% of people aged under 50 who do not believe in God are married; 50% of under 50s who believe are married; 50% of over 50s who believe are married; and 75% of over 50s who believe are married)
egen married = fill(0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 ///	
	0 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1 ///	
	0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 ///	
	0 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1)
tab married
table (married) (age belief), nototals
table (married) (age belief), statistic(percent, across(married)) nototals

** Calculate 'true' odds ratio (OR = 9) of age on RSBB (excluding 'married', as we know that 'married' is a collider, as it is caused by both 'age' and 'religious belief') - Will calculate both manually and via logistic regression
tab belief age, col
display (0.75/0.25) / (0.25/0.75) 
logistic belief age

* Now, if we include the collider 'married' in this model, the odds ratio is biased to 6.75
logistic belief age married

* Can also show this manually within each strata of 'married'
tab belief age if married == 1, col
display (0.6/0.4) / (0.1818/0.8182)
logistic belief age if married == 1

tab belief age if married == 0, col
display (0.8182/0.1818) / (0.4/0.6)
logistic belief age if married == 0


*** Can also bias away from the null if age and belief go in opposite directions for the collider (say, age increases the probability of being married but belief in God decreases it)
egen married2 = fill(0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ///	
	1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 ///	
	0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ///	
	1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1)
tab married2
logistic belief age
logistic belief age married2

clear
	