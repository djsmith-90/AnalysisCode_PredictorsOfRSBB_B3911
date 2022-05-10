*** Predictors of RSBB (B3911) - G0 partner/father analysis script
*** Created 23/11/2021 by Dan Smith
*** Stata v17.0

*** This script reads in the cleaned G0 partner/father data, explores associations between exposures, and then conducts an analysis exploring how all the demographic and SES variables are associated with various facets of RSBB.


**********************************************************************************
**** Set working directory, start a log file, and read in cleaned dataset

cd "X:\Groups\ARC\DanS\Descriptive_PredictorsOfRSBB_B3911"

capture log close
log using ".\G0Partner_Results\Desc_RSBB_B3911_G0PartnerAnalysis_log", replace text

use "G0Partner_PredictorsOfRSBB_B3911.dta", clear


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
order aln pb150 pb153 pb153_grp pb155 FC3000 FC3040 FC3040_grp FC3080 FC3080_OccNever FC3080_OccYr FC3153 FC3153_cat FC3160 FC3170 FC3155 FC3155_cat

* Keep just pregnancy RSBB variables (belief in God, religious affiliation and frequency of church attendance)
drop FC3000 FC3040 FC3040_grp FC3080 FC3080_OccNever FC3080_OccYr FC3153 FC3153_cat FC3160 FC3170 FC3155 FC3155_cat

* Also drop the cognitive/psychological variables (as will focus on these in another paper), and also the adverse childhood experiences variables
drop pb546-pb551 pa782 esteem_prorated
drop pb481 pb482 

** Using the 'distinct' command (see above), for each variable inspect the number of unque values; if < 10 then display table, while if >= 10 then displays means/SDs/IQRs

* Outcomes
foreach var of varlist pb150-pb155 {
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
foreach var of varlist a053-neighbour_qual {
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
misstable sum pb150-pb155, all

* Exposures
misstable sum a053-neighbour_qual, all



**********************************************************************************
*** Correlations between exposures

** Explore correlations between exposures to see how inter-related these factors are
desc a053-neighbour_qual

* Most variables are either continuous, ordered categories, or binary variables, so will just use standard pearson correlations for these. Only unordered categories are home ownership (owned/mortaged vs renting vs counci/housing association vs other) and marital status (never married vs married vs widowed/divorced/separated). So will exclude these two variables from the correlation matrix and calculate their associations using these variables as outcomes in a multinomial regression and square-rooting the pseudo-R2 value (similar to how the 'rexposome' package in R for ExWAS analyses does)

* First, order variables into categories of 'demographic' and 'socioeconomic/material insecurity' (FC9992 is age at @28 RSBB questionnaire, so will pop at end as not needed here)
order aln-pb155 ///
pb910 c801_grp pa065_grp pa005_grp jan1993ur01ind_grp b032_grp ///
c666a pb359a pb376a c765_grp pb_sc_pgm_grp pb_sc_pgf_grp logavinceq jan1993imd2010q5_M jan1993Townsendq5_M a006_grp pb184_grp pd685 pb479a a053 a551 neighbour_qual ///
FC9992

* Next, rename all these exposures so they are more intuitive and can be read easier on the correlation heatmaps below
rename pb910 ageInPreg
rename c801_grp nonWhiteEthnic
rename pa065_grp maritalStatus
rename pa005_grp mobility 
rename jan1993ur01ind_grp rural
rename b032_grp parity
rename c666a education
rename pb359a maternalEdu
rename pb376a paternalEdu
rename c765_grp highSocClass
rename pb_sc_pgm_grp highSocClass_mat
rename pb_sc_pgf_grp highSocClass_pat
rename logavinceq income
rename jan1993imd2010q5_M IMD
rename jan1993Townsendq5_M townsendDep
rename pb479a poorerChildhood
rename a053 accessToCar
rename a006_grp housing 
rename pb184_grp financeDiffs
rename pd685 financeDiffsScore
rename a551 crowding
rename neighbour_qual neighPercept
rename FC9992 ageAt28

* Associations between demographic factors (excluding marital status; pa065_grp) - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
corr ageInPreg nonWhiteEthnic mobility rural parity

matrix cor_demo = r(C)
matrix list cor_demo

heatplot cor_demo, values(format(%9.3f)) color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) legend(off) xlabel(, angle(45))

* Save heatmap
graph export ".\G0Partner_Results\corr_heatplot_demoOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Partner_Results\corrMatrix_demoOnly.xlsx", replace
putexcel A1=matrix(cor_demo), names nformat(number_d2)

* And now for associations of demographic variables with unordered categorical variable marital status (and save to a CSV file)
capture postclose marital_corrs_demoOnly
postfile marital_corrs_demoOnly str30 variable corr ///
	using ".\G0Partner_Results\marital_corrs_demoOnly.dta", replace

foreach var of varlist ageInPreg nonWhiteEthnic mobility rural parity {
	quietly mlogit maritalStatus `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on marital status is: " round((sqrt(`r2')), .001)
	
	post marital_corrs_demoOnly ("`var'") (sqrt(`r2'))
}

postclose marital_corrs_demoOnly

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G0Partner_Results\marital_corrs_demoOnly.dta", clear
list, clean
outsheet using ".\G0Partner_Results\marital_corrs_demoOnly.csv", comma replace
restore


** Now repeat for socioeconomic/material insecurity variables (to exclude housing status [housing], as is an unordered categorical variable)
corr education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept

matrix cor_socio = r(C)

* As lots of entries, the correlation coefficients are hard to read, so will drop these values and include the legend
heatplot cor_socio, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))
	
* Save heatmap
graph export ".\G0Partner_Results\corr_heatplot_socioOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Partner_Results\corrMatrix_socioOnly.xlsx", replace
putexcel A1=matrix(cor_socio), names nformat(number_d2)


* And now for associations of socioeconomic variables with unordered categorical variable home ownership status (and save to a CSV file)
capture postclose housing_corrs_socioOnly
postfile housing_corrs_socioOnly str30 variable corr ///
	using ".\G0Partner_Results\housing_corrs_socioOnly.dta", replace
	
foreach var of varlist education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept {
	quietly mlogit housing `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on home ownership status is: " round((sqrt(`r2')), .001)

	post housing_corrs_socioOnly ("`var'") (sqrt(`r2'))
}

postclose housing_corrs_socioOnly

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G0Partner_Results\housing_corrs_socioOnly.dta", clear
list, clean
outsheet using ".\G0Partner_Results\housing_corrs_socioOnly.csv", comma replace
restore

	
*** Finally, repeat this on all of the exposures together (excluding unordered cateogorical variables housing status and marital status)
corr ageInPreg nonWhiteEthnic mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept

matrix cor_all = r(C)

heatplot cor_all, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle(""))

* Save heatmap
graph export ".\G0Partner_Results\corr_heatplot_all.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Partner_Results\corrMatrix_all.xlsx", replace
putexcel A1=matrix(cor_all), names nformat(number_d2)

	
* And now for associations of all other exposures variables with unordered categorical variables marital status and home ownership status

* Marital status
capture postclose marital_corrs_all
postfile marital_corrs_all str30 variable corr ///
	using ".\G0Partner_Results\marital_corrs_all.dta", replace
	
foreach var of varlist ageInPreg nonWhiteEthnic mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep housing financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept {
	quietly mlogit maritalStatus `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on marital status is: " round((sqrt(`r2')), .001)

	post marital_corrs_all ("`var'") (sqrt(`r2'))
}

postclose marital_corrs_all

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G0Partner_Results\marital_corrs_all.dta", clear
list, clean
outsheet using ".\G0Partner_Results\marital_corrs_all.csv", comma replace
restore


* Housing status
capture postclose housing_corrs_all
postfile housing_corrs_all str30 variable corr ///
	using ".\G0Partner_Results\housing_corrs_all.dta", replace
	
foreach var of varlist ageInPreg nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept {
	quietly mlogit housing `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on home ownership status is: " round((sqrt(`r2')), .001)


	post housing_corrs_all ("`var'") (sqrt(`r2'))
}

postclose housing_corrs_all

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G0Partner_Results\housing_corrs_all.dta", clear
list, clean
outsheet using ".\G0Partner_Results\housing_corrs_all.csv", comma replace
restore



************************************************************************************
*** Next, we want to run the actual analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab pb150, m

** This will be quite complicated, as want to post results to file, but as exposures differ extracting the results will be variable specific. To adjust for multiple corrections will use conservative bonferroni adjustment when constructing confidence intervals and interpreting p-values - As 22 exposures, a Bonferroni p-value threshold of 0.05/22 = 0.0023.
display 0.05/22

** We also want to store both estimates adjusting for age (other than for the age-only model), and also the interaction between age and the exposure, to see whether it's moderated by age. Again, this makes the set-up a bit more complicated.

** Create a postfile to post results to, then start the loop - Will create three postfiles; one for coefficients and CIs, another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model), and then a third for the (pseduo-)R2 values to assess model fit - NOTE: Have to store pvalues as 'double' format, else really tiny p-values get coded as '0' (as default if float format, which has minimum value of -3.40282346639e+38 [https://blog.stata.com/2012/04/02/the-penultimate-guide-to-precision/]).
capture postclose partner_belief
postfile partner_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_belief_results.dta", replace

capture postclose partner_belief_lr
postfile partner_belief_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_belief_results_lr.dta", replace
	
capture postclose partner_belief_r2
postfile partner_belief_r2 str30 exposure r2_main r2_int ///
	using ".\G0Partner_Results\partner_belief_results_r2.dta", replace

foreach var of varlist ageInPreg nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep housing financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageInPreg' separately first as will be adjusted for in all other models
	if "`var'" == "ageInPreg" {
		mlogit pb150 `var', baseoutcome(3) rrr
		
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
		
		post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit pb150 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb150 `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb150 `var', baseoutcome(3) rrr
		local r2_main = e(r2_p)
		
		// As no interaction model for age, will just fill with missing values
		local lr_p_int = .
		local r2_int = .
		
		post partner_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_belief_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "highSocClass_mat" | "`var'" == "highSocClass_pat" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "financeDiffsScore" | "`var'" == "poorerChildhood" | "`var'" == "accessToCar" | "`var'" == "neighPercept" {
		
		mlogit pb150 ageInPreg `var', baseoutcome(3) rrr
		
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
		mlogit pb150 c.ageInPreg##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/not sure)
		mlogit pb150 ageInPreg `var', baseoutcome(3) rrr
		
		local outcome_level = "Not sure (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
				
		// Now for interaction model
		mlogit pb150 c.ageInPreg##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit pb150 ageInPreg if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb150 ageInPreg `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb150 ageInPreg if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit pb150 ageInPreg `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit pb150 c.ageInPreg##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit pb150 c.ageInPreg##c.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post partner_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_belief_r2 ("`exp'") (`r2_main') (`r2_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Move to the next category of the exposure (category 3)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,17]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures and R2 values
		mlogit pb150 ageInPreg if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb150 ageInPreg if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post partner_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_belief_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose partner_belief
postclose partner_belief_lr
postclose partner_belief_r2
	

***********************************************************************************
*** Now to the next RSBB outcome: Religious affiliation

*** As this is an unordered categorical variable, will again use multinomial logistic model (with 'none' as reference)
tab pb153_grp

* Do need to re-order the categories so can simply copy and paste the script from above without having to faff around with editing the cells to take the statistics from
recode pb153_grp (1 = 3) (2 = 1) (3 = 2)
label define relig_lb 1 "Christian" 2 "Other" 3 "None", modify
numlabel relig_lb, add
tab pb153_grp


*** Now run the loop to save all the results
capture postclose partner_relig
postfile partner_relig str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_relig_results.dta", replace

capture postclose partner_relig_lr
postfile partner_relig_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_relig_results_lr.dta", replace
	
capture postclose partner_relig_r2
postfile partner_relig_r2 str30 exposure r2_main r2_int ///
	using ".\G0Partner_Results\partner_relig_results_r2.dta", replace

foreach var of varlist ageInPreg nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep housing financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageInPreg' separately first as will be adjusted for in all other models
	if "`var'" == "ageInPreg" {
		mlogit pb153_grp `var', baseoutcome(3) rrr
		
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
		
		post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit pb153_grp if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb153_grp `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb153_grp `var', baseoutcome(3) rrr
		local r2_main = e(r2_p)
		
		// As no interaction model for age, will just fill with missing values
		local lr_p_int = .
		local r2_int = .
		
		post partner_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_relig_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "highSocClass_mat" | "`var'" == "highSocClass_pat" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "financeDiffsScore" | "`var'" == "poorerChildhood" | "`var'" == "accessToCar" | "`var'" == "neighPercept" {
		
		mlogit pb153_grp ageInPreg `var', baseoutcome(3) rrr
		
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
		mlogit pb153_grp c.ageInPreg##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/Other)
		mlogit pb153_grp ageInPreg `var', baseoutcome(3) rrr
		
		local outcome_level = "Other (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
				
		// Now for interaction model
		mlogit pb153_grp c.ageInPreg##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit pb153_grp ageInPreg if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb153_grp ageInPreg `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb153_grp ageInPreg if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit pb153_grp ageInPreg `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit pb153_grp c.ageInPreg##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit pb153_grp c.ageInPreg##c.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post partner_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_relig_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Move to the next category of the exposure (category 3)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,17]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
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
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures nd R2 values
		mlogit pb153_grp ageInPreg if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb153_grp ageInPreg if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post partner_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_relig_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose partner_relig
postclose partner_relig_lr
postclose partner_relig_r2


************************************************************************************
*** Now to the next RSBB outcome: Attendance at church/place of worship

*** As this is an ordered categorical variable, originally planned to use ordinal regression models. However, as the G0 mother data violated the proportional odds assumption (as does the G0 partner/father data; see below), will just use multinomial models regardless so that results are comparable.

** For consistency with previous models, will recode so that higher values indicate greater RSBB
tab pb155

recode pb155 (1 = 3) (2 = 2) (3 = 1) (4 = 0), gen(pb155_rev)
label define attend_rev_lb 0 "Not at all" 1 "Min once a year" 2 "Min once a month" 3 "Min once a week"
numlabel attend_rev_lb, add
label values pb155_rev attend_rev_lb
tab pb155_rev


** Quick test of whether proportional odds assumption been violated in most basic model (with just age at birth). Ah, it has been violated, and gets worse if add in additional predictors. 
ologit pb155_rev ageInPreg, or
brant, detail

ologit pb155_rev ageInPreg i.IMD, or
brant, detail


** So instead will just run multinomial models with 'not at all' as the baseline/reference category (this also means all outcomes are using the same model, making them easier to compare).


*** Now run the loop to save all the results
capture postclose partner_attend
postfile partner_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_attend_results.dta", replace

capture postclose partner_attend_lr
postfile partner_attend_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_attend_results_lr.dta", replace
	
capture postclose partner_attend_r2
postfile partner_attend_r2 str30 exposure r2_main r2_int ///
	using ".\G0Partner_Results\partner_attend_results_r2.dta", replace

foreach var of varlist ageInPreg nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu paternalEdu highSocClass highSocClass_mat highSocClass_pat income IMD townsendDep housing financeDiffs financeDiffsScore poorerChildhood accessToCar crowding neighPercept {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageInPreg' separately first as will be adjusted for in all other models
	if "`var'" == "ageInPreg" {
		mlogit pb155_rev `var', baseoutcome(0) rrr
		
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
		
		post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit pb155_rev if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit pb155_rev `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb155_rev `var', baseoutcome(0) rrr
		local r2_main = e(r2_p)
		
		// As no interaction model for age, will just fill with missing values
		local lr_p_int = .
		loca r2_int = .
		
		post partner_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_attend_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "highSocClass_mat" | "`var'" == "highSocClass_pat" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "financeDiffsScore" | "`var'" == "poorerChildhood" | "`var'" == "accessToCar" | "`var'" == "neighPercept" {
		
		mlogit pb155_rev ageInPreg `var', baseoutcome(0) rrr
		
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
		mlogit pb155_rev c.ageInPreg##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/once month)
		mlogit pb155_rev ageInPreg `var', baseoutcome(0) rrr
		
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit pb155_rev c.ageInPreg##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/once week)
		mlogit pb155_rev ageInPreg `var', baseoutcome(0) rrr
		
		local outcome_level = "Min once week (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit pb155_rev c.ageInPreg##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit pb155_rev ageInPreg if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit pb155_rev ageInPreg `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb155_rev ageInPreg if `var' != ., baseoutcome(0) rrr
		local r2_base = e(r2_p)
		mlogit pb155_rev ageInPreg `var', baseoutcome(0) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit pb155_rev c.ageInPreg##c.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit pb155_rev c.ageInPreg##c.`var', baseoutcome(0) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post partner_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_attend_r2 ("`exp'") (`r2_main') (`r2_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (2/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "maternal_edu" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "paternal_edu" {
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
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
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post partner_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures and R2 values
		mlogit pb155_rev ageInPreg if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit pb155_rev ageInPreg if `var' != ., baseoutcome(0) rrr
		local r2_base = e(r2_p)
		mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post partner_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post partner_attend_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose partner_attend
postclose partner_attend_lr
postclose partner_attend_r2



* And close the log file
log close

* Save this file if we want to use it later
save ".\G0Partner_Results\G0Partner_PredictorsOfRSBB_B3911_postAnalysis.dta", replace



***********************************************************************************
***********************************************************************************
***** Next step is to make some nice plots of these results

**** Start with outcome of belief in God/divine power 
use ".\G0Partner_Results\partner_belief_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 22 exposures, will do 0.05/22)
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Main effect") ///
	name(belief_main, replace)
	
graph export ".\G0Partner_Results\belief_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/21)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Age interaction") ///
	name(belief_int, replace)
	
graph export ".\G0Partner_Results\belief_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(belief_both, replace)

graph export ".\G0Partner_Results\belief_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\G0Partner_Results\belief_pvalues.dta", replace


**** Now read in the next outcome - religious affiliation
use ".\G0Partner_Results\partner_relig_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 22 exposures, will do 0.05/22)
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Main effect") ///
	name(relig_main, replace)
	
graph export ".\G0Partner_Results\relig_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/21)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Age interaction") ///
	name(relig_int, replace)
	
graph export ".\G0Partner_Results\relig_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(relig_both, replace)

graph export ".\G0Partner_Results\relig_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\G0Partner_Results\relig_pvalues.dta", replace


**** Now read in the next outcome - church attendance
use ".\G0Partner_Results\partner_attend_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 22 exposures, will do 0.05/22)
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Main effect") ///
	name(attend_main, replace)
	
graph export ".\G0Partner_Results\attend_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/21)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Age interaction") ///
	name(attend_int, replace)
	
graph export ".\G0Partner_Results\attend_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/22)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(attend_both, replace)

graph export ".\G0Partner_Results\attend_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'church attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\G0Partner_Results\attend_pvalues.dta", replace


*** Combine all these datasets together
use ".\G0Partner_Results\belief_pvalues.dta", clear
append using ".\G0Partner_Results\relig_pvalues.dta"
append using ".\G0Partner_Results\attend_pvalues.dta"


** Now look at combined results

* Pregnancy vars main effects
local bon_thresh = -log10(0.05/22)
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
	ylabel(1(1)22, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(preg_main, replace)

graph export ".\G0Partner_Results\beliefReligAttend_mainEffects_pvalues.pdf", replace

* Pregnancy vars interaction effects
local bon_thresh = -log10(0.05/21)
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
	ylabel(2(1)22, valuelabel labsize(small) angle(0)) ///
	title("Age interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(preg_int, replace)

graph export ".\G0Partner_Results\beliefReligAttend_ageInt_pvalues.pdf", replace


** Combine all these graphs together
graph combine preg_main preg_int, ysize(3) xsize(6)

graph export ".\G0Partner_Results\allData_pvalues.pdf", replace

* And save as .EPS format
graph combine preg_main preg_int, ysize(10) xsize(20)
graph export ".\G0Partner_Results\allData_pvalues.eps", replace

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

outsheet exp_num-lr_p_intAttend using ".\G0Partner_Results\pvalue_results.csv", comma replace


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
use ".\G0Partner_Results\partner_belief_results_r2.dta", clear

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
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Main effect") ///
	name(belief_main, replace)
	
graph export ".\G0Partner_Results\belief_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Age interaction") ///
	name(belief_int, replace)
	
graph export ".\G0Partner_Results\belief_ageInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(belief_both, replace)

graph export ".\G0Partner_Results\belief_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\G0Partner_Results\belief_r2.dta", replace


**** Now read in the next outcome - religious affiliation
use ".\G0Partner_Results\partner_relig_results_r2.dta", clear

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
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Main effect") ///
	name(relig_main, replace)
	
graph export ".\G0Partner_Results\relig_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Age interaction") ///
	name(relig_int, replace)
	
graph export ".\G0Partner_Results\relig_ageInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(relig_both, replace)

graph export ".\G0Partner_Results\relig_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\G0Partner_Results\relig_r2.dta", replace


**** Now read in the next outcome - church attendance
use ".\G0Partner_Results\partner_attend_results_r2.dta", clear

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
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Main effect") ///
	name(attend_main, replace)
	
graph export ".\G0Partner_Results\attend_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Age interaction") ///
	name(attend_int, replace)
	
graph export ".\G0Partner_Results\attend_ageInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(attend_both, replace)

graph export ".\G0Partner_Results\attend_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'church attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\G0Partner_Results\attend_r2.dta", replace


*** Combine all these datasets together
use ".\G0Partner_Results\belief_r2.dta", clear
append using ".\G0Partner_Results\relig_r2.dta"
append using ".\G0Partner_Results\attend_r2.dta"


** Now look at combined results

* Pregnancy vars main effects
twoway (scatter exp_num r2_main if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)22, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(preg_main, replace)

graph export ".\G0Partner_Results\beliefReligAttend_mainEffects_r2.pdf", replace

* Pregnancy vars interaction effects
twoway (scatter exp_num r2_int if outcome == "Belief" & exp_num != 1, ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Religious affil." & exp_num != 1, ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Church attendance" & exp_num != 1, ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)22, valuelabel labsize(small) angle(0)) ///
	title("Age interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(preg_int, replace)

graph export ".\G0Partner_Results\beliefReligAttend_ageInt_r2.pdf", replace


** Combine all these graphs together
graph combine preg_main preg_int, ysize(3) xsize(6)

graph export ".\G0Partner_Results\allData_r2.pdf", replace

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

outsheet exp_num-r2_intAttend using ".\G0Partner_Results\r2_results.csv", comma replace



********************************************************************************
*** Next, want to plot some of the actual coefficient results

** Will read the datasets in then combine together into one single dataset
use ".\G0Partner_Results\partner_belief_results.dta", clear
gen outcome = "Belief"

append using ".\G0Partner_Results\partner_relig_results.dta"
replace outcome = "Relig" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_attend_results.dta"
replace outcome = "Attend" if outcome == ""
tab outcome, m


** Save these results as CSV files to add to the SI

* Save each result in turn
format coef lci uci coef_int lci_int uci_int %9.3f
format p p_int %9.4f

outsheet exposure-p_int using ".\G0Partner_Results\belief_coefs.csv" if outcome == "Belief", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\relig_coefs.csv" if outcome == "Relig", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\attend_coefs.csv" if outcome == "Attend", comma replace

* Convert format back to default format (so axis on plots display correctly)
format coef lci uci coef_int lci_int uci_int %9.0g
format p p_int %10.0g


** First, make a plot for the age results - As some outcomes are on different scales, will just use the results from multinomial regression for all outcomes (inc. intrinsic and total/DUREL religiosity, even though also ran linear regressions on these as well) - Having all variables on the same plot on the same scale makes things easier to visualise.
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

sum lci uci if (exposure == "ageInPreg" | exposure == "ageAt28") & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "ageInPreg", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "ageInPreg", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "ageInPreg", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "ageInPreg", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "ageInPreg", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "ageInPreg", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Age and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(1 1.02 1.04 1.06 1.08, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(age_cat, replace)
		
graph export ".\G0Partner_Results\ageResults.pdf", replace


** Create plot for ethnicity (ref = white)
sum lci uci if exposure == "nonWhiteEthnic" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "nonWhiteEthnic", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "nonWhiteEthnic", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "nonWhiteEthnic", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "nonWhiteEthnic", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "nonWhiteEthnic", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "nonWhiteEthnic", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = White)") ///
		title("Other than White ethnicity and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(1 2 5 10 20, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(ethnic_cat, replace)
		
graph export ".\G0Partner_Results\ethnicityResults.pdf", replace


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
		xlabel(0.5 1 2 3 5 10, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Married" 7 "Widowed/Divorced/Separated")) ///
		name(marital_cat, replace)
		
graph export ".\G0Partner_Results\maritalStatusResults.pdf", replace


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
		xlabel(0.5 1 2 3 5 10, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Vocational" 7 "O-levels" 13 "A-levels" ///
			19 "Degree") rows(1)) ///
		name(edu_cat, replace)
		
graph export ".\G0Partner_Results\eduResults.pdf", replace
	

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
		xlabel(0.8 1 1.5 2 3 4, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(socClass_cat, replace)
		
graph export ".\G0Partner_Results\socClassResults.pdf", replace


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
		xlabel(1 1.5 2 3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(income_cat, replace)
		
graph export ".\G0Partner_Results\incomeResults.pdf", replace


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
		xlabel(0.35 0.5 0.7 1 1.5 2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "2" 7 "3" 13 "4" 19 "5/Most dep.") rows(1)) ///
		name(imd_cat, replace)
		
graph export ".\G0Partner_Results\imdResults.pdf", replace


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
		xlabel(0.25 0.5 1 1.5 2 3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Rented" 7 "Council/HA" 13 "Other") ///
			rows(1)) ///
		name(housing_cat, replace)
		
graph export ".\G0Partner_Results\housingResults.pdf", replace


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
		xlabel(0.5 1 2 3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "1" 7 "2" 13 "3" 19 "4" 25 "5+") rows(1)) ///
		name(mobility_cat, replace)

graph export ".\G0Partner_Results\mobilityResults.pdf", replace


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
		xlabel(0.8 0.85 0.9 0.95 1 1.05 1.1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Vocational" 7 "O-levels" 13 "A-levels" ///
			19 "Degree") rows(1)) ///
		name(edu_int, replace)
		
graph export ".\G0Partner_Results\eduResults_int.pdf", replace
	

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
		xlabel(0.94 0.96 0.98 1 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(socClass_int, replace)
		
graph export ".\G0Partner_Results\socClassResults_int.pdf", replace


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
		xlabel(0.94 0.96 0.98 1 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(income_int, replace)
		
graph export ".\G0Partner_Results\incomeResults_int.pdf", replace


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
		xlabel(0.95 1 1.05 1.1 1.15, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Rented" 7 "Council/HA" 13 "Other") ///
			rows(1)) ///
		name(housing_int, replace)
		
graph export ".\G0Partner_Results\housingResults_int.pdf", replace


graph close _all



*********************************************************************************
** For the multinomial regression results, as interpretation not intuitive, could convert to predicted probabilities using the 'margins' command? (see: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/). Some - pretty rough - predicted probability plots are below.

* Age only model (plot with 95% CIs)
use ".\G0Partner_Results\G0Partner_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit pb150 ageInPreg, rrr baseoutcome(3)
margins, at(ageInPreg = (15(1)70))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageInPreg = fill(15 16)
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

twoway (line prob_yes ageInPreg, col(black)) ///
	(rarea lci_yes uci_yes ageInPreg, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure ageInPreg, col(red)) ///
	(rarea lci_notSure uci_notSure ageInPreg, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no ageInPreg, col(blue)) ///
	(rarea lci_no uci_no ageInPreg, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Yes" 3 "Not sure" 5 "No") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(age_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religion
use ".\G0Partner_Results\G0Partner_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit pb153_grp ageInPreg, rrr baseoutcome(3)
margins, at(ageInPreg = (15(1)70))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageInPreg = fill(15 16)
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

twoway (line prob_xian ageInPreg, col(black)) ///
	(rarea lci_xian uci_xian ageInPreg, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other ageInPreg, col(red)) ///
	(rarea lci_other uci_other ageInPreg, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none ageInPreg, col(blue)) ///
	(rarea lci_none uci_none ageInPreg, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(age_relig, replace)
	
	
* Attend church
use ".\G0Partner_Results\G0Partner_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit pb155_rev ageInPreg, rrr baseoutcome(0)
margins, at(ageInPreg = (15(1)70))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen ageInPreg = fill(15 16)
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

twoway (line prob_no ageInPreg, col(black)) ///
	(rarea lci_no uci_no ageInPreg, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr ageInPreg, col(red)) ///
	(rarea lci_yr uci_yr ageInPreg, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth ageInPreg, col(blue)) ///
	(rarea lci_mth uci_mth ageInPreg, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk ageInPreg, col(green)) ///
	(rarea lci_wk uci_wk ageInPreg, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(13 72)) xlabel(15(5)72, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Not at all" 3 "1/yr" 5 "1/mth" 7 "1/wk") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(age_attend, replace)
	
	
* Combine these plots together
graph combine age_bel age_relig age_attend, iscale(0.5) rows(2)

graph export ".\G0Partner_Results\agePredProbs_combined.pdf", replace

graph close _all


** And now interaction between age and education
use ".\G0Partner_Results\G0Partner_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

* Belief in God
mlogit pb150 c.ageInPreg##i.education, rrr baseoutcome(3)
margins i.education, at(ageInPreg = (20 30 40))
capture drop p1 p2 p3
predict p1 p2 p3
sort ageInPreg
twoway (line p1 ageInPreg if education == 1) ///
	(line p1 ageInPreg if education == 2) ///
	(line p1 ageInPreg if education == 3) ///
	(line p1 ageInPreg if education == 4) ///
	(line p1 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious belief - Yes", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(yes, replace)
	
twoway (line p2 ageInPreg if education == 1) ///
	(line p2 ageInPreg if education == 2) ///
	(line p2 ageInPreg if education == 3) ///
	(line p2 ageInPreg if education == 4) ///
	(line p2 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious belief - Not sure", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(notSure, replace)
	
twoway (line p3 ageInPreg if education == 1) ///
	(line p3 ageInPreg if education == 2) ///
	(line p3 ageInPreg if education == 3) ///
	(line p3 ageInPreg if education == 4) ///
	(line p3 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious belief - No", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(no, replace)
	
grc1leg yes notSure no, ycommon

graph export ".\G0Partner_Results\eduAgeIntPredProbs_belief.pdf", replace

* Religious affiliation
mlogit pb153_grp c.ageInPreg##i.education, rrr baseoutcome(3)
margins i.education, at(ageInPreg = (20 30 40))
capture drop p1 p2 p3
predict p1 p2 p3
sort ageInPreg
twoway (line p1 ageInPreg if education == 1) ///
	(line p1 ageInPreg if education == 2) ///
	(line p1 ageInPreg if education == 3) ///
	(line p1 ageInPreg if education == 4) ///
	(line p1 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious affiliation - Christian", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(Xian, replace)
	
twoway (line p2 ageInPreg if education == 1) ///
	(line p2 ageInPreg if education == 2) ///
	(line p2 ageInPreg if education == 3) ///
	(line p2 ageInPreg if education == 4) ///
	(line p2 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious affiliation - Other", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(other, replace)
	
twoway (line p3 ageInPreg if education == 1) ///
	(line p3 ageInPreg if education == 2) ///
	(line p3 ageInPreg if education == 3) ///
	(line p3 ageInPreg if education == 4) ///
	(line p3 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious affiliation - None", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(none, replace)
	
grc1leg Xian other none, ycommon

graph export ".\G0Partner_Results\eduAgeIntPredProbs_relig.pdf", replace

* Church attendance
mlogit pb155_rev c.ageInPreg##i.education, rrr baseoutcome(0)
margins i.education, at(ageInPreg = (20 30 40))
capture drop p1 p2 p3
predict p1 p2 p3 p4
sort ageInPreg
twoway (line p1 ageInPreg if education == 1) ///
	(line p1 ageInPreg if education == 2) ///
	(line p1 ageInPreg if education == 3) ///
	(line p1 ageInPreg if education == 4) ///
	(line p1 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious attendance - Not at all", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(never, replace)
	
twoway (line p2 ageInPreg if education == 1) ///
	(line p2 ageInPreg if education == 2) ///
	(line p2 ageInPreg if education == 3) ///
	(line p2 ageInPreg if education == 4) ///
	(line p2 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious attendance - Min 1/year", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(year, replace)
	
twoway (line p3 ageInPreg if education == 1) ///
	(line p3 ageInPreg if education == 2) ///
	(line p3 ageInPreg if education == 3) ///
	(line p3 ageInPreg if education == 4) ///
	(line p3 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious attendance - Min 1/month", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(month, replace)
	
twoway (line p4 ageInPreg if education == 1) ///
	(line p4 ageInPreg if education == 2) ///
	(line p4 ageInPreg if education == 3) ///
	(line p4 ageInPreg if education == 4) ///
	(line p4 ageInPreg if education == 5), ///
	xscale(range(13 72)) xlabel(15(5)70, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Age in pregnancy") ytitle("Predicted probability") ///
	title("Religious attendance - Min 1/week", size(medium)) ///
	legend(order(1 "CSE/None" 2 "Vocational" 3 "O-level" 4 "A-level" 5 "Degree") ///
	cols(5) size(small)) ///
	name(week, replace)
	
grc1leg never year month week, ycommon

graph export ".\G0Partner_Results\eduAgeIntPredProbs_attend.pdf", replace

graph close _all




