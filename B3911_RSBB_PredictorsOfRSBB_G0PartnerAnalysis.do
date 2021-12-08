*** Predictors of RSBB (B3911) - G0 partner/father analysis script
*** Created 23/11/2021 by Dan Smith
*** Stata v16.0

*** This script reads in the cleaned G0 partner/father data, explores associations between exposures, and then conducts an 'exposome-wide association analysis' (ExWAS) on all these variables to examine how they are associated with various facets of RSBB.


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


**********************************************************************************
*** Descriptive statistics

** Put the RSBB variables at the start of the dataset
order aln pb150 pb153 pb153_grp pb155 FC3000 FC3040 FC3040_grp FC3080 FC3080_OccNever FC3080_OccYr FC3153 FC3153_cat FC3160 FC3170 FC3155 FC3155_cat

** Using the 'distinct' command (see above), for each variable inspect the number of unque values; if < 10 then display table, while if >= 10 then displays means/SDs/IQRs

* Outcomes
foreach var of varlist pb150-FC3155_cat {
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
foreach var of varlist a551-esteem_prorated {
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


** Amount of missing data in outcomes and exposures

* Outcomes
misstable sum pb150-FC3155_cat, all

* Exposures
misstable sum a551-esteem_prorated, all



**********************************************************************************
*** Correlations between exposures

** Explore correlations between exposures to see how inter-related these factors are
desc a551-esteem_prorated

* Most variables are either continuous, ordered categories, or binary variables, so will just use standard pearson correlations for these. Only unordered categories are home ownership (owned/mortaged vs renting vs counci/housing association vs other) and marital status (never married vs married vs widowed/divorced/separated). So will exclude these two variables from the correlation matrix and calculate their associations using these variables as outcomes in a multinomial regression and square-rooting the pseudo-R2 value (similar to how the 'rexposome' package in R for ExWAS analyses does)

* First, order variables into categories of 'demographic', 'socioeconomic/material insecurity' and 'cognitive/psychological' (FC9992 is age at @28 RSBB questionnaire, so will pop at end as not needed here)
order aln-FC3155_cat ///
pb910 c801_grp pa065_grp pa005_grp jan1993ur01ind_grp b032_grp ///
c666a partner_mum_edu c765_grp logavinceq jan1993imd2010q5_M jan1993Townsendq5_M a006_grp pb184_grp pb481 pb482 a551 neighbour_qual ///
pb546-pb551 pa782 esteem_prorated ///
FC9992

* Next, rename all these exposures so they are more intuitive and can be read easier on the correlation heatmaps below
rename pb910 ageInPreg
rename c801_grp nonWhiteEthnic
rename pa065_grp maritalStatus
rename pa005_grp mobility 
rename jan1993ur01ind_grp rural
rename b032_grp parity
rename c666a education
rename partner_mum_edu maternalEdu
rename c765_grp highSocClass
rename logavinceq income
rename jan1993imd2010q5_M IMD
rename jan1993Townsendq5_M townsendDep
rename a006_grp housing 
rename pb184_grp financeDiffs
rename pb481 chLifeEvents_wgt
rename pb482 chLifeEvents_total
rename a551 crowding
rename neighbour_qual neighPercept
rename pb546 IPSM_interpAware
rename pb547 IPSM_approval
rename pb548 IPSM_sepAnx
rename pb549 IPSM_timidity
rename pb550 IPSM_fragility
rename pb551 IPSM_total
rename pa782 LoC_external
rename esteem_prorated selfEsteem
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
putexcel A1=matrix(cor_demo), names

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
corr education maternalEdu highSocClass income IMD townsendDep financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept

matrix cor_socio = r(C)

* As lots of entries, the correlation coefficients are hard to read, so will drop these values and include the legend
heatplot cor_socio, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))
	
* Save heatmap
graph export ".\G0Partner_Results\corr_heatplot_socioOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Partner_Results\corrMatrix_socioOnly.xlsx", replace
putexcel A1=matrix(cor_socio), names


* And now for associations of socioeconomic variables with unordered categorical variable home ownership status (and save to a CSV file)
capture postclose housing_corrs_socioOnly
postfile housing_corrs_socioOnly str30 variable corr ///
	using ".\G0Partner_Results\housing_corrs_socioOnly.dta", replace
	
foreach var of varlist education maternalEdu highSocClass income IMD townsendDep financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept {
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


** Now repeat for cognitive/psychological variables (no unordered categorical variables, so can include all variables)
corr IPSM_interpAware-selfEsteem

matrix cor_cog = r(C)

heatplot cor_cog, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))

* Save heatmap
graph export ".\G0Partner_Results\corr_heatplot_cogOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Partner_Results\corrMatrix_cogOnly.xlsx", replace
putexcel A1=matrix(cor_cog), names

	
*** Finally, repeat this on all of the exposures together (excluding unordered cateogorical variables housing status and marital status)
corr ageInPreg nonWhiteEthnic mobility rural parity education maternalEdu highSocClass income IMD townsendDep financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem

matrix cor_all = r(C)

heatplot cor_all, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle(""))

* Save heatmap
graph export ".\G0Partner_Results\corr_heatplot_all.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Partner_Results\corrMatrix_all.xlsx", replace
putexcel A1=matrix(cor_all), names

	
* And now for associations of all other exposures variables with unordered categorical variables marital status and home ownership status

* Marital status
capture postclose marital_corrs_all
postfile marital_corrs_all str30 variable corr ///
	using ".\G0Partner_Results\marital_corrs_all.dta", replace
	
foreach var of varlist ageInPreg nonWhiteEthnic mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
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
	
foreach var of varlist ageInPreg nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
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



**************************************************************************************
*** Next, we want to run the actual ExWAS analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab pb150, m

** This will be quite complicated, as want to post results to file, but as exposures differ extracting the results will be variable specific. To adjust for multiple corrections will use conservative bonferroni adjustment when constructing confidence intervals and interpreting p-values - As 26 exposures, a Bonferroni p-value threshold of 0.05/26 = 0.0019.
display round(100 - 5/26, 0.01)
display 0.05/26

** We also want to store both estimates adjusting for age (other than for the age-only model), and also the interaction between age and the exposure, to see whether it's moderated by age. Again, this makes the set-up a bit more complicated.

** Create a postfile to post results to, then start the loop - Will create two postfiles; one for coefficients and CIs, and another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model) - NOTE: Have to store pvalues as 'double' format, else really tiny p-values get coded as '0' (as default if float format, which has minimum value of -3.40282346639e+38 [https://blog.stata.com/2012/04/02/the-penultimate-guide-to-precision/]).
capture postclose partner_belief
postfile partner_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_belief_results.dta", replace

capture postclose partner_belief_lr
postfile partner_belief_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_belief_results_lr.dta", replace

foreach var of varlist ageInPreg nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
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
		
		// And finally run the likelihood ratio tests
		mlogit pb150 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb150 `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
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
		
		// And finally run the likelihood ratio tests
		mlogit pb150 ageInPreg if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb150 ageInPreg `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit pb150 c.ageInPreg##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit pb150 ageInPreg if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb150 ageInPreg i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit pb150 c.ageInPreg##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_belief
postclose partner_belief_lr
	

****************************************************************************************
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

foreach var of varlist ageInPreg nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
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
		
		// And finally run the likelihood ratio tests
		mlogit pb153_grp if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb153_grp `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
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
		
		// And finally run the likelihood ratio tests
		mlogit pb153_grp ageInPreg if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb153_grp ageInPreg `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit pb153_grp c.ageInPreg##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit pb153_grp ageInPreg if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit pb153_grp ageInPreg i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit pb153_grp c.ageInPreg##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_relig
postclose partner_relig_lr



****************************************************************************************
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


** So instead will just run multinomial models with 'not at all' as the baseline/reference category.


*** Now run the loop to save all the results
capture postclose partner_attend
postfile partner_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_attend_results.dta", replace

capture postclose partner_attend_lr
postfile partner_attend_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_attend_results_lr.dta", replace

foreach var of varlist ageInPreg nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
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
		
		// And finally run the likelihood ratio tests
		mlogit pb155_rev if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit pb155_rev `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
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
		
		// And finally run the likelihood ratio tests
		mlogit pb155_rev ageInPreg if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit pb155_rev ageInPreg `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit pb155_rev c.ageInPreg##c.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit pb155_rev ageInPreg if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit pb155_rev ageInPreg i.`var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit pb155_rev c.ageInPreg##i.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_attend
postclose partner_attend_lr



****************************************************************************************
*** Now to the next RSBB outcome: Intrinsic religiosity

** As this was collected 28 years post-birth and not all partners who completed this will have an age at birth (as enrolled in later phases), will use age at questionnaire completion instead.

* Intrinsic religiosity outcome is very non-normal (spike at lowesst value '3', then uniform), so will also run multinomial model using categories of this variable (with '3' as baseline category).
tab1 FC3153 FC3153_cat
sum FC3153
hist FC3153, freq width(1)


*** Start with continuous measure and linear regression

*** Now run the loop to save all the results
capture postclose partner_intrinsic
postfile partner_intrinsic str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_intrinsic_results.dta", replace

capture postclose partner_intrinsic_lr
postfile partner_intrinsic_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_intrinsic_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		regress FC3153 `var',
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
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
		
		post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// And finally run the likelihood ratio tests
		regress FC3153 if `var' != .,
		est store base
		regress FC3153 `var',
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		regress FC3153 ageAt28 `var',
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// Now for interaction model
		regress FC3153 c.ageAt28##c.`var',
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
		// And finally run the likelihood ratio tests
		regress FC3153 ageAt28 if `var' != .,
		est store base
		regress FC3153 ageAt28 `var',
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress FC3153 c.ageAt28##c.`var',
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
		
			// No reference level of outcome, so set as NA
			local outcome_level = "NA"
			
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
			
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
					
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 4)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
						
			
			// Move to the next category of the exposure (category 4)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
							
			// Move to the next category of the exposure (category 5)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
			
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
			
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 5)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)
		
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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
							
			
			// Move to the next category of the exposure (category 6)
			regress FC3153 ageAt28 i.`var',
		
			local n = e(N)

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
			regress FC3153 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post partner_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		regress FC3153 ageAt28 if `var' != .,
		est store base
		regress FC3153 ageAt28 i.`var',
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress FC3153 c.ageAt28##i.`var',
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_intrinsic
postclose partner_intrinsic_lr


**** And now repeat for the categorical intrinsic religiosity variable as a sensitivity analysis to ensure results robust and not due to odd distribution of outcome
tab FC3153_cat

*** Now run the loop to save all the results
capture postclose partner_intrinsic_cat
postfile partner_intrinsic_cat str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_intrinsic_cat_results.dta", replace

capture postclose partner_intrinsic_cat_lr
postfile partner_intrinsic_cat_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_intrinsic_cat_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit FC3153_cat `var', baseoutcome(1) rrr
		
		local n = e(N)
		
		// Start with the first reference category (2/moderate IR)
		local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
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
		
		post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/high IR)
		local outcome_level = "High IR/8-11 (ref = lowest/3)"
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
		
		post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// Now onto the next reference category (4/highest IR)
		local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
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
		
		post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3153_cat if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3153_cat `var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit FC3153_cat ageAt28 `var', baseoutcome(1) rrr
		
		local n = e(N)
		
		// Start with the first reference category (2/moderate)
		local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Now for interaction model
		mlogit FC3153_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/high)
		mlogit FC3153_cat ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "High IR/8-11 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit FC3153_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (4/highest)
		mlogit FC3153_cat ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit FC3153_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3153_cat ageAt28 if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3153_cat ageAt28 `var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3153_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post partner_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit FC3153_cat ageAt28 if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3153_cat ageAt28 i.`var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_intrinsic_cat
postclose partner_intrinsic_cat_lr



**********************************************************************************
*** Next, want to explore extrinsic religiosity - As this is two separate questions (make friends and pray for protection), will have two run two separate loops for each outcome
tab1 FC3160 FC3170

** To reduce number of categories and combine low cell counts, will combine 'mildly' and 'strongly' agree and disagree together, to give 4 categories: Agree, not sure, disagree and not applicable. Will have 'agree' as the baseline, which is the most 'extrinsic' response.
recode FC3160 (1 2 = 1) (3 = 2) (4 5 = 3) (6 = 4), gen(FC3160_new)

label define extrinsic_lb 1 "Agree" 2 "Not sure" 3 "Disagree" 4 "Not applicable"
numlabel extrinsic_lb, add
label values FC3160_new extrinsic_lb
tab FC3160_new

recode FC3170 (1 2 = 1) (3 = 2) (4 5 = 3) (6 = 4), gen(FC3170_new)
label values FC3170_new extrinsic_lb
tab FC3170_new


*** Start with FC3160 - Attends church to help them make friends

*** Now run the loop to save all the results
capture postclose partner_extrinsic_friends
postfile partner_extrinsic_friends str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_extrinsic_friends_results.dta", replace

capture postclose partner_extrinsic_friends_lr
postfile partner_extrinsic_friends_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_extrinsic_friends_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit FC3160_new `var', baseoutcome(1) rrr
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
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
		
		post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		local outcome_level = "Disagree ER (ref = Agree)"
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
		
		post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// Now onto the next reference category (4/NA)
		local outcome_level = "Not applicable ER (ref = Agree)"
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
		
		post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3160_new if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3160_new `var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit FC3160_new ageAt28 `var', baseoutcome(1) rrr
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Now for interaction model
		mlogit FC3160_new c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		mlogit FC3160_new ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "Disagree ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit FC3160_new c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (4/NA)
		mlogit FC3160_new ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "Not applicable ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit FC3160_new c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3160_new ageAt28 if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3160_new ageAt28 `var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3160_new c.ageAt28##c.`var', baseoutcome(1) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post partner_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit FC3160_new ageAt28 if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3160_new ageAt28 i.`var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3160_new c.ageAt28##i.`var', baseoutcome(1) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_extrinsic_friends
postclose partner_extrinsic_friends_lr



*** And now repeat with FC3170 - Prays for relief and protection

*** Now run the loop to save all the results
capture postclose partner_extrinsic_prayer
postfile partner_extrinsic_prayer str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_extrinsic_prayer_results.dta", replace

capture postclose partner_extrinsic_prayer_lr
postfile partner_extrinsic_prayer_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_extrinsic_prayer_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit FC3170_new `var', baseoutcome(1) rrr
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
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
		
		post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		local outcome_level = "Disagree ER (ref = Agree)"
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
		
		post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// Now onto the next reference category (4/NA)
		local outcome_level = "Not applicable ER (ref = Agree)"
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
		
		post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3170_new if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3170_new `var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit FC3170_new ageAt28 `var', baseoutcome(1) rrr
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Now for interaction model
		mlogit FC3170_new c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		mlogit FC3170_new ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "Disagree ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit FC3170_new c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (4/NA)
		mlogit FC3170_new ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "Not applicable ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit FC3170_new c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3170_new ageAt28 if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3170_new ageAt28 `var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3170_new c.ageAt28##c.`var', baseoutcome(1) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post partner_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit FC3170_new ageAt28 if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3170_new ageAt28 i.`var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3170_new c.ageAt28##i.`var', baseoutcome(1) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_extrinsic_prayer
postclose partner_extrinsic_prayer_lr



**************************************************************************************
*** Now finally to the last RSBB outcome: Total DUREL religiosity score

* Total DUREL religiosity outcome is very non-normal (spike at lowesst value '5', then uniform), so will also run multinomial model using categories of this variable (with '5' as baseline category).
tab1 FC3155 FC3155_cat
sum FC3155
hist FC3155, freq width(1)


*** Start with continuous measure and linear regression

*** Now run the loop to save all the results
capture postclose partner_DUREL
postfile partner_DUREL str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_DUREL_results.dta", replace

capture postclose partner_DUREL_lr
postfile partner_DUREL_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_DUREL_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		regress FC3155 `var',
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
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
		
		post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// And finally run the likelihood ratio tests
		regress FC3155 if `var' != .,
		est store base
		regress FC3155 `var',
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		regress FC3155 ageAt28 `var',
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// Now for interaction model
		regress FC3155 c.ageAt28##c.`var',
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
		// And finally run the likelihood ratio tests
		regress FC3155 ageAt28 if `var' != .,
		est store base
		regress FC3155 ageAt28 `var',
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress FC3155 c.ageAt28##c.`var',
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
		
			// No reference level of outcome, so set as NA
			local outcome_level = "NA"
			
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
			
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
					
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 4)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
						
			
			// Move to the next category of the exposure (category 4)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
							
			// Move to the next category of the exposure (category 5)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
			
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
			
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 5)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)
		
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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
							
			
			// Move to the next category of the exposure (category 6)
			regress FC3155 ageAt28 i.`var',
		
			local n = e(N)

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
			regress FC3155 c.ageAt28##i.`var',
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post partner_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		regress FC3155 ageAt28 if `var' != .,
		est store base
		regress FC3155 ageAt28 i.`var',
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress FC3155 c.ageAt28##i.`var',
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_DUREL
postclose partner_DUREL_lr


**** And now repeat for the categorical total DUREL religiosity variable as a sensitivity analysis to ensure results robust and not due to odd distribution of outcome
tab FC3155_cat

*** Now run the loop to save all the results
capture postclose partner_DUREL_cat
postfile partner_DUREL_cat str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_DUREL_cat_results.dta", replace

capture postclose partner_DUREL_cat_lr
postfile partner_DUREL_cat_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_DUREL_cat_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit FC3155_cat `var', baseoutcome(1) rrr
		
		local n = e(N)
		
		// Start with the first reference category (2/6-10 DUREL)
		local outcome_level = "DUREL 6-10 (ref = lowest/5)"
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
		
		post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/11-15 DUREL)
		local outcome_level = "11-15 DUREL (ref = lowest/5)"
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
		
		post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// Now onto the next reference category (4/16-20 DUREL)
		local outcome_level = "16-20 DUREL (ref = lowest/5)"
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
		
		post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (5/21-26 DUREL)
		local outcome_level = "21-26 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,9]
		local lci = res[5,9]
		local uci = res[6,9]
		local p = res[4,9]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local age_main = .
		local exp_main = .
		
		post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3155_cat if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3155_cat `var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit FC3155_cat ageAt28 `var', baseoutcome(1) rrr
		
		local n = e(N)
		
		// Start with the first reference category (2/6-10 DUREL)
		local outcome_level = "DUREL 6-10 (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Now for interaction model
		mlogit FC3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/11-15 DUREL)
		mlogit FC3155_cat ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "11-15 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit FC3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (4/16-20 DUREL)
		mlogit FC3155_cat ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "16-20 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit FC3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (5/21-26 DUREL)
		mlogit FC3155_cat ageAt28 `var', baseoutcome(1) rrr
		
		local outcome_level = "21-26 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,14]
		local lci = res[5,14]
		local uci = res[6,14]
		local p = res[4,14]
				
		// Now for interaction model
		mlogit FC3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		
		matrix res = r(table)
		local coef_int = res[1,19]
		local lci_int = res[5,19]
		local uci_int = res[6,19]
		local p_int = res[4,19]
		local age_main = res[1,17]
		local exp_main = res[1,18]
		
		post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3155_cat ageAt28 if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3155_cat ageAt28 `var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,33]
			local exp_main = res[1,35]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,33]
			local exp_main = res[1,36]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
						
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,41]
			local exp_main = res[1,43]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,48]
			local lci_int = res[5,48]
			local uci_int = res[6,48]
			local p_int = res[4,48]
			local age_main = res[1,41]
			local exp_main = res[1,44]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,49]
			local lci_int = res[5,49]
			local uci_int = res[6,49]
			local p_int = res[4,49]
			local age_main = res[1,41]
			local exp_main = res[1,45]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,56]
			local lci_int = res[5,56]
			local uci_int = res[6,56]
			local p_int = res[4,56]
			local age_main = res[1,49]
			local exp_main = res[1,51]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,57]
			local lci_int = res[5,57]
			local uci_int = res[6,57]
			local p_int = res[4,57]
			local age_main = res[1,49]
			local exp_main = res[1,52]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,58]
			local lci_int = res[5,58]
			local uci_int = res[6,58]
			local p_int = res[4,58]
			local age_main = res[1,49]
			local exp_main = res[1,53]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,59]
			local lci_int = res[5,59]
			local uci_int = res[6,59]
			local p_int = res[4,59]
			local age_main = res[1,49]
			local exp_main = res[1,54]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,35]
			local lci = res[5,35]
			local uci = res[6,35]
			local p = res[4,35]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,65]
			local lci_int = res[5,65]
			local uci_int = res[6,65]
			local p_int = res[4,65]
			local age_main = res[1,57]
			local exp_main = res[1,59]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,36]
			local lci = res[5,36]
			local uci = res[6,36]
			local p = res[4,36]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,66]
			local lci_int = res[5,66]
			local uci_int = res[6,66]
			local p_int = res[4,66]
			local age_main = res[1,57]
			local exp_main = res[1,60]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,37]
			local lci = res[5,37]
			local uci = res[6,37]
			local p = res[4,37]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,67]
			local lci_int = res[5,67]
			local uci_int = res[6,67]
			local p_int = res[4,67]
			local age_main = res[1,57]
			local exp_main = res[1,61]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,38]
			local lci = res[5,38]
			local uci = res[6,38]
			local p = res[4,38]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,68]
			local lci_int = res[5,68]
			local uci_int = res[6,68]
			local p_int = res[4,68]
			local age_main = res[1,57]
			local exp_main = res[1,62]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,39]
			local lci = res[5,39]
			local uci = res[6,39]
			local p = res[4,39]
				
			// Now for interaction model
			mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		
			matrix res = r(table)
			local coef_int = res[1,69]
			local lci_int = res[5,69]
			local uci_int = res[6,69]
			local p_int = res[4,69]
			local age_main = res[1,57]
			local exp_main = res[1,63]
		
			post partner_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit FC3155_cat ageAt28 if `var' != ., baseoutcome(1) rrr
		est store base
		mlogit FC3155_cat ageAt28 i.`var', baseoutcome(1) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_DUREL_cat
postclose partner_DUREL_cat_lr


********************************************************************************
**** Also want to get the likelihood ratio p-values for belief in God, religious affiliation and church attendance for the 2019 data as well, as more directly comparable with the intrinsic/extrinsic/total religiosity data (will also save coefficients as well, in case decide to use this too)
tab1 pb150 FC3000, m

*** Belief in God (2019)
capture postclose partner_belief_2019
postfile partner_belief_2019 str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_belief_2019_results.dta", replace

capture postclose partner_belief_2019_lr
postfile partner_belief_2019_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_belief_2019_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit FC3000 `var', baseoutcome(3) rrr
		
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
		
		post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3000 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit FC3000 `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_belief_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit FC3000 ageAt28 `var', baseoutcome(3) rrr
		
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
		mlogit FC3000 c.ageAt28##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/not sure)
		mlogit FC3000 ageAt28 `var', baseoutcome(3) rrr
		
		local outcome_level = "Not sure (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
				
		// Now for interaction model
		mlogit FC3000 c.ageAt28##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3000 ageAt28 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit FC3000 ageAt28 `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3000 c.ageAt28##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_belief_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,17]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_belief_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit FC3000 ageAt28 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit FC3000 ageAt28 i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3000 c.ageAt28##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_belief_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_belief_2019
postclose partner_belief_2019_lr


*** Religious affiliation (2019)
*** As this is an unordered categorical variable, will again use multinomial logistic model (with 'none' as reference)
tab1 pb153_grp FC3040_grp

* Do need to re-order the categories so can simply copy and paste the script from above without having to faff around with editing the cells to take the statistics from
recode FC3040_grp (1 = 3) (2 = 1) (3 = 2)
tab FC3040_grp


*** Now run the loop to save all the results
capture postclose partner_relig_2019
postfile partner_relig_2019 str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_relig_2019_results.dta", replace

capture postclose partner_relig_2019_lr
postfile partner_relig_2019_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_relig_2019_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit FC3040_grp `var', baseoutcome(3) rrr
		
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
		
		post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3040_grp if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit FC3040_grp `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_relig_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit FC3040_grp ageAt28 `var', baseoutcome(3) rrr
		
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
		mlogit FC3040_grp c.ageAt28##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/Other)
		mlogit FC3040_grp ageAt28 `var', baseoutcome(3) rrr
		
		local outcome_level = "Other (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
				
		// Now for interaction model
		mlogit FC3040_grp c.ageAt28##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3040_grp ageAt28 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit FC3040_grp ageAt28 `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3040_grp c.ageAt28##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_relig_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,17]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
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
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_relig_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit FC3040_grp ageAt28 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit FC3040_grp ageAt28 i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3040_grp c.ageAt28##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_relig_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_relig_2019
postclose partner_relig_2019_lr


**** Attendance at church/place of worship (2019)

*** As this is an ordered categorical variable, will use ordinal regression model. For consistency with previous models, will recode so that higher values indicate greater RSBB. For the 2019 data, this question was asked different and included an 'occasionally' option; here, I will combine this this 'at least once a year' so is the same coding as the pregnancy data
tab1 pb155 FC3080_OccYr

recode FC3080_OccYr (1 = 3) (2 = 2) (3 = 1) (4 = 0), gen(FC3080_OccYr_rev)
label values FC3080_OccYr_rev attend_rev_lb
tab FC3080_OccYr_rev


** Quick test of whether proportional odds assumption been violated in most basic model (with just age at birth). Not violated here, but for consistency with other analyses will just use multinomial models.
ologit FC3080_OccYr_rev ageAt28, or
brant, detail

ologit FC3080_OccYr_rev ageAt28 i.IMD, or
brant, detail


** Will just run multinomial models with 'not at all' as the baseline/reference category.


*** Now run the loop to save all the results
capture postclose partner_attend_2019
postfile partner_attend_2019 str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\G0Partner_Results\partner_attend_2019_results.dta", replace

capture postclose partner_attend_2019_lr
postfile partner_attend_2019_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G0Partner_Results\partner_attend_2019_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit FC3080_OccYr_rev `var', baseoutcome(0) rrr
		
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
		
		post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3080_OccYr_rev if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit FC3080_OccYr_rev `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post partner_attend_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit FC3080_OccYr_rev ageAt28 `var', baseoutcome(0) rrr
		
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
		mlogit FC3080_OccYr_rev c.ageAt28##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (2/once month)
		mlogit FC3080_OccYr_rev ageAt28 `var', baseoutcome(0) rrr
		
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit FC3080_OccYr_rev c.ageAt28##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/once week)
		mlogit FC3080_OccYr_rev ageAt28 `var', baseoutcome(0) rrr
		
		local outcome_level = "Min once week (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit FC3080_OccYr_rev c.ageAt28##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit FC3080_OccYr_rev ageAt28 if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit FC3080_OccYr_rev ageAt28 `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3080_OccYr_rev c.ageAt28##c.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_attend_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Vocational (ref = None/CSE/GCSE)"
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "A-level (ref = None/CSE/GCSE)"
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "maternalEdu" {
				local exp_level = "Degree (ref = None/CSE/GCSE)"
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (2/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,19]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
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
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/once week)
			mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post partner_attend_2019 ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit FC3080_OccYr_rev ageAt28 if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit FC3080_OccYr_rev ageAt28 i.`var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit FC3080_OccYr_rev c.ageAt28##i.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post partner_attend_2019_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose partner_attend_2019
postclose partner_attend_2019_lr


* Save this file if we want to use it later
save ".\G0Partner_Results\G0Partner_PredictorsOfRSBB_B3911_postAnalysis.dta", replace

* And close the log file
log close


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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 26 exposures, will do 0.05/26)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Main effect") ///
	name(belief_main, replace)
	
graph export ".\G0Partner_Results\belief_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief in God/divine power - Age interaction") ///
	name(belief_int, replace)
	
graph export ".\G0Partner_Results\belief_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 26 exposures, will do 0.05/26)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Main effect") ///
	name(relig_main, replace)
	
graph export ".\G0Partner_Results\relig_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Age interaction") ///
	name(relig_int, replace)
	
graph export ".\G0Partner_Results\relig_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 26 exposures, will do 0.05/26)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Main effect") ///
	name(attend_main, replace)
	
graph export ".\G0Partner_Results\attend_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Church attendance - Age interaction") ///
	name(attend_int, replace)
	
graph export ".\G0Partner_Results\attend_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
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


**** Now read in the next outcome - intrinsic religiosity (will use the multinomial results, as the distribution for the linear model is very non-normal, and using multinomial results makes these p-value plots consistent with the coefficient plots below, which also use the multinomial results)
use ".\G0Partner_Results\partner_intrinsic_cat_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 26 exposures, will do 0.05/26)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Intrinsic religiosity - Main effect") ///
	name(intrin_main, replace)
	
graph export ".\G0Partner_Results\intrin_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Intrinsic religiosity - Age interaction") ///
	name(intrin_int, replace)
	
graph export ".\G0Partner_Results\intrin_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Intrinsic religiosity") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(intrin_both, replace)

graph export ".\G0Partner_Results\intrin_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'intrinsic relig' as a variable, then save this file
gen outcome = "Intrinsic relig."
recast str30 outcome
order outcome

save ".\G0Partner_Results\intrin_pvalues.dta", replace


**** Now read in the next outcome - extrinsic religiosity (friends)
use ".\G0Partner_Results\partner_extrinsic_friends_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 26 exposures, will do 0.05/26)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Extrinsic religiosity (friends) - Main effect") ///
	name(extrinFriends_main, replace)
	
graph export ".\G0Partner_Results\extrinFriends_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Extrinsic religiosity (friends) - Age interaction") ///
	name(extrinFriends_int, replace)
	
graph export ".\G0Partner_Results\extrinFriends_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Extrinsic religiosity (friends)") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(extrinFriends_both, replace)

graph export ".\G0Partner_Results\extrinFriends_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'extrinsic relig (friends)' as a variable, then save this file
gen outcome = "Extrinsic relig. (friends)"
recast str30 outcome
order outcome

save ".\G0Partner_Results\extrinFriends_pvalues.dta", replace


**** Now read in the next outcome - extrinsic religiosity (prayer)
use ".\G0Partner_Results\partner_extrinsic_prayer_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 26 exposures, will do 0.05/26)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Extrinsic religiosity (prayer) - Main effect") ///
	name(extrinPray_main, replace)
	
graph export ".\G0Partner_Results\extrinPray_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Extrinsic religiosity (prayer) - Age interaction") ///
	name(extrinPray_int, replace)
	
graph export ".\G0Partner_Results\extrinPray_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Extrinsic religiosity (prayer)") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(extrinPray_both, replace)

graph export ".\G0Partner_Results\extrinPray_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'extrinsic relig (prayer)' as a variable, then save this file
gen outcome = "Extrinsic relig. (prayer)"
recast str30 outcome
order outcome

save ".\G0Partner_Results\extrinPrayer_pvalues.dta", replace


**** Now read in the final outcome - Total DUREL religiosity score (will use the multinomial results, as the distribution for the linear model is very non-normal, and using multinomial results makes these p-value plots consistent with the coefficient plots below, which also use the multinomial results)
use ".\G0Partner_Results\partner_DUREL_cat_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 26 exposures, will do 0.05/26)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Total DUREL religiosity score - Main effect") ///
	name(durel_main, replace)
	
graph export ".\G0Partner_Results\durel_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Total DUREL religiosity score - Age interaction") ///
	name(durel_int, replace)
	
graph export ".\G0Partner_Results\durel_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Total DUREL religiosity score") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
	name(durel_both, replace)

graph export ".\G0Partner_Results\durel_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'DUREL' as a variable, then save this file
gen outcome = "DUREL"
recast str30 outcome
order outcome

save ".\G0Partner_Results\durel_pvalues.dta", replace


*** Combine all these datasets together
use ".\G0Partner_Results\belief_pvalues.dta", clear
append using ".\G0Partner_Results\relig_pvalues.dta"
append using ".\G0Partner_Results\attend_pvalues.dta"
append using ".\G0Partner_Results\intrin_pvalues.dta"
append using ".\G0Partner_Results\extrinFriends_pvalues.dta"
append using ".\G0Partner_Results\extrinPrayer_pvalues.dta"
append using ".\G0Partner_Results\durel_pvalues.dta"


** Now look at combined results - As data come from different sources with different sample sizes the p-values wont be directly comparable, so will compare belief/religious affiliation/church attendance (from pregnancy) in one plot, and intrinsic/extrinsic/total religiosity (from the recent RSBB data collection @28) in another plot.

* Pregnancy vars main effects
local bon_thresh = -log10(0.05/26)
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
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief, religion and attendance - Main effects") ///
	legend(order(1 "Belief in God" 2 "Religious affiliation" ///
		3 "Church attendance") rows(1) size(small)) ///
	name(preg_main, replace)

graph export ".\G0Partner_Results\beliefReligAttend_mainEffects_pvalues.pdf", replace

* Pregnancy vars interaction effects
local bon_thresh = -log10(0.05/25)
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
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Belief, religion and attendance - Age interaction") ///
	legend(order(1 "Belief in God" 2 "Religious affiliation" ///
		3 "Church attendance") rows(1) size(small)) ///
	name(preg_int, replace)

graph export ".\G0Partner_Results\beliefReligAttend_ageInt_pvalues.pdf", replace

* Age 28 vars main effects
local bon_thresh = -log10(0.05/26)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main if outcome == "Intrinsic relig.", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_main if outcome == "Extrinsic relig. (friends)", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num logp_main if outcome == "Extrinsic relig. (prayer)", ///
		col(blue) msize(small) msym(D)) ///
	(scatter exp_num logp_main if outcome == "DUREL", ///
		col(green) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Intrinsic, extrinsic and DUREL - Main effects") ///
	legend(order(1 "Intrinsic religiosity" 2 "Extrinsic religiosity (friends)" ///
		3 "Extrinsic religiosity (prayer)" 4 "DUREL religiosity") ///
		rows(2) size(small)) ///
	name(age28_main, replace)

graph export ".\G0Partner_Results\intExtDUREL_mainEffects_pvalues.pdf", replace

* Age 28 vars interaction effects
local bon_thresh = -log10(0.05/25)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if outcome == "Intrinsic relig." & exp_num != 1, ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Extrinsic relig. (friends)" & exp_num != 1, ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Extrinsic relig. (prayer)" & exp_num != 1, ///
		col(blue) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "DUREL" & exp_num != 1, ///
		col(green) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)26, valuelabel labsize(vsmall) angle(0)) ///
	title("Intrinsic, extrinsic and DUREL - Age interaction") ///
	legend(order(1 "Intrinsic religiosity" 2 "Extrinsic religiosity (friends)" ///
		3 "Extrinsic religiosity (prayer)" 4 "DUREL religiosity") ///
		rows(2) size(small)) ///
	name(age28_int, replace)

graph export ".\G0Partner_Results\intExtDUREL_ageInt_pvalues.pdf", replace


** Combine all these graphs together
graph combine preg_main preg_int age28_main age28_int, imargin(0 0 0 0) iscale(0.5)

graph export ".\G0Partner_Results\allData_pvalues.pdf", replace

graph close _all


** Save these p-values as CSV files
list outcome exposure lr_p_main lr_p_int in 1/10

keep outcome exp_num lr_p_main lr_p_int
order outcome exp_num lr_p_main lr_p_int

* Convert to wide format (need to edit outcomes strings so have no spaces)
replace outcome = "Attend" if outcome == "Church attendance"
replace outcome = "Extrin_fr" if outcome == "Extrinsic relig. (friends)"
replace outcome = "Extrin_pr" if outcome == "Extrinsic relig. (prayer)"
replace outcome = "Intrins" if outcome == "Intrinsic relig."
replace outcome = "Religion" if outcome == "Religious affil."
tab outcome

reshape wide lr_p_main lr_p_int, i(exp_num) j (outcome) string

order exp_num lr_p_mainBelief lr_p_intBelief lr_p_mainReligion lr_p_intReligion lr_p_mainAttend lr_p_intAttend lr_p_mainIntrins lr_p_intIntrins lr_p_mainExtrin_fr lr_p_intExtrin_fr lr_p_mainExtrin_pr lr_p_intExtrin_pr lr_p_mainDUREL lr_p_intDUREL

format %9.4f lr_p_mainBelief-lr_p_intDUREL

outsheet exp_num-lr_p_intDUREL using ".\G0Partner_Results\pvalue_results.csv", comma replace


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

* Intrinsic religiosity - main effect
count if lr_p_mainIntrins < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainIntrins < 0.05
display (r(N) / _N) * 100

* Intrinsic religiosity - interaction
count if lr_p_intIntrins < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_intIntrins < 0.05
display (r(N) / (_N - 1)) * 100

* Extrinsic religiosity (friends) - main effect
count if lr_p_mainExtrin_fr < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainExtrin_fr < 0.05
display (r(N) / _N) * 100

* Extrinsic religiosity (friends) - interaction
count if lr_p_intExtrin_fr < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_intExtrin_fr < 0.05
display (r(N) / (_N - 1)) * 100

* Extrinsic religiosity (prayer) - main effect
count if lr_p_mainExtrin_pr < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainExtrin_pr < 0.05
display (r(N) / _N) * 100

* Extrinsic religiosity (prayer) - interaction
count if lr_p_intExtrin_pr < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_intExtrin_pr < 0.05
display (r(N) / (_N - 1)) * 100

* Total religiosity - main effect
count if lr_p_mainDUREL < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainDUREL < 0.05
display (r(N) / _N) * 100

* Total religiosity - interaction
count if lr_p_intDUREL < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_intDUREL < 0.05
display (r(N) / (_N - 1)) * 100


*** Now get these for the 2019 data

** Belief in God
use ".\G0Partner_Results\partner_belief_2019_results_lr.dta", clear

* Main effect
count if lr_p_main < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_main < 0.05
display (r(N) / _N) * 100

* Interaction
count if lr_p_int < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_int < 0.05
display (r(N) / (_N - 1)) * 100

** Religious affiliation
use ".\G0Partner_Results\partner_relig_2019_results_lr.dta", clear

* Main effect
count if lr_p_main < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_main < 0.05
display (r(N) / _N) * 100

* Interaction
count if lr_p_int < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_int < 0.05
display (r(N) / (_N - 1)) * 100

** Church attendance
use ".\G0Partner_Results\partner_attend_2019_results_lr.dta", clear

* Main effect
count if lr_p_main < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_main < 0.05
display (r(N) / _N) * 100

* Interaction
count if lr_p_int < 0.05/(_N - 1)
display (r(N) / (_N - 1)) * 100

count if lr_p_int < 0.05
display (r(N) / (_N - 1)) * 100


*** Next, want to plot some of the actual results

** Will read the datasets in then combine together into one single dataset
use ".\G0Partner_Results\partner_belief_results.dta", clear
gen outcome = "Belief"

append using ".\G0Partner_Results\partner_relig_results.dta"
replace outcome = "Relig" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_attend_results.dta"
replace outcome = "Attend" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_intrinsic_results.dta"
replace outcome = "Intrinsic" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_intrinsic_cat_results.dta"
replace outcome = "Intrinsic (cat)" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_extrinsic_friends_results.dta"
replace outcome = "Extrinsic - friends" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_extrinsic_prayer_results.dta"
replace outcome = "Extrinsic - prayer" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_DUREL_results.dta"
replace outcome = "DUREL" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_DUREL_cat_results.dta"
replace outcome = "DUREL (cat)" if outcome == ""
tab outcome, m


** Save these results as CSV files to add to the SI

* Add in the 2019 data (then drop after)
append using ".\G0Partner_Results\partner_belief_2019_results.dta"
replace outcome = "Belief_2019" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_relig_2019_results.dta"
replace outcome = "Relig_2019" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_attend_2019_results.dta"
replace outcome = "Attend_2019" if outcome == ""
tab outcome, m

* Save each result in turn
format coef lci uci coef_int lci_int uci_int %9.3f
format p p_int %9.4f

outsheet exposure-p_int using ".\G0Partner_Results\belief_coefs.csv" if outcome == "Belief", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\relig_coefs.csv" if outcome == "Relig", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\attend_coefs.csv" if outcome == "Attend", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\IntrinsicCat_coefs.csv" if outcome == "Intrinsic (cat)", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\ExtrinFriends_coefs.csv" if outcome == "Extrinsic - friends", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\ExtrinPray_coefs.csv" if outcome == "Extrinsic - prayer", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\DUREL_coefs.csv" if outcome == "DUREL (cat)", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\belief_2019_coefs.csv" if outcome == "Belief_2019", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\relig_2019_coefs.csv" if outcome == "Relig_2019", comma replace

outsheet exposure-p_int using ".\G0Partner_Results\attend_2019_coefs.csv" if outcome == "Attend_2019", comma replace

* And now drop the 2019 data (as don't want to include these in the plots below)
drop if outcome == "Belief_2019" | outcome == "Relig_2019" | outcome == "Attend_2019"
tab outcome, m


** First, make a plot for the age results - As some outcomes are on different scales, will just use the results from multinomial regression for all outcomes (inc. intrinsic and total/DUREL religiosity, even though also ran linear regressions on these as well) - Having all variables on the same plot on the same scale makes things easier to visualise.
capture drop level_num
gen level_num = 0
replace level_num = 1 if outcome_level == "Not sure (ref = No)"
replace level_num = 3 if outcome_level == "Christian (ref = None)"
replace level_num = 4 if outcome_level == "Other (ref = None)"
replace level_num = 6 if outcome_level == "Min once year (ref = Not at al"
replace level_num = 7 if outcome_level == "Min once month (ref = Not at a"
replace level_num = 8 if outcome_level == "Min once week (ref = Not at al"
replace level_num = 10 if outcome_level == "Moderate IR/4-7 (ref = lowest/"
replace level_num = 11 if outcome_level == "High IR/8-11 (ref = lowest/3)"
replace level_num = 12 if outcome_level == "Highest IR/12-15 (ref = lowest"
replace level_num = 14 if outcome_level == "Not sure ER (ref = Agree)" & outcome == "Extrinsic - friends"
replace level_num = 15 if outcome_level == "Disagree ER (ref = Agree)" & outcome == "Extrinsic - friends"
replace level_num = 16 if outcome_level == "Not applicable ER (ref = Agree" & outcome == "Extrinsic - friends"
replace level_num = 18 if outcome_level == "Not sure ER (ref = Agree)" & outcome == "Extrinsic - prayer"
replace level_num = 19 if outcome_level == "Disagree ER (ref = Agree)" & outcome == "Extrinsic - prayer"
replace level_num = 20 if outcome_level == "Not applicable ER (ref = Agree" & outcome == "Extrinsic - prayer"
replace level_num = 22 if outcome_level == "DUREL 6-10 (ref = lowest/5)"
replace level_num = 23 if outcome_level == "11-15 DUREL (ref = lowest/5)"
replace level_num = 24 if outcome_level == "16-20 DUREL (ref = lowest/5)"
replace level_num = 25 if outcome_level == "21-26 DUREL (ref = lowest/5)"

label define level_lb 0 "Belief in God - Yes (ref = No)" 1 "Belief in God - Not sure (ref = No)" 3 "Religious affiliation - Christian (ref = None)" 4 "Religious affiliation - Other (ref = None)" 6 "Church attendance - Min once a year (ref = Not at all)" 7 "Church attendance - Min once a month (ref = Not at all)" 8 "Church attendance - Min once a week (ref = Not at all)" 10 "Intrinsic religiosity - Moderate (4-7 score; ref = lowest/3)" 11 "Intrinsic religiosity - High (8-11 score; ref = lowest/3)" 12 "Intrinsic religiosity - Highest (12-15 score; ref = lowest/3)" 14 "Extrinsic religiosity (friends) - Not sure (ref = Agree)" 15 "Extrinsic religiosity (friends) - Disagree (ref = Agree)" 16 "Extrinsic religiosity (friends) - Not applicable (ref = Agree)" 18 "Extrinsic religiosity (prayer) - Not sure (ref = Agree)" 19 "Extrinsic religiosity (prayer) - Disagree (ref = Agree)" 20 "Extrinsic religiosity (prayer) - Not applicable (ref = Agree)" 22 "DUREL total religiosity - 6-10 score (ref = lowest/5)" 23 "DUREL total religiosity - 11-15 score (ref = lowest/5)" 24 "DUREL total religiosity - 16-20 score (ref = lowest/5)" 25 "DUREL total religiosity - 21-26 score (ref = lowest/5)", replace
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
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "ageAt28", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "ageAt28", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "ageAt28", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "ageAt28", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "ageAt28", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "ageAt28", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "ageAt28", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "ageAt28", horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Age and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(age_cat, replace)
		
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
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "nonWhiteEthnic", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "nonWhiteEthnic", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "nonWhiteEthnic", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "nonWhiteEthnic", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "nonWhiteEthnic", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "nonWhiteEthnic", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "nonWhiteEthnic", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "nonWhiteEthnic", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = White)") ///
		title("Other than White ethnicity and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.1 0.2 0.5 1 2 5 10 20 30, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(ethnic_cat, replace)
		
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
			"Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Married (ref = Never married)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Wid/Div/Sep (ref = Never married)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Wid/Div/Sep (ref = Never married)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal col(red)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = Never married)") ///
		title("Marital status and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.01 0.03 0.1 0.3 1 3 10 30 100, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Married" 15 "Widowed/Divorced/Separated")) ///
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
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = CSE/None)") ///
		title("Education and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.05 0.1 0.2 0.3 0.5 1 2 3 5 10, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Vocational" 15 "O-levels" 29 "A-levels" ///
			43 "Degree") rows(1)) ///
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
			"highSocClass", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "highSocClass", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "highSocClass", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "highSocClass", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "highSocClass", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = lower [III manual/IV/V])") ///
		title("Occupational Social Class and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.3 0.5 0.7 1 1.5 2 3 5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(socClass_cat, replace)
		
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
			"income", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "income", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "income", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "income", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "income", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "income", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "income", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "income", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "income", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in log income)") ///
		title("Household income and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.3 0.5 0.7 1 1.5 2 3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(income_cat, replace)
		
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
			"2 (ref = 1/Least dep.)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"2 (ref = 1/Least dep.)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"2 (ref = 1/Least dep.)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"2 (ref = 1/Least dep.)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"2 (ref = 1/Least dep.)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"2 (ref = 1/Least dep.)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"2 (ref = 1/Least dep.)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"2 (ref = 1/Least dep.)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "2 (ref = 1/Least dep.)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "2 (ref = 1/Least dep.)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "2 (ref = 1/Least dep.)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "2 (ref = 1/Least dep.)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"2 (ref = 1/Least dep.)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"2 (ref = 1/Least dep.)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"3 (ref = 1/Least dep.)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"3 (ref = 1/Least dep.)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"3 (ref = 1/Least dep.)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"3 (ref = 1/Least dep.)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"3 (ref = 1/Least dep.)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"3 (ref = 1/Least dep.)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"3 (ref = 1/Least dep.)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"3 (ref = 1/Least dep.)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "3 (ref = 1/Least dep.)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "3 (ref = 1/Least dep.)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "3 (ref = 1/Least dep.)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "3 (ref = 1/Least dep.)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"3 (ref = 1/Least dep.)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"3 (ref = 1/Least dep.)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"4 (ref = 1/Least dep.)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"4 (ref = 1/Least dep.)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"4 (ref = 1/Least dep.)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"4 (ref = 1/Least dep.)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"4 (ref = 1/Least dep.)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"4 (ref = 1/Least dep.)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"4 (ref = 1/Least dep.)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"4 (ref = 1/Least dep.)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "4 (ref = 1/Least dep.)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "4 (ref = 1/Least dep.)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "4 (ref = 1/Least dep.)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "4 (ref = 1/Least dep.)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"4 (ref = 1/Least dep.)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"4 (ref = 1/Least dep.)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "5/Most dep. (ref = 1/Least dep.)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "5/Most dep. (ref = 1/Least dep.)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "5/Most dep. (ref = 1/Least dep.)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "5/Most dep. (ref = 1/Least dep.)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"5/Most dep. (ref = 1/Least dep.)", horizontal col(green) msize(vtiny) lwidth(thin)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = 1/Least Deprived)") ///
		title("IMD and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.06 0.1 0.2 0.5 1 2 3 5 7, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "2" 15 "3" 29 "4" 43 "5/Most dep.") rows(1)) ///
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
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = Owned/Mortgaged)") ///
		title("Housing status and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.03 0.1 0.2 0.5 1 2 3 5 10, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Rented" 15 "Council/HA" 29 "Other") ///
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
			"1 move (ref = 0 moves)", col(black) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"1 move (ref = 0 moves)", horizontal col(black) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"1 move (ref = 0 moves)", col(black) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"1 move (ref = 0 moves)", horizontal col(black) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"1 move (ref = 0 moves)", col(black) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"1 move (ref = 0 moves)", horizontal col(black) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"1 move (ref = 0 moves)", col(black) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"1 move (ref = 0 moves)", horizontal col(black) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "1 move (ref = 0 moves)", col(black) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "1 move (ref = 0 moves)", horizontal col(black) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "1 move (ref = 0 moves)", col(black) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "1 move (ref = 0 moves)", horizontal col(black) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"1 move (ref = 0 moves)", col(black) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"1 move (ref = 0 moves)", horizontal col(black) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"2 moves (ref = 0 moves)", col(red) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"2 moves (ref = 0 moves)", horizontal col(red) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"2 moves (ref = 0 moves)", col(red) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"2 moves (ref = 0 moves)", horizontal col(red) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"2 moves (ref = 0 moves)", col(red) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"2 moves (ref = 0 moves)", horizontal col(red) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"2 moves (ref = 0 moves)", col(red) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"2 moves (ref = 0 moves)", horizontal col(red) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "2 moves (ref = 0 moves)", col(red) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "2 moves (ref = 0 moves)", horizontal col(red) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "2 moves (ref = 0 moves)", col(red) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "2 moves (ref = 0 moves)", horizontal col(red) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"2 moves (ref = 0 moves)", col(red) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"2 moves (ref = 0 moves)", horizontal col(red) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"3 moves (ref = 0 moves)", col(blue) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"3 moves (ref = 0 moves)", horizontal col(blue) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"3 moves (ref = 0 moves)", col(blue) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"3 moves (ref = 0 moves)", horizontal col(blue) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"3 moves (ref = 0 moves)", col(blue) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"3 moves (ref = 0 moves)", horizontal col(blue) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"3 moves (ref = 0 moves)", col(blue) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"3 moves (ref = 0 moves)", horizontal col(blue) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "3 moves (ref = 0 moves)", col(blue) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "3 moves (ref = 0 moves)", horizontal col(blue) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "3 moves (ref = 0 moves)", col(blue) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "3 moves (ref = 0 moves)", horizontal col(blue) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"3 moves (ref = 0 moves)", col(blue) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"3 moves (ref = 0 moves)", horizontal col(blue) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"4 moves (ref = 0 moves)", col(green) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"4 moves (ref = 0 moves)", horizontal col(green) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"4 moves (ref = 0 moves)", col(green) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"4 moves (ref = 0 moves)", horizontal col(green) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"4 moves (ref = 0 moves)", col(green) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"4 moves (ref = 0 moves)", horizontal col(green) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"4 moves (ref = 0 moves)", col(green) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"4 moves (ref = 0 moves)", horizontal col(green) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "4 moves (ref = 0 moves)", col(green) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "4 moves (ref = 0 moves)", horizontal col(green) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "4 moves (ref = 0 moves)", col(green) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "4 moves (ref = 0 moves)", horizontal col(green) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"4 moves (ref = 0 moves)", col(green) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"4 moves (ref = 0 moves)", horizontal col(green) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"5 + moves (ref = 0 moves)", col(orange) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"5 + moves (ref = 0 moves)", horizontal col(orange) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"5 + moves (ref = 0 moves)", col(orange) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"5 + moves (ref = 0 moves)", horizontal col(orange) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"5 + moves (ref = 0 moves)", col(orange) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"5 + moves (ref = 0 moves)", horizontal col(orange) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"5 + moves (ref = 0 moves)", col(orange) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"5 + moves (ref = 0 moves)", horizontal col(orange) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "5 + moves (ref = 0 moves)", col(orange) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "5 + moves (ref = 0 moves)", horizontal col(orange) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "5 + moves (ref = 0 moves)", col(orange) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "5 + moves (ref = 0 moves)", horizontal col(orange) msize(vtiny) lwidth(vthin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"5 + moves (ref = 0 moves)", col(orange) msize(vtiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"5 + moves (ref = 0 moves)", horizontal col(orange) msize(vtiny) lwidth(vthin)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = 0 moves)") ///
		title("Residential mobility and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.02 0.1 0.2 0.5 1 2 5 10 15, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "1" 15 "2" 29 "3" 43 "4" 57 "5+") rows(1)) ///
		name(mobility_cat, replace)

graph export ".\G0Partner_Results\mobilityResults.pdf", replace


** Create plot for IPSM total score
sum lci uci if exposure == "IPSM_total" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"IPSM_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"IPSM_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"IPSM_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"IPSM_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"IPSM_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"IPSM_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "IPSM_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "IPSM_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "IPSM_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "IPSM_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "IPSM_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "IPSM_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "IPSM_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "IPSM_total", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in IPSM)") ///
		title("Inter-personal sensitivity and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.96 0.98 1 1.02 1.04, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(ipsm_cat, replace)
		
graph export ".\G0Partner_Results\ipsmResults.pdf", replace


** Create plot for locus of control score
sum lci uci if exposure == "LoC_external" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"LoC_external", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"LoC_external", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"LoC_external", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "LoC_external", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "LoC_external", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "LoC_external", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "LoC_external", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in LoC)") ///
		title("External Locus of Control and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(loc_cat, replace)
		
graph export ".\G0Partner_Results\locResults.pdf", replace


** Create plot for number of negative life events aged < 16
sum lci uci if exposure == "chLifeEvents_total" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"chLifeEvents_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"chLifeEvents_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"chLifeEvents_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"chLifeEvents_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"chLifeEvents_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"chLifeEvents_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "chLifeEvents_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "chLifeEvents_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "chLifeEvents_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "chLifeEvents_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "chLifeEvents_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "chLifeEvents_total", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "chLifeEvents_total", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "chLifeEvents_total", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per additional life event)") ///
		title("Adverse childhood experiences and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.85 0.9 0.95 1 1.05 1.1 1.15, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(ace_cat, replace)
		
graph export ".\G0Partner_Results\aceResults.pdf", replace


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
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Intrinsic (cat)" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - friends" & exp_level ///
			== "Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - prayer" & exp_level ///
			== "Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "DUREL (cat)" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Intrinsic (cat)" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - friends" & exp_level ///
			== "O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - prayer" & exp_level ///
			== "O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "DUREL (cat)" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "DUREL (cat)" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Intrinsic (cat)" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - friends" & exp_level ///
			== "A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - prayer" & exp_level ///
			== "A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "DUREL (cat)" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "DUREL (cat)" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Intrinsic (cat)" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - friends" & exp_level ///
			== "Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - prayer" & exp_level ///
			== "Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "DUREL (cat)" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal col(green) msize(vtiny) lwidth(thin)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = CSE/None)") ///
		title("Education*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.8 0.9 1 1.1 1.2 1.3 1.4, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Vocational" 15 "O-levels" 29 "A-levels" ///
			43 "Degree") rows(1)) ///
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
			"highSocClass", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Intrinsic (cat)" & ///
			exposure == "highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "highSocClass", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - friends" & ///
			exposure == "highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - friends" & ///
			exposure == "highSocClass", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - prayer" & ///
			exposure == "highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "highSocClass", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "DUREL (cat)" & ///
			exposure == "highSocClass", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "DUREL (cat)" & ///
			exposure == "highSocClass", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = lower [III manual/IV/V])") ///
		title("Occupational Class*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.8 0.9 1 1.1 1.2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(socClass_int, replace)
		
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
			"income", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Intrinsic (cat)" & ///
			exposure == "income", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "income", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - friends" & ///
			exposure == "income", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - friends" & ///
			exposure == "income", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - prayer" & ///
			exposure == "income", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "income", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "DUREL (cat)" & ///
			exposure == "income", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "DUREL (cat)" & ///
			exposure == "income", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in log income)") ///
		title("Household income*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.85 0.9 0.95 1 1.05 1.1 , labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(income_int, replace)
		
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
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Intrinsic (cat)" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - friends" & exp_level ///
			== "Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - prayer" & exp_level ///
			== "Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "DUREL (cat)" & exp_level == ///
			"Rent (ref = Own/Mortgage)", col(black) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Rent (ref = Own/Mortgage)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Intrinsic (cat)" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - friends" & exp_level ///
			== "Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - prayer" & exp_level ///
			== "Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "DUREL (cat)" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", col(red) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Council/HA (ref = Own/Mortgage)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Belief" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Belief" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Relig" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Relig" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Attend" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Attend" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Intrinsic (cat)" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - friends" & exp_level ///
			== "Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "Extrinsic - prayer" & exp_level ///
			== "Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef_int if outcome == "DUREL (cat)" & exp_level == ///
			"Other (ref = Own/Mortgage)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci_int uci_int level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Other (ref = Own/Mortgage)", horizontal col(blue) msize(vtiny) lwidth(thin)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = Owned/Mortgaged)") ///
		title("Housing status*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.6 0.7 0.8 0.9 1 1.2 1.4 1.7, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Rented" 15 "Council/HA" 29 "Other") ///
			rows(1)) ///
		name(housing_int, replace)
		
graph export ".\G0Partner_Results\housingResults_int.pdf", replace


** Create plot for locus of control by age interaction
sum lci_int uci_int if exposure == "LoC_external" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == ///
			"LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == ///
			"LoC_external", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == ///
			"LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == ///
			"LoC_external", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == ///
			"LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == ///
			"LoC_external", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Intrinsic (cat)" & ///
			exposure == "LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "LoC_external", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - friends" & ///
			exposure == "LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - friends" & ///
			exposure == "LoC_external", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - prayer" & ///
			exposure == "LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "LoC_external", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "DUREL (cat)" & ///
			exposure == "LoC_external", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "DUREL (cat)" & ///
			exposure == "LoC_external", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in LoC)") ///
		title("External LoC*Age Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.96 0.98 1 1.02 1.04, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(loc_int, replace)
		
graph export ".\G0Partner_Results\locResults_int.pdf", replace


graph close _all



*********************************************************************************
*** Comparing pregnancy vs 2019 data side-by-side (for belief in God, religious affiliation and church attendance) to see whether effect sizes (and not just p-values) differ between these different RSBB waves

** Wont do this for all exposures (as too many), but will focus on a few key ones
use ".\G0Partner_Results\partner_belief_results.dta", clear
gen outcome = "Belief_preg"

append using ".\G0Partner_Results\partner_belief_2019_results.dta"
replace outcome = "Belief_2019" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_relig_results.dta"
replace outcome = "Relig_preg" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_relig_2019_results.dta"
replace outcome = "Relig_2019" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_attend_results.dta"
replace outcome = "Attend_preg" if outcome == ""
tab outcome, m

append using ".\G0Partner_Results\partner_attend_2019_results.dta"
replace outcome = "Attend_2019" if outcome == ""
tab outcome, m

* Keep just the exposures we're interested in
keep if exposure == "highSocClass" | exposure == "IPSM_total" | exposure == "LoC_external" | exposure == "ageAt28" | exposure == "ageInPreg" | exposure == "education" | exposure == "income" | exposure == "intel_factor" | exposure == "maritalStatus" | exposure == "nonWhiteEthnic"

* Recode the outcomes
capture drop level_num
gen level_num = 0
replace level_num = 1 if outcome_level == "Not sure (ref = No)"
replace level_num = 3 if outcome_level == "Christian (ref = None)"
replace level_num = 4 if outcome_level == "Other (ref = None)"
replace level_num = 6 if outcome_level == "Min once year (ref = Not at al"
replace level_num = 7 if outcome_level == "Min once month (ref = Not at a"
replace level_num = 8 if outcome_level == "Min once week (ref = Not at al"

label define level_lb 0 "Belief in God - Yes (ref = No)" 1 "Belief in God - Not sure (ref = No)" 3 "Religious affiliation - Christian (ref = None)" 4 "Religious affiliation - Other (ref = None)" 6 "Church attendance - Min 1/Yr (ref = Not at all)" 7 "Church attendance - Min 1/Mth (ref = Not at all)" 8 "Church attendance - Min 1/Wk (ref = Not at all)", replace
label values level_num level_lb
tab level_num

* And split by whether from pregnancy or 2019
gen time = "preg" if outcome == "Belief_preg" | outcome == "Relig_preg" | outcome == "Attend_preg"
replace time = "2019" if outcome == "Belief_2019" | outcome == "Relig_2019" | outcome == "Attend_2019"
tab time

capture drop level_split
gen level_split = level_num - 0.2 if time == "preg"
replace level_split = level_num + 0.2 if time == "2019"
label values level_split level_lb
tab level_split


** Start with age plot
replace exposure = "Age" if exposure == "ageInPreg" | exposure == "ageAt28"
tab exposure

* Min and max x-axis values
sum lci uci if level_split < . & exposure == "Age"

twoway (scatter level_split coef if time == "preg" & exposure == "Age", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exposure == "Age", ///
			horizontal lwidth(thin) col(black)) ///
		(scatter level_split coef if time == "2019" & exposure == "Age", ///
			col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exposure == "Age", ///
			horizontal lwidth(thin) col(red)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio") ///
		title("Age and RSBB - Pregnancy vs 2019", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.98 1 1.02 1.04 1.06 1.08, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(order(1 "Pregnancy data" 3 "2019 data")) ///
		name(age, replace)
		
graph export ".\G0Partner_Results\Age_pregVs2019.pdf", replace


** Now ethnicity plot
sum lci uci if level_split < . & exposure == "nonWhiteEthnic"

* Drop some cases with no data
replace coef = . if exposure == "nonWhiteEthnic" & outcome_level == "Min once month (ref = Not at a" & time == "2019"
replace lci = . if exposure == "nonWhiteEthnic" & outcome_level == "Min once month (ref = Not at a" & time == "2019"

sum lci uci if level_split < . & exposure == "nonWhiteEthnic"

twoway (scatter level_split coef if time == "preg" & exposure == "nonWhiteEthnic", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exposure == "nonWhiteEthnic", ///
			horizontal lwidth(thin) col(black)) ///
		(scatter level_split coef if time == "2019" & exposure == "nonWhiteEthnic", ///
			col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exposure == "nonWhiteEthnic", ///
			horizontal lwidth(thin) col(red)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio") ///
		title("Ethnicity and RSBB - Pregnancy vs 2019", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.3 0.5 1 2 3 5 10 25, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(order(1 "Pregnancy data" 3 "2019 data")) ///
		name(ethnic, replace)
		
graph export ".\G0Partner_Results\Ethnic_pregVs2019.pdf", replace


** Next marital status

* Need to re-arrange level splitting here, as more than one marital status outcome
capture drop level_split
gen level_split = level_num - 0.2 - 0.08 if time == "preg" & exp_level == "Married (ref = Never married)"
replace level_split = level_num - 0.2 + 0.08 if time == "2019" & exp_level == "Married (ref = Never married)"
replace level_split = level_num + 0.2 - 0.08 if time == "preg" & exp_level == "Wid/Div/Sep (ref = Never married)"
replace level_split = level_num + 0.2 + 0.08 if time == "2019" & exp_level == "Wid/Div/Sep (ref = Never married)"
label values level_split level_lb
tab level_split

sum lci uci if level_split < . & exposure == "maritalStatus"

twoway (scatter level_split coef if time == "preg" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exp_level == ///
			"Married (ref = Never married)", horizontal lwidth(thin) col(black)) ///
		(scatter level_split coef if time == "2019" & exp_level == ///
			"Married (ref = Never married)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exp_level == ///
			"Married (ref = Never married)", horizontal lwidth(thin) ///
			col(black) lpattern(shortdash)) ///
		(scatter level_split coef if time == "preg" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal lwidth(thin) col(red)) ///
		(scatter level_split coef if time == "2019" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exp_level == ///
			"Wid/Div/Sep (ref = Never married)", horizontal lwidth(thin) ///#
			col(red) lpattern(shortdash)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio") ///
		title("Marital status and RSBB - Pregnancy vs 2019", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.25 0.5 1 2 5 10 20 50, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(2 "Married - preg" 4 "Married - 2019" ///
			6 "Wid/Div/Sep - preg" 8 "Wid/Div/Sep - 2019") size(small)) ///
		name(marital, replace)
		
graph export ".\G0Partner_Results\Marital_pregVs2019.pdf", replace


** Next education

* Need to re-arrange level splitting here, as more than one education outcome
capture drop level_split
gen level_split = level_num - 0.3 - 0.04 if time == "preg" & exp_level == "Vocational (ref = CSE/None)"
replace level_split = level_num - 0.3 + 0.04 if time == "2019" & exp_level == "Vocational (ref = CSE/None)"
replace level_split = level_num - 0.1 - 0.04 if time == "preg" & exp_level == "O-level (ref = CSE/None)"
replace level_split = level_num - 0.1 + 0.04 if time == "2019" & exp_level == "O-level (ref = CSE/None)"
replace level_split = level_num + 0.1 - 0.04 if time == "preg" & exp_level == "A-level (ref = CSE/None)"
replace level_split = level_num + 0.1 + 0.04 if time == "2019" & exp_level == "A-level (ref = CSE/None)"
replace level_split = level_num + 0.3 - 0.04 if time == "preg" & exp_level == "Degree (ref = CSE/None)"
replace level_split = level_num + 0.3 + 0.04 if time == "2019" & exp_level == "Degree (ref = CSE/None)"
label values level_split level_lb
tab level_split

sum lci uci if level_split < . & exposure == "education"

twoway (scatter level_split coef if time == "preg" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal lwidth(vthin) col(black)) ///
		(scatter level_split coef if time == "2019" & exp_level == ///
			"Vocational (ref = CSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exp_level == ///
			"Vocational (ref = CSE/None)", horizontal lwidth(vthin) ///
			col(black) lpattern(shortdash)) ///
		(scatter level_split coef if time == "preg" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal lwidth(vthin) col(red)) ///
		(scatter level_split coef if time == "2019" & exp_level == ///
			"O-level (ref = CSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exp_level == ///
			"O-level (ref = CSE/None)", horizontal lwidth(vthin) ///#
			col(red) lpattern(shortdash)) ///
		(scatter level_split coef if time == "preg" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal lwidth(vthin) col(blue)) ///
		(scatter level_split coef if time == "2019" & exp_level == ///
			"A-level (ref = CSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exp_level == ///
			"A-level (ref = CSE/None)", horizontal lwidth(vthin) ///
			col(blue) lpattern(shortdash)) ///
		(scatter level_split coef if time == "preg" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal lwidth(vthin) col(green)) ///
		(scatter level_split coef if time == "2019" & exp_level == ///
			"Degree (ref = CSE/None)", col(green) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exp_level == ///
			"Degree (ref = CSE/None)", horizontal lwidth(vthin) ///#
			col(green) lpattern(shortdash)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio") ///
		title("Education and RSBB - Pregnancy vs 2019", size(small)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.05 0.1 0.2 0.5 1 2 3 5 10, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(2 "Vocational - preg" 4 "Vocational - 2019" ///
			6 "O-level - preg" 8 "O-level - 2019" 10 "A-level - preg" ///
			12 "A-level - 2019" 14 "Degree - preg" 16 "Degree - 2019") ///
			size(vsmall) cols(4)) ///
		name(edu, replace) xsize(8)
		
graph export ".\G0Partner_Results\Education_pregVs2019.pdf", replace


** Next to occupational social class
capture drop level_split
gen level_split = level_num - 0.2 if time == "preg"
replace level_split = level_num + 0.2 if time == "2019"
label values level_split level_lb
tab level_split

sum lci uci if level_split < . & exposure == "highSocClass"

twoway (scatter level_split coef if time == "preg" & exposure == "highSocClass", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exposure == "highSocClass", ///
			horizontal lwidth(thin) col(black)) ///
		(scatter level_split coef if time == "2019" & exposure == "highSocClass", ///
			col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exposure == "highSocClass", ///
			horizontal lwidth(thin) col(red)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio") ///
		title("Soc. Class and RSBB - Pregnancy vs 2019", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.6 0.8 1 1.5 2 3 5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(order(1 "Pregnancy data" 3 "2019 data")) ///
		name(socClass, replace)
		
graph export ".\G0Partner_Results\socClass_pregVs2019.pdf", replace


** Next is income
sum lci uci if level_split < . & exposure == "income"

twoway (scatter level_split coef if time == "preg" & exposure == "income", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exposure == "income", ///
			horizontal lwidth(thin) col(black)) ///
		(scatter level_split coef if time == "2019" & exposure == "income", ///
			col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exposure == "income", ///
			horizontal lwidth(thin) col(red)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio") ///
		title("Income and RSBB - Pregnancy vs 2019", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.5 0.7 1 1.5 2 3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(order(1 "Pregnancy data" 3 "2019 data")) ///
		name(income, replace)
		
graph export ".\G0Partner_Results\income_pregVs2019.pdf", replace


** Total IPSM
sum lci uci if level_split < . & exposure == "IPSM_total"

twoway (scatter level_split coef if time == "preg" & exposure == "IPSM_total", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exposure == "IPSM_total", ///
			horizontal lwidth(thin) col(black)) ///
		(scatter level_split coef if time == "2019" & exposure == "IPSM_total", ///
			col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exposure == "IPSM_total", ///
			horizontal lwidth(thin) col(red)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio") ///
		title("IPSM and RSBB - Pregnancy vs 2019", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.98 0.99 1 1.01 1.02 1.03, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(order(1 "Pregnancy data" 3 "2019 data")) ///
		name(ipsm, replace)
		
graph export ".\G0Partner_Results\ipsm_pregVs2019.pdf", replace


** External locus of control
sum lci uci if level_split < . & exposure == "LoC_external"

twoway (scatter level_split coef if time == "preg" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "preg" & exposure == "LoC_external", ///
			horizontal lwidth(thin) col(black)) ///
		(scatter level_split coef if time == "2019" & exposure == "LoC_external", ///
			col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if time == "2019" & exposure == "LoC_external", ///
			horizontal lwidth(thin) col(red)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio") ///
		title("External LoC and RSBB - Pregnancy vs 2019", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.7 0.8 0.9 1 1.1 1.2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) ///
		angle(0)) legend(order(1 "Pregnancy data" 3 "2019 data")) ///
		name(LoC, replace)
		
graph export ".\G0Partner_Results\loc_pregVs2019.pdf", replace

graph close _all


** For the multinomial regression results, as interpretation not intuitive, could convert to predicted probabilities using the 'margins' command? (see: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/) - Have started this in the G0 mothers file; see there for example code.




