*** Predictors of RSBB (B3911) - G0 mother analysis script
*** Created 16/11/2021 by Dan Smith
*** Stata v16.0

*** This script reads in the cleaned G0 mother data, explores associations between exposures, and then conducts an 'exposome-wide association analysis' (ExWAS) on all these variables to examine how they are associated with various facets of RSBB.


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

** And also install 'spost' commands for testing proportional odds assumption of ordinal regression models (to install 'spost', type 'search spost13' and install the 'spost13_ado' package)


**********************************************************************************
*** Descriptive statistics

** Put the RSBB variables at the start of the dataset
order aln d810 d813 d813_grp d816 Y3000 Y3040 Y3040_grp Y3080 Y3080_OccNever Y3080_OccYr Y3153 Y3153_cat Y3160 Y3170 Y3155 Y3155_cat

** Using the 'distinct' command (see above), for each variable inspect the number of unque values; if < 10 then display table, while if >= 10 then displays means/SDs/IQRs

* Outcomes
foreach var of varlist d810-Y3155_cat {
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
foreach var of varlist mz028b-fom_cog_factor1 {
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
misstable sum d810-Y3155_cat, all

* Exposures
misstable sum mz028b-fom_cog_factor1, all



**********************************************************************************
*** Correlations between exposures

** Explore correlations between exposures to see how inter-related these factors are
desc mz028b-fom_cog_factor1

* Most variables are either continuous, ordered categories, or binary variables, so will just use standard pearson correlations for these. Only unordered categories are home ownership (owned/mortaged vs renting vs counci/housing association vs other) and marital status (never married vs married vs widowed/divorced/separated). So will exclude these two variables from the correlation matrix and calculate their associations using these variables as outcomes in a multinomial regression and square-rooting the pseudo-R2 value (similar to how the 'rexposome' package in R for ExWAS analyses does)

* First, order variables into categories of 'demographic', 'socioeconomic/material insecurity' and 'cognitive/psychological' (Y9992 is age at @28 RSBB questionnaire, so will pop at end as not needed here)
order aln-Y3155_cat ///
mz028b c800_grp a525_grp a005_grp jan1993ur01ind_grp b032_grp ///
c645a c686a c755_grp logavinceq jan1993imd2010q5_M jan1993Townsendq5_M a006_grp b594_grp c432 c433 a551 a636 partner_ab ///
logic_mem-fom_cog_factor1 b916-b921 d842 h151b ///
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
rename c755_grp highSocClass
rename logavinceq income
rename jan1993imd2010q5_M IMD
rename jan1993Townsendq5_M townsendDep
rename a006_grp housing 
rename b594_grp financeDiffs
rename c432 chLifeEvents_wgt
rename c433 chLifeEvents_total
rename a551 crowding
rename a636 neighPercept
rename partner_ab partnerAbsence
rename logic_mem logicMemory
rename digit_back digitBack
rename spot_word spotWord
rename digit_symbol digitSymbol
rename logic_mem_delay logicMemory_delay 
rename fom_cog_factor1 intel_factor
rename b916 IPSM_interpAware
rename b917 IPSM_approval
rename b918 IPSM_sepAnx
rename b919 IPSM_timidity
rename b920 IPSM_fragility
rename b921 IPSM_total
rename d842 LoC_external
rename h151b selfEsteem
rename Y9992 ageAt28

* Associations between demographic factors (excluding marital status; a5252_grp) - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
corr ageAtBirth nonWhiteEthnic mobility rural parity

matrix cor_demo = r(C)
matrix list cor_demo

heatplot cor_demo, values(format(%9.3f)) color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) legend(off) xlabel(, angle(45))

* Save heatmap
graph export ".\G0Mother_Results\corr_heatplot_demoOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Mother_Results\corrMatrix_demoOnly.xlsx", replace
putexcel A1=matrix(cor_demo), names

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
corr education maternalEdu highSocClass income IMD townsendDep financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence

matrix cor_socio = r(C)

* As lots of entries, the correlation coefficients are hard to read, so will drop these values and include the legend
heatplot cor_socio, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))
	
* Save heatmap
graph export ".\G0Mother_Results\corr_heatplot_socioOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Mother_Results\corrMatrix_socioOnly.xlsx", replace
putexcel A1=matrix(cor_socio), names


* And now for associations of socioeconomic variables with unordered categorical variable home ownership status (and save to a CSV file)
capture postclose housing_corrs_socioOnly
postfile housing_corrs_socioOnly str30 variable corr ///
	using ".\G0Mother_Results\housing_corrs_socioOnly.dta", replace
	
foreach var of varlist education maternalEdu highSocClass income IMD townsendDep financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence {
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


** Now repeat for cognitive/psychological variables (no unordered categorical variables, so can include all variables)
corr logicMemory-selfEsteem

matrix cor_cog = r(C)

heatplot cor_cog, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))

* Save heatmap
graph export ".\G0Mother_Results\corr_heatplot_cogOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Mother_Results\corrMatrix_cogOnly.xlsx", replace
putexcel A1=matrix(cor_cog), names

	
*** Finally, repeat this on all of the exposures together (excluding unordered cateogorical variables housing status and marital status)
corr ageAtBirth nonWhiteEthnic mobility rural parity education maternalEdu highSocClass income IMD townsendDep financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem

matrix cor_all = r(C)

heatplot cor_all, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle(""))

* Save heatmap
graph export ".\G0Mother_Results\corr_heatplot_all.pdf", replace

* Save matrix as Excel file
putexcel set ".\G0Mother_Results\corrMatrix_all.xlsx", replace
putexcel A1=matrix(cor_all), names

	
* And now for associations of all other exposures variables with unordered categorical variables marital status and home ownership status

* Marital status
capture postclose marital_corrs_all
postfile marital_corrs_all str30 variable corr ///
	using ".\G0Mother_Results\marital_corrs_all.dta", replace
	
foreach var of varlist ageAtBirth nonWhiteEthnic mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
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
	
foreach var of varlist ageAtBirth nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
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



**************************************************************************************
*** Next, we want to run the actual ExWAS analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab d810, m

** This will be quite complicated, as want to post results to file, but as exposures differ extracting the results will be variable specific. To adjust for multiple corrections will use conservative bonferroni adjustment when constructing confidence intervals and interpreting p-values - As 34 exposures, will use 99.85% confidence intervals (as this is 100 - 5/34) and a p-value threshold of 0.05/34 = 0.0015.
display round(100 - 5/34, 0.01)
display 0.05/34

** We also want to store both estimates adjusting for age (other than for the age-only model), and also the interaction between age and the exposure, to see whether it's moderated by age. Again, this makes the set-up a bit more complicated.

** Create a postfile to post results to, then start the loop - Will create two postfiles; one for coefficients and CIs, and another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model)
capture postclose mother_belief
postfile mother_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_belief_results.dta", replace

capture postclose mother_belief_lr
postfile mother_belief_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_belief_results_lr.dta", replace

foreach var of varlist ageAtBirth nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAtBirth' separately first as will be adjusted for in all other models
	if "`var'" == "ageAtBirth" {
		mlogit d810 `var', baseoutcome(3) rrr level(99.85)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit d810 if `var' != ., baseoutcome(3) rrr level(99.85)
		est store base
		mlogit d810 `var', baseoutcome(3) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit d810 ageAtBirth `var', baseoutcome(3) rrr level(99.85)
		
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
		mlogit d810 c.ageAtBirth##c.`var', baseoutcome(3) rrr level(99.85)
		
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
		mlogit d810 ageAtBirth `var', baseoutcome(3) rrr level(99.85)
		
		local outcome_level = "Not sure (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,4]
		local lci = res[5,4]
		local uci = res[6,4]
		local p = res[4,4]
				
		// Now for interaction model
		mlogit d810 c.ageAtBirth##c.`var', baseoutcome(3) rrr level(99.85)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit d810 ageAtBirth if `var' != ., baseoutcome(3) rrr level(99.85)
		est store base
		mlogit d810 ageAtBirth `var', baseoutcome(3) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit d810 c.ageAtBirth##c.`var', baseoutcome(3) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
		
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
		
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
		
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit d810 ageAtBirth if `var' != ., baseoutcome(3) rrr level(99.85)
		est store base
		mlogit d810 ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit d810 c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_belief
postclose mother_belief_lr
	

****************************************************************************************
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
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_relig_results.dta", replace

capture postclose mother_relig_lr
postfile mother_relig_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_relig_results_lr.dta", replace

foreach var of varlist ageAtBirth nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAtBirth' separately first as will be adjusted for in all other models
	if "`var'" == "ageAtBirth" {
		mlogit d813_grp `var', baseoutcome(3) rrr level(99.85)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit d813_grp if `var' != ., baseoutcome(3) rrr level(99.85)
		est store base
		mlogit d813_grp `var', baseoutcome(3) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit d813_grp ageAtBirth `var', baseoutcome(3) rrr level(99.85)
		
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
		mlogit d813_grp c.ageAtBirth##c.`var', baseoutcome(3) rrr level(99.85)
		
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
		mlogit d813_grp ageAtBirth `var', baseoutcome(3) rrr level(99.85)
		
		local outcome_level = "Other (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,4]
		local lci = res[5,4]
		local uci = res[6,4]
		local p = res[4,4]
				
		// Now for interaction model
		mlogit d813_grp c.ageAtBirth##c.`var', baseoutcome(3) rrr level(99.85)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit d813_grp ageAtBirth if `var' != ., baseoutcome(3) rrr level(99.85)
		est store base
		mlogit d813_grp ageAtBirth `var', baseoutcome(3) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit d813_grp c.ageAtBirth##c.`var', baseoutcome(3) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
		
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,9]
			local lci = res[5,9]
			local uci = res[6,9]
			local p = res[4,9]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
		
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
		
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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
			mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit d813_grp ageAtBirth if `var' != ., baseoutcome(3) rrr level(99.85)
		est store base
		mlogit d813_grp ageAtBirth i.`var', baseoutcome(3) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit d813_grp c.ageAtBirth##i.`var', baseoutcome(3) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_relig
postclose mother_relig_lr



****************************************************************************************
*** Now to the next RSBB outcome: Attendance at church/place of worship

*** As this is an ordered categorical variable, will use ordinal regression model. For consistency with previous models, will recode so that higher values indicate greater RSBB
tab d816

recode d816 (1 = 3) (2 = 2) (3 = 1) (4 = 0), gen(d816_rev)
label define attend_rev_lb 0 "Not at all" 1 "Min once a year" 2 "Min once a month" 3 "Min once a week"
numlabel attend_rev_lb, add
label values d816_rev attend_rev_lb
tab d816_rev


** Quick test of whether proportional odds assumption been violated in most basic model (with just age at birth). Ah, it has been violated, and gets worse if add in additional predictors. 
ologit d816_rev ageAtBirth, or level(99.85)
brant, detail

ologit d816_rev ageAtBirth i.IMD, or level(99.85)
brant, detail


** So instead will just run multinomial models with 'not at all' as the baseline/reference category.


*** Now run the loop to save all the results
capture postclose mother_attend
postfile mother_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_attend_results.dta", replace

capture postclose mother_attend_lr
postfile mother_attend_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_attend_results_lr.dta", replace

foreach var of varlist ageAtBirth nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAtBirth' separately first as will be adjusted for in all other models
	if "`var'" == "ageAtBirth" {
		mlogit d816_rev `var', baseoutcome(0) rrr level(99.85)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit d816_rev if `var' != ., baseoutcome(0) rrr level(99.85)
		est store base
		mlogit d816_rev `var', baseoutcome(0) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr level(99.85)
		
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
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr level(99.85)
		
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
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr level(99.85)
		
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr level(99.85)
		
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
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr level(99.85)
		
		local outcome_level = "Min once week (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr level(99.85)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit d816_rev ageAtBirth if `var' != ., baseoutcome(0) rrr level(99.85)
		est store base
		mlogit d816_rev ageAtBirth `var', baseoutcome(0) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit d816_rev c.ageAtBirth##c.`var', baseoutcome(0) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
		
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
				
			// Now onto the next reference category (3/once year)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
		
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
		
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,29]
		
			post mother_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (2/once month)
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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
			mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		
			local outcome_level = "Min once week (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit d816_rev ageAtBirth if `var' != ., baseoutcome(0) rrr level(99.85)
		est store base
		mlogit d816_rev ageAtBirth i.`var', baseoutcome(0) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit d816_rev c.ageAtBirth##i.`var', baseoutcome(0) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_attend
postclose mother_attend_lr



****************************************************************************************
*** Now to the next RSBB outcome: Intrinsic religiosity

** As this was collected 28 years post-birth and not all mothers who completed this will have an age at birth (as enrolled in later phases), will use age at questionnaire completion instead.

* Intrinsic religiosity outcome is very non-normal (spike at lowesst value '3', then uniform), so will also run multinomial model using categories of this variable (with '3' as baseline category).
tab1 Y3153 Y3153_cat
sum Y3153
hist Y3153, freq width(1)


*** Start with continuous measure and linear regression

*** Now run the loop to save all the results
capture postclose mother_intrinsic
postfile mother_intrinsic str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_intrinsic_results.dta", replace

capture postclose mother_intrinsic_lr
postfile mother_intrinsic_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_intrinsic_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		regress Y3153 `var', level(99.85)
		
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
		
		post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// And finally run the likelihood ratio tests
		regress Y3153 if `var' != ., level(99.85)
		est store base
		regress Y3153 `var', level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		regress Y3153 ageAt28 `var', level(99.85)
		
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
		regress Y3153 c.ageAt28##c.`var', level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
		// And finally run the likelihood ratio tests
		regress Y3153 ageAt28 if `var' != ., level(99.85)
		est store base
		regress Y3153 ageAt28 `var', level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress Y3153 c.ageAt28##c.`var', level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			regress Y3153 ageAt28 i.`var', level(99.85)
		
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
					
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			regress Y3153 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 4)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			regress Y3153 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
						
			
			// Move to the next category of the exposure (category 4)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
							
			// Move to the next category of the exposure (category 5)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			regress Y3153 ageAt28 i.`var', level(99.85)
		
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 5)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
							
			
			// Move to the next category of the exposure (category 6)
			regress Y3153 ageAt28 i.`var', level(99.85)
		
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
			regress Y3153 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post mother_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		regress Y3153 ageAt28 if `var' != ., level(99.85)
		est store base
		regress Y3153 ageAt28 i.`var', level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress Y3153 c.ageAt28##i.`var', level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_intrinsic
postclose mother_intrinsic_lr


**** And now repeat for the categorical intrinsic religiosity variable as a sensitivity analysis to ensure results robust and not due to odd distribution of outcome
tab Y3153_cat

*** Now run the loop to save all the results
capture postclose mother_intrinsic_cat
postfile mother_intrinsic_cat str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_intrinsic_cat_results.dta", replace

capture postclose mother_intrinsic_cat_lr
postfile mother_intrinsic_cat_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_intrinsic_cat_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit Y3153_cat `var', baseoutcome(1) rrr level(99.85)
		
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
		
		post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit Y3153_cat if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3153_cat `var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit Y3153_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
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
		mlogit Y3153_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/high)
		mlogit Y3153_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "High IR/8-11 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit Y3153_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (4/highest)
		mlogit Y3153_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit Y3153_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit Y3153_cat ageAt28 if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3153_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit Y3153_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,29]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post mother_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit Y3153_cat ageAt28 if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3153_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit Y3153_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_intrinsic_cat
postclose mother_intrinsic_cat_lr



**********************************************************************************
*** Next, want to explore extrinsic religiosity - As this is two separate questions (make friends and pray for protection), will have two run two separate loops for each outcome
tab1 Y3160 Y3170

** To reduce number of categories and combine low cell counts, will combine 'mildly' and 'strongly' agree and disagree together, to give 4 categories: Agree, not sure, disagree and not applicable. Will have 'agree' as the baseline, which is the most 'extrinsic' response.
recode Y3160 (1 2 = 1) (3 = 2) (4 5 = 3) (6 = 4), gen(Y3160_new)

label define extrinsic_lb 1 "Agree" 2 "Not sure" 3 "Disagree" 4 "Not applicable"
numlabel extrinsic_lb, add
label values Y3160_new extrinsic_lb
tab Y3160_new

recode Y3170 (1 2 = 1) (3 = 2) (4 5 = 3) (6 = 4), gen(Y3170_new)
label values Y3170_new extrinsic_lb
tab Y3170_new


*** Start with Y3160 - Attends church to help them make friends

*** Now run the loop to save all the results
capture postclose mother_extrinsic_friends
postfile mother_extrinsic_friends str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_extrinsic_friends_results.dta", replace

capture postclose mother_extrinsic_friends_lr
postfile mother_extrinsic_friends_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_extrinsic_friends_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit Y3160_new `var', baseoutcome(1) rrr level(99.85)
		
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
		
		post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit Y3160_new if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3160_new `var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit Y3160_new ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
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
		mlogit Y3160_new c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		mlogit Y3160_new ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "Disagree ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit Y3160_new c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (4/NA)
		mlogit Y3160_new ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "Not applicable ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit Y3160_new c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit Y3160_new ageAt28 if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3160_new ageAt28 `var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit Y3160_new c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,29]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post mother_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit Y3160_new ageAt28 if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3160_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit Y3160_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_extrinsic_friends
postclose mother_extrinsic_friends_lr



*** And now repeat with Y3170 - Prays for relief and protection

*** Now run the loop to save all the results
capture postclose mother_extrinsic_prayer
postfile mother_extrinsic_prayer str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_extrinsic_prayer_results.dta", replace

capture postclose mother_extrinsic_prayer_lr
postfile mother_extrinsic_prayer_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_extrinsic_prayer_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit Y3170_new `var', baseoutcome(1) rrr level(99.85)
		
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
		
		post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit Y3170_new if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3170_new `var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit Y3170_new ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
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
		mlogit Y3170_new c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		mlogit Y3170_new ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "Disagree ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit Y3170_new c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (4/NA)
		mlogit Y3170_new ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "Not applicable ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit Y3170_new c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit Y3170_new ageAt28 if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3170_new ageAt28 `var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit Y3170_new c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,29]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post mother_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit Y3170_new ageAt28 if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3170_new ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit Y3170_new c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_extrinsic_prayer
postclose mother_extrinsic_prayer_lr



**************************************************************************************
*** Now finally to the last RSBB outcome: Total DUREL religiosity score

* Total DUREL religiosity outcome is very non-normal (spike at lowesst value '5', then uniform), so will also run multinomial model using categories of this variable (with '5' as baseline category).
tab1 Y3155 Y3155_cat
sum Y3155
hist Y3155, freq width(1)


*** Start with continuous measure and linear regression

*** Now run the loop to save all the results
capture postclose mother_DUREL
postfile mother_DUREL str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_DUREL_results.dta", replace

capture postclose mother_DUREL_lr
postfile mother_DUREL_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_DUREL_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		regress Y3155 `var', level(99.85)
		
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
		
		post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// And finally run the likelihood ratio tests
		regress Y3155 if `var' != ., level(99.85)
		est store base
		regress Y3155 `var', level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		regress Y3155 ageAt28 `var', level(99.85)
		
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
		regress Y3155 c.ageAt28##c.`var', level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local age_main = res[1,1]
		local exp_main = res[1,2]
		
		post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
		// And finally run the likelihood ratio tests
		regress Y3155 ageAt28 if `var' != ., level(99.85)
		est store base
		regress Y3155 ageAt28 `var', level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress Y3155 c.ageAt28##c.`var', level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			regress Y3155 ageAt28 i.`var', level(99.85)
		
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,6]
			local lci_int = res[5,6]
			local uci_int = res[6,6]
			local p_int = res[4,6]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
					
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			regress Y3155 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 4)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			regress Y3155 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
						
			
			// Move to the next category of the exposure (category 4)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
							
			// Move to the next category of the exposure (category 5)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			regress Y3155 ageAt28 i.`var', level(99.85)
		
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local age_main = res[1,1]
			local exp_main = res[1,3]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local age_main = res[1,1]
			local exp_main = res[1,4]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (1/once year)
			local outcome_level = "Min once year (ref = Not at all)"
			
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local age_main = res[1,1]
			local exp_main = res[1,5]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 5)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local age_main = res[1,1]
			local exp_main = res[1,6]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
							
			
			// Move to the next category of the exposure (category 6)
			regress Y3155 ageAt28 i.`var', level(99.85)
		
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
			regress Y3155 c.ageAt28##i.`var', level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local age_main = res[1,1]
			local exp_main = res[1,7]
		
			post mother_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		regress Y3155 ageAt28 if `var' != ., level(99.85)
		est store base
		regress Y3155 ageAt28 i.`var', level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress Y3155 c.ageAt28##i.`var', level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_DUREL
postclose mother_DUREL_lr


**** And now repeat for the categorical total DUREL religiosity variable as a sensitivity analysis to ensure results robust and not due to odd distribution of outcome
tab Y3155_cat

*** Now run the loop to save all the results
capture postclose mother_DUREL_cat
postfile mother_DUREL_cat str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int age_main exp_main ///
	using ".\G0Mother_Results\mother_DUREL_cat_results.dta", replace

capture postclose mother_DUREL_cat_lr
postfile mother_DUREL_cat_lr str30 exposure lr_p_main lr_p_int ///
	using ".\G0Mother_Results\mother_DUREL_cat_results_lr.dta", replace

foreach var of varlist ageAt28 nonWhiteEthnic maritalStatus mobility rural parity education maternalEdu highSocClass income IMD townsendDep housing financeDiffs chLifeEvents_wgt chLifeEvents_total crowding neighPercept partnerAbsence logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit Y3155_cat `var', baseoutcome(1) rrr level(99.85)
		
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
		
		post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
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
		
		post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit Y3155_cat if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3155_cat `var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for age, will just fill with missing value
		local lr_p_int = .
		
		post mother_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "chLifeEvents_wgt" | "`var'" == "chLifeEvents_total" | "`var'" == "neighPercept" | "`var'" == "partnerAbsence" | "`var'" == "logicMemory" | "`var'" == "digitBack" | "`var'" == "spotWord" | "`var'" == "digitSymbol" | "`var'" == "verbal" | "`var'" == "logicMemory_delay" | "`var'" == "intel_factor" | "`var'" == "IPSM_interpAware" | "`var'" == "IPSM_approval" | "`var'" == "IPSM_sepAnx" | "`var'" == "IPSM_timidity" | "`var'" == "IPSM_fragility" | "`var'" == "IPSM_total" | "`var'" == "LoC_external" | "`var'" == "selfEsteem" {
		
		mlogit Y3155_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
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
		mlogit Y3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local age_main = res[1,5]
		local exp_main = res[1,6]
		
		post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (3/11-15 DUREL)
		mlogit Y3155_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "11-15 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
				
		// Now for interaction model
		mlogit Y3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local age_main = res[1,9]
		local exp_main = res[1,10]
		
		post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (4/16-20 DUREL)
		mlogit Y3155_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "16-20 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit Y3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local age_main = res[1,13]
		local exp_main = res[1,14]
		
		post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		// Now onto the next reference category (5/21-26 DUREL)
		mlogit Y3155_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		
		local outcome_level = "21-26 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,14]
		local lci = res[5,14]
		local uci = res[6,14]
		local p = res[4,14]
				
		// Now for interaction model
		mlogit Y3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		
		matrix res = r(table)
		local coef_int = res[1,19]
		local lci_int = res[5,19]
		local uci_int = res[6,19]
		local p_int = res[4,19]
		local age_main = res[1,17]
		local exp_main = res[1,18]
		
		post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit Y3155_cat ageAt28 if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3155_cat ageAt28 `var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit Y3155_cat c.ageAt28##c.`var', baseoutcome(1) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local age_main = res[1,9]
			local exp_main = res[1,11]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,17]
			local exp_main = res[1,19]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,33]
			local exp_main = res[1,35]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,15]
			local lci_int = res[5,15]
			local uci_int = res[6,15]
			local p_int = res[4,15]
			local age_main = res[1,9]
			local exp_main = res[1,12]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,17]
			local exp_main = res[1,20]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,33]
			local exp_main = res[1,36]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
						
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local age_main = res[1,11]
			local exp_main = res[1,13]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,21]
			local exp_main = res[1,23]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,31]
			local exp_main = res[1,33]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,41]
			local exp_main = res[1,43]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,18]
			local lci_int = res[5,18]
			local uci_int = res[6,18]
			local p_int = res[4,18]
			local age_main = res[1,11]
			local exp_main = res[1,14]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local age_main = res[1,21]
			local exp_main = res[1,24]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,31]
			local exp_main = res[1,34]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,48]
			local lci_int = res[5,48]
			local uci_int = res[6,48]
			local p_int = res[4,48]
			local age_main = res[1,41]
			local exp_main = res[1,44]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local age_main = res[1,11]
			local exp_main = res[1,15]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local age_main = res[1,21]
			local exp_main = res[1,25]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,31]
			local exp_main = res[1,35]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,49]
			local lci_int = res[5,49]
			local uci_int = res[6,49]
			local p_int = res[4,49]
			local age_main = res[1,41]
			local exp_main = res[1,45]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local age_main = res[1,25]
			local exp_main = res[1,27]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local age_main = res[1,37]
			local exp_main = res[1,39]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,56]
			local lci_int = res[5,56]
			local uci_int = res[6,56]
			local p_int = res[4,56]
			local age_main = res[1,49]
			local exp_main = res[1,51]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local age_main = res[1,13]
			local exp_main = res[1,16]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,33]
			local lci_int = res[5,33]
			local uci_int = res[6,33]
			local p_int = res[4,33]
			local age_main = res[1,25]
			local exp_main = res[1,28]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,45]
			local lci_int = res[5,45]
			local uci_int = res[6,45]
			local p_int = res[4,45]
			local age_main = res[1,37]
			local exp_main = res[1,40]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,57]
			local lci_int = res[5,57]
			local uci_int = res[6,57]
			local p_int = res[4,57]
			local age_main = res[1,49]
			local exp_main = res[1,52]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local age_main = res[1,13]
			local exp_main = res[1,15]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local age_main = res[1,25]
			local exp_main = res[1,29]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,46]
			local lci_int = res[5,46]
			local uci_int = res[6,46]
			local p_int = res[4,46]
			local age_main = res[1,37]
			local exp_main = res[1,41]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,58]
			local lci_int = res[5,58]
			local uci_int = res[6,58]
			local p_int = res[4,58]
			local age_main = res[1,49]
			local exp_main = res[1,53]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "maternalEdu" {
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,13]
			local exp_main = res[1,18]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local age_main = res[1,25]
			local exp_main = res[1,30]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,47]
			local lci_int = res[5,47]
			local uci_int = res[6,47]
			local p_int = res[4,47]
			local age_main = res[1,37]
			local exp_main = res[1,42]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,59]
			local lci_int = res[5,59]
			local uci_int = res[6,59]
			local p_int = res[4,59]
			local age_main = res[1,49]
			local exp_main = res[1,54]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local age_main = res[1,15]
			local exp_main = res[1,17]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local age_main = res[1,29]
			local exp_main = res[1,31]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local age_main = res[1,43]
			local exp_main = res[1,45]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,35]
			local lci = res[5,35]
			local uci = res[6,35]
			local p = res[4,35]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,65]
			local lci_int = res[5,65]
			local uci_int = res[6,65]
			local p_int = res[4,65]
			local age_main = res[1,57]
			local exp_main = res[1,59]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local age_main = res[1,15]
			local exp_main = res[1,18]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local age_main = res[1,29]
			local exp_main = res[1,32]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local age_main = res[1,43]
			local exp_main = res[1,46]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,36]
			local lci = res[5,36]
			local uci = res[6,36]
			local p = res[4,36]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,66]
			local lci_int = res[5,66]
			local uci_int = res[6,66]
			local p_int = res[4,66]
			local age_main = res[1,57]
			local exp_main = res[1,60]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local age_main = res[1,15]
			local exp_main = res[1,29]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,39]
			local lci_int = res[5,39]
			local uci_int = res[6,39]
			local p_int = res[4,39]
			local age_main = res[1,29]
			local exp_main = res[1,33]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local age_main = res[1,43]
			local exp_main = res[1,47]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,37]
			local lci = res[5,37]
			local uci = res[6,37]
			local p = res[4,37]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,67]
			local lci_int = res[5,67]
			local uci_int = res[6,67]
			local p_int = res[4,67]
			local age_main = res[1,57]
			local exp_main = res[1,61]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local age_main = res[1,15]
			local exp_main = res[1,20]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local age_main = res[1,29]
			local exp_main = res[1,34]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local age_main = res[1,43]
			local exp_main = res[1,48]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,38]
			local lci = res[5,38]
			local uci = res[6,38]
			local p = res[4,38]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,68]
			local lci_int = res[5,68]
			local uci_int = res[6,68]
			local p_int = res[4,68]
			local age_main = res[1,57]
			local exp_main = res[1,62]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
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
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local age_main = res[1,15]
			local exp_main = res[1,21]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local age_main = res[1,29]
			local exp_main = res[1,35]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local age_main = res[1,43]
			local exp_main = res[1,49]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,39]
			local lci = res[5,39]
			local uci = res[6,39]
			local p = res[4,39]
				
			// Now for interaction model
			mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		
			matrix res = r(table)
			local coef_int = res[1,69]
			local lci_int = res[5,69]
			local uci_int = res[6,69]
			local p_int = res[4,69]
			local age_main = res[1,57]
			local exp_main = res[1,63]
		
			post mother_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`age_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit Y3155_cat ageAt28 if `var' != ., baseoutcome(1) rrr level(99.85)
		est store base
		mlogit Y3155_cat ageAt28 i.`var', baseoutcome(1) rrr level(99.85)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit Y3155_cat c.ageAt28##i.`var', baseoutcome(1) rrr level(99.85)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post mother_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose mother_DUREL_cat
postclose mother_DUREL_cat_lr



***********************************************************************************
***********************************************************************************
*** Next step is to make some nice plots of these results


** For the multinomial regression results, as interpretation not intuitive, could convert to predicted probabilities using the 'margins' command? (see: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/)



