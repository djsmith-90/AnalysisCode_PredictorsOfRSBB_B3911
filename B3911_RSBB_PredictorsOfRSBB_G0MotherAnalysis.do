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

** Create a postfile to post results to, then start the loop
capture postclose mother_belief
postfile mother_belief str30 exposure str30 outcome_level str30 exp_level /// 
	n coef lci uci p coef_int lci_int uci_int p_int ///
	using ".\G0Mother_Results\mother_belief_results.dta", replace


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
		
		post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int')
			
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
		
		post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int')
		
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
		
		post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int')
			
		// Now onto the next reference category (2/not sure) - Need to run the models again
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
		
		post mother_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed
		
}

postclose mother_belief
	
	
mlogit d810 ageAtBirth, baseoutcome(3)
mlogit d810 ageAtBirth, baseoutcome(3) rrr

matrix res = r(table)
matrix list res

local bon_level = round(100 - 5/34, 0.01)
display `bon_level'

mlogit d810 ageAtBirth parity, baseoutcome(3) rrr level(99.85)

matrix res = r(table)
matrix list res

mlogit d810 ageAtBirth financeDiffs, baseoutcome(3) rrr level(99.85)
mlogit d810 c.ageAtBirth##c.financeDiffs, baseoutcome(3) rrr level(99.85)

matrix res = r(table)
matrix list res
	