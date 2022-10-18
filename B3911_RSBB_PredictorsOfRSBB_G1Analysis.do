*** Predictors of RSBB (B3911) - G1 analysis script
*** Created 16/11/2021 by Dan Smith
*** Stata v17.0

*** This script reads in the cleaned G1 data, explores associations between exposures, and then conducts an analysis exploring how all the demographic and SES variables are associated with various facets of RSBB.


**********************************************************************************
**** Set working directory, start a log file, and read in cleaned dataset

cd "X:\Groups\ARC\DanS\Descriptive_PredictorsOfRSBB_B3911"

capture log close
log using ".\G1_Results\Desc_RSBB_B3911_G1Analysis_log", replace text

use "G1_PredictorsOfRSBB_B3911.dta", clear


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
order aln qlet YPG3000 YPG3040 YPG3040_grp YPG3080 YPG3080_OccNever YPG3080_OccYr YPG3153 YPG3153_cat YPG3160 YPG3170 YPG3155 YPG3155_cat

* Keep just RSBB variables asked to G0's in pregnancy (belief in God, religious affiliation and frequency of church attendance)
drop YPG3153 YPG3153_cat YPG3160 YPG3170 YPG3155 YPG3155_cat

* Also drop the cognitive/psychological variables (as will focus on these in another paper), and also the adverse childhood experiences variables
drop f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126
drop clon100-clon109 clon114-clon123

** Using the 'distinct' command (see above), for each variable inspect the number of unque values; if < 10 then display table, while if >= 10 then displays means/SDs/IQRs

* Recode the attendance var first, so that categories are 'at least once a month', 'at least once a year', 'occasionally' and 'not at all'
tab1 YPG3080*, m
recode YPG3080 (1 2 = 1) (3 = 2) (4 = 3) (5 = 4), gen(YPG3080_alt)
label define attend2_lb 1 "Min once a month" 2 "Min once a year" 3 "Occasionally" 4 "Not at all"
numlabel attend2_lb, add
label values YPG3080_alt attend2_lb
tab YPG3080_alt

* Outcomes
foreach var of varlist YPG3000-YPG3080_OccYr YPG3080_alt {
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
foreach var of varlist kz021-age_FA {
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

* Change sex to a binary marker where 1 = male
replace kz021 = 0 if kz021 == 2
rename kz021 male
label define male_lb 0 "Female" 1 "Male"
numlabel male_lb, add
label values male male_lb

tab male

** Amount of missing data in outcomes and exposures
* Outcomes
misstable sum YPG3000-YPG3080_OccYr, all

* Exposures
misstable sum male-age_FA, all


** Also want to get descriptive stats for each exposure split by each RSBB outcome category
foreach var of varlist male-age_FA {
	quietly distinct `var'
	local unique = r(ndistinct)
	
	display ""
	display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	display "Variable " "`var'" " has " `unique' " values."
	
	if `unique' < 10 {
		tab YPG3000 `var', col
		tab YPG3040_grp `var', col
		tab YPG3080_alt `var', col
	}
	else {
		tab YPG3000, sum(`var')
		tab YPG3040_grp, sum(`var')
		tab YPG3080_alt, sum(`var')
	}
}


**********************************************************************************
*** Correlations between exposures

** Explore correlations between exposures to see how inter-related these factors are
desc male-age_FA

* Most variables are either continuous, ordered categories, or binary variables, so will just use standard pearson correlations for these. Only unordered categories are home ownership status (of YP and mother; owned/mortaged vs renting vs counci/housing association vs other) and mother;s marital status (never married vs married vs widowed/divorced/separated). So will exclude these three variables from the correlation matrix and calculate their associations using these variables as outcomes in a multinomial regression and square-rooting the pseudo-R2 value (similar to how the 'rexposome' package in R for ExWAS analyses does)

* First, order variables into categories of 'demographic' and 'socioeconomic/material insecurity'
order aln-YPG3080_OccYr YPG3080_alt ///
YPG8000 mz028b male c804 YPG1052 a525_grp a005_grp jan2020ur01ind_grp jan1993ur01ind_grp parent b032_grp ///
yp_edu c645a c666a yp_employed c755_grp c765_grp clon111 yp_income logavinceq jan2020imd2010q5 jan1993Townsendq5_M jan2020Townsendq5 jan1993imd2010q5_M yp_housing a006_grp yp_finDiffs clon112 c525 a053 a551 clon113 a636 father_ab age_FA

* Next, rename all these exposures so they are more intuitive and can be read easier on the correlation heatmaps below
rename YPG8000 ageAt28
rename mz028b mother_ageAtBirth
rename c804 nonWhiteEthnic
rename YPG1052 livePartner
rename a525_grp maritalStatus_mat
rename a005_grp mobility_mat
rename jan2020ur01ind_grp rural
rename jan1993ur01ind_grp rural_mat
rename b032_grp parity_mat
rename yp_edu education
rename c645a education_mat
rename c666a education_pat
rename yp_employed employed
rename c755_grp highSocClass_mat
rename c765_grp highSocClass_pat
rename clon111 lowSocClass_0_16
rename yp_income income
rename logavinceq income_parents
rename jan2020imd2010q5 IMD
rename jan1993Townsendq5_M IMD_mat
rename jan2020Townsendq5 townsendDep
rename jan1993imd2010q5_M townsendDep_mat
rename yp_housing housing
rename a006_grp housing_mat
rename yp_finDiffs financeDiffs
rename clon112 financeDiffs_0_16
rename c525 financeDiffsScore_mat
rename a053 accessToCar_parents
rename a551 crowding_birth
rename clon113 badNeigh_0_16
rename a636 neighPercept_mat
rename father_ab fatherAbsence
rename age_FA age_fatherAbsence


* Associations between demographic factors (excluding marital status; a525_grp) - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
spearman ageAt28 mother_ageAtBirth male nonWhiteEthnic livePartner mobility_mat rural rural_mat parent parity_mat, pw

matrix cor_demo = r(Rho)
matrix list cor_demo

heatplot cor_demo, values(format(%9.3f)) color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) legend(off) xlabel(, angle(45))

* Save heatmap
graph export ".\G1_Results\corr_heatplot_demoOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G1_Results\corrMatrix_demoOnly.xlsx", replace
putexcel A1=matrix(cor_demo), names

* And now for associations of demographic variables with unordered categorical variable marital status (and save to a CSV file)
capture postclose marital_corrs_demoOnly
postfile marital_corrs_demoOnly str30 variable corr ///
	using ".\G1_Results\marital_corrs_demoOnly.dta", replace

foreach var of varlist ageAt28 mother_ageAtBirth male nonWhiteEthnic livePartner mobility_mat rural rural_mat parent parity_mat {
	quietly mlogit maritalStatus_mat `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on mother's marital status is: " round((sqrt(`r2')), .001)
	
	post marital_corrs_demoOnly ("`var'") (sqrt(`r2'))
}

postclose marital_corrs_demoOnly

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G1_Results\marital_corrs_demoOnly.dta", clear
list, clean
outsheet using ".\G1_Results\marital_corrs_demoOnly.csv", comma replace
restore


** Now repeat for socioeconomic/material insecurity variables (to exclude housing status [housing and housing_mat], as are unordered categorical variables)
spearman education-townsendDep_mat financeDiffs-age_fatherAbsence, pw

matrix cor_socio = r(Rho)

* As lots of entries, the correlation coefficients are hard to read, so will drop these values and include the legend
heatplot cor_socio, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))
	
* Save heatmap
graph export ".\G1_Results\corr_heatplot_socioOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G1_Results\corrMatrix_socioOnly.xlsx", replace
putexcel A1=matrix(cor_socio), names


* And now for associations of socioeconomic variables with unordered categorical variable home ownership status (and save to a CSV file)
capture postclose housing_corrs_socioOnly
postfile housing_corrs_socioOnly str30 variable corr ///
	using ".\G1_Results\housing_corrs_socioOnly.dta", replace
	
foreach var of varlist education-townsendDep_mat housing_mat financeDiffs-age_fatherAbsence {
	quietly mlogit housing `var'
	local r2 = e(r2_p)
	
	* Using categorical notation for unordered variables
	if "`var'" == "housing_mat" {
		quietly mlogit housing i.`var'
		local r2 = e(r2_p)
	}
	
	display "Estimated correlation for " "`var'" " on home ownership status is: " round((sqrt(`r2')), .001)

	post housing_corrs_socioOnly ("`var'") (sqrt(`r2'))
}

postclose housing_corrs_socioOnly

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G1_Results\housing_corrs_socioOnly.dta", clear
list, clean
outsheet using ".\G1_Results\housing_corrs_socioOnly.csv", comma replace
restore


** And repeat for mother's home ownership status
capture postclose housingMat_corrs_socioOnly
postfile housingMat_corrs_socioOnly str30 variable corr ///
	using ".\G1_Results\housingMat_corrs_socioOnly.dta", replace
	
foreach var of varlist education-townsendDep_mat housing financeDiffs-age_fatherAbsence {
	quietly mlogit housing_mat `var'
	local r2 = e(r2_p)
	
	* Using categorical notation for unordered variables
	if "`var'" == "housing" {
		quietly mlogit housing_mat i.`var'
		local r2 = e(r2_p)
	}
	
	display "Estimated correlation for " "`var'" " on mother's home ownership status is: " round((sqrt(`r2')), .001)

	post housingMat_corrs_socioOnly ("`var'") (sqrt(`r2'))
}

postclose housingMat_corrs_socioOnly

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G1_Results\housingMat_corrs_socioOnly.dta", clear
list, clean
outsheet using ".\G1_Results\housingMat_corrs_socioOnly.csv", comma replace
restore


	
*** Finally, repeat this on all of the exposures together (excluding unordered cateogorical variables housing status and marital status)
spearman ageAt28 mother_ageAtBirth male nonWhiteEthnic livePartner mobility_mat rural rural_mat parent parity_mat education-townsendDep_mat financeDiffs-age_fatherAbsence, pw

matrix cor_all = r(Rho)

heatplot cor_all, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle(""))

* Save heatmap
graph export ".\G1_Results\corr_heatplot_all.pdf", replace

* And save in .EPS format
heatplot cor_all, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle("")) ysize(10) xsize(14)
	
graph export ".\G1_Results\corr_heatplot_all.eps", replace

* Save matrix as Excel file
putexcel set ".\G1_Results\corrMatrix_all.xlsx", replace
putexcel A1=matrix(cor_all), names

	
* And now for associations of all other exposures variables with unordered categorical variables marital status, YPs home ownership status and mother's home ownership status

* Marital status
capture postclose marital_corrs_all
postfile marital_corrs_all str30 variable corr ///
	using ".\G1_Results\marital_corrs_all.dta", replace
	
foreach var of varlist ageAt28 mother_ageAtBirth male nonWhiteEthnic livePartner mobility_mat rural rural_mat parent parity_mat education-townsendDep_mat housing housing_mat financeDiffs-age_fatherAbsence {
	quietly mlogit maritalStatus `var'
	local r2 = e(r2_p)
	
	* Using categorical notation for unordered variables
	if "`var'" == "housing" | "`var'" == "housing_mat" {
		quietly mlogit maritalStatus_mat i.`var'
		local r2 = e(r2_p)
	}
	
	display "Estimated correlation for " "`var'" " on marital status is: " round((sqrt(`r2')), .001)

	post marital_corrs_all ("`var'") (sqrt(`r2'))
}

postclose marital_corrs_all

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G1_Results\marital_corrs_all.dta", clear
list, clean
outsheet using ".\G1_Results\marital_corrs_all.csv", comma replace
restore


* Housing status
capture postclose housing_corrs_all
postfile housing_corrs_all str30 variable corr ///
	using ".\G1_Results\housing_corrs_all.dta", replace
	
foreach var of varlist ageAt28 mother_ageAtBirth male nonWhiteEthnic livePartner maritalStatus_mat mobility_mat rural rural_mat parent parity_mat education-townsendDep_mat housing_mat financeDiffs-age_fatherAbsence {
	quietly mlogit housing `var'
	
	* Using categorical notation for unordered variables
	if "`var'" == "maritalStatus_mat" | "`var'" == "housing_mat" {
		quietly mlogit housing i.`var'
		local r2 = e(r2_p)
	}
	
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on home ownership status is: " round((sqrt(`r2')), .001)


	post housing_corrs_all ("`var'") (sqrt(`r2'))
}

postclose housing_corrs_all

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G1_Results\housing_corrs_all.dta", clear
list, clean
outsheet using ".\G1_Results\housing_corrs_all.csv", comma replace
restore


* Mother's housing status
capture postclose housingMat_corrs_all
postfile housingMat_corrs_all str30 variable corr ///
	using ".\G1_Results\housingMat_corrs_all.dta", replace
	
foreach var of varlist ageAt28 mother_ageAtBirth male nonWhiteEthnic livePartner maritalStatus_mat mobility_mat rural rural_mat parent parity_mat education-townsendDep_mat housing financeDiffs-age_fatherAbsence {
	quietly mlogit housing_mat `var'
	
	* Using categorical notation for unordered variables
	if "`var'" == "maritalStatus_mat" | "`var'" == "housing" {
		quietly mlogit housing_mat i.`var'
		local r2 = e(r2_p)
	}
	
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on mother's home ownership status is: " round((sqrt(`r2')), .001)


	post housingMat_corrs_all ("`var'") (sqrt(`r2'))
}

postclose housingMat_corrs_all

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G1_Results\housingMat_corrs_all.dta", clear
list, clean
outsheet using ".\G1_Results\housingMat_corrs_all.csv", comma replace
restore



************************************************************************************
*** Next, we want to run the actual analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab YPG3000, m

** This will be quite complicated, as want to post results to file, but as exposures differ extracting the results will be variable specific. To adjust for multiple corrections will use conservative bonferroni adjustment when interpreting p-values - As 35 exposures, will a Bonferroni p-value threshold of 0.05/35 = 0.0014.
display 0.05/35

** We also want to store both estimates adjusting for age and sex, and also the interaction between sex and the exposure (as all G1s are of a similar age, will only explore interactions with sex, not age), to see whether it's moderated by sex. Again, this makes the set-up a bit more complicated.

** Create a postfile to post results to, then start the loop - Will create three postfiles; one for coefficients and CIs, another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model), and then a third for the (pseduo-)R2 values to assess model fit - NOTE: Have to store pvalues as 'double' format, else really tiny p-values get coded as '0' (as default if float format, which has minimum value of -3.40282346639e+38 [https://blog.stata.com/2012/04/02/the-penultimate-guide-to-precision/]).
capture postclose G1_belief
postfile G1_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_belief_results.dta", replace

capture postclose G1_belief_lr
postfile G1_belief_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_belief_results_lr.dta", replace
	
capture postclose G1_belief_r2
postfile G1_belief_r2 str30 exposure r2_main r2_int ///
	using ".\G1_Results\G1_belief_results_r2.dta", replace

foreach var of varlist ageAt28-age_fatherAbsence {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' and 'sex' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3000 male `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/yes)
		local outcome_level = "Yes (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// Interaction between age and sex
		mlogit YPG3000 c.male##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local sex_main = res[1,1]
		local exp_main = res[1,2]
		
		post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/not sure)
		local outcome_level = "Not sure (ref = No)"
		local exp_level = "NA"
		
		mlogit YPG3000 male `var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3000 c.male##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local sex_main = res[1,5]
		local exp_main = res[1,6]
		
		post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit YPG3000 male if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit YPG3000 male `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3000 male if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit YPG3000 male `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3000 c.male##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3000 c.male##c.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_belief_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3000 ageAt28 `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/yes)
		local outcome_level = "Yes (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local sex_main = .
		local exp_main = .
		
		post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/not sure)
		local outcome_level = "Not sure (ref = No)"
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
		local sex_main = .
		local exp_main = .
		
		post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit YPG3000 ageAt28 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit YPG3000 ageAt28 `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3000 ageAt28 if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit YPG3000 ageAt28 `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// As no interaction model for sex, will just fill with missing values
		local lr_p_int = .
		local r2_int = .
		
		post G1_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_belief_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "livePartner" | "`var'" == "rural" | "`var'" == "rural_mat" | "`var'" == "parent" | "`var'" == "employed" | "`var'" == "highSocClass_mat" |  "`var'" == "highSocClass_pat" |  "`var'" == "lowSocClass_0_16" | "`var'" == "income_parents" | "`var'" == "financeDiffs" | "`var'" == "financeDiffs_0_16" | "`var'" == " financeDiffsScore_mat" | "`var'" == "accessToCar_parents" | "`var'" == "badNeigh_0_16" | "`var'" == "neighPercept_mat" | "`var'" == "fatherAbsence" {
		
		mlogit YPG3000 ageAt28 male `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/yes)
		local outcome_level = "Yes (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,3]
		local lci = res[5,3]
		local uci = res[6,3]
		local p = res[4,3]
		
		// Now for interaction model
		mlogit YPG3000 ageAt28 c.male##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,4]
		local lci_int = res[5,4]
		local uci_int = res[6,4]
		local p_int = res[4,4]
		local sex_main = res[1,2]
		local exp_main = res[1,3]
		
		post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/not sure)
		mlogit YPG3000 ageAt28 male `var', baseoutcome(3) rrr
		
		local outcome_level = "Not sure (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
				
		// Now for interaction model
		mlogit YPG3000 ageAt28 c.male##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,9]
		local lci_int = res[5,9]
		local uci_int = res[6,9]
		local p_int = res[4,9]
		local sex_main = res[1,7]
		local exp_main = res[1,8]
		
		post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit YPG3000 ageAt28 male if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit YPG3000 ageAt28 male `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3000 ageAt28 male if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit YPG3000 ageAt28 male `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3000 ageAt28 c.male##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3000 ageAt28 c.male##c.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_belief_r2 ("`exp'") (`r2_main') (`r2_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus_mat" {
				local exp_level = "Married (ref = Never married)"
			}
			else if "`var'" == "parity_mat" {
				local exp_level = "1 (ref = 0)"
			}
			else if "`var'" == "age_fatherAbsence" {
				local exp_level = "FA age < 5 (ref = no FA)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,16]
			local lci_int = res[5,16]
			local uci_int = res[6,16]
			local p_int = res[4,16]
			local sex_main = res[1,11]
			local exp_main = res[1,13]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus_mat" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity_mat" {
				local exp_level = "2 or more (ref = 0)"
			}
			else if "`var'" == "age_fatherAbsence" {
				local exp_level = "FA age 5+ (ref = no FA)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local sex_main = res[1,11]
			local exp_main = res[1,14]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local sex_main = res[1,13]
			local exp_main = res[1,15]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local sex_main = res[1,1]
			local exp_main = res[1,5]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local sex_main = res[1,13]
			local exp_main = res[1,16]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local sex_main = res[1,13]
			local exp_main = res[1,17]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£500-£999 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£1000-£1499 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local sex_main = res[1,15]
			local exp_main = res[1,18]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£1500-£1999 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local sex_main = res[1,15]
			local exp_main = res[1,19]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£2000 and above (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local sex_main = res[1,1]
			local exp_main = res[1,7]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,15]
			local exp_main = res[1,20]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,17]
			local exp_main = res[1,19]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,17]
			local exp_main = res[1,20]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local sex_main = res[1,17]
			local exp_main = res[1,21]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local sex_main = res[1,2]
			local exp_main = res[1,7]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local sex_main = res[1,17]
			local exp_main = res[1,22]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local sex_main = res[1,2]
			local exp_main = res[1,8]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/not sure)
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local sex_main = res[1,17]
			local exp_main = res[1,23]
		
			post G1_belief ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures and R2 values
		mlogit YPG3000 ageAt28 male if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3000 ageAt28 male if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_belief_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose G1_belief
postclose G1_belief_lr
postclose G1_belief_r2
	

	
************************************************************************************
*** Now to the next RSBB outcome: Religious affiliation

*** As this is an unordered categorical variable, will again use multinomial logistic model (with 'none' as reference)
tab YPG3040_grp

* Do need to re-order the categories so can simply copy and paste the script from above without having to faff around with editing the cells to take the statistics from
recode YPG3040_grp (1 = 3) (2 = 1) (3 = 2)
label define relig_lb 1 "Christian" 2 "Other" 3 "None", modify
numlabel relig_lb, add
tab YPG3040_grp


*** Now run the loop to save all the results
capture postclose G1_relig
postfile G1_relig str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_relig_results.dta", replace

capture postclose G1_relig_lr
postfile G1_relig_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_relig_results_lr.dta", replace
	
capture postclose G1_relig_r2
postfile G1_relig_r2 str30 exposure r2_main r2_int ///
	using ".\G1_Results\G1_relig_results_r2.dta", replace

foreach var of varlist ageAt28-age_fatherAbsence  {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' and 'sex' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3040_grp male `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/Chrstian)
		local outcome_level = "Christian (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// Interaction between age and sex
		mlogit YPG3040_grp c.male##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local sex_main = res[1,1]
		local exp_main = res[1,2]
		
		post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/Other)
		local outcome_level = "Other (ref = None)"
		local exp_level = "NA"
		
		mlogit YPG3040_grp male `var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3040_grp c.male##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local sex_main = res[1,5]
		local exp_main = res[1,6]
		
		post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit YPG3040_grp male if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit YPG3040_grp male `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3040_grp male if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit YPG3040_grp male `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3040_grp c.male##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3040_grp c.male##c.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_relig_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
		// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3040_grp ageAt28 `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/Christian)
		local outcome_level = "Christian (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local sex_main = .
		local exp_main = .
		
		post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/other)
		local outcome_level = "Other (ref = None)"
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
		local sex_main = .
		local exp_main = .
		
		post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit YPG3040_grp ageAt28 if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit YPG3040_grp ageAt28 `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3040_grp ageAt28 if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit YPG3040_grp ageAt28 `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// As no interaction model for sex, will just fill with missing values
		local lr_p_int = .
		local r2_int = .
		
		post G1_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_relig_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "livePartner" | "`var'" == "rural" | "`var'" == "rural_mat" | "`var'" == "parent" | "`var'" == "employed" | "`var'" == "highSocClass_mat" |  "`var'" == "highSocClass_pat" |  "`var'" == "lowSocClass_0_16" | "`var'" == "income_parents" | "`var'" == "financeDiffs" | "`var'" == "financeDiffs_0_16" | "`var'" == " financeDiffsScore_mat" | "`var'" == "accessToCar_parents" | "`var'" == "badNeigh_0_16" | "`var'" == "neighPercept_mat" | "`var'" == "fatherAbsence" {
		
		mlogit YPG3040_grp ageAt28 male `var', baseoutcome(3) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/Chrstian)
		local outcome_level = "Christian (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,3]
		local lci = res[5,3]
		local uci = res[6,3]
		local p = res[4,3]
		
		// Now for interaction model
		mlogit YPG3040_grp ageAt28 c.male##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,4]
		local lci_int = res[5,4]
		local uci_int = res[6,4]
		local p_int = res[4,4]
		local sex_main = res[1,2]
		local exp_main = res[1,3]
		
		post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/Other)
		mlogit YPG3040_grp ageAt28 male `var', baseoutcome(3) rrr
		
		local outcome_level = "Other (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
				
		// Now for interaction model
		mlogit YPG3040_grp ageAt28 c.male##c.`var', baseoutcome(3) rrr
		
		matrix res = r(table)
		local coef_int = res[1,9]
		local lci_int = res[5,9]
		local uci_int = res[6,9]
		local p_int = res[4,9]
		local sex_main = res[1,7]
		local exp_main = res[1,8]
		
		post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit YPG3040_grp ageAt28 male if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit YPG3040_grp ageAt28 male `var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3040_grp ageAt28 male if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit YPG3040_grp ageAt28 male `var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3040_grp ageAt28 c.male##c.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3040_grp ageAt28 c.male##c.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_relig_r2 ("`exp'") (`r2_main') (`r2_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus_mat" {
				local exp_level = "Married (ref = Never married)"
			}
			else if "`var'" == "parity_mat" {
				local exp_level = "1 (ref = 0)"
			}
			else if "`var'" == "age_fatherAbsence" {
				local exp_level = "FA age < 5 (ref = no FA)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,16]
			local lci_int = res[5,16]
			local uci_int = res[6,16]
			local p_int = res[4,16]
			local sex_main = res[1,11]
			local exp_main = res[1,13]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus_mat" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity_mat" {
				local exp_level = "2 or more (ref = 0)"
			}
			else if "`var'" == "age_fatherAbsence" {
				local exp_level = "FA age 5+ (ref = no FA)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local sex_main = res[1,11]
			local exp_main = res[1,14]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local sex_main = res[1,13]
			local exp_main = res[1,15]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local sex_main = res[1,13]
			local exp_main = res[1,16]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local sex_main = res[1,13]
			local exp_main = res[1,17]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£500-£999 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£1000-£1499 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local sex_main = res[1,15]
			local exp_main = res[1,18]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£1500-£1999 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local sex_main = res[1,15]
			local exp_main = res[1,19]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£2000 and above (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local sex_main = res[1,2]
			local exp_main = res[1,7]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,15]
			local exp_main = res[1,20]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,17]
			local exp_main = res[1,19]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,17]
			local exp_main = res[1,20]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local sex_main = res[1,17]
			local exp_main = res[1,21]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local sex_main = res[1,2]
			local exp_main = res[1,7]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local sex_main = res[1,17]
			local exp_main = res[1,22]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local sex_main = res[1,2]
			local exp_main = res[1,8]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/Other)
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local sex_main = res[1,17]
			local exp_main = res[1,23]
		
			post G1_relig ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures and R2 values
		mlogit YPG3040_grp ageAt28 male if `var' != ., baseoutcome(3) rrr
		est store base
		mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3040_grp ageAt28 male if `var' != ., baseoutcome(3) rrr
		local r2_base = e(r2_p)
		mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_relig_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose G1_relig
postclose G1_relig_lr
postclose G1_relig_r2



************************************************************************************
*** Now to the next RSBB outcome: Attendance at church/place of worship

*** As this is an ordered categorical variable, originally planned to use ordinal regression models. However, as the G0 mother and G0 partner/father data violated the proportional odds assumption (as does the G1 data; see below), will just use multinomial models regardless so that results are comparable.

** For consistency with previous models, will recode so that higher values indicate greater RSBB - Due to small sample sizes of those attending at least once a month and different way question was asked to G1s (it includes the option 'occasionally' between 'at least once a year' and 'not at all'), will code into categories of: min once a month, min once a year, occasionally, and not at all (in contrast to G0 data which is: min once a week, min once a month, min once a year and not at all).
tab YPG3080

recode YPG3080 (1 2 = 3) (3 = 2) (4 = 1) (5 = 0), gen(YPG3080_rev)
label define attend_rev_lb 0 "Not at all" 1 "Occasionally" 2 "Min once a year" 3 "Min once a month"
numlabel attend_rev_lb, add
label values YPG3080_rev attend_rev_lb
tab YPG3080_rev


** Quick test of whether proportional odds assumption been violated in most basic model (with just sex and age at 28). Ah, it has been violated. 
ologit YPG3080_rev ageAt28 male, or
brant, detail

ologit YPG3080_rev ageAt28 male i.IMD, or
brant, detail


** So instead will just run multinomial models with 'not at all' as the baseline/reference category (this also means all outcomes are using the same model, making them easier to compare).


*** Now run the loop to save all the results
capture postclose G1_attend
postfile G1_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_attend_results.dta", replace

capture postclose G1_attend_lr
postfile G1_attend_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_attend_results_lr.dta", replace
	
capture postclose G1_attend_r2
postfile G1_attend_r2 str30 exposure r2_main r2_int ///
	using ".\G1_Results\G1_attend_results_r2.dta", replace

foreach var of varlist ageAt28-age_fatherAbsence {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' and 'sex' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/Occasionally)
		local outcome_level = "Occasionally (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local sex_main = res[1,5]
		local exp_main = res[1,6]
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/once year)
		local outcome_level = "Min once year (ref = Not at all)"
		local exp_level = "NA"
		
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
		
		// Interaction between age and sex
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local sex_main = res[1,9]
		local exp_main = res[1,10]
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// Now onto the next reference category (3/once month)
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
		
		// Interaction between age and sex
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local sex_main = res[1,13]
		local exp_main = res[1,14]
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and store R2 values
		mlogit YPG3080_rev male if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3080_rev male if `var' != ., baseoutcome(0) rrr
		local r2_base = e(r2_p)
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_attend_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3080_rev ageAt28 `var', baseoutcome(0) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/Occasionally)
		local outcome_level = "Occasionally (ref = Not at all)"
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
		local sex_main = .
		local exp_main = .
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/once year)
		local outcome_level = "Min once year (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local sex_main = .
		local exp_main = .
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/once month)
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local sex_main = .
		local exp_main = .
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit YPG3080_rev ageAt28 if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit YPG3080_rev ageAt28 `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3080_rev ageAt28 if `var' != ., baseoutcome(0) rrr
		local r2_base = e(r2_p)
		mlogit YPG3080_rev ageAt28 `var', baseoutcome(0) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// As no interaction model for sex, will just fill with missing values
		local lr_p_int = .
		loca r2_int = .
		
		post G1_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_attend_r2 ("`exp'") (`r2_main') (`r2_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "livePartner" | "`var'" == "rural" | "`var'" == "rural_mat" | "`var'" == "parent" | "`var'" == "employed" | "`var'" == "highSocClass_mat" |  "`var'" == "highSocClass_pat" |  "`var'" == "lowSocClass_0_16" | "`var'" == "income_parents" | "`var'" == "financeDiffs" | "`var'" == "financeDiffs_0_16" | "`var'" == " financeDiffsScore_mat" | "`var'" == "accessToCar_parents" | "`var'" == "badNeigh_0_16" | "`var'" == "neighPercept_mat" | "`var'" == "fatherAbsence" {
		
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr
		
		local n = e(N)
		
		// Start with the first reference category (1/Occasionally)
		local outcome_level = "Occasionally (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
		
		// Now for interaction model
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,9]
		local lci_int = res[5,9]
		local uci_int = res[6,9]
		local p_int = res[4,9]
		local sex_main = res[1,7]
		local exp_main = res[1,8]
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (2/once year)
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr
		
		local outcome_level = "Min once year (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,14]
		local lci_int = res[5,14]
		local uci_int = res[6,14]
		local p_int = res[4,14]
		local sex_main = res[1,12]
		local exp_main = res[1,13]
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/once month)
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr
		
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,15]
		local lci = res[5,15]
		local uci = res[6,15]
		local p = res[4,15]
				
		// Now for interaction model
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr
		
		matrix res = r(table)
		local coef_int = res[1,19]
		local lci_int = res[5,19]
		local uci_int = res[6,19]
		local p_int = res[4,19]
		local sex_main = res[1,17]
		local exp_main = res[1,18]
		
		post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests and R2 values
		mlogit YPG3080_rev ageAt28 male if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3080_rev ageAt28 male if `var' != ., baseoutcome(0) rrr
		local r2_base = e(r2_p)
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_attend_r2 ("`exp'") (`r2_main') (`r2_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus_mat" {
				local exp_level = "Married (ref = Never married)"
			}
			else if "`var'" == "parity_mat" {
				local exp_level = "1 (ref = 0)"
			}
			else if "`var'" == "age_fatherAbsence" {
				local exp_level = "FA age < 5 (ref = no FA)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,16]
			local lci_int = res[5,16]
			local uci_int = res[6,16]
			local p_int = res[4,16]
			local sex_main = res[1,11]
			local exp_main = res[1,13]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,20]
			local exp_main = res[1,22]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local sex_main = res[1,29]
			local exp_main = res[1,31]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus_mat" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity_mat" {
				local exp_level = "2 or more (ref = 0)"
			}
			else if "`var'" == "age_fatherAbsence" {
				local exp_level = "FA age 5+ (ref = no FA)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local sex_main = res[1,11]
			local exp_main = res[1,14]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,20]
			local exp_main = res[1,23]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,29]
			local exp_main = res[1,32]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local sex_main = res[1,13]
			local exp_main = res[1,15]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local sex_main = res[1,24]
			local exp_main = res[1,26]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,35]
			local exp_main = res[1,37]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local sex_main = res[1,13]
			local exp_main = res[1,16]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local sex_main = res[1,24]
			local exp_main = res[1,27]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,35]
			local exp_main = res[1,38]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "housing_mat" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding_birth" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local sex_main = res[1,13]
			local exp_main = res[1,17]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local sex_main = res[1,24]
			local exp_main = res[1,28]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,35]
			local exp_main = res[1,39]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "Vocational (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "2 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£500-£999 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,28]
			local exp_main = res[1,30]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (2/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,48]
			local lci_int = res[5,48]
			local uci_int = res[6,48]
			local p_int = res[4,48]
			local sex_main = res[1,41]
			local exp_main = res[1,43]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "O-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "3 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£1000-£1499 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local sex_main = res[1,15]
			local exp_main = res[1,18]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,36]
			local lci_int = res[5,36]
			local uci_int = res[6,36]
			local p_int = res[4,36]
			local sex_main = res[1,28]
			local exp_main = res[1,31]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,49]
			local lci_int = res[5,49]
			local uci_int = res[6,49]
			local p_int = res[4,49]
			local sex_main = res[1,41]
			local exp_main = res[1,44]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "A-level (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "4 (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£1500-£1999 (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local sex_main = res[1,28]
			local exp_main = res[1,32]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,50]
			local lci_int = res[5,50]
			local uci_int = res[6,50]
			local p_int = res[4,50]
			local sex_main = res[1,41]
			local exp_main = res[1,45]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education_mat" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			if "`var'" == "education_pat" {
				local exp_level = "Degree (ref = CSE/None)"
			}
			else if "`var'" == "IMD" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "IMD_mat" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "townsendDep_mat" {
				local exp_level = "5/Most dep. (ref = 1/Least dep.)"
			}
			else if "`var'" == "income" {
				local exp_level = "£2000 and above (ref = £0-£499)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,15]
			local exp_main = res[1,20]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local sex_main = res[1,28]
			local exp_main = res[1,33]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local sex_main = res[1,41]
			local exp_main = res[1,46]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,17]
			local exp_main = res[1,19]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local sex_main = res[1,32]
			local exp_main = res[1,34]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local sex_main = res[1,47]
			local exp_main = res[1,49]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,17]
			local exp_main = res[1,20]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,32]
			local exp_main = res[1,35]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,56]
			local lci_int = res[5,56]
			local uci_int = res[6,56]
			local p_int = res[4,56]
			local sex_main = res[1,47]
			local exp_main = res[1,50]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local sex_main = res[1,17]
			local exp_main = res[1,21]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,32]
			local exp_main = res[1,36]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,57]
			local lci_int = res[5,57]
			local uci_int = res[6,57]
			local p_int = res[4,57]
			local sex_main = res[1,47]
			local exp_main = res[1,51]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local sex_main = res[1,17]
			local exp_main = res[1,22]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,32]
			local exp_main = res[1,37]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,58]
			local lci_int = res[5,58]
			local uci_int = res[6,58]
			local p_int = res[4,58]
			local sex_main = res[1,47]
			local exp_main = res[1,52]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility_mat" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local sex_main = res[1,17]
			local exp_main = res[1,23]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (2/once year)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local sex_main = res[1,32]
			local exp_main = res[1,38]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/once month)
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,35]
			local lci = res[5,35]
			local uci = res[6,35]
			local p = res[4,35]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		
			matrix res = r(table)
			local coef_int = res[1,59]
			local lci_int = res[5,59]
			local uci_int = res[6,59]
			local p_int = res[4,59]
			local sex_main = res[1,47]
			local exp_main = res[1,53]
		
			post G1_attend ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures and R2 values
		mlogit YPG3080_rev ageAt28 male if `var' != ., baseoutcome(0) rrr
		est store base
		mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		mlogit YPG3080_rev ageAt28 male if `var' != ., baseoutcome(0) rrr
		local r2_base = e(r2_p)
		mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr
		local r2_main = e(r2_p) - `r2_base'
		
		// And the interaction model
		mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr
		local r2_int = e(r2_p) - (`r2_main' + `r2_base')
		
		post G1_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		post G1_attend_r2 ("`exp'") (`r2_main') (`r2_int')
				
	}
		
}

postclose G1_attend
postclose G1_attend_lr
postclose G1_attend_r2


* And close the log file
log close

* Save this file if we want to use it later
save ".\G1_Results\G1_PredictorsOfRSBB_B3911_postAnalysis.dta", replace




***********************************************************************************
***********************************************************************************
***** Next step is to make some nice plots of these results

**** Start with outcome of belief in God/divine power 
use ".\G1_Results\G1_belief_results_lr.dta", clear

* Put 'male' as the first variable, so that interaction plots exclude this variable without needing to re-write the code.
replace exposure = "aamale" if exposure == "male"
sort exposure in 1/3
replace exposure = "male" if exposure == "aamale"

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 35 exposures, will do 0.05/35)
local bon_thresh = -log10(0.05/35)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious belief - Main effect") ///
	name(belief_main, replace)
	
graph export ".\G1_Results\belief_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/34)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)35, valuelabel labsize(tiny) angle(0)) ///
	title("Religious belief - Sex interaction") ///
	name(belief_int, replace)
	
graph export ".\G1_Results\belief_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/35)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(tiny) angle(0)) ///
	title("Religious belief") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(belief_both, replace)

graph export ".\G1_Results\belief_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\G1_Results\belief_pvalues.dta", replace


**** Now read in the next outcome - religious affiliation
use ".\G1_Results\G1_relig_results_lr.dta", clear

* Put 'male' as the first variable, so that interaction plots exclude this variable without needing to re-write the code.
replace exposure = "aamale" if exposure == "male"
sort exposure in 1/3
replace exposure = "male" if exposure == "aamale"

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 35 exposures, will do 0.05/35)
local bon_thresh = -log10(0.05/35)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(tiny) angle(0)) ///
	title("Religious affiliation - Main effect") ///
	name(relig_main, replace)
	
graph export ".\G1_Results\relig_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/34)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)35, valuelabel labsize(tiny) angle(0)) ///
	title("Religious affiliation - Sex interaction") ///
	name(relig_int, replace)
	
graph export ".\G1_Results\relig_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/35)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(tiny) angle(0)) ///
	title("Religious affiliation") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(relig_both, replace)

graph export ".\G1_Results\relig_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\G1_Results\relig_pvalues.dta", replace


**** Now read in the next outcome - church attendance
use ".\G1_Results\G1_attend_results_lr.dta", clear

* Put 'male' as the first variable, so that interaction plots exclude this variable without needing to re-write the code.
replace exposure = "aamale" if exposure == "male"
sort exposure in 1/3
replace exposure = "male" if exposure == "aamale"

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 35 exposures, will do 0.05/35)
local bon_thresh = -log10(0.05/35)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(tiny) angle(0)) ///
	title("Religious attendance - Main effect") ///
	name(attend_main, replace)
	
graph export ".\G1_Results\attend_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/34)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)35, valuelabel labsize(tiny) angle(0)) ///
	title("Religious attendance - Sex interaction") ///
	name(attend_int, replace)
	
graph export ".\G1_Results\attend_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/35)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(tiny) angle(0)) ///
	title("Religious attendance") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(attend_both, replace)

graph export ".\G1_Results\attend_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'church attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\G1_Results\attend_pvalues.dta", replace


*** Combine all these datasets together
use ".\G1_Results\belief_pvalues.dta", clear
append using ".\G1_Results\relig_pvalues.dta"
append using ".\G1_Results\attend_pvalues.dta"


** Now look at combined results

* Belief/religion/church vars main effects
local bon_thresh = -log10(0.05/35)
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
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(belRelCh_main, replace)

graph export ".\G1_Results\beliefReligAttend_mainEffects_pvalues.pdf", replace

* Belief/religion/church vars interaction effects
local bon_thresh = -log10(0.05/34)
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
	ylabel(2(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Sex interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(belRelCh_int, replace)

graph export ".\G1_Results\beliefReligAttend_sexInt_pvalues.pdf", replace


** Combine all these graphs together
graph combine belRelCh_main belRelCh_int, ysize(3) xsize(6)

graph export ".\G1_Results\allData_pvalues.pdf", replace

* And convert to EPS format for Wellcome Open Research formatting (also have to enlarge size, else resolution is terrible)
graph combine belRelCh_main belRelCh_int, ysize(10) xsize(20)
graph export ".\G1_Results\allData_pvalues.eps", replace

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

outsheet exp_num-lr_p_intAttend using ".\G1_Results\pvalue_results.csv", comma replace


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
use ".\G1_Results\G1_belief_results_r2.dta", clear

* Put 'male' as the first variable, so that interaction plots exclude this variable without needing to re-write the code.
replace exposure = "aamale" if exposure == "male"
sort exposure in 1/3
replace exposure = "male" if exposure == "aamale"

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
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious belief - Main effect") ///
	name(belief_main, replace)
	
graph export ".\G1_Results\belief_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious belief - Sex interaction") ///
	name(belief_int, replace)
	
graph export ".\G1_Results\belief_sexInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious belief") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(belief_both, replace)

graph export ".\G1_Results\belief_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\G1_Results\belief_r2.dta", replace


**** Now read in the next outcome - religious affiliation
use ".\G1_Results\G1_relig_results_r2.dta", clear

* Put 'male' as the first variable, so that interaction plots exclude this variable without needing to re-write the code.
replace exposure = "aamale" if exposure == "male"
sort exposure in 1/3
replace exposure = "male" if exposure == "aamale"

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
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Main effect") ///
	name(relig_main, replace)
	
graph export ".\G1_Results\relig_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation - Sex interaction") ///
	name(relig_int, replace)
	
graph export ".\G1_Results\relig_sexInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious affiliation") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(relig_both, replace)

graph export ".\G1_Results\relig_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\G1_Results\relig_r2.dta", replace


**** Now read in the next outcome - church attendance
use ".\G1_Results\G1_attend_results_r2.dta", clear

* Put 'male' as the first variable, so that interaction plots exclude this variable without needing to re-write the code.
replace exposure = "aamale" if exposure == "male"
sort exposure in 1/3
replace exposure = "male" if exposure == "aamale"

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
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious attendance - Main effect") ///
	name(attend_main, replace)
	
graph export ".\G1_Results\attend_mainEffect_r2.pdf", replace
	
* And repeat for interaction effect (exclude 'age' here, as can't interact with itself!)
twoway (scatter exp_num r2_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious attendance - Sex interaction") ///
	name(attend_int, replace)
	
graph export ".\G1_Results\attend_sexInteraction_r2.pdf", replace
	
** Combine these results on the same plot
twoway (scatter exp_num r2_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Religious attendance") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(attend_both, replace)

graph export ".\G1_Results\attend_mainAndInt_r2.pdf", replace

graph close _all
	
** Add 'church attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\G1_Results\attend_r2.dta", replace


*** Combine all these datasets together
use ".\G1_Results\belief_r2.dta", clear
append using ".\G1_Results\relig_r2.dta"
append using ".\G1_Results\attend_r2.dta"


** Now look at combined results

* Main effects
twoway (scatter exp_num r2_main if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(r2_main, replace)

graph export ".\G1_Results\beliefReligAttend_mainEffects_r2.pdf", replace

* Interaction effects
twoway (scatter exp_num r2_int if outcome == "Belief" & exp_num != 1, ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Religious affil." & exp_num != 1, ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Church attendance" & exp_num != 1, ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(2(1)35, valuelabel labsize(vsmall) angle(0)) ///
	title("Sex interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(r2_int, replace)

graph export ".\G1_Results\beliefReligAttend_sexInt_r2.pdf", replace


** Combine all these graphs together
graph combine r2_main r2_int, ysize(3) xsize(6)

graph export ".\G1_Results\allData_r2.pdf", replace

* And convert to EPS format
graph combine r2_main r2_int, ysize(10) xsize(20)
graph export ".\G1_Results\allData_r2.eps", replace

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

outsheet exp_num-r2_intAttend using ".\G1_Results\r2_results.csv", comma replace


********************************************************************************
*** Next, want to plot some of the actual coefficient results

** Will read the datasets in then combine together into one single dataset
use ".\G1_Results\G1_belief_results.dta", clear
gen outcome = "Belief"

append using ".\G1_Results\G1_relig_results.dta"
replace outcome = "Relig" if outcome == ""
tab outcome, m

append using ".\G1_Results\G1_attend_results.dta"
replace outcome = "Attend" if outcome == ""
tab outcome, m


** Save these results as CSV files to add to the SI

* Save each result in turn
format coef lci uci coef_int lci_int uci_int %9.3f
format p p_int %9.4f

outsheet exposure-p_int using ".\G1_Results\belief_coefs.csv" if outcome == "Belief", comma replace

outsheet exposure-p_int using ".\G1_Results\relig_coefs.csv" if outcome == "Relig", comma replace

outsheet exposure-p_int using ".\G1_Results\attend_coefs.csv" if outcome == "Attend", comma replace


* Convert format back to default format (so axis on plots display correctly)
format coef lci uci coef_int lci_int uci_int %9.0g
format p p_int %10.0g


** First, make a plot for the age results
capture drop level_num
gen level_num = 0
replace level_num = 1 if outcome_level == "Not sure (ref = No)"
replace level_num = 3 if outcome_level == "Christian (ref = None)"
replace level_num = 4 if outcome_level == "Other (ref = None)"
replace level_num = 6 if outcome_level == "Occasionally (ref = Not at all"
replace level_num = 7 if outcome_level == "Min once year (ref = Not at al"
replace level_num = 8 if outcome_level == "Min once month (ref = Not at a"

label define level_lb 0 "Belief - Yes (ref = No)" 1 "Belief - Not sure (ref = No)" 3 "Affiliation - Christian (ref = None)" 4 "Affiliation - Other (ref = None)" 6 "Attendance - Occasionally (ref = Not at all)" 7 "Attendance - Min 1/year (ref = Not at all)" 8 "Attendance - Min 1/month (ref = Not at all)", replace
label values level_num level_lb
tab level_num

twoway (scatter level_num coef if outcome == "Belief" & exposure == "ageAt28", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "ageAt28", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "ageAt28", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "ageAt28", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "ageAt28", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "ageAt28", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Age and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(age_cat, replace)
		
graph export ".\G1_Results\ageResults.pdf", replace


** Create plot for ethnicity (ref = white)

* Min and max x-axis values
sum lci uci if exposure == "nonWhiteEthnic" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "nonWhiteEthnic", ///
			col(black) msize(small) msym(D)) ///
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
		xlabel(0.7 1 2 3 5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(ethnic_cat, replace)
		
graph export ".\G1_Results\ethnicityResults.pdf", replace


** Create plot for sex (ref = female)

* Min and max x-axis values
sum lci uci if exposure == "male" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "male", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "male", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "male", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "male", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "male", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "male", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = female)") ///
		title("Male and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.6 0.7 0.8 0.9 1 1.1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(sex_cat, replace)
		
graph export ".\G1_Results\sexResults.pdf", replace


** Create plot for living with a partner (ref = No)

* Min and max x-axis values
sum lci uci if exposure == "livePartner" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "livePartner", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "livePartner", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "livePartner", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "livePartner", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "livePartner", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "livePartner", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = female)") ///
		title("Living with a partner and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.5 0.6 0.8 1 1.2 1.4, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(partner_cat, replace)
		
graph export ".\G1_Results\livePartnerResults.pdf", replace


** Create plot for mother's marital status (ref = never married)

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
		title("Mother's marital status and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.5 0.7 1 2 5 10, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Married" 7 "Widowed/Divorced/Separated")) ///
		name(marital_cat, replace)
		
graph export ".\G1_Results\maritalStatusResults.pdf", replace


** Create plot for maternal education (ref = CSE/None)

* As four exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "education_mat" & exp_level == "Vocational (ref = CSE/None)"
replace level_split = level_num - 0.1 if exposure == "education_mat" & exp_level == "O-level (ref = CSE/None)"
replace level_split = level_num + 0.1 if exposure == "education_mat" & exp_level == "A-level (ref = CSE/None)"
replace level_split = level_num + 0.3 if exposure == "education_mat" & exp_level == "Degree (ref = CSE/None)"
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
		title("Maternal education and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.7 1 1.5 2 3 5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Vocational" 7 "O-levels" 13 "A-levels" ///
			19 "Degree") rows(1)) ///
		name(matedu_cat, replace)
		
graph export ".\G1_Results\matEduResults.pdf", replace


** Create plot for YP education (ref = GCSE/None)

* As four exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "education" & exp_level == "Vocational (ref = GCSE/None)"
replace level_split = level_num if exposure == "education" & exp_level == "AS/A level (ref = GCSE/None)"
replace level_split = level_num + 0.3 if exposure == "education" & exp_level == "Degree (ref = GCSE/None)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci uci if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = GCSE/None)") ///
		title("Education and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.7 1 1.5 2 3 5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Vocational" 7 "AS/A levels" 13 "Degree") rows(1)) ///
		name(edu_cat, replace)
		
graph export ".\G1_Results\eduResults.pdf", replace

* And convert to EPS format
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = GCSE/None)") ///
		title("Education and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.7 1 1.5 2 3 5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Vocational" 7 "AS/A levels" 13 "Degree") rows(1)) ///
		name(edu_cat, replace) ysize(10) xsize(14)
		
graph export ".\G1_Results\eduResults.eps", replace
	

** Create plot for being a parent (ref = not a parent)
sum lci uci if exposure == "parent" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"parent", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"parent", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"parent", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"parent", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"parent", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"parent", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = Not a parent)") ///
		title("Parental status and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.5 0.7 1 1.5 2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(parent_cat, replace)
		
graph export ".\G1_Results\parentResults.pdf", replace


** Create plot for father absence (ref = no father absence)
sum lci uci if exposure == "fatherAbsence" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"fatherAbsence", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"fatherAbsence", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"fatherAbsence", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"fatherAbsence", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"fatherAbsence", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"fatherAbsence", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = No FA)") ///
		title("Father absence and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.4 0.6 0.8 1 1.2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(FA_cat, replace)
		
graph export ".\G1_Results\fatherAbResults.pdf", replace


** Create plot for maternal age
sum lci uci if exposure == "mother_ageAtBirth" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"mother_ageAtBirth", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"mother_ageAtBirth", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"mother_ageAtBirth", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"mother_ageAtBirth", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"mother_ageAtBirth", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"mother_ageAtBirth", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase)") ///
		title("Maternal age at birth and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.98 1 1.02 1.04 1.06, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(mumAge_cat, replace)
		
graph export ".\G1_Results\mumAgeResults.pdf", replace


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
		xlabel(0.5 0.7 1 1.5 2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "2" 7 "3" 13 "4" 19 "5/Most dep.") rows(1)) ///
		name(imd_cat, replace)
		
graph export ".\G1_Results\imdResults.pdf", replace


** Create plot for income (ref = £0-£499)

* As four exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "income" & exp_level == "£500-£999 (ref = £0-£499)"
replace level_split = level_num - 0.1 if exposure == "income" & exp_level == "£1000-£1499 (ref = £0-£499)"
replace level_split = level_num + 0.1 if exposure == "income" & exp_level == "£1500-£1999 (ref = £0-£499)"
replace level_split = level_num + 0.3 if exposure == "income" & exp_level == "£2000 and above (ref = £0-£499)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci uci if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"£500-£999 (ref = £0-£499)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"£500-£999 (ref = £0-£499)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"£500-£999 (ref = £0-£499)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"£500-£999 (ref = £0-£499)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"£500-£999 (ref = £0-£499)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"£500-£999 (ref = £0-£499)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"£2000 and above (ref = £0-£499)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"£2000 and above (ref = £0-£499)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"£2000 and above (ref = £0-£499)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"£2000 and above (ref = £0-£499)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"£2000 and above (ref = £0-£499)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"£2000 and above (ref = £0-£499)", horizontal col(green)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = £0-£499)") ///
		title("Take-home income and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.2 0.3 0.5 0.7 1 1.5 2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "£500-£999" 7 "£1000-£1499" 13 "£1500-£1999" 19 "£2000 and above") rows(2)) ///
		name(income_cat, replace)
		
graph export ".\G1_Results\incomeResults.pdf", replace

* And convert to EPS format
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"£500-£999 (ref = £0-£499)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"£500-£999 (ref = £0-£499)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"£500-£999 (ref = £0-£499)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"£500-£999 (ref = £0-£499)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"£500-£999 (ref = £0-£499)", col(black) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"£500-£999 (ref = £0-£499)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", col(red) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"£1000-£1499 (ref = £0-£499)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", col(blue) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"£1500-£1999 (ref = £0-£499)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"£2000 and above (ref = £0-£499)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"£2000 and above (ref = £0-£499)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"£2000 and above (ref = £0-£499)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"£2000 and above (ref = £0-£499)", horizontal col(green)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"£2000 and above (ref = £0-£499)", col(green) msize(vsmall) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"£2000 and above (ref = £0-£499)", horizontal col(green)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = £0-£499)") ///
		title("Take-home income and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.2 0.3 0.5 0.7 1 1.5 2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "£500-£999" 7 "£1000-£1499" 13 "£1500-£1999" 19 "£2000 and above") rows(2)) ///
		name(income_cat, replace) ysize(10) xsize(14)
		
graph export ".\G1_Results\incomeResults.eps", replace


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
		xlabel(0.35 0.5 0.7 1 2 4, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "Rented" 7 "Council/HA" 13 "Other") rows(1)) ///
		name(housing_cat, replace)
		
graph export ".\G1_Results\housingResults.pdf", replace


** Create plot for crowding at birth (ref = <= 0.5)

* As three exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.25 if exposure == "crowding_birth" & exp_level == "> 0.5 to 0.75 (ref = <= 0.5)"
replace level_split = level_num - 0 if exposure == "crowding_birth" & exp_level == "> 0.75 to 1 (ref = <= 0.5)"
replace level_split = level_num + 0.25 if exposure == "crowding_birth" & exp_level == "> 1 (ref = <= 0.5)"
label values level_split level_lb
tab level_split

* Min and max x-axis values
sum lci uci if level_split < . & outcome_level != "NA"

* Now make the graph
twoway (scatter level_split coef if outcome == "Belief" & exp_level == ///
			"> 0.5 to 0.75 (ref = <= 0.5)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"> 0.5 to 0.75 (ref = <= 0.5)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"> 0.5 to 0.75 (ref = <= 0.5)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"> 0.5 to 0.75 (ref = <= 0.5)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"> 0.5 to 0.75 (ref = <= 0.5)", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"> 0.5 to 0.75 (ref = <= 0.5)", horizontal col(black)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"> 0.75 to 1 (ref = <= 0.5)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"> 0.75 to 1 (ref = <= 0.5)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"> 0.75 to 1 (ref = <= 0.5)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"> 0.75 to 1 (ref = <= 0.5)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"> 0.75 to 1 (ref = <= 0.5)", col(red) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"> 0.75 to 1 (ref = <= 0.5)", horizontal col(red)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"> 1 (ref = <= 0.5)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"> 1 (ref = <= 0.5)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"> 1 (ref = <= 0.5)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"> 1 (ref = <= 0.5)", horizontal col(blue)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"> 1 (ref = <= 0.5)", col(blue) msize(small) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"> 1 (ref = <= 0.5)", horizontal col(blue)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = <= 0.5)") ///
		title("Household crowding and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.2 0.3 0.5 0.7 1 2 4, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "> 0.5 to 0.75" 7 "> 0.75 to 1" 13 "> 1") rows(1)) ///
		name(crowding_cat, replace)
		
graph export ".\G1_Results\crowdingResults.pdf", replace


** Create plot for mobility (ref = 1/least deprived)

* As five exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "mobility_mat" & exp_level == "1 move (ref = 0 moves)"
replace level_split = level_num - 0.15 if exposure == "mobility_mat" & exp_level == "2 moves (ref = 0 moves)"
replace level_split = level_num if exposure == "mobility_mat" & exp_level == "3 moves (ref = 0 moves)"
replace level_split = level_num + 0.15 if exposure == "mobility_mat" & exp_level == "4 moves (ref = 0 moves)"
replace level_split = level_num + 0.3 if exposure == "mobility_mat" & exp_level == "5 + moves (ref = 0 moves)"
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
		xlabel(0.35 0.5 0.7 1 1.5 2 2.5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(order(1 "1" 7 "2" 13 "3" 19 "4" 25 "5+") rows(1)) ///
		name(mobility_cat, replace)

graph export ".\G1_Results\mobilityResults.pdf", replace

graph close _all


****** And now for some interaction plots (although really no interactions were found...)

** Create plot for parental income by age interaction
sum lci_int uci_int if exposure == "income_parents" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == ///
			"income_parents", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == ///
			"income_parents", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == ///
			"income_parents", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == ///
			"income_parents", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == ///
			"income_parents", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == ///
			"income_parents", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in log income)") ///
		title("Parental income*Sex Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.3 0.5 0.7 1 1.5, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(income_int, replace)
		
graph export ".\G1_Results\incomeResults_int.pdf", replace

graph close _all



*********************************************************************************
** For the multinomial regression results, as interpretation not intuitive. Will convert some results to predicted probabilities using the 'margins' command to provide more context/understanding of the effect sizes involved (see: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/)

use ".\G1_Results\G1_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

* YP education and attendance
mlogit YPG3080_rev ageAt28 male i.education, baseoutcome(0) rrr
margins education

* YP income and RSBB
mlogit YPG3000 ageAt28 male i.income, baseoutcome(3) rrr
margins income

mlogit YPG3040_grp ageAt28 male i.income, baseoutcome(3) rrr
margins income

mlogit YPG3080_rev ageAt28 male i.income, baseoutcome(0) rrr
margins income




********************************************************************************
*** Demonstrating in some simple simulations how differences in selection of the exposure and outcome determines the direction and magnitude of selection bias

* Set-up a log, so can store these results
capture log close
log using ".\G1_Results\G1_selection_sims_log", replace text

** Will use 10 million observations in each simulation, to remove issues of random variability. Will also focus on scenario with 2 binary variables for simplicity, but same principles extent to other types of data as well (with some exceptions, such as there being no bias if both the exposure and outcome are associated with selection independently; as we're simulating missing data under a logistic model here, there is an implicit dependence between x and y, even without an explicit interaction term, hence why bias is observed here - for more information, and a simple simulated example demonstraing this, see supplementary section S2 of Millard et al. 'Exploring selection bias in COVID-19 research: Simulations and prospective analyses of two UK cohort studies' [https://www.medrxiv.org/content/10.1101/2021.12.10.21267363v1]).

* Start with no causal effect of exposure on outcome
clear
set obs 10000000
set seed 3911

gen x = rbinomial(1, 0.5)
tab x

gen y_prob = invlogit(log(1) + log(1) * x)
sum y_prob
gen y = rbinomial(1, y_prob)
tab y

* Check these are independent and that odds ratio is as we simulated
logistic y x


*** Now create missingness mechanisms and explore extent and direction of bias

** x and y both positively predict having data (weak effects)
gen data_bothPos_weak_p = invlogit(log(0.5) + (log(2) * x) + (log(2) * y))
gen data_bothPos_weak = rbinomial(1, data_bothPos_weak_p)
tab data_bothPos_weak

logistic y x if data_bothPos_weak == 1

** x and y both positively predict having data (strong effects)
gen data_bothPos_strong_p = invlogit(log(0.1) + (log(10) * x) + (log(10) * y))
gen data_bothPos_strong = rbinomial(1, data_bothPos_strong_p)
tab data_bothPos_strong

logistic y x if data_bothPos_strong == 1

** x positively predicts having data, while y negatively predicts is (weak effects)
gen data_xPos_yNeg_weak_p = invlogit(log(1) + (log(2) * x) + (log(0.5) * y))
gen data_xPos_yNeg_weak = rbinomial(1, data_xPos_yNeg_weak_p)
tab data_xPos_yNeg_weak

logistic y x if data_xPos_yNeg_weak == 1

** x positively predicts having data, while y negatively predicts is (strong effects)
gen data_xPos_yNeg_strong_p = invlogit(log(1) + (log(10) * x) + (log(0.1) * y))
gen data_xPos_yNeg_strong = rbinomial(1, data_xPos_yNeg_strong_p)
tab data_xPos_yNeg_strong

logistic y x if data_xPos_yNeg_strong == 1

** x negatively predicts having data, while y positively predicts is (weak effects)
gen data_xNeg_yPos_weak_p = invlogit(log(1) + (log(0.5) * x) + (log(2) * y))
gen data_xNeg_yPos_weak = rbinomial(1, data_xNeg_yPos_weak_p)
tab data_xNeg_yPos_weak

logistic y x if data_xNeg_yPos_weak == 1

** x negatively predicts having data, while y positively predicts is (strong effects)
gen data_xNeg_yPos_strong_p = invlogit(log(1) + (log(0.1) * x) + (log(10) * y))
gen data_xNeg_yPos_strong = rbinomial(1, data_xNeg_yPos_strong_p)
tab data_xNeg_yPos_strong

logistic y x if data_xNeg_yPos_strong == 1

** x and y both negatively predict having data (weak effects)
gen data_bothNeg_weak_p = invlogit(log(2) + (log(0.5) * x) + (log(0.5) * y))
gen data_bothNeg_weak = rbinomial(1, data_bothNeg_weak_p)
tab data_bothNeg_weak

logistic y x if data_bothNeg_weak == 1

** x and y both negatively predict having data (strong effects)
gen data_bothNeg_strong_p = invlogit(log(10) + (log(0.1) * x) + (log(0.1) * y))
gen data_bothNeg_strong = rbinomial(1, data_bothNeg_strong_p)
tab data_bothNeg_strong

logistic y x if data_bothNeg_strong == 1


** Quick example to show that if selection of x and y are independent then there will be no selection bias. Here, the probability of having data for x or y is 0.5 if x or y take the value 1, meaning the probability of having data when both x and y are 1 is 0.25 (0.5 * 0.5).
gen data_xyIndep = .
replace data_xyIndep = 1 if x == 0 & y == 0
replace data_xyIndep = rbinomial(1, 0.5) if x == 0 & y == 1
replace data_xyIndep = rbinomial(1, 0.5) if x == 1 & y == 0
replace data_xyIndep = rbinomial(1, 0.25) if x == 1 & y == 1
tab data_xyIndep

logistic y x if data_xyIndep == 1


*** Now repeat when x causes an increase in y
clear
set obs 10000000
set seed 3912

gen x = rbinomial(1, 0.5)
tab x

gen y_prob = invlogit(log(0.448) + log(5) * x)
sum y_prob
gen y = rbinomial(1, y_prob)
tab y

* Check that odds ratio is as we simulated
logistic y x


*** Now create missingness mechanisms and explore extent and direction of bias

** x and y both positively predict having data (weak effects)
gen data_bothPos_weak_p = invlogit(log(0.5) + (log(2) * x) + (log(2) * y))
gen data_bothPos_weak = rbinomial(1, data_bothPos_weak_p)
tab data_bothPos_weak

logistic y x if data_bothPos_weak == 1

** x and y both positively predict having data (strong effects)
gen data_bothPos_strong_p = invlogit(log(0.1) + (log(10) * x) + (log(10) * y))
gen data_bothPos_strong = rbinomial(1, data_bothPos_strong_p)
tab data_bothPos_strong

logistic y x if data_bothPos_strong == 1

** x positively predicts having data, while y negatively predicts is (weak effects)
gen data_xPos_yNeg_weak_p = invlogit(log(1) + (log(2) * x) + (log(0.5) * y))
gen data_xPos_yNeg_weak = rbinomial(1, data_xPos_yNeg_weak_p)
tab data_xPos_yNeg_weak

logistic y x if data_xPos_yNeg_weak == 1

** x positively predicts having data, while y negatively predicts is (strong effects)
gen data_xPos_yNeg_strong_p = invlogit(log(1) + (log(10) * x) + (log(0.1) * y))
gen data_xPos_yNeg_strong = rbinomial(1, data_xPos_yNeg_strong_p)
tab data_xPos_yNeg_strong

logistic y x if data_xPos_yNeg_strong == 1

** x negatively predicts having data, while y positively predicts is (weak effects)
gen data_xNeg_yPos_weak_p = invlogit(log(1) + (log(0.5) * x) + (log(2) * y))
gen data_xNeg_yPos_weak = rbinomial(1, data_xNeg_yPos_weak_p)
tab data_xNeg_yPos_weak

logistic y x if data_xNeg_yPos_weak == 1

** x negatively predicts having data, while y positively predicts is (strong effects)
gen data_xNeg_yPos_strong_p = invlogit(log(1) + (log(0.1) * x) + (log(10) * y))
gen data_xNeg_yPos_strong = rbinomial(1, data_xNeg_yPos_strong_p)
tab data_xNeg_yPos_strong

logistic y x if data_xNeg_yPos_strong == 1

** x and y both negatively predict having data (weak effects)
gen data_bothNeg_weak_p = invlogit(log(2) + (log(0.5) * x) + (log(0.5) * y))
gen data_bothNeg_weak = rbinomial(1, data_bothNeg_weak_p)
tab data_bothNeg_weak

logistic y x if data_bothNeg_weak == 1

** x and y both negatively predict having data (strong effects)
gen data_bothNeg_strong_p = invlogit(log(10) + (log(0.1) * x) + (log(0.1) * y))
gen data_bothNeg_strong = rbinomial(1, data_bothNeg_strong_p)
tab data_bothNeg_strong

logistic y x if data_bothNeg_strong == 1


** Quick example to show that if selection of x and y are independent then there will be no selection bias. Here, the probability of having data for x or y is 0.5 if x or y take the value 1, meaning the probability of having data when both x and y are 1 is 0.25 (0.5 * 0.5).
gen data_xyIndep = .
replace data_xyIndep = 1 if x == 0 & y == 0
replace data_xyIndep = rbinomial(1, 0.5) if x == 0 & y == 1
replace data_xyIndep = rbinomial(1, 0.5) if x == 1 & y == 0
replace data_xyIndep = rbinomial(1, 0.25) if x == 1 & y == 1
tab data_xyIndep

logistic y x if data_xyIndep == 1


*** Now repeat when x causes a decrease in y
clear
set obs 10000000
set seed 3913

gen x = rbinomial(1, 0.5)
tab x

gen y_prob = invlogit(log(2.24) + log(0.2) * x)
sum y_prob
gen y = rbinomial(1, y_prob)
tab y

* Check that odds ratio is as we simulated
logistic y x


*** Now create missingness mechanisms and explore extent and direction of bias

** x and y both positively predict having data (weak effects)
gen data_bothPos_weak_p = invlogit(log(0.5) + (log(2) * x) + (log(2) * y))
gen data_bothPos_weak = rbinomial(1, data_bothPos_weak_p)
tab data_bothPos_weak

logistic y x if data_bothPos_weak == 1

** x and y both positively predict having data (strong effects)
gen data_bothPos_strong_p = invlogit(log(0.1) + (log(10) * x) + (log(10) * y))
gen data_bothPos_strong = rbinomial(1, data_bothPos_strong_p)
tab data_bothPos_strong

logistic y x if data_bothPos_strong == 1

** x positively predicts having data, while y negatively predicts is (weak effects)
gen data_xPos_yNeg_weak_p = invlogit(log(1) + (log(2) * x) + (log(0.5) * y))
gen data_xPos_yNeg_weak = rbinomial(1, data_xPos_yNeg_weak_p)
tab data_xPos_yNeg_weak

logistic y x if data_xPos_yNeg_weak == 1

** x positively predicts having data, while y negatively predicts is (strong effects)
gen data_xPos_yNeg_strong_p = invlogit(log(1) + (log(10) * x) + (log(0.1) * y))
gen data_xPos_yNeg_strong = rbinomial(1, data_xPos_yNeg_strong_p)
tab data_xPos_yNeg_strong

logistic y x if data_xPos_yNeg_strong == 1

** x negatively predicts having data, while y positively predicts is (weak effects)
gen data_xNeg_yPos_weak_p = invlogit(log(1) + (log(0.5) * x) + (log(2) * y))
gen data_xNeg_yPos_weak = rbinomial(1, data_xNeg_yPos_weak_p)
tab data_xNeg_yPos_weak

logistic y x if data_xNeg_yPos_weak == 1

** x negatively predicts having data, while y positively predicts is (strong effects)
gen data_xNeg_yPos_strong_p = invlogit(log(1) + (log(0.1) * x) + (log(10) * y))
gen data_xNeg_yPos_strong = rbinomial(1, data_xNeg_yPos_strong_p)
tab data_xNeg_yPos_strong

logistic y x if data_xNeg_yPos_strong == 1

** x and y both negatively predict having data (weak effects)
gen data_bothNeg_weak_p = invlogit(log(2) + (log(0.5) * x) + (log(0.5) * y))
gen data_bothNeg_weak = rbinomial(1, data_bothNeg_weak_p)
tab data_bothNeg_weak

logistic y x if data_bothNeg_weak == 1

** x and y both negatively predict having data (strong effects)
gen data_bothNeg_strong_p = invlogit(log(10) + (log(0.1) * x) + (log(0.1) * y))
gen data_bothNeg_strong = rbinomial(1, data_bothNeg_strong_p)
tab data_bothNeg_strong

logistic y x if data_bothNeg_strong == 1


** Quick example to show that if selection of x and y are independent then there will be no selection bias. Here, the probability of having data for x or y is 0.5 if x or y take the value 1, meaning the probability of having data when both x and y are 1 is 0.25 (0.5 * 0.5).
gen data_xyIndep = .
replace data_xyIndep = 1 if x == 0 & y == 0
replace data_xyIndep = rbinomial(1, 0.5) if x == 0 & y == 1
replace data_xyIndep = rbinomial(1, 0.5) if x == 1 & y == 0
replace data_xyIndep = rbinomial(1, 0.25) if x == 1 & y == 1
tab data_xyIndep

logistic y x if data_xyIndep == 1


log close




********************************************************************************
*** Next, will show an example of using multiple imputation to try and explore/address potential selection bias

** Read in the 'post-analysis' dataset, and merge with the mother's and father's data to get maternal and paternal RSBB variables

* First, open raw dataset, and just keep parental RSBB variables, and drop if QLET == B
use "Desc_RSBB_B3911.dta", clear

* Keep just relevant admin and RSBB vars
keep aln qlet d810 d813 d816 pb150 pb153 pb155

* Drop if second-born twin
tab qlet, m
drop if qlet == "B"

* Check no duplicates
duplicates report

* Drop QLET variable
drop qlet

* Check the RSBB variables and tidy
tab1 d810 d813 d816 pb150 pb153 pb155, m

replace d810 = . if d810 < 0
tab d810, m
tab d810

replace d813 = . if d813 < 0
tab d813, m
recode d813 (0 = 1) (1/6 = 2) (7/13 = 3), gen(d813_grp)
label define relig_lb 1 "None" 2 "Christian" 3 "Other"
numlabel relig_lb, add
label values d813_grp relig_lb
tab d813_grp, m
tab d813_grp

replace d816 = . if d816 < 0
tab d816, m
tab d816

replace pb150 = . if pb150 < 0
tab pb150, m
tab pb150

replace pb153 = . if pb153 < 0
tab pb153, m
recode pb153 (0 = 1) (1/6 = 2) (7/13 = 3), gen(pb153_grp)
label values pb153_grp relig_lb
tab pb153_grp, m
tab pb153_grp

replace pb155 = . if pb155 < 0
tab pb155, m
tab pb155

* Drop some old vars and add num labels
drop d813 pb153

numlabel, add

* Drop any WoCs
tab1 d810 d813_grp d816 pb150 pb153_grp pb155, m

drop if d810 == .a | pb150 == .c

* Save this file
save ".\G1_Results\G0_RSBB_forMerging.dta", replace


** Now load the G1 post-analysis dataset and merge together
use ".\G1_Results\G1_PredictorsOfRSBB_B3911_postAnalysis.dta", clear

merge m:1 aln using ".\G1_Results\G0_RSBB_forMerging.dta", keepusing(d810 d813_grp d816 pb150 pb153_grp pb155)

* Keep just those in the G1 file
tab _merge

drop if _merge == 2
drop _merge

tab1 d810 d813_grp d816 pb150 pb153_grp pb155, m


** For this imputation example, we are going to impute a handful of sociodemographic exposures and see how these imputed results differ from the complete case results. The exposures we're going to focus on are:
*  - Sex (no missing data in exposure, but definitely causes selection [with males being less likely to continue participating in ALSPAC])
*  - Ethnicity (little missing data in exposure, and likely to cause selection [with ethnicities other than White less likely to continue participating])
*  - YP education (lots of missing data [as only measured at age 27]), and likely to cause selection [with lower education less likely to continue participating])
*  - YP income (lots of missing data [as only measured at age 26]), and likely to cause selection [with lower income less likely to continue participating])

** As the YP RSBB outcomes also have high levels of missingness (as only measured at age 28), will have to impute these as well. Will use maternal and paternal RSBB from pregnancy as auxiliary variables, as both have less missing data, and likely predict child RSBB.

** For YP education and income, will use a range of parental SEP and demographic variables as auxiliarly variables (maternal and paternal education, maternal and paternal occupational social class, maternal IMD, household income, maternal age at birth, mother's marital status, parity, mother's urban/rural location, home ownership status and parental access to a car), as have less missing data and predict child education and income.

** This multiple imputation approach rests on various assumptions, the majority of which are untestable and may not be met. For instance, using parental RSBB and SEP as auxiliary variables assumes that these strongly predict child RSBB and SEP; while some association is evident, whether this is sufficient to satisfy the MAR assumption is unclear. This approach also assumes that the imputation model is correctly specified, and that there are no additional variables which also relate to missing data in the exposure and/or outcome which have not been included here.


*** Keep just variables of interest, check amounts of missing data, then set up the imputation model
keep aln qlet YPG3000 YPG3040_grp YPG3080_rev ageAt28 ///
male nonWhiteEthnic education income ///
mother_ageAtBirth maritalStatus_mat rural_mat parity_mat education_mat education_pat highSocClass_mat highSocClass_pat income_parents IMD_mat housing_mat accessToCar_parents ///
d810 d813_grp d816 pb150 pb153_grp pb155

order aln qlet YPG3000 YPG3040_grp YPG3080_rev ageAt28 ///
male nonWhiteEthnic education income ///
mother_ageAtBirth maritalStatus_mat rural_mat parity_mat education_mat education_pat highSocClass_mat highSocClass_pat income_parents IMD_mat housing_mat accessToCar_parents ///
d810 d813_grp d816 pb150 pb153_grp pb155

* Check missing data (only sex has no missing data)
misstable summarize YPG3000-pb155, all

* Set up the imputation model
mi set flong
mi register regular male
mi register imputed YPG3000 YPG3040_grp YPG3080_rev ageAt28 ///
nonWhiteEthnic education income ///
mother_ageAtBirth maritalStatus_mat rural_mat parity_mat education_mat education_pat highSocClass_mat highSocClass_pat income_parents IMD_mat housing_mat accessToCar_parents ///
d810 d813_grp d816 pb150 pb153_grp pb155

* Set-up a dry-run, just to make sure imputation works and models specified correctly
mi impute chained ///
	(regress) mother_ageAtBirth income_parents accessToCar_parents ///
	(pmm, knn(5)) ageAt28 ///
	(logit) nonWhiteEthnic rural_mat highSocClass_mat highSocClass_pat ///
	(mlogit) YPG3000 YPG3040_grp maritalStatus_mat housing_mat d810 d813_grp pb150 pb153_grp ///
	(ologit) YPG3080_rev education income parity_mat education_mat education_pat IMD_mat d816 pb155 ///
	= male, ///
	add(50) burnin(10) rseed(12345) dryrun
	
* Looks okay, so now run the actual imputations. Will perform 50 imputations with a burn-in period of 10 (and saving the tracplot to check that convergence reached). Use 'dots' option to show progress, and 'augment' to avoid perfect prediction during imputations. As lots of variables in imputation model, and relatively large dataset, this takes a good few hours.
mi impute chained ///
	(regress) mother_ageAtBirth income_parents accessToCar_parents ///
	(pmm, knn(5)) ageAt28 ///
	(logit) nonWhiteEthnic rural_mat highSocClass_mat highSocClass_pat ///
	(mlogit) YPG3000 YPG3040_grp maritalStatus_mat housing_mat d810 d813_grp pb150 pb153_grp ///
	(ologit) YPG3080_rev education income parity_mat education_mat education_pat IMD_mat d816 pb155 ///
	= male, ///
	add(50) burnin(10) rseed(12345) dots augment ///
	savetrace(".\G1_Results\imp_trace.dta", replace)
	
	
** Save this imputed dataset, so not have to run whole imputation again to access results
save ".\G1_Results\imp_G1_RSBB.dta", replace


** Check convergence and that imputation chains are well-mixed
* Read in the trace dataset
use ".\G1_Results\imp_trace.dta", clear

sum 

* Save the mean value to add as a line in the plot - Do this for all outcomes and exposures
sum nonWhiteEthnic_mean
local mean_nonWhiteEthnic = r(mean)
display `mean_nonWhiteEthnic'

sum education_mean
local mean_education = r(mean)
display `mean_education'

sum income_mean
local mean_income = r(mean)
display `mean_income'

sum YPG3000_mean
local mean_belief = r(mean)
display `mean_belief'

sum YPG3040_grp_mean
local mean_denom = r(mean)
display `mean_denom'

sum YPG3080_rev_mean
local mean_attend = r(mean)
display `mean_attend'


* Convert the data from long to wide format (is necessary to create the plots)
reshape wide *mean *sd, i(iter) j(m)

* Set the iteration variable as the 'time' variable
tsset iter

* Make the plots - These all look relatively well-mixed and converged
tsline nonWhiteEthnic_mean*, yline(`mean_nonWhiteEthnic') legend(off) name(ethnic, replace)
tsline education_mean*, yline(`mean_education') legend(off) name(edu, replace)
tsline income_mean*, yline(`mean_income') legend(off) name(income, replace)
tsline YPG3000_mean*, yline(`mean_belief') legend(off) name(belief, replace)
tsline YPG3040_grp_mean*, yline(`mean_denom') legend(off) name(relig, replace)
tsline YPG3080_rev_mean*, yline(`mean_attend') legend(off) name(attend, replace)

graph close _all


**** Now run the models on the imputed data and combine using Rubin's Rules
use ".\G1_Results\imp_G1_RSBB.dta", clear


* Set-up a log, so can store these results
capture log close
log using ".\G1_Results\G1_imp_analysis_log", replace text


** Quick check of imputed vs observed results (all look reasonable and no obvious outliers/causes for concern)
tab nonWhiteEthnic if _mi_m == 0
tab nonWhiteEthnic if _mi_m == 1

tab education if _mi_m == 0
tab education if _mi_m == 1

tab income if _mi_m == 0
tab income if _mi_m == 1

tab YPG3000 if _mi_m == 0
tab YPG3000 if _mi_m == 1

tab YPG3040_grp if _mi_m == 0
tab YPG3040_grp if _mi_m == 1

tab YPG3080_rev if _mi_m == 0
tab YPG3080_rev if _mi_m == 1


**** Now run the actual models comparing complete-cases analysis to imputed analysis

*** Sex and RSBB

** Religious belief

* CCA
mlogit YPG3000 ageAt28 male if _mi_m == 0, baseoutcome(3) rrr
matrix res = r(table)
matrix list res

* MI (can't display results in exponentiated format from mi estimate, so have to exponentiate manually)
mi estimate: mlogit YPG3000 ageAt28 male, baseoutcome(3)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Y: male"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: male"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: male"])
di "p-value = " %9.4f res["pvalue", "1__Y: male"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: male"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: male"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: male"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: male"]


** Religious affiliation

* CCA
mlogit YPG3040_grp ageAt28 male if _mi_m == 0, baseoutcome(3) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3040_grp ageAt28 male, baseoutcome(3)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Christian: male"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: male"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: male"])
di "p-value = " %9.4f res["pvalue", "1__Christian: male"]

di "RRR = " %9.3f exp(res["b", "2__Other: male"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: male"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: male"])
di "p-value = " %9.4f res["pvalue", "2__Other: male"]


** Religious attendance

* CCA
mlogit YPG3080_rev ageAt28 male if _mi_m == 0, baseoutcome(0) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3080_rev ageAt28 male, baseoutcome(0)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Occasionally: male"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: male"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: male"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: male"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: male"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: male"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: male"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: male"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: male"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: male"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: male"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: male"]


*** Ethnicity and RSBB

** Religious belief

* CCA
mlogit YPG3000 ageAt28 male nonWhiteEthnic if _mi_m == 0, baseoutcome(3) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3000 ageAt28 male nonWhiteEthnic, baseoutcome(3)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Y: nonWhiteEthnic"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: nonWhiteEthnic"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: nonWhiteEthnic"])
di "p-value = " %9.4f res["pvalue", "1__Y: nonWhiteEthnic"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: nonWhiteEthnic"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: nonWhiteEthnic"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: nonWhiteEthnic"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: nonWhiteEthnic"]


** Religious affiliation

* CCA
mlogit YPG3040_grp ageAt28 male nonWhiteEthnic if _mi_m == 0, baseoutcome(3) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3040_grp ageAt28 male nonWhiteEthnic, baseoutcome(3)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Christian: nonWhiteEthnic"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: nonWhiteEthnic"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: nonWhiteEthnic"])
di "p-value = " %9.4f res["pvalue", "1__Christian: nonWhiteEthnic"]

di "RRR = " %9.3f exp(res["b", "2__Other: nonWhiteEthnic"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: nonWhiteEthnic"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: nonWhiteEthnic"])
di "p-value = " %9.4f res["pvalue", "2__Other: nonWhiteEthnic"]


** Religious attendance

* CCA
mlogit YPG3080_rev ageAt28 male nonWhiteEthnic if _mi_m == 0, baseoutcome(0) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3080_rev ageAt28 male nonWhiteEthnic, baseoutcome(0)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Occasionally: nonWhiteEthnic"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: nonWhiteEthnic"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: nonWhiteEthnic"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: nonWhiteEthnic"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: nonWhiteEthnic"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: nonWhiteEthnic"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: nonWhiteEthnic"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: nonWhiteEthnic"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: nonWhiteEthnic"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: nonWhiteEthnic"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: nonWhiteEthnic"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: nonWhiteEthnic"]


*** YP education and RSBB

** Religious belief

* CCA
mlogit YPG3000 ageAt28 male i.education if _mi_m == 0, baseoutcome(3) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3000 ageAt28 male i.education, baseoutcome(3)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Y: 2.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: 2.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: 2.education"])
di "p-value = " %9.4f res["pvalue", "1__Y: 2.education"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: 2.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: 2.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: 2.education"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: 2.education"]

di "RRR = " %9.3f exp(res["b", "1__Y: 3.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: 3.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: 3.education"])
di "p-value = " %9.4f res["pvalue", "1__Y: 3.education"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: 3.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: 3.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: 3.education"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: 3.education"]

di "RRR = " %9.3f exp(res["b", "1__Y: 4.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: 4.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: 4.education"])
di "p-value = " %9.4f res["pvalue", "1__Y: 4.education"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: 4.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: 4.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: 4.education"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: 4.education"]


** Religious affiliation

* CCA
mlogit YPG3040_grp ageAt28 male i.education if _mi_m == 0, baseoutcome(3) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3040_grp ageAt28 male i.education, baseoutcome(3)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Christian: 2.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: 2.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: 2.education"])
di "p-value = " %9.4f res["pvalue", "1__Christian: 2.education"]

di "RRR = " %9.3f exp(res["b", "2__Other: 2.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: 2.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: 2.education"])
di "p-value = " %9.4f res["pvalue", "2__Other: 2.education"]

di "RRR = " %9.3f exp(res["b", "1__Christian: 3.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: 3.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: 3.education"])
di "p-value = " %9.4f res["pvalue", "1__Christian: 3.education"]

di "RRR = " %9.3f exp(res["b", "2__Other: 3.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: 3.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: 3.education"])
di "p-value = " %9.4f res["pvalue", "2__Other: 3.education"]

di "RRR = " %9.3f exp(res["b", "1__Christian: 4.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: 4.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: 4.education"])
di "p-value = " %9.4f res["pvalue", "1__Christian: 4.education"]

di "RRR = " %9.3f exp(res["b", "2__Other: 4.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: 4.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: 4.education"])
di "p-value = " %9.4f res["pvalue", "2__Other: 4.education"]


** Religious attendance

* CCA
mlogit YPG3080_rev ageAt28 male i.education if _mi_m == 0, baseoutcome(0) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3080_rev ageAt28 male i.education, baseoutcome(0)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Occasionally: 2.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: 2.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: 2.education"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: 2.education"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: 2.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: 2.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: 2.education"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: 2.education"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: 2.education"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: 2.education"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: 2.education"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: 2.education"]

di "RRR = " %9.3f exp(res["b", "1__Occasionally: 3.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: 3.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: 3.education"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: 3.education"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: 3.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: 3.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: 3.education"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: 3.education"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: 3.education"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: 3.education"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: 3.education"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: 3.education"]

di "RRR = " %9.3f exp(res["b", "1__Occasionally: 4.education"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: 4.education"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: 4.education"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: 4.education"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: 4.education"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: 4.education"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: 4.education"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: 4.education"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: 4.education"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: 4.education"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: 4.education"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: 4.education"]



*** YP income and RSBB

** Religious belief

* CCA
mlogit YPG3000 ageAt28 male i.income if _mi_m == 0, baseoutcome(3) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3000 ageAt28 male i.income, baseoutcome(3)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Y: 2.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: 2.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: 2.income"])
di "p-value = " %9.4f res["pvalue", "1__Y: 2.income"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: 2.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: 2.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: 2.income"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: 2.income"]

di "RRR = " %9.3f exp(res["b", "1__Y: 3.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: 3.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: 3.income"])
di "p-value = " %9.4f res["pvalue", "1__Y: 3.income"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: 3.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: 3.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: 3.income"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: 3.income"]

di "RRR = " %9.3f exp(res["b", "1__Y: 4.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: 4.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: 4.income"])
di "p-value = " %9.4f res["pvalue", "1__Y: 4.income"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: 4.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: 4.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: 4.income"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: 4.income"]

di "RRR = " %9.3f exp(res["b", "1__Y: 5.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Y: 5.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Y: 5.income"])
di "p-value = " %9.4f res["pvalue", "1__Y: 5.income"]

di "RRR = " %9.3f exp(res["b", "2__Not_sure: 5.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Not_sure: 5.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Not_sure: 5.income"])
di "p-value = " %9.4f res["pvalue", "2__Not_sure: 5.income"]


** Religious affiliation

* CCA
mlogit YPG3040_grp ageAt28 male i.income if _mi_m == 0, baseoutcome(3) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3040_grp ageAt28 male i.income, baseoutcome(3)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Christian: 2.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: 2.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: 2.income"])
di "p-value = " %9.4f res["pvalue", "1__Christian: 2.income"]

di "RRR = " %9.3f exp(res["b", "2__Other: 2.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: 2.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: 2.income"])
di "p-value = " %9.4f res["pvalue", "2__Other: 2.income"]

di "RRR = " %9.3f exp(res["b", "1__Christian: 3.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: 3.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: 3.income"])
di "p-value = " %9.4f res["pvalue", "1__Christian: 3.income"]

di "RRR = " %9.3f exp(res["b", "2__Other: 3.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: 3.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: 3.income"])
di "p-value = " %9.4f res["pvalue", "2__Other: 3.income"]

di "RRR = " %9.3f exp(res["b", "1__Christian: 4.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: 4.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: 4.income"])
di "p-value = " %9.4f res["pvalue", "1__Christian: 4.income"]

di "RRR = " %9.3f exp(res["b", "2__Other: 4.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: 4.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: 4.income"])
di "p-value = " %9.4f res["pvalue", "2__Other: 4.income"]

di "RRR = " %9.3f exp(res["b", "1__Christian: 5.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Christian: 5.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Christian: 5.income"])
di "p-value = " %9.4f res["pvalue", "1__Christian: 5.income"]

di "RRR = " %9.3f exp(res["b", "2__Other: 5.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Other: 5.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Other: 5.income"])
di "p-value = " %9.4f res["pvalue", "2__Other: 5.income"]


** Religious attendance

* CCA
mlogit YPG3080_rev ageAt28 male i.income if _mi_m == 0, baseoutcome(0) rrr
matrix res = r(table)
matrix list res

* MI
mi estimate: mlogit YPG3080_rev ageAt28 male i.income, baseoutcome(0)
matrix res = r(table)
matrix list res

di "RRR = " %9.3f exp(res["b", "1__Occasionally: 2.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: 2.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: 2.income"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: 2.income"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: 2.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: 2.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: 2.income"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: 2.income"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: 2.income"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: 2.income"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: 2.income"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: 2.income"]

di "RRR = " %9.3f exp(res["b", "1__Occasionally: 3.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: 3.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: 3.income"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: 3.income"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: 3.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: 3.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: 3.income"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: 3.income"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: 3.income"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: 3.income"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: 3.income"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: 3.income"]

di "RRR = " %9.3f exp(res["b", "1__Occasionally: 4.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: 4.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: 4.income"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: 4.income"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: 4.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: 4.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: 4.income"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: 4.income"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: 4.income"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: 4.income"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: 4.income"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: 4.income"]

di "RRR = " %9.3f exp(res["b", "1__Occasionally: 5.income"])
di "Lower CI = " %9.3f exp(res["ll", "1__Occasionally: 5.income"])
di "Upper CI =" %9.3f exp(res["ul", "1__Occasionally: 5.income"])
di "p-value = " %9.4f res["pvalue", "1__Occasionally: 5.income"]

di "RRR = " %9.3f exp(res["b", "2__Min_once_a_year: 5.income"])
di "Lower CI = " %9.3f exp(res["ll", "2__Min_once_a_year: 5.income"])
di "Upper CI =" %9.3f exp(res["ul", "2__Min_once_a_year: 5.income"])
di "p-value = " %9.4f res["pvalue", "2__Min_once_a_year: 5.income"]

di "RRR = " %9.3f exp(res["b", "3__Min_once_a_month: 5.income"])
di "Lower CI = " %9.3f exp(res["ll", "3__Min_once_a_month: 5.income"])
di "Upper CI =" %9.3f exp(res["ul", "3__Min_once_a_month: 5.income"])
di "p-value = " %9.4f res["pvalue", "3__Min_once_a_month: 5.income"]


* Close the log file and clear data
log close

clear

