*** Predictors of RSBB (B3911) - G1 analysis script
*** Created 16/11/2021 by Dan Smith
*** Stata v16.0

*** This script reads in the cleaned G1 data, explores associations between exposures, and then conducts an 'exposome-wide association analysis' (ExWAS) on all these variables to examine how they are associated with various facets of RSBB.


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

** And also install 'spost' commands for testing proportional odds assumption of ordinal regression models (to install 'spost', type 'search spost13' and install the 'spost13_ado' package)


**********************************************************************************
*** Descriptive statistics

** Put the RSBB variables at the start of the dataset
order aln qlet YPG3000 YPG3040 YPG3040_grp YPG3080 YPG3080_OccNever YPG3080_OccYr YPG3153 YPG3153_cat YPG3160 YPG3170 YPG3155 YPG3155_cat

** Using the 'distinct' command (see above), for each variable inspect the number of unque values; if < 10 then display table, while if >= 10 then displays means/SDs/IQRs

* Outcomes
foreach var of varlist YPG3000-YPG3155_cat {
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
foreach var of varlist kz021-prosocial25 {
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
misstable sum YPG3000-YPG3155_cat, all

* Exposures
misstable sum male-prosocial25, all



**********************************************************************************
*** Correlations between exposures

** Explore correlations between exposures to see how inter-related these factors are
desc male-prosocial25

* Most variables are either continuous, ordered categories, or binary variables, so will just use standard pearson correlations for these. Only unordered categories are mother's home ownership (owned/mortaged vs renting vs counci/housing association vs other) and marital status (never married vs married vs widowed/divorced/separated). So will exclude these two variables from the correlation matrix and calculate their associations using these variables as outcomes in a multinomial regression and square-rooting the pseudo-R2 value (similar to how the 'rexposome' package in R for ExWAS analyses does)

* First, order variables into categories of 'demographic', 'socioeconomic/material insecurity' and 'cognitive/psychological' (Y9992 is age at @28 RSBB questionnaire, so will pop at end as not needed here)
order aln-YPG3155_cat ///
YPG8000 mz028b male c804 a525_grp a005_grp jan2020ur01ind_grp b032_grp parent ///
yp_edu c645a c755_grp logavinceq jan2020imd2010q5 jan2020Townsendq5 a006_grp yp_finDiffs physical_abuse_0_16yrs-ACEcat_classic_0_16yrs a551 a636 father_ab age_FA ///
f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126

* Next, rename all these exposures so they are more intuitive and can be read easier on the correlation heatmaps below
rename YPG8000 ageAt28
rename mz028b mother_ageAtBirth
rename c804 nonWhiteEthnic
rename a525_grp maritalStatus
rename a005_grp mobility 
rename jan2020ur01ind_grp rural
rename b032_grp parity
rename yp_edu education
rename c645a maternalEdu
rename c755_grp highSocClass
rename logavinceq income
rename jan2020imd2010q5 IMD
rename jan2020Townsendq5 townsendDep
rename a006_grp housing 
rename yp_finDiffs financeDiffs

drop physical_abuse_0_16yrs-parent_child_bond_0_16yrs ACEcat_extended_0_16yrs ACEcat_classic_0_16yrs
rename ACEscore_extended_0_16yrs ACEscore_13items
rename ACEscore_classic_0_16yrs ACEscore_10items

rename a551 crowding
rename a636 neighPercept
rename father_ab fatherAbsence
rename age_FA age_fatherAbsence
rename f8ws110 verbalIQ_age8
rename f8ws111 performanceIQ_age8
rename f8ws112 totalIQ_age8
rename fh6280 totalIQ_age15
rename FKWI1030 digitSymbol_age24
rename FKWI1050 vocab_age24
rename fg7360 extraversion_age13
rename fg7361 agreeableness_age13
rename fg7362 conscientiousness_age13
rename fg7363 emotionalStab_age13
rename fg7364 Openess_age13
rename f8lc125 loc_age8
rename FJCQ1001 negCogStyles_age17
rename f8dv440a emoRec_faces_age8
rename kr554a skuseSocCog_age8
rename kq348a SDQ_prosocial_age8
rename tc4025e SDQ_prosocial_age13
rename CCXD860a esteem_bachman_age17
rename f8se125 scholasticEsteem_age8
rename f8se126 globalEsteem_age8
rename triangles_total emoRec_triangles_age13
rename skuse16 skuseSocCog_age16
rename autism25 autismSpec_age25
rename prosocial25 SDQ_prosocial_age25

* Associations between demographic factors (excluding marital status; a525_grp) - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
corr ageAt28 mother_ageAtBirth nonWhiteEthnic mobility rural parity parent

matrix cor_demo = r(C)
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

foreach var of varlist ageAt28 mother_ageAtBirth nonWhiteEthnic mobility rural parity parent {
	quietly mlogit maritalStatus `var'
	local r2 = e(r2_p)
	display "Estimated correlation for " "`var'" " on marital status is: " round((sqrt(`r2')), .001)
	
	post marital_corrs_demoOnly ("`var'") (sqrt(`r2'))
}

postclose marital_corrs_demoOnly

* Read in this file and save as CSV file, so easier to access (use 'preserve' and 'restore' to keep main file in memory)
preserve
use ".\G1_Results\marital_corrs_demoOnly.dta", clear
list, clean
outsheet using ".\G1_Results\marital_corrs_demoOnly.csv", comma replace
restore


** Now repeat for socioeconomic/material insecurity variables (to exclude housing status [housing], as is an unordered categorical variable)
corr education-townsendDep financeDiffs-age_fatherAbsence

matrix cor_socio = r(C)

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
	
foreach var of varlist education-townsendDep financeDiffs-age_fatherAbsence {
	quietly mlogit housing `var'
	local r2 = e(r2_p)
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


** Now repeat for cognitive/psychological variables (no unordered categorical variables, so can include all variables)
corr verbalIQ_age8-globalEsteem_age8

matrix cor_cog = r(C)

heatplot cor_cog, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle(""))

* Save heatmap
graph export ".\G1_Results\corr_heatplot_cogOnly.pdf", replace

* Save matrix as Excel file
putexcel set ".\G1_Results\corrMatrix_cogOnly.xlsx", replace
putexcel A1=matrix(cor_cog), names

	
*** Finally, repeat this on all of the exposures together (excluding unordered cateogorical variables housing status and marital status)
corr ageAt28 mother_ageAtBirth nonWhiteEthnic mobility rural parity parent education-townsendDep financeDiffs-age_fatherAbsence verbalIQ_age8-globalEsteem_age8

matrix cor_all = r(C)

heatplot cor_all, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(tiny)) ///
	ylabel(, labsize(tiny)) legend(subtitle(""))

* Save heatmap
graph export ".\G1_Results\corr_heatplot_all.pdf", replace

* Save matrix as Excel file
putexcel set ".\G1_Results\corrMatrix_all.xlsx", replace
putexcel A1=matrix(cor_all), names

	
* And now for associations of all other exposures variables with unordered categorical variables marital status and home ownership status

* Marital status
capture postclose marital_corrs_all
postfile marital_corrs_all str30 variable corr ///
	using ".\G1_Results\marital_corrs_all.dta", replace
	
foreach var of varlist ageAt28 mother_ageAtBirth nonWhiteEthnic mobility rural parity parent education-townsendDep financeDiffs-age_fatherAbsence verbalIQ_age8-globalEsteem_age8 {
	quietly mlogit maritalStatus `var'
	local r2 = e(r2_p)
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
	
foreach var of varlist ageAt28 mother_ageAtBirth nonWhiteEthnic mobility rural parity parent education-townsendDep financeDiffs-age_fatherAbsence verbalIQ_age8-globalEsteem_age8 {
	quietly mlogit housing `var'
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



**************************************************************************************
*** Next, we want to run the actual ExWAS analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab YPG3000, m

** This will be quite complicated, as want to post results to file, but as exposures differ extracting the results will be variable specific. To adjust for multiple corrections will use conservative bonferroni adjustment when constructing confidence intervals and interpreting p-values - As 48 exposures, will use 99.9% confidence intervals (as this is 100 - 5/48) and a p-value threshold of 0.05/48 = 0.0010.
display round(100 - 5/48, 0.01)
display 0.05/48

** We also want to store both estimates adjusting for age and sex, and also the interaction between sex and the exposure (as all G1s are of a similar age, will only explore interactions with sex, not age), to see whether it's moderated by sex. Again, this makes the set-up a bit more complicated.

** Create a postfile to post results to, then start the loop - Will create two postfiles; one for coefficients and CIs, and another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model) - NOTE: Have to store pvalues as 'double' format, else really tiny p-values get coded as '0' (as default if float format, which has minimum value of -3.40282346639e+38 [https://blog.stata.com/2012/04/02/the-penultimate-guide-to-precision/]).
capture postclose G1_belief
postfile G1_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_belief_results.dta", replace

capture postclose G1_belief_lr
postfile G1_belief_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_belief_results_lr.dta", replace

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' and 'sex' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3000 male `var', baseoutcome(3) rrr level(99.9)
		
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
		mlogit YPG3000 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		
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
		
		mlogit YPG3000 male `var', baseoutcome(3) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3000 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3000 male if `var' != ., baseoutcome(3) rrr level(99.9)
		est store base
		mlogit YPG3000 male `var', baseoutcome(3) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3000 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3000 ageAt28 `var', baseoutcome(3) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3000 ageAt28 if `var' != ., baseoutcome(3) rrr level(99.9)
		est store base
		mlogit YPG3000 ageAt28 `var', baseoutcome(3) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		mlogit YPG3000 ageAt28 male `var', baseoutcome(3) rrr level(99.9)
		
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
		mlogit YPG3000 ageAt28 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		
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
		mlogit YPG3000 ageAt28 male `var', baseoutcome(3) rrr level(99.9)
		
		local outcome_level = "Not sure (ref = No)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
				
		// Now for interaction model
		mlogit YPG3000 ageAt28 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3000 ageAt28 male if `var' != ., baseoutcome(3) rrr level(99.9)
		est store base
		mlogit YPG3000 ageAt28 male `var', baseoutcome(3) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3000 ageAt28 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
		
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
		
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
		
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/yes)
			local outcome_level = "Yes (ref = No)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
		
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Not sure (ref = No)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit YPG3000 ageAt28 male if `var' != ., baseoutcome(3) rrr level(99.9)
		est store base
		mlogit YPG3000 ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3000 ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_belief_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_belief
postclose G1_belief_lr
	

**************************************************************************************
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

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' and 'sex' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3040_grp male `var', baseoutcome(3) rrr level(99.9)
		
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
		mlogit YPG3040_grp c.male##c.`var', baseoutcome(3) rrr level(99.9)
		
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
		
		mlogit YPG3040_grp male `var', baseoutcome(3) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3040_grp c.male##c.`var', baseoutcome(3) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3040_grp male if `var' != ., baseoutcome(3) rrr level(99.9)
		est store base
		mlogit YPG3040_grp male `var', baseoutcome(3) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3040_grp c.male##c.`var', baseoutcome(3) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
		// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3040_grp ageAt28 `var', baseoutcome(3) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3000 ageAt28 if `var' != ., baseoutcome(3) rrr level(99.9)
		est store base
		mlogit YPG3000 ageAt28 `var', baseoutcome(3) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		mlogit YPG3040_grp ageAt28 male `var', baseoutcome(3) rrr level(99.9)
		
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
		mlogit YPG3040_grp ageAt28 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		
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
		mlogit YPG3040_grp ageAt28 male `var', baseoutcome(3) rrr level(99.9)
		
		local outcome_level = "Other (ref = None)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
				
		// Now for interaction model
		mlogit YPG3040_grp ageAt28 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3040_grp ageAt28 male if `var' != ., baseoutcome(3) rrr level(99.9)
		est store base
		mlogit YPG3040_grp ageAt28 male `var', baseoutcome(3) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3040_grp ageAt28 c.male##c.`var', baseoutcome(3) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,10]
			local lci = res[5,10]
			local uci = res[6,10]
			local p = res[4,10]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
		
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
		
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
		
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Chrstian)
			local outcome_level = "Christian (ref = None)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
		
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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
			mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		
			local outcome_level = "Other (ref = None)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit YPG3040_grp ageAt28 male if `var' != ., baseoutcome(3) rrr level(99.9)
		est store base
		mlogit YPG3040_grp ageAt28 male i.`var', baseoutcome(3) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3040_grp ageAt28 c.male##i.`var', baseoutcome(3) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_relig_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_relig
postclose G1_relig_lr



**************************************************************************************
*** Now to the next RSBB outcome: Attendance at church/place of worship

*** As this is an ordered categorical variable, originally planned to use ordinal regression models. However, as the G0 mother and G0 partner/father data violated the proportional odds assumption (as does the G1 data; see below), will just use multinomial models regardless so that results are comparable.

** For consistency with previous models, will recode so that higher values indicate greater RSBB - Due to small sample sizes of those attending at least once a month and different way question was asked to G1s (it includes the option 'occasionally' between 'at least once a year' and 'not at all'), will code into categories of: min once a month, min once a year, occasionally, and not at all (in contrast to G0 data which is: min once a week, min once a month, min once a year and not at all).
tab YPG3080

recode YPG3080 (1 2 = 3) (3 = 2) (4 = 1) (5 = 0), gen(YPG3080_rev)
label define attend_rev_lb 0 "Not at all" 1 "Occasionally" 2 "Min once a year" 3 "Min once a month"
numlabel attend_rev_lb, add
label values YPG3080_rev attend_rev_lb
tab YPG3080_rev


** Quick test of whether proportional odds assumption been violated in most basic model (with just age at birth). Ah, it has been violated. 
ologit YPG3080_rev ageAt28, or level(99.9)
brant, detail

ologit YPG3080_rev ageAt28 i.IMD, or level(99.9)
brant, detail


** So instead will just run multinomial models with 'not at all' as the baseline/reference category.


*** Now run the loop to save all the results
capture postclose G1_attend
postfile G1_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_attend_results.dta", replace

capture postclose G1_attend_lr
postfile G1_attend_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_attend_results_lr.dta", replace

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'ageAt28' and 'sex' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr level(99.9)
		
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
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr level(99.9)
		
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
		
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
		
		// Interaction between age and sex
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr level(99.9)
		
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
		
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
		
		// Interaction between age and sex
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3080_rev male if `var' != ., baseoutcome(0) rrr level(99.9)
		est store base
		mlogit YPG3080_rev male `var', baseoutcome(0) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3080_rev c.male##c.`var', baseoutcome(0) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3080_rev ageAt28 `var', baseoutcome(0) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3080_rev ageAt28 if `var' != ., baseoutcome(0) rrr level(99.9)
		est store base
		mlogit YPG3080_rev ageAt28 `var', baseoutcome(0) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr level(99.9)
		
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
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr level(99.9)
		
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
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr level(99.9)
		
		local outcome_level = "Min once year (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr level(99.9)
		
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
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr level(99.9)
		
		local outcome_level = "Min once month (ref = Not at all)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,15]
		local lci = res[5,15]
		local uci = res[6,15]
		local p = res[4,15]
				
		// Now for interaction model
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr level(99.9)
		
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
		
		// And finally run the likelihood ratio tests
		mlogit YPG3080_rev ageAt28 male if `var' != ., baseoutcome(0) rrr level(99.9)
		est store base
		mlogit YPG3080_rev ageAt28 male `var', baseoutcome(0) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3080_rev ageAt28 c.male##c.`var', baseoutcome(0) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Married (ref = Never married)"
			}
			else if "`var'" == "parity" {
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
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity" {
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
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
		
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
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
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
		
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
		
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (1/Occasionally)
			local outcome_level = "Occasionally (ref = Not at all)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
		
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once year (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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
			mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		
			local outcome_level = "Min once month (ref = Not at all)"
		
			matrix res = r(table)
			local coef = res[1,35]
			local lci = res[5,35]
			local uci = res[6,35]
			local p = res[4,35]
				
			// Now for interaction model
			mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		
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

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit YPG3080_rev ageAt28 male if `var' != ., baseoutcome(0) rrr level(99.9)
		est store base
		mlogit YPG3080_rev ageAt28 male i.`var', baseoutcome(0) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3080_rev ageAt28 c.male##i.`var', baseoutcome(0) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_attend_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_attend
postclose G1_attend_lr



**************************************************************************************
*** Now to the next RSBB outcome: Intrinsic religiosity

* Intrinsic religiosity outcome is very non-normal (spike at lowesst value '3', then uniform), so will also run multinomial model using categories of this variable (with '3' as baseline category).
tab1 YPG3153 YPG3153_cat
sum YPG3153
hist YPG3153, freq width(1)


*** Start with continuous measure and linear regression

*** Now run the loop to save all the results
capture postclose G1_intrinsic
postfile G1_intrinsic str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_intrinsic_results.dta", replace

capture postclose G1_intrinsic_lr
postfile G1_intrinsic_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_intrinsic_results_lr.dta", replace

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		regress YPG3153 male `var', level(99.9)
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// Interaction between age and sex
		regress YPG3153 c.male##c.`var', level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local sex_main = res[1,1]
		local exp_main = res[1,2]
		
		post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// And finally run the likelihood ratio tests
		regress YPG3153 male if `var' != ., level(99.9)
		est store base
		regress YPG3153 male `var', level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress YPG3153 c.male##c.`var', level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		regress YPG3153 ageAt28 `var', level(99.9)
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
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
		
		post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// And finally run the likelihood ratio tests
		regress YPG3153 ageAt28 if `var' != ., level(99.9)
		est store base
		regress YPG3153 ageAt28 `var', level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		regress YPG3153 ageAt28 male `var', level(99.9)
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,3]
		local lci = res[5,3]
		local uci = res[6,3]
		local p = res[4,3]
		
		// Now for interaction model
		regress YPG3153 ageAt28 c.male##c.`var', level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,4]
		local lci_int = res[5,4]
		local uci_int = res[6,4]
		local p_int = res[4,4]
		local sex_main = res[1,2]
		local exp_main = res[1,3]
		
		post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
		// And finally run the likelihood ratio tests
		regress YPG3153 ageAt28 male if `var' != ., level(99.9)
		est store base
		regress YPG3153 ageAt28 male `var', level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress YPG3153 ageAt28 c.male##c.`var', level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
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
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity" {
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
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
					
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 4)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
						
			
			// Move to the next category of the exposure (category 4)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
							
			// Move to the next category of the exposure (category 5)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local sex_main = res[1,2]
			local exp_main = res[1,7]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 5)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
		
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local sex_main = res[1,2]
			local exp_main = res[1,7]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
							
			
			// Move to the next category of the exposure (category 6)
			regress YPG3153 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)

			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
		
			// Now for interaction model
			regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local sex_main = res[1,2]
			local exp_main = res[1,8]
		
			post G1_intrinsic ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		regress YPG3153 ageAt28 male if `var' != ., level(99.9)
		est store base
		regress YPG3153 ageAt28 male i.`var', level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress YPG3153 ageAt28 c.male##i.`var', level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_intrinsic_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_intrinsic
postclose G1_intrinsic_lr


**** And now repeat for the categorical intrinsic religiosity variable as a sensitivity analysis to ensure results robust and not due to odd distribution of outcome
tab YPG3153_cat

*** Now run the loop to save all the results
capture postclose G1_intrinsic_cat
postfile G1_intrinsic_cat str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_intrinsic_cat_results.dta", replace

capture postclose G1_intrinsic_cat_lr
postfile G1_intrinsic_cat_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_intrinsic_cat_results_lr.dta", replace

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3153_cat male `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/moderate IR)
		local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3153_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local sex_main = res[1,5]
		local exp_main = res[1,6]
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/high IR)
		local outcome_level = "High IR/8-11 (ref = lowest/3)"
		local exp_level = "NA"
		
		mlogit YPG3153_cat male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
		
		// Interaction between age and sex
		mlogit YPG3153_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local sex_main = res[1,9]
		local exp_main = res[1,10]
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// Now onto the next reference category (4/highest IR)
		local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		local exp_level = "NA"
		
		mlogit YPG3153_cat male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
		
		// Interaction between age and sex
		mlogit YPG3153_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local sex_main = res[1,13]
		local exp_main = res[1,14]		
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3153_cat male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3153_cat male `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3153_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3153_cat ageAt28 `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/moderate IR)
		local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
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
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/high IR)
		local outcome_level = "High IR/8-11 (ref = lowest/3)"
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
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (4/highest IR)
		local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
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
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3153_cat ageAt28 if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3153_cat ageAt28 `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		mlogit YPG3153_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/moderate)
		local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
		
		// Now for interaction model
		mlogit YPG3153_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,9]
		local lci_int = res[5,9]
		local uci_int = res[6,9]
		local p_int = res[4,9]
		local sex_main = res[1,7]
		local exp_main = res[1,8]
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/high)
		mlogit YPG3153_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "High IR/8-11 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit YPG3153_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,14]
		local lci_int = res[5,14]
		local uci_int = res[6,14]
		local p_int = res[4,14]
		local sex_main = res[1,12]
		local exp_main = res[1,13]
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (4/highest)
		mlogit YPG3153_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,15]
		local lci = res[5,15]
		local uci = res[6,15]
		local p = res[4,15]
				
		// Now for interaction model
		mlogit YPG3153_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,19]
		local lci_int = res[5,19]
		local uci_int = res[6,19]
		local p_int = res[4,19]
		local sex_main = res[1,17]
		local exp_main = res[1,18]
		
		post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3153_cat ageAt28 male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3153_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3153_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
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
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,16]
			local lci_int = res[5,16]
			local uci_int = res[6,16]
			local p_int = res[4,16]
			local sex_main = res[1,11]
			local exp_main = res[1,13]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,20]
			local exp_main = res[1,22]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local sex_main = res[1,29]
			local exp_main = res[1,31]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
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
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local sex_main = res[1,11]
			local exp_main = res[1,14]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,20]
			local exp_main = res[1,23]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,29]
			local exp_main = res[1,32]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
		
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local sex_main = res[1,13]
			local exp_main = res[1,15]
			
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local sex_main = res[1,24]
			local exp_main = res[1,26]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,35]
			local exp_main = res[1,37]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local sex_main = res[1,13]
			local exp_main = res[1,16]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local sex_main = res[1,24]
			local exp_main = res[1,27]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,35]
			local exp_main = res[1,38]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local sex_main = res[1,13]
			local exp_main = res[1,17]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local sex_main = res[1,24]
			local exp_main = res[1,28]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,35]
			local exp_main = res[1,39]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,28]
			local exp_main = res[1,30]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
			
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,48]
			local lci_int = res[5,48]
			local uci_int = res[6,48]
			local p_int = res[4,48]
			local sex_main = res[1,41]
			local exp_main = res[1,43]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local sex_main = res[1,15]
			local exp_main = res[1,18]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,36]
			local lci_int = res[5,36]
			local uci_int = res[6,36]
			local p_int = res[4,36]
			local sex_main = res[1,28]
			local exp_main = res[1,31]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,49]
			local lci_int = res[5,49]
			local uci_int = res[6,49]
			local p_int = res[4,49]
			local sex_main = res[1,41]
			local exp_main = res[1,44]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local sex_main = res[1,28]
			local exp_main = res[1,32]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,50]
			local lci_int = res[5,50]
			local uci_int = res[6,50]
			local p_int = res[4,50]
			local sex_main = res[1,41]
			local exp_main = res[1,45]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,15]
			local exp_main = res[1,20]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local sex_main = res[1,28]
			local exp_main = res[1,33]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local sex_main = res[1,41]
			local exp_main = res[1,46]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,17]
			local exp_main = res[1,19]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local sex_main = res[1,32]
			local exp_main = res[1,34]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local sex_main = res[1,47]
			local exp_main = res[1,49]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,17]
			local exp_main = res[1,20]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,32]
			local exp_main = res[1,35]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,56]
			local lci_int = res[5,56]
			local uci_int = res[6,56]
			local p_int = res[4,56]
			local sex_main = res[1,47]
			local exp_main = res[1,50]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local sex_main = res[1,17]
			local exp_main = res[1,21]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,32]
			local exp_main = res[1,36]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,57]
			local lci_int = res[5,57]
			local uci_int = res[6,57]
			local p_int = res[4,57]
			local sex_main = res[1,47]
			local exp_main = res[1,51]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local sex_main = res[1,17]
			local exp_main = res[1,22]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,32]
			local exp_main = res[1,37]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,58]
			local lci_int = res[5,58]
			local uci_int = res[6,58]
			local p_int = res[4,58]
			local sex_main = res[1,47]
			local exp_main = res[1,52]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/moderate)
			local outcome_level = "Moderate IR/4-7 (ref = lowest/3)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
		
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local sex_main = res[1,17]
			local exp_main = res[1,23]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/high)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "High IR/8-11 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local sex_main = res[1,32]
			local exp_main = res[1,38]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/highest)
			mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Highest IR/12-15 (ref = lowest/3)"
		
			matrix res = r(table)
			local coef = res[1,35]
			local lci = res[5,35]
			local uci = res[6,35]
			local p = res[4,35]
				
			// Now for interaction model
			mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,59]
			local lci_int = res[5,59]
			local uci_int = res[6,59]
			local p_int = res[4,59]
			local sex_main = res[1,47]
			local exp_main = res[1,53]
		
			post G1_intrinsic_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit YPG3153_cat ageAt28 male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3153_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3153_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_intrinsic_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_intrinsic_cat
postclose G1_intrinsic_cat_lr



**********************************************************************************
*** Next, want to explore extrinsic religiosity - As this is two separate questions (make friends and pray for protection), will have two run two separate loops for each outcome
tab1 YPG3160 YPG3170

** To reduce number of categories and combine low cell counts, will combine 'mildly' and 'strongly' agree and disagree together, to give 4 categories: Agree, not sure, disagree and not applicable. Will have 'agree' as the baseline, which is the most 'extrinsic' response.
recode YPG3160 (1 2 = 1) (3 = 2) (4 5 = 3) (6 = 4), gen(YPG3160_new)

label define extrinsic_lb 1 "Agree" 2 "Not sure" 3 "Disagree" 4 "Not applicable"
numlabel extrinsic_lb, add
label values YPG3160_new extrinsic_lb
tab YPG3160_new

recode YPG3170 (1 2 = 1) (3 = 2) (4 5 = 3) (6 = 4), gen(YPG3170_new)
label values YPG3170_new extrinsic_lb
tab YPG3170_new


*** Start with YPG3160 - Attends church to help them make friends

*** Now run the loop to save all the results
capture postclose G1_extrinsic_friends
postfile G1_extrinsic_friends str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_extrinsic_friends_results.dta", replace

capture postclose G1_extrinsic_friends_lr
postfile G1_extrinsic_friends_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_extrinsic_friends_results_lr.dta", replace

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3160_new male `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3160_new c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local sex_main = res[1,5]
		local exp_main = res[1,6]
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		local outcome_level = "Disagree ER (ref = Agree)"
		local exp_level = "NA"
		
		mlogit YPG3160_new male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
		
		// Interaction between age and sex
		mlogit YPG3160_new c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local sex_main = res[1,9]
		local exp_main = res[1,10]
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// Now onto the next reference category (4/NA)
		local outcome_level = "Not applicable ER (ref = Agree)"
		local exp_level = "NA"
		
		mlogit YPG3160_new male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
		
		// Interaction between age and sex
		mlogit YPG3160_new c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local sex_main = res[1,13]
		local exp_main = res[1,14]
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3160_new male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3160_new male `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3160_new c.male##c.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3160_new ageAt28 `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
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
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		local outcome_level = "Disagree ER (ref = Agree)"
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
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (4/NA)
		local outcome_level = "Not applicable ER (ref = Agree)"
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
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3160_new ageAt28 if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3160_new ageAt28 `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		mlogit YPG3160_new ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
		
		// Now for interaction model
		mlogit YPG3160_new ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,9]
		local lci_int = res[5,9]
		local uci_int = res[6,9]
		local p_int = res[4,9]
		local sex_main = res[1,7]
		local exp_main = res[1,8]
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		mlogit YPG3160_new ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "Disagree ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit YPG3160_new ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,14]
		local lci_int = res[5,14]
		local uci_int = res[6,14]
		local p_int = res[4,14]
		local sex_main = res[1,12]
		local exp_main = res[1,13]
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (4/NA)
		mlogit YPG3160_new ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "Not applicable ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,15]
		local lci = res[5,15]
		local uci = res[6,15]
		local p = res[4,15]
				
		// Now for interaction model
		mlogit YPG3160_new ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,19]
		local lci_int = res[5,19]
		local uci_int = res[6,19]
		local p_int = res[4,19]
		local sex_main = res[1,17]
		local exp_main = res[1,18]
		
		post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3160_new ageAt28 male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3160_new ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3160_new ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
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
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,16]
			local lci_int = res[5,16]
			local uci_int = res[6,16]
			local p_int = res[4,16]
			local sex_main = res[1,11]
			local exp_main = res[1,13]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,20]
			local exp_main = res[1,22]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local sex_main = res[1,29]
			local exp_main = res[1,31]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
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
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local sex_main = res[1,11]
			local exp_main = res[1,14]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,20]
			local exp_main = res[1,23]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,29]
			local exp_main = res[1,32]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
			
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local sex_main = res[1,13]
			local exp_main = res[1,15]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local sex_main = res[1,24]
			local exp_main = res[1,26]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,35]
			local exp_main = res[1,37]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local sex_main = res[1,13]
			local exp_main = res[1,16]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local sex_main = res[1,24]
			local exp_main = res[1,27]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
			
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,35]
			local exp_main = res[1,38]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local sex_main = res[1,13]
			local exp_main = res[1,17]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local sex_main = res[1,24]
			local exp_main = res[1,28]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,35]
			local exp_main = res[1,39]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,28]
			local exp_main = res[1,30]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,48]
			local lci_int = res[5,48]
			local uci_int = res[6,48]
			local p_int = res[4,48]
			local sex_main = res[1,41]
			local exp_main = res[1,43]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local sex_main = res[1,15]
			local exp_main = res[1,18]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,36]
			local lci_int = res[5,36]
			local uci_int = res[6,36]
			local p_int = res[4,36]
			local sex_main = res[1,28]
			local exp_main = res[1,31]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,49]
			local lci_int = res[5,49]
			local uci_int = res[6,49]
			local p_int = res[4,49]
			local sex_main = res[1,41]
			local exp_main = res[1,44]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local sex_main = res[1,28]
			local exp_main = res[1,32]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,50]
			local lci_int = res[5,50]
			local uci_int = res[6,50]
			local p_int = res[4,50]
			local sex_main = res[1,41]
			local exp_main = res[1,45]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,15]
			local exp_main = res[1,20]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local sex_main = res[1,28]
			local exp_main = res[1,33]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local sex_main = res[1,41]
			local exp_main = res[1,46]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,17]
			local exp_main = res[1,19]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local sex_main = res[1,32]
			local exp_main = res[1,34]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local sex_main = res[1,47]
			local exp_main = res[1,49]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,17]
			local exp_main = res[1,20]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,32]
			local exp_main = res[1,35]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,56]
			local lci_int = res[5,56]
			local uci_int = res[6,56]
			local p_int = res[4,56]
			local sex_main = res[1,47]
			local exp_main = res[1,50]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local sex_main = res[1,17]
			local exp_main = res[1,21]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,32]
			local exp_main = res[1,36]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,57]
			local lci_int = res[5,57]
			local uci_int = res[6,57]
			local p_int = res[4,57]
			local sex_main = res[1,47]
			local exp_main = res[1,51]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local sex_main = res[1,17]
			local exp_main = res[1,22]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,32]
			local exp_main = res[1,37]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,58]
			local lci_int = res[5,58]
			local uci_int = res[6,58]
			local p_int = res[4,58]
			local sex_main = res[1,47]
			local exp_main = res[1,52]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
		
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local sex_main = res[1,17]
			local exp_main = res[1,23]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local sex_main = res[1,32]
			local exp_main = res[1,38]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,35]
			local lci = res[5,35]
			local uci = res[6,35]
			local p = res[4,35]
				
			// Now for interaction model
			mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,59]
			local lci_int = res[5,59]
			local uci_int = res[6,59]
			local p_int = res[4,59]
			local sex_main = res[1,47]
			local exp_main = res[1,53]
		
			post G1_extrinsic_friends ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit YPG3160_new ageAt28 male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3160_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3160_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_extrinsic_friends_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_extrinsic_friends
postclose G1_extrinsic_friends_lr



*** And now repeat with YPG3170 - Prays for relief and protection

*** Now run the loop to save all the results
capture postclose G1_extrinsic_prayer
postfile G1_extrinsic_prayer str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_extrinsic_prayer_results.dta", replace

capture postclose G1_extrinsic_prayer_lr
postfile G1_extrinsic_prayer_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_extrinsic_prayer_results_lr.dta", replace

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3170_new male `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3170_new c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local sex_main = res[1,5]
		local exp_main = res[1,6]
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		local outcome_level = "Disagree ER (ref = Agree)"
		local exp_level = "NA"
		
		mlogit YPG3170_new male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
		
		// Interaction between age and sex
		mlogit YPG3170_new c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local sex_main = res[1,9]
		local exp_main = res[1,10]
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// Now onto the next reference category (4/NA)
		local outcome_level = "Not applicable ER (ref = Agree)"
		local exp_level = "NA"
		
		mlogit YPG3170_new male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
		
		// Interaction between age and sex
		mlogit YPG3170_new c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local sex_main = res[1,13]
		local exp_main = res[1,14]
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3170_new male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3170_new male `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3170_new c.male##c.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3170_new ageAt28 `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
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
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		local outcome_level = "Disagree ER (ref = Agree)"
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
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (4/NA)
		local outcome_level = "Not applicable ER (ref = Agree)"
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
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3170_new ageAt28 if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3170_new ageAt28 `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		mlogit YPG3170_new ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/not sure)
		local outcome_level = "Not sure ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
		
		// Now for interaction model
		mlogit YPG3170_new ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,9]
		local lci_int = res[5,9]
		local uci_int = res[6,9]
		local p_int = res[4,9]
		local sex_main = res[1,7]
		local exp_main = res[1,8]
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/disagree)
		mlogit YPG3170_new ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "Disagree ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit YPG3170_new ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,14]
		local lci_int = res[5,14]
		local uci_int = res[6,14]
		local p_int = res[4,14]
		local sex_main = res[1,12]
		local exp_main = res[1,13]
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (4/NA)
		mlogit YPG3170_new ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "Not applicable ER (ref = Agree)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,15]
		local lci = res[5,15]
		local uci = res[6,15]
		local p = res[4,15]
				
		// Now for interaction model
		mlogit YPG3170_new ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,19]
		local lci_int = res[5,19]
		local uci_int = res[6,19]
		local p_int = res[4,19]
		local sex_main = res[1,17]
		local exp_main = res[1,18]
		
		post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3170_new ageAt28 male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3170_new ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3170_new ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
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
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,16]
			local lci_int = res[5,16]
			local uci_int = res[6,16]
			local p_int = res[4,16]
			local sex_main = res[1,11]
			local exp_main = res[1,13]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,20]
			local exp_main = res[1,22]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local sex_main = res[1,29]
			local exp_main = res[1,31]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
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
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local sex_main = res[1,11]
			local exp_main = res[1,14]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,20]
			local exp_main = res[1,23]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,29]
			local exp_main = res[1,32]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
			
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local sex_main = res[1,13]
			local exp_main = res[1,15]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local sex_main = res[1,24]
			local exp_main = res[1,26]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,35]
			local exp_main = res[1,37]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local sex_main = res[1,13]
			local exp_main = res[1,16]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local sex_main = res[1,24]
			local exp_main = res[1,27]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
			
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,35]
			local exp_main = res[1,38]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local sex_main = res[1,13]
			local exp_main = res[1,17]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local sex_main = res[1,24]
			local exp_main = res[1,28]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,35]
			local exp_main = res[1,39]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,28]
			local exp_main = res[1,30]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,48]
			local lci_int = res[5,48]
			local uci_int = res[6,48]
			local p_int = res[4,48]
			local sex_main = res[1,41]
			local exp_main = res[1,43]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local sex_main = res[1,15]
			local exp_main = res[1,18]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,36]
			local lci_int = res[5,36]
			local uci_int = res[6,36]
			local p_int = res[4,36]
			local sex_main = res[1,28]
			local exp_main = res[1,31]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,49]
			local lci_int = res[5,49]
			local uci_int = res[6,49]
			local p_int = res[4,49]
			local sex_main = res[1,41]
			local exp_main = res[1,44]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local sex_main = res[1,28]
			local exp_main = res[1,32]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,50]
			local lci_int = res[5,50]
			local uci_int = res[6,50]
			local p_int = res[4,50]
			local sex_main = res[1,41]
			local exp_main = res[1,45]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,15]
			local exp_main = res[1,20]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local sex_main = res[1,28]
			local exp_main = res[1,33]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local sex_main = res[1,41]
			local exp_main = res[1,46]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,17]
			local exp_main = res[1,19]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local sex_main = res[1,32]
			local exp_main = res[1,34]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local sex_main = res[1,47]
			local exp_main = res[1,49]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,17]
			local exp_main = res[1,20]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,32]
			local exp_main = res[1,35]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,56]
			local lci_int = res[5,56]
			local uci_int = res[6,56]
			local p_int = res[4,56]
			local sex_main = res[1,47]
			local exp_main = res[1,50]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local sex_main = res[1,17]
			local exp_main = res[1,21]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,32]
			local exp_main = res[1,36]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,57]
			local lci_int = res[5,57]
			local uci_int = res[6,57]
			local p_int = res[4,57]
			local sex_main = res[1,47]
			local exp_main = res[1,51]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local sex_main = res[1,17]
			local exp_main = res[1,22]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,32]
			local exp_main = res[1,37]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,58]
			local lci_int = res[5,58]
			local uci_int = res[6,58]
			local p_int = res[4,58]
			local sex_main = res[1,47]
			local exp_main = res[1,52]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/not sure)
			local outcome_level = "Not sure ER (ref = Agree)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
		
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local sex_main = res[1,17]
			local exp_main = res[1,23]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/disagree)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Disagree ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local sex_main = res[1,32]
			local exp_main = res[1,38]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/NA)
			mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "Not applicable ER (ref = Agree)"
		
			matrix res = r(table)
			local coef = res[1,35]
			local lci = res[5,35]
			local uci = res[6,35]
			local p = res[4,35]
				
			// Now for interaction model
			mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,59]
			local lci_int = res[5,59]
			local uci_int = res[6,59]
			local p_int = res[4,59]
			local sex_main = res[1,47]
			local exp_main = res[1,53]
		
			post G1_extrinsic_prayer ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit YPG3170_new ageAt28 male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3170_new ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3170_new ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_extrinsic_prayer_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_extrinsic_prayer
postclose G1_extrinsic_prayer_lr



**************************************************************************************
*** Now finally to the last RSBB outcome: Total DUREL religiosity score

* Total DUREL religiosity outcome is very non-normal (spike at lowesst value '5', then uniform), so will also run multinomial model using categories of this variable (with '5' as baseline category).
tab1 YPG3155 YPG3155_cat
sum YPG3155
hist YPG3155, freq width(1)


*** Start with continuous measure and linear regression

*** Now run the loop to save all the results
capture postclose G1_DUREL
postfile G1_DUREL str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_DUREL_results.dta", replace

capture postclose G1_DUREL_lr
postfile G1_DUREL_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_DUREL_results_lr.dta", replace

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		regress YPG3155 male `var', level(99.9)
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,2]
		local lci = res[5,2]
		local uci = res[6,2]
		local p = res[4,2]
		
		// Interaction between age and sex
		regress YPG3155 c.male##c.`var', level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,3]
		local lci_int = res[5,3]
		local uci_int = res[6,3]
		local p_int = res[4,3]
		local sex_main = res[1,1]
		local exp_main = res[1,2]
		
		post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// And finally run the likelihood ratio tests
		regress YPG3155 male if `var' != ., level(99.9)
		est store base
		regress YPG3155 male `var', level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress YPG3155 c.male##c.`var', level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		regress YPG3155 ageAt28 `var', level(99.9)
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
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
		
		post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// And finally run the likelihood ratio tests
		regress YPG3155 ageAt28 if `var' != ., level(99.9)
		est store base
		regress YPG3155 ageAt28 `var', level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		regress YPG3155 ageAt28 male `var', level(99.9)
		
		local n = e(N)
		
		// Save estimates
		local outcome_level = "NA"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,3]
		local lci = res[5,3]
		local uci = res[6,3]
		local p = res[4,3]
		
		// Now for interaction model
		regress YPG3155 ageAt28 c.male##c.`var', level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,4]
		local lci_int = res[5,4]
		local uci_int = res[6,4]
		local p_int = res[4,4]
		local sex_main = res[1,2]
		local exp_main = res[1,3]
		
		post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
		// And finally run the likelihood ratio tests
		regress YPG3155 ageAt28 male if `var' != ., level(99.9)
		est store base
		regress YPG3155 ageAt28 male `var', level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress YPG3155 ageAt28 c.male##c.`var', level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
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
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,7]
			local lci_int = res[5,7]
			local uci_int = res[6,7]
			local p_int = res[4,7]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maritalStatus" {
				local exp_level = "Wid/Div/Sep (ref = Never married)"
			}
			else if "`var'" == "parity" {
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
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
					
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,8]
			local lci_int = res[5,8]
			local uci_int = res[6,8]
			local p_int = res[4,8]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 4)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,9]
			local lci_int = res[5,9]
			local uci_int = res[6,9]
			local p_int = res[4,9]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
						
			// Move to the next category of the exposure (category 3)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
						
			
			// Move to the next category of the exposure (category 4)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
			
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
							
			// Move to the next category of the exposure (category 5)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local sex_main = res[1,2]
			local exp_main = res[1,7]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
		
			// No reference level for outcome, so set as NA
			local outcome_level = "NA"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,4]
			local lci = res[5,4]
			local uci = res[6,4]
			local p = res[4,4]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,10]
			local lci_int = res[5,10]
			local uci_int = res[6,10]
			local p_int = res[4,10]
			local sex_main = res[1,2]
			local exp_main = res[1,4]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,5]
			local lci = res[5,5]
			local uci = res[6,5]
			local p = res[4,5]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,11]
			local lci_int = res[5,11]
			local uci_int = res[6,11]
			local p_int = res[4,11]
			local sex_main = res[1,2]
			local exp_main = res[1,5]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,6]
			local lci = res[5,6]
			local uci = res[6,6]
			local p = res[4,6]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,12]
			local lci_int = res[5,12]
			local uci_int = res[6,12]
			local p_int = res[4,12]
			local sex_main = res[1,2]
			local exp_main = res[1,6]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 5)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)
		
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,7]
			local lci = res[5,7]
			local uci = res[6,7]
			local p = res[4,7]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,13]
			local lci_int = res[5,13]
			local uci_int = res[6,13]
			local p_int = res[4,13]
			local sex_main = res[1,2]
			local exp_main = res[1,7]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
							
			
			// Move to the next category of the exposure (category 6)
			regress YPG3155 ageAt28 male i.`var', level(99.9)
		
			local n = e(N)

			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,8]
			local lci = res[5,8]
			local uci = res[6,8]
			local p = res[4,8]
		
			// Now for interaction model
			regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,14]
			local lci_int = res[5,14]
			local uci_int = res[6,14]
			local p_int = res[4,14]
			local sex_main = res[1,2]
			local exp_main = res[1,8]
		
			post G1_DUREL ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		regress YPG3155 ageAt28 male if `var' != ., level(99.9)
		est store base
		regress YPG3155 ageAt28 male i.`var', level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		regress YPG3155 ageAt28 c.male##i.`var', level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_DUREL_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_DUREL
postclose G1_DUREL_lr


**** And now repeat for the categorical total DUREL religiosity variable as a sensitivity analysis to ensure results robust and not due to odd distribution of outcome
tab YPG3155_cat

*** Now run the loop to save all the results
capture postclose G1_DUREL_cat
postfile G1_DUREL_cat str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) sex_main exp_main ///
	using ".\G1_Results\G1_DUREL_cat_results.dta", replace

capture postclose G1_DUREL_cat_lr
postfile G1_DUREL_cat_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\G1_Results\G1_DUREL_cat_results_lr.dta", replace

foreach var of varlist ageAt28-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
	
	// Next, how we run the analysis will depend on the type of variable - So need to specify whether variable is continuous/binary (as these can be treated the same), or categorical. Will start with cont/binary variables - Although need to analyse 'age' separately first as will be adjusted for in all other models
	if "`var'" == "ageAt28" {
		mlogit YPG3155_cat male `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/6-10 DUREL)
		local outcome_level = "DUREL 6-10 (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,5]
		local lci = res[5,5]
		local uci = res[6,5]
		local p = res[4,5]
		
		// Interaction between age and sex
		mlogit YPG3155_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,7]
		local lci_int = res[5,7]
		local uci_int = res[6,7]
		local p_int = res[4,7]
		local sex_main = res[1,5]
		local exp_main = res[1,6]
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/11-15 DUREL)
		local outcome_level = "11-15 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		mlogit YPG3155_cat male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,8]
		local lci = res[5,8]
		local uci = res[6,8]
		local p = res[4,8]
		
		// Interaction between age and sex
		mlogit YPG3155_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,11]
		local lci_int = res[5,11]
		local uci_int = res[6,11]
		local p_int = res[4,11]
		local sex_main = res[1,9]
		local exp_main = res[1,10]
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// Now onto the next reference category (4/16-20 DUREL)
		local outcome_level = "16-20 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		mlogit YPG3155_cat male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
		
		// Interaction between age and sex
		mlogit YPG3155_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,15]
		local lci_int = res[5,15]
		local uci_int = res[6,15]
		local p_int = res[4,15]
		local sex_main = res[1,13]
		local exp_main = res[1,14]
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (5/21-26 DUREL)
		local outcome_level = "21-26 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		mlogit YPG3155_cat male `var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef = res[1,14]
		local lci = res[5,14]
		local uci = res[6,14]
		local p = res[4,14]
		
		// Interaction between age and sex
		mlogit YPG3155_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,19]
		local lci_int = res[5,19]
		local uci_int = res[6,19]
		local p_int = res[4,19]
		local sex_main = res[1,17]
		local exp_main = res[1,18]
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3155_cat male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3155_cat male `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3155_cat c.male##c.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the 'sex' variable
	else if "`var'" == "male" {
		mlogit YPG3155_cat ageAt28 `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/6-10 DUREL)
		local outcome_level = "DUREL 6-10 (ref = lowest/5)"
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
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/11-15 DUREL)
		local outcome_level = "11-15 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
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
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// Now onto the next reference category (4/16-20 DUREL)
		local outcome_level = "16-20 DUREL (ref = lowest/5)"
		local exp_level = "NA"

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
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (5/21-26 DUREL)
		local outcome_level = "21-26 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		local coef = res[1,14]
		local lci = res[5,14]
		local uci = res[6,14]
		local p = res[4,14]
		
		// As no interaction model, will fill with blank values
		local coef_int = .
		local lci_int = .
		local uci_int = .
		local p_int = .
		local sex_main = .
		local exp_main = .
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3155_cat ageAt28 if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3155_cat ageAt28 `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// As no interaction model for sex, will just fill with missing value
		local lr_p_int = .
		
		post G1_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
		
	}
	
	// Next, analyse the rest of the continuous/binary variables
	else if "`var'" == "mother_ageAtBirth" | "`var'" == "nonWhiteEthnic" | "`var'" == "rural" | "`var'" == "parent" | "`var'" == "highSocClass" | "`var'" == "income" | "`var'" == "financeDiffs" | "`var'" == "ACEscore_13items" | "`var'" == "ACEscore_10items" | "`var'" == "neighPercept" | "`var'" == "fatherAbsence" | "`var'" == "verbalIQ_age8" | "`var'" == "performanceIQ_age8" | "`var'" == "totalIQ_age8" | "`var'" == "totalIQ_age15" | "`var'" == "digitSymbol_age24" | "`var'" == "vocaba_age24" | "`var'" == "extraversion_age13" | "`var'" == "agreeableness_age13" | "`var'" == "conscientiousness_age13" | "`var'" == "emotionalStab_age13" | "`var'" == "Openess_age13" | "`var'" == "loc_age8" | "`var'" == "loc_age16" | "`var'" == "negCogStyles_age17" | "`var'" == "emoRec_faces_age8" | "`var'" == "emoRec_triangles_age13" | "`var'" == "skuseSocCog_age8" | "`var'" == "skuseSocCog_age16" | "`var'" == "autismSpec_age25" | "`var'" == "SDQ_prosocial_age8" | "`var'" == "SDQ_prosocial_age13" | "`var'" == "SDQ_prosocial_age25" | "`var'" == "esteem_bachman_age17" | "`var'" == "scholasticEsteem_age8" | "`var'" == "globalEsteem_age8" {
		
		mlogit YPG3155_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local n = e(N)
		
		// Start with the first reference category (2/6-10 DUREL)
		local outcome_level = "DUREL 6-10 (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,7]
		local lci = res[5,7]
		local uci = res[6,7]
		local p = res[4,7]
		
		// Now for interaction model
		mlogit YPG3155_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,9]
		local lci_int = res[5,9]
		local uci_int = res[6,9]
		local p_int = res[4,9]
		local sex_main = res[1,7]
		local exp_main = res[1,8]
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (3/11-15 DUREL)
		mlogit YPG3155_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "11-15 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,11]
		local lci = res[5,11]
		local uci = res[6,11]
		local p = res[4,11]
				
		// Now for interaction model
		mlogit YPG3155_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,14]
		local lci_int = res[5,14]
		local uci_int = res[6,14]
		local p_int = res[4,14]
		local sex_main = res[1,12]
		local exp_main = res[1,13]
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (4/16-20 DUREL)
		mlogit YPG3155_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "16-20 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,15]
		local lci = res[5,15]
		local uci = res[6,15]
		local p = res[4,15]
				
		// Now for interaction model
		mlogit YPG3155_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,19]
		local lci_int = res[5,19]
		local uci_int = res[6,19]
		local p_int = res[4,19]
		local sex_main = res[1,17]
		local exp_main = res[1,18]
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		// Now onto the next reference category (5/21-26 DUREL)
		mlogit YPG3155_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		
		local outcome_level = "21-26 DUREL (ref = lowest/5)"
		local exp_level = "NA"
		
		matrix res = r(table)
		local coef = res[1,19]
		local lci = res[5,19]
		local uci = res[6,19]
		local p = res[4,19]
				
		// Now for interaction model
		mlogit YPG3155_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		
		matrix res = r(table)
		local coef_int = res[1,24]
		local lci_int = res[5,24]
		local uci_int = res[6,24]
		local p_int = res[4,24]
		local sex_main = res[1,22]
		local exp_main = res[1,23]
		
		post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
			(`n') (`coef') (`lci') (`uci') (`p') ///
			(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
		
		// And finally run the likelihood ratio tests
		mlogit YPG3155_cat ageAt28 male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3155_cat ageAt28 male `var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3155_cat ageAt28 c.male##c.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
			
	}
	
	// Next, go through the remaining categorical variables and code as needed - In all cases will treat lowest category as reference
	else {
	
		// First, need to know how many categories these vars have, and edit the number of cycles depending on this - As the number of categories shifts the number of columns in the results matrix, need to take variables in turn, depending on number of categories
		quietly distinct `var'
		local cats = r(ndistinct) - 1
		
		// Start with variables that have 2 categories (exc. reference)
		if `cats' == 2 {
		
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
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
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,16]
			local lci_int = res[5,16]
			local uci_int = res[6,16]
			local p_int = res[4,16]
			local sex_main = res[1,11]
			local exp_main = res[1,13]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,20]
			local exp_main = res[1,22]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,34]
			local lci_int = res[5,34]
			local uci_int = res[6,34]
			local p_int = res[4,34]
			local sex_main = res[1,29]
			local exp_main = res[1,31]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,38]
			local exp_main = res[1,40]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
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
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,17]
			local lci_int = res[5,17]
			local uci_int = res[6,17]
			local p_int = res[4,17]
			local sex_main = res[1,11]
			local exp_main = res[1,14]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,20]
			local exp_main = res[1,23]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,29]
			local exp_main = res[1,32]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local sex_main = res[1,38]
			local exp_main = res[1,41]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
						
		}
		
		
		// Now to variables that have 3 categories (exc. reference)
		if `cats' == 3 {
		
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Vocational (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Rent (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.5 to 0.75 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,11]
			local lci = res[5,11]
			local uci = res[6,11]
			local p = res[4,11]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,19]
			local lci_int = res[5,19]
			local uci_int = res[6,19]
			local p_int = res[4,19]
			local sex_main = res[1,13]
			local exp_main = res[1,15]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,18]
			local lci = res[5,18]
			local uci = res[6,18]
			local p = res[4,18]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,30]
			local lci_int = res[5,30]
			local uci_int = res[6,30]
			local p_int = res[4,30]
			local sex_main = res[1,24]
			local exp_main = res[1,26]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,35]
			local exp_main = res[1,37]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,52]
			local lci_int = res[5,52]
			local uci_int = res[6,52]
			local p_int = res[4,52]
			local sex_main = res[1,46]
			local exp_main = res[1,48]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "AS/A level (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Council/HA (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 0.75 to 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,20]
			local lci_int = res[5,20]
			local uci_int = res[6,20]
			local p_int = res[4,20]
			local sex_main = res[1,13]
			local exp_main = res[1,16]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,19]
			local lci = res[5,19]
			local uci = res[6,19]
			local p = res[4,19]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,31]
			local lci_int = res[5,31]
			local uci_int = res[6,31]
			local p_int = res[4,31]
			local sex_main = res[1,24]
			local exp_main = res[1,27]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,35]
			local exp_main = res[1,38]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,53]
			local lci_int = res[5,53]
			local uci_int = res[6,53]
			local p_int = res[4,53]
			local sex_main = res[1,46]
			local exp_main = res[1,49]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "education" {
				local exp_level = "Degree (ref = GCSE/None)"
			}
			else if "`var'" == "housing" {
				local exp_level = "Other (ref = Own/Mortgage)"
			}
			else if "`var'" == "crowding" {
				local exp_level = "> 1 (ref = <= 0.5)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,21]
			local lci_int = res[5,21]
			local uci_int = res[6,21]
			local p_int = res[4,21]
			local sex_main = res[1,13]
			local exp_main = res[1,17]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,32]
			local lci_int = res[5,32]
			local uci_int = res[6,32]
			local p_int = res[4,32]
			local sex_main = res[1,24]
			local exp_main = res[1,28]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,27]
			local lci = res[5,27]
			local uci = res[6,27]
			local p = res[4,27]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,35]
			local exp_main = res[1,39]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,54]
			local lci_int = res[5,54]
			local uci_int = res[6,54]
			local p_int = res[4,54]
			local sex_main = res[1,46]
			local exp_main = res[1,50]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
		}
				
			
		// Now to variables that have 4 categories (exc. reference)
		if `cats' == 4 {
		
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,12]
			local lci = res[5,12]
			local uci = res[6,12]
			local p = res[4,12]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,22]
			local lci_int = res[5,22]
			local uci_int = res[6,22]
			local p_int = res[4,22]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,20]
			local lci = res[5,20]
			local uci = res[6,20]
			local p = res[4,20]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,35]
			local lci_int = res[5,35]
			local uci_int = res[6,35]
			local p_int = res[4,35]
			local sex_main = res[1,28]
			local exp_main = res[1,30]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,28]
			local lci = res[5,28]
			local uci = res[6,28]
			local p = res[4,28]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,48]
			local lci_int = res[5,48]
			local uci_int = res[6,48]
			local p_int = res[4,48]
			local sex_main = res[1,41]
			local exp_main = res[1,43]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,36]
			local lci = res[5,36]
			local uci = res[6,36]
			local p = res[4,36]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,61]
			local lci_int = res[5,61]
			local uci_int = res[6,61]
			local p_int = res[4,61]
			local sex_main = res[1,54]
			local exp_main = res[1,56]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,23]
			local lci_int = res[5,23]
			local uci_int = res[6,23]
			local p_int = res[4,23]
			local sex_main = res[1,15]
			local exp_main = res[1,18]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,21]
			local lci = res[5,21]
			local uci = res[6,21]
			local p = res[4,21]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,36]
			local lci_int = res[5,36]
			local uci_int = res[6,36]
			local p_int = res[4,36]
			local sex_main = res[1,28]
			local exp_main = res[1,31]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,29]
			local lci = res[5,29]
			local uci = res[6,29]
			local p = res[4,29]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,49]
			local lci_int = res[5,49]
			local uci_int = res[6,49]
			local p_int = res[4,49]
			local sex_main = res[1,41]
			local exp_main = res[1,44]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,37]
			local lci = res[5,37]
			local uci = res[6,37]
			local p = res[4,37]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,62]
			local lci_int = res[5,62]
			local uci_int = res[6,62]
			local p_int = res[4,62]
			local sex_main = res[1,54]
			local exp_main = res[1,57]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,24]
			local lci_int = res[5,24]
			local uci_int = res[6,24]
			local p_int = res[4,24]
			local sex_main = res[1,15]
			local exp_main = res[1,17]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,37]
			local lci_int = res[5,37]
			local uci_int = res[6,37]
			local p_int = res[4,37]
			local sex_main = res[1,28]
			local exp_main = res[1,32]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,30]
			local lci = res[5,30]
			local uci = res[6,30]
			local p = res[4,30]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,50]
			local lci_int = res[5,50]
			local uci_int = res[6,50]
			local p_int = res[4,50]
			local sex_main = res[1,41]
			local exp_main = res[1,45]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,38]
			local lci = res[5,38]
			local uci = res[6,38]
			local p = res[4,38]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,63]
			local lci_int = res[5,63]
			local uci_int = res[6,63]
			local p_int = res[4,63]
			local sex_main = res[1,54]
			local exp_main = res[1,58]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "maternalEdu" {
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
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,15]
			local exp_main = res[1,20]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,38]
			local lci_int = res[5,38]
			local uci_int = res[6,38]
			local p_int = res[4,38]
			local sex_main = res[1,28]
			local exp_main = res[1,33]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,51]
			local lci_int = res[5,51]
			local uci_int = res[6,51]
			local p_int = res[4,51]
			local sex_main = res[1,41]
			local exp_main = res[1,46]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,39]
			local lci = res[5,39]
			local uci = res[6,39]
			local p = res[4,39]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,64]
			local lci_int = res[5,64]
			local uci_int = res[6,64]
			local p_int = res[4,64]
			local sex_main = res[1,54]
			local exp_main = res[1,59]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}
		
		
		// Now to variables that have 5 categories (exc. reference)
		if `cats' == 5 {
		
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "1 move (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,13]
			local lci = res[5,13]
			local uci = res[6,13]
			local p = res[4,13]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,25]
			local lci_int = res[5,25]
			local uci_int = res[6,25]
			local p_int = res[4,25]
			local sex_main = res[1,17]
			local exp_main = res[1,19]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,22]
			local lci = res[5,22]
			local uci = res[6,22]
			local p = res[4,22]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,40]
			local lci_int = res[5,40]
			local uci_int = res[6,40]
			local p_int = res[4,40]
			local sex_main = res[1,32]
			local exp_main = res[1,34]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,31]
			local lci = res[5,31]
			local uci = res[6,31]
			local p = res[4,31]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,55]
			local lci_int = res[5,55]
			local uci_int = res[6,55]
			local p_int = res[4,55]
			local sex_main = res[1,47]
			local exp_main = res[1,49]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,40]
			local lci = res[5,40]
			local uci = res[6,40]
			local p = res[4,40]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,70]
			local lci_int = res[5,70]
			local uci_int = res[6,70]
			local p_int = res[4,70]
			local sex_main = res[1,62]
			local exp_main = res[1,64]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 3)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "2 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,14]
			local lci = res[5,14]
			local uci = res[6,14]
			local p = res[4,14]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,26]
			local lci_int = res[5,26]
			local uci_int = res[6,26]
			local p_int = res[4,26]
			local sex_main = res[1,17]
			local exp_main = res[1,20]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,23]
			local lci = res[5,23]
			local uci = res[6,23]
			local p = res[4,23]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,41]
			local lci_int = res[5,41]
			local uci_int = res[6,41]
			local p_int = res[4,41]
			local sex_main = res[1,32]
			local exp_main = res[1,35]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,32]
			local lci = res[5,32]
			local uci = res[6,32]
			local p = res[4,32]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,56]
			local lci_int = res[5,56]
			local uci_int = res[6,56]
			local p_int = res[4,56]
			local sex_main = res[1,47]
			local exp_main = res[1,50]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,41]
			local lci = res[5,41]
			local uci = res[6,41]
			local p = res[4,41]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,71]
			local lci_int = res[5,71]
			local uci_int = res[6,71]
			local p_int = res[4,71]
			local sex_main = res[1,62]
			local exp_main = res[1,65]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			
			// Move to the next category of the exposure (category 4)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "3 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,15]
			local lci = res[5,15]
			local uci = res[6,15]
			local p = res[4,15]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,27]
			local lci_int = res[5,27]
			local uci_int = res[6,27]
			local p_int = res[4,27]
			local sex_main = res[1,17]
			local exp_main = res[1,21]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,24]
			local lci = res[5,24]
			local uci = res[6,24]
			local p = res[4,24]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,42]
			local lci_int = res[5,42]
			local uci_int = res[6,42]
			local p_int = res[4,42]
			local sex_main = res[1,32]
			local exp_main = res[1,36]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,33]
			local lci = res[5,33]
			local uci = res[6,33]
			local p = res[4,33]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,57]
			local lci_int = res[5,57]
			local uci_int = res[6,57]
			local p_int = res[4,57]
			local sex_main = res[1,47]
			local exp_main = res[1,51]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,42]
			local lci = res[5,42]
			local uci = res[6,42]
			local p = res[4,42]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,72]
			local lci_int = res[5,72]
			local uci_int = res[6,72]
			local p_int = res[4,72]
			local sex_main = res[1,62]
			local exp_main = res[1,66]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
				
			// Move to the next category of the exposure (category 5)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "4 moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,16]
			local lci = res[5,16]
			local uci = res[6,16]
			local p = res[4,16]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,28]
			local lci_int = res[5,28]
			local uci_int = res[6,28]
			local p_int = res[4,28]
			local sex_main = res[1,17]
			local exp_main = res[1,22]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,25]
			local lci = res[5,25]
			local uci = res[6,25]
			local p = res[4,25]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,43]
			local lci_int = res[5,43]
			local uci_int = res[6,43]
			local p_int = res[4,43]
			local sex_main = res[1,32]
			local exp_main = res[1,37]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,34]
			local lci = res[5,34]
			local uci = res[6,34]
			local p = res[4,34]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,58]
			local lci_int = res[5,58]
			local uci_int = res[6,58]
			local p_int = res[4,58]
			local sex_main = res[1,47]
			local exp_main = res[1,52]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,43]
			local lci = res[5,43]
			local uci = res[6,43]
			local p = res[4,43]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,73]
			local lci_int = res[5,73]
			local uci_int = res[6,73]
			local p_int = res[4,73]
			local sex_main = res[1,62]
			local exp_main = res[1,67]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			
			// Move to the next category of the exposure (category 6)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local n = e(N)
		
			// Start with the first reference category (2/6-10 DUREL)
			local outcome_level = "DUREL 6-10 (ref = lowest/5)"
			
			// Specify the level of the categorical exposure variable
			if "`var'" == "mobility" {
				local exp_level = "5 + moves (ref = 0 moves)"
			}
			
			// Now extract the relevant statistics		
			matrix res = r(table)
			local coef = res[1,17]
			local lci = res[5,17]
			local uci = res[6,17]
			local p = res[4,17]
		
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,29]
			local lci_int = res[5,29]
			local uci_int = res[6,29]
			local p_int = res[4,29]
			local sex_main = res[1,17]
			local exp_main = res[1,23]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (3/11-15 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "11-15 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,26]
			local lci = res[5,26]
			local uci = res[6,26]
			local p = res[4,26]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,44]
			local lci_int = res[5,44]
			local uci_int = res[6,44]
			local p_int = res[4,44]
			local sex_main = res[1,32]
			local exp_main = res[1,38]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
			// Now onto the next reference category (4/16-20 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "16-20 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,35]
			local lci = res[5,35]
			local uci = res[6,35]
			local p = res[4,35]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,59]
			local lci_int = res[5,59]
			local uci_int = res[6,59]
			local p_int = res[4,59]
			local sex_main = res[1,47]
			local exp_main = res[1,53]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
				
			// Now onto the next reference category (5/21-26 DUREL)
			mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		
			local outcome_level = "21-26 DUREL (ref = lowest/5)"
		
			matrix res = r(table)
			local coef = res[1,44]
			local lci = res[5,44]
			local uci = res[6,44]
			local p = res[4,44]
				
			// Now for interaction model
			mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		
			matrix res = r(table)
			local coef_int = res[1,74]
			local lci_int = res[5,74]
			local uci_int = res[6,74]
			local p_int = res[4,74]
			local sex_main = res[1,62]
			local exp_main = res[1,68]
		
			post G1_DUREL_cat ("`exp'") ("`outcome_level'") ("`exp_level'") ///
				(`n') (`coef') (`lci') (`uci') (`p') ///
				(`coef_int') (`lci_int') (`uci_int') (`p_int') (`sex_main') (`exp_main')
			
		}

		
		// And finally run the likelihood ratio tests for all these categorical exposures
		mlogit YPG3155_cat ageAt28 male if `var' != ., baseoutcome(1) rrr level(99.9)
		est store base
		mlogit YPG3155_cat ageAt28 male i.`var', baseoutcome(1) rrr level(99.9)
		est store main
		
		lrtest base main
		local lr_p_main = r(p)
		
		// And the interaction model
		mlogit YPG3155_cat ageAt28 c.male##i.`var', baseoutcome(1) rrr level(99.9)
		est store inter
		
		lrtest main inter
		local lr_p_int = r(p)
		
		post G1_DUREL_cat_lr ("`exp'") (`lr_p_main') (`lr_p_int')
				
	}
		
}

postclose G1_DUREL_cat
postclose G1_DUREL_cat_lr

* And close the log file
log close



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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 48 exposures, will do 0.05/48)
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Belief in God/divine power - Main effect") ///
	name(belief_main, replace)
	
graph export ".\G1_Results\belief_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Belief in God/divine power - Age interaction") ///
	name(belief_int, replace)
	
graph export ".\G1_Results\belief_ageInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Belief in God/divine power") ///
	legend(order(1 "Main effect" 2 "Age interaction") size(small)) ///
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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 48 exposures, will do 0.05/48)
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Religious affiliation - Main effect") ///
	name(relig_main, replace)
	
graph export ".\G1_Results\relig_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Religious affiliation - Sex interaction") ///
	name(relig_int, replace)
	
graph export ".\G1_Results\relig_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 48 exposures, will do 0.05/48)
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Church attendance - Main effect") ///
	name(attend_main, replace)
	
graph export ".\G1_Results\attend_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Church attendance - Sex interaction") ///
	name(attend_int, replace)
	
graph export ".\G1_Results\attend_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Church attendance") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(attend_both, replace)

graph export ".\G1_Results\attend_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'church attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\G1_Results\attend_pvalues.dta", replace


**** Now read in the next outcome - intrinsic religiosity (will use the multinomial results, as the distribution for the linear model is very non-normal, and using multinomial results makes these p-value plots consistent with the coefficient plots below, which also use the multinomial results)
use ".\G1_Results\G1_intrinsic_cat_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 48 exposures, will do 0.05/48)
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Intrinsic religiosity - Main effect") ///
	name(intrin_main, replace)
	
graph export ".\G1_Results\intrin_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Intrinsic religiosity - Sex interaction") ///
	name(intrin_int, replace)
	
graph export ".\G1_Results\intrin_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Intrinsic religiosity") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(intrin_both, replace)

graph export ".\G1_Results\intrin_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'intrinsic relig' as a variable, then save this file
gen outcome = "Intrinsic relig."
recast str30 outcome
order outcome

save ".\G1_Results\intrin_pvalues.dta", replace


**** Now read in the next outcome - extrinsic religiosity (friends)
use ".\G1_Results\G1_extrinsic_friends_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 48 exposures, will do 0.05/48)
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Extrinsic religiosity (friends) - Main effect") ///
	name(extrinFriends_main, replace)
	
graph export ".\G1_Results\extrinFriends_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Extrinsic religiosity (friends) - Sex interaction") ///
	name(extrinFriends_int, replace)
	
graph export ".\G1_Results\extrinFriends_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Extrinsic religiosity (friends)") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(extrinFriends_both, replace)

graph export ".\G1_Results\extrinFriends_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'extrinsic relig (friends)' as a variable, then save this file
gen outcome = "Extrinsic relig. (friends)"
recast str30 outcome
order outcome

save ".\G1_Results\extrinFriends_pvalues.dta", replace


**** Now read in the next outcome - extrinsic religiosity (prayer)
use ".\G1_Results\G1_extrinsic_prayer_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 48 exposures, will do 0.05/48)
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Extrinsic religiosity (prayer) - Main effect") ///
	name(extrinPray_main, replace)
	
graph export ".\G1_Results\extrinPray_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Extrinsic religiosity (prayer) - Sex interaction") ///
	name(extrinPray_int, replace)
	
graph export ".\G1_Results\extrinPray_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///, ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Extrinsic religiosity (prayer)") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(extrinPray_both, replace)

graph export ".\G1_Results\extrinPray_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'extrinsic relig (prayer)' as a variable, then save this file
gen outcome = "Extrinsic relig. (prayer)"
recast str30 outcome
order outcome

save ".\G1_Results\extrinPrayer_pvalues.dta", replace


**** Now read in the final outcome - Total DUREL religiosity score (will use the multinomial results, as the distribution for the linear model is very non-normal, and using multinomial results makes these p-value plots consistent with the coefficient plots below, which also use the multinomial results)
use ".\G1_Results\G1_DUREL_cat_results_lr.dta", clear

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

* Display two thresholds; standard 0.05 and Bonferroni-corrected one (as 48 exposures, will do 0.05/48)
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Total DUREL religiosity score - Main effect") ///
	name(durel_main, replace)
	
graph export ".\G1_Results\durel_mainEffect_pvalues.pdf", replace
	
* And repeat for interaction effect (exclude 'sex' here, as can't interact with itself!)
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if exp_num != 1, ///
		col(black) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Total DUREL religiosity score - Sex interaction") ///
	name(durel_int, replace)
	
graph export ".\G1_Results\durel_sexInteraction_pvalues.pdf", replace
	
** Combine these results on the same plot
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main, col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if exp_num != 1, ///
		col(red) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Total DUREL religiosity score") ///
	legend(order(1 "Main effect" 2 "Sex interaction") size(small)) ///
	name(durel_both, replace)

graph export ".\G1_Results\durel_mainAndInt_pvalues.pdf", replace

graph close _all
	
** Add 'DUREL' as a variable, then save this file
gen outcome = "DUREL"
recast str30 outcome
order outcome

save ".\G1_Results\durel_pvalues.dta", replace


*** Combine all these datasets together
use ".\G1_Results\belief_pvalues.dta", clear
append using ".\G1_Results\relig_pvalues.dta"
append using ".\G1_Results\attend_pvalues.dta"
append using ".\G1_Results\intrin_pvalues.dta"
append using ".\G1_Results\extrinFriends_pvalues.dta"
append using ".\G1_Results\extrinPrayer_pvalues.dta"
append using ".\G1_Results\extrinPrayer_pvalues.dta"
append using ".\G1_Results\durel_pvalues.dta"


** Now look at combined results - Even though here the data come from the same time-point, for readability reasons (and to match parental data), will compare belief/religious affiliation/church attendance in one plot, and intrinsic/extrinsic/total religiosity in another plot.

* Belief/religion/church vars main effects
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main if outcome == "Belief", ///
		col(black) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_main if outcome == "Religious affil.", ///
		col(red) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_main if outcome == "Church attendance", ///
		col(blue) msize(vsmall) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Belief, religion and attendance - Main effects") ///
	legend(order(1 "Belief in God" 2 "Religious affiliation" ///
		3 "Church attendance") rows(1) size(small)) ///
	name(belRelCh_main, replace)

graph export ".\G1_Results\beliefReligAttend_mainEffects_pvalues.pdf", replace

* Belief/religion/church vars interaction effects
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if outcome == "Belief" & exp_num != 1, ///
		col(black) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Religious affil." & exp_num != 1, ///
		col(red) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Church attendance" & exp_num != 1, ///
		col(blue) msize(vsmall) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Belief, religion and attendance - Sex interaction") ///
	legend(order(1 "Belief in God" 2 "Religious affiliation" ///
		3 "Church attendance") rows(1) size(small)) ///
	name(belRelCh_int, replace)

graph export ".\G1_Results\beliefReligAttend_sexInt_pvalues.pdf", replace

* intrinsic, extrincic and total religiosity vars main effects
local bon_thresh = -log10(0.05/48)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_main if outcome == "Intrinsic relig.", ///
		col(black) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_main if outcome == "Extrinsic relig. (friends)", ///
		col(red) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_main if outcome == "Extrinsic relig. (prayer)", ///
		col(blue) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_main if outcome == "DUREL", ///
		col(green) msize(vsmall) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Intrinsic, extrinsic and DUREL - Main effects") ///
	legend(order(1 "Intrinsic religiosity" 2 "Extrinsic religiosity (friends)" ///
		3 "Extrinsic religiosity (prayer)" 4 "DUREL religiosity") ///
		rows(2) size(small)) ///
	name(intExtTot_main, replace)

graph export ".\G1_Results\intExtDUREL_mainEffects_pvalues.pdf", replace

* intrinsic, extrincic and total religiosity vars interaction effects
local bon_thresh = -log10(0.05/47)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if outcome == "Intrinsic relig." & exp_num != 1, ///
		col(black) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Extrinsic relig. (friends)" & exp_num != 1, ///
		col(red) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Extrinsic relig. (prayer)" & exp_num != 1, ///
		col(blue) msize(vsmall) msym(D)) ///
	(scatter exp_num logp_int if outcome == "DUREL" & exp_num != 1, ///
		col(green) msize(vsmall) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(2(1)48, valuelabel labsize(tiny) angle(0)) ///
	title("Intrinsic, extrinsic and DUREL - Sex interaction") ///
	legend(order(1 "Intrinsic religiosity" 2 "Extrinsic religiosity (friends)" ///
		3 "Extrinsic religiosity (prayer)" 4 "DUREL religiosity") ///
		rows(2) size(small)) ///
	name(intExtTot_int, replace)

graph export ".\G1_Results\intExtDUREL_sexInt_pvalues.pdf", replace


** Combine all these graphs together
graph combine belRelCh_main belRelCh_int intExtTot_main intExtTot_int, imargin(0 0 0 0) iscale(0.5)

graph export ".\G1_Results\allData_pvalues.pdf", replace

graph close _all


*** Next, want to plot some of the actual results

** Will read the datasets in then combine together into one single dataset
use ".\G1_Results\G1_belief_results.dta", clear
gen outcome = "Belief"

append using ".\G1_Results\G1_relig_results.dta"
replace outcome = "Relig" if outcome == ""
tab outcome, m

append using ".\G1_Results\G1_attend_results.dta"
replace outcome = "Attend" if outcome == ""
tab outcome, m

append using ".\G1_Results\G1_intrinsic_results.dta"
replace outcome = "Intrinsic" if outcome == ""
tab outcome, m

append using ".\G1_Results\G1_intrinsic_cat_results.dta"
replace outcome = "Intrinsic (cat)" if outcome == ""
tab outcome, m

append using ".\G1_Results\G1_extrinsic_friends_results.dta"
replace outcome = "Extrinsic - friends" if outcome == ""
tab outcome, m

append using ".\G1_Results\G1_extrinsic_prayer_results.dta"
replace outcome = "Extrinsic - prayer" if outcome == ""
tab outcome, m

append using ".\G1_Results\G1_DUREL_results.dta"
replace outcome = "DUREL" if outcome == ""
tab outcome, m

append using ".\G1_Results\G1_DUREL_cat_results.dta"
replace outcome = "DUREL (cat)" if outcome == ""
tab outcome, m


** First, make a plot for the age results - As some outcomes are on different scales, will just use the results from multinomial regression for all outcomes (inc. intrinsic and total/DUREL religiosity, even though also ran linear regressions on these as well) - Having all variables on the same plot on the same scale makes things easier to visualise.
capture drop level_num
gen level_num = 0
replace level_num = 1 if outcome_level == "Not sure (ref = No)"
replace level_num = 3 if outcome_level == "Christian (ref = None)"
replace level_num = 4 if outcome_level == "Other (ref = None)"
replace level_num = 6 if outcome_level == "Occasionally (ref = Not at all"
replace level_num = 7 if outcome_level == "Min once year (ref = Not at al"
replace level_num = 8 if outcome_level == "Min once month (ref = Not at a"
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

label define level_lb 0 "Belief in God - Yes (ref = No)" 1 "Belief in God - Not sure (ref = No)" 3 "Religious affiliation - Christian (ref = None)" 4 "Religious affiliation - Other (ref = None)" 6 "Church attendance - Occasionally (ref = Not at all)" 7 "Church attendance - Min once a year (ref = Not at all)" 8 "Church attendance - Min once a month (ref = Not at all)" 10 "Intrinsic religiosity - Moderate (4-7 score; ref = lowest/3)" 11 "Intrinsic religiosity - High (8-11 score; ref = lowest/3)" 12 "Intrinsic religiosity - Highest (12-15 score; ref = lowest/3)" 14 "Extrinsic religiosity (friends) - Not sure (ref = Agree)" 15 "Extrinsic religiosity (friends) - Disagree (ref = Agree)" 16 "Extrinsic religiosity (friends) - Not applicable (ref = Agree)" 18 "Extrinsic religiosity (prayer) - Not sure (ref = Agree)" 19 "Extrinsic religiosity (prayer) - Disagree (ref = Agree)" 20 "Extrinsic religiosity (prayer) - Not applicable (ref = Agree)" 22 "DUREL total religiosity - 6-10 score (ref = lowest/5)" 23 "DUREL total religiosity - 11-15 score (ref = lowest/5)" 24 "DUREL total religiosity - 16-20 score (ref = lowest/5)" 25 "DUREL total religiosity - 21-26 score (ref = lowest/5)", replace
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
		
graph export ".\G1_Results\ageResults.pdf", replace


** Create plot for ethnicity (ref = white)

* Min and max x-axis values
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
		xlabel(0.2 0.3 0.5 1 2 3 5 10, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(ethnic_cat, replace)
		
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
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "male", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "male", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "male", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "male", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "male", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "male", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "male", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "male", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = female)") ///
		title("Male sex and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.4 0.6 0.8 1 1.5 2 3 4, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(sex_cat, replace)
		
graph export ".\G1_Results\sexResults.pdf", replace


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
		xlabel(0.02 0.05 0.1 0.2 0.5 1 2 5 10 20, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Married" 15 "Widowed/Divorced/Separated")) ///
		name(marital_cat, replace)
		
graph export ".\G1_Results\maritalStatusResults.pdf", replace


** Create plot for maternal education (ref = CSE/None)

* As four exposure levels, need to split these up
capture drop level_split
gen level_split = level_num - 0.3 if exposure == "maternalEdu" & exp_level == "Vocational (ref = CSE/None)"
replace level_split = level_num - 0.1 if exposure == "maternalEdu" & exp_level == "O-level (ref = CSE/None)"
replace level_split = level_num + 0.1 if exposure == "maternalEdu" & exp_level == "A-level (ref = CSE/None)"
replace level_split = level_num + 0.3 if exposure == "maternalEdu" & exp_level == "Degree (ref = CSE/None)"
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
		title("Maternal education and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 40, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Vocational" 15 "O-levels" 29 "A-levels" ///
			43 "Degree") rows(1)) ///
		name(matedu_cat, replace)
		
graph export ".\G1_Results\matEduResults.pdf", replace


** Create plot for maternal education (ref = CSE/None)

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
			"Vocational (ref = GCSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Vocational (ref = GCSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Vocational (ref = GCSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Vocational (ref = GCSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Vocational (ref = GCSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Vocational (ref = GCSE/None)", col(black) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Vocational (ref = GCSE/None)", horizontal col(black) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "AS/A level (ref = GCSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "AS/A level (ref = GCSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "AS/A level (ref = GCSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "AS/A level (ref = GCSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"AS/A level (ref = GCSE/None)", col(red) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"AS/A level (ref = GCSE/None)", horizontal col(red) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Belief" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Belief" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Relig" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Relig" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Attend" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Attend" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Intrinsic (cat)" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Intrinsic (cat)" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - friends" & exp_level ///
			== "Degree (ref = GCSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - friends" & exp_level ///
			== "Degree (ref = GCSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "Extrinsic - prayer" & exp_level ///
			== "Degree (ref = GCSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "Extrinsic - prayer" & exp_level ///
			== "Degree (ref = GCSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)) ///
		(scatter level_split coef if outcome == "DUREL (cat)" & exp_level == ///
			"Degree (ref = GCSE/None)", col(blue) msize(tiny) msym(D)) ///
		(rspike lci uci level_split if outcome == "DUREL (cat)" & exp_level == ///
			"Degree (ref = GCSE/None)", horizontal col(blue) msize(vtiny) lwidth(thin)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = GCSE/None)") ///
		title("Education and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash) lwidth(thin)) xscale(log) ///
		xlabel(0.06 0.1 0.2 0.5 1 2 5 10, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Vocational" 15 "AS/A levels" 29 "Degree") rows(1)) ///
		name(edu_cat, replace)
		
graph export ".\G1_Results\eduResults.pdf", replace
	

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
			"parent", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "parent", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "parent", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "parent", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "parent", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "parent", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "parent", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "parent", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "parent", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = Not a parent") ///
		title("Parental status and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.35 0.5 0.7 1 1.5 2 3 4, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(parent_cat, replace)
		
graph export ".\G1_Results\parentResults.pdf", replace


** Create plot for being father absence (ref = no father absence)
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
			"fatherAbsence", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "fatherAbsence", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "fatherAbsence", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "fatherAbsence", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "fatherAbsence", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "fatherAbsence", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "fatherAbsence", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "fatherAbsence", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "fatherAbsence", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (ref = No FA") ///
		title("Father absence and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.25 0.4 0.7 1 1.5 2 3 4 6, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(FA_cat, replace)
		
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
			"mother_ageAtBirth", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "mother_ageAtBirth", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "mother_ageAtBirth", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "mother_ageAtBirth", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "mother_ageAtBirth", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "mother_ageAtBirth", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "mother_ageAtBirth", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "mother_ageAtBirth", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "mother_ageAtBirth", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase)") ///
		title("Maternal age at birth and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.88 0.92 0.96 1 1.04 1.08 1.12, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(mumAge_cat, replace)
		
graph export ".\G1_Results\mumAgeResults.pdf", replace


** Create plot for total IQ at age 8
sum lci uci if exposure == "totalIQ_age8" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"totalIQ_age8", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"totalIQ_age8", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"totalIQ_age8", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"totalIQ_age8", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"totalIQ_age8", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"totalIQ_age8", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "totalIQ_age8", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "totalIQ_age8", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "totalIQ_age8", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "totalIQ_age8", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "totalIQ_age8", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "totalIQ_age8", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "totalIQ_age8", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "totalIQ_age8", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in IQ)") ///
		title("Total IQ Age 8 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.92 0.94 0.96 0.98 1 1.02 1.04, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(iq8_cat, replace)
		
graph export ".\G1_Results\iq8Results.pdf", replace


** Create plot for total IQ at age 15
sum lci uci if exposure == "totalIQ_age15" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "totalIQ_age15", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in IQ)") ///
		title("Total IQ Age 15 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.92 0.94 0.96 0.98 1 1.02 1.04 1.06, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(iq15_cat, replace)
		
graph export ".\G1_Results\iq15Results.pdf", replace


** Create plot for conscientiousness at age 13
sum lci uci if exposure == "conscientiousness_age13" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"conscientiousness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"conscientiousness_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"conscientiousness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"conscientiousness_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"conscientiousness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"conscientiousness_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "conscientiousness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "conscientiousness_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "conscientiousness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "conscientiousness_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "conscientiousness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "conscientiousness_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "conscientiousness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "conscientiousness_age13", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase)") ///
		title("Conscientiousness Age 13 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.88 0.9 0.92 0.94 0.96 0.98 1 1.02 1.04 1.06 1.08 1.1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(cons13_cat, replace)
		
graph export ".\G1_Results\cons13Results.pdf", replace


** Create plot for emotional recognition triangles task at age 13
sum lci uci if exposure == "emoRec_triangles_age13" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == ///
			"emoRec_triangles_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == ///
			"emoRec_triangles_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == ///
			"emoRec_triangles_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == ///
			"emoRec_triangles_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == ///
			"emoRec_triangles_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == ///
			"emoRec_triangles_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Intrinsic (cat)" & ///
			exposure == "emoRec_triangles_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "emoRec_triangles_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - friends" & ///
			exposure == "emoRec_triangles_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - friends" & ///
			exposure == "emoRec_triangles_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "Extrinsic - prayer" & ///
			exposure == "emoRec_triangles_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "emoRec_triangles_age13", horizontal col(black)) ///
		(scatter level_num coef if outcome == "DUREL (cat)" & ///
			exposure == "emoRec_triangles_age13", col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "DUREL (cat)" & ///
			exposure == "emoRec_triangles_age13", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase)") ///
		title("Emo. Recog. (triangles task) Age 13 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.88 0.9 0.92 0.94 0.96 0.98 1 1.02 1.04 1.06 1.08, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(tri13_cat, replace)
		
graph export ".\G1_Results\tri13Results.pdf", replace


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
		xlabel(0.2 0.3 0.5 1 2 3 5 8, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "2" 15 "3" 29 "4" 43 "5/Most dep.") rows(1)) ///
		name(imd_cat, replace)
		
graph export ".\G1_Results\imdResults.pdf", replace


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
		xlabel(0.02 0.05 0.2 0.5 1 2 5 10 20, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "Rented" 15 "Council/HA" 29 "Other") ///
			rows(1)) ///
		name(housing_cat, replace)
		
graph export ".\G1_Results\housingResults.pdf", replace


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
		xlabel(0.005 0.03 0.1 0.3 1 2 5 10 15, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(order(1 "1" 15 "2" 29 "3" 43 "4" 57 "5+") rows(1)) ///
		name(mobility_cat, replace)

graph export ".\G1_Results\mobilityResults.pdf", replace

graph close _all


****** And now for some interaction plots

** Create plot for IQ 15 by age interaction
sum lci_int uci_int if exposure == "totalIQ_age15" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == ///
			"totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == ///
			"totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == ///
			"totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == ///
			"totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == ///
			"totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == ///
			"totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Intrinsic (cat)" & ///
			exposure == "totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - friends" & ///
			exposure == "totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - friends" & ///
			exposure == "totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - prayer" & ///
			exposure == "totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "totalIQ_age15", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "DUREL (cat)" & ///
			exposure == "totalIQ_age15", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "DUREL (cat)" & ///
			exposure == "totalIQ_age15", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in IQ)") ///
		title("IQ age 15*Sex Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.88 0.92 0.96 1 1.04 1.08, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(iq15_int, replace)
		
graph export ".\G1_Results\iq15Results_int.pdf", replace


** Create plot for agreeableness age 13 by age interaction
sum lci_int uci_int if exposure == "agreeableness_age13" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == ///
			"agreeableness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == ///
			"agreeableness_age13", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == ///
			"agreeableness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == ///
			"agreeableness_age13", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == ///
			"agreeableness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == ///
			"agreeableness_age13", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Intrinsic (cat)" & ///
			exposure == "agreeableness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Intrinsic (cat)" & ///
			exposure == "agreeableness_age13", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - friends" & ///
			exposure == "agreeableness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - friends" & ///
			exposure == "agreeableness_age13", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Extrinsic - prayer" & ///
			exposure == "agreeableness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Extrinsic - prayer" & ///
			exposure == "agreeableness_age13", horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "DUREL (cat)" & ///
			exposure == "agreeableness_age13", col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "DUREL (cat)" & ///
			exposure == "agreeableness_age13", horizontal col(black)), ///
		yscale(reverse)	ytitle("") ///
		xtitle("Relative risk ratio (per unit increase in agree.)") ///
		title("Agreeableness*Sex Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.7 0.8 0.9 1 1.1 1.2 1.3 1.4, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(agree13_int, replace)
		
graph export ".\G1_Results\agree13Results_int.pdf", replace


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
		title("Household income*Sex Interaction and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.03 0.1 0.3 0.5 1 2 4 6, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8 10 11 12 14 15 16 18 19 20 22 23 24 25 ///
		, valuelabel labsize(vsmall) ///
		angle(0)) legend(off) name(income_int, replace)
		
graph export ".\G1_Results\incomeResults_int.pdf", replace

graph close _all


** For the multinomial regression results, as interpretation not intuitive, could convert to predicted probabilities using the 'margins' command? (see: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/) - Have started this in the G0 mothers file; see there for example code.



