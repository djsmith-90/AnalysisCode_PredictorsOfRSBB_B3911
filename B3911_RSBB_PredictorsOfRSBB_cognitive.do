*** Predictors of RSBB (B3911) - Cognitive/psychological variables analysis script
*** Created 8/2/2022 by Dan Smith
*** Stata v17.0

*** This script reads in the cleaned G0 mother, G0 partner/father and G1 child data, explores associations between exposures, and then conducts an analysis exploring how all the cognitive/psychological variables are associated with various facets of RSBB for each cohort.


**********************************************************************************
**** Set working directory and start a log file

cd "X:\Groups\ARC\DanS\Descriptive_PredictorsOfRSBB_B3911"

capture log close
log using ".\Cognitive_Results\Desc_RSBB_B3911_CognitiveAnalysis_log", replace text


** As want to use the number of distinct values in a variable later in the script, need to install the user-written 'distinct' package (type 'search distinct' and install package 'dm0042_2')

** Also want to install the 'heatplot' package for displaying correlation plots (plus some dependencies)
*ssc install heatplot, replace
*ssc install palettes, replace
*ssc install colrspace, replace

** Install 'grc1leg' to merge plots together with a single legend
*ssc install grc1leg, replace



**********************************************************************************
**** Start with the G0 mother's analysis
use "G0Mother_PredictorsOfRSBB_B3911.dta", clear

* Drop if not enrolled during pregnancy
drop if in_core == 2

*** Descriptive statistics

** Put the RSBB variables at the start of the dataset
order aln d810 d813 d813_grp d816 Y3000 Y3040 Y3040_grp Y3080 Y3080_OccNever Y3080_OccYr Y3153 Y3153_cat Y3160 Y3170 Y3155 Y3155_cat

* Keep just pregnancy RSBB variables (belief in God, religious affiliation and frequency of church attendance)
drop Y3000 Y3040 Y3040_grp Y3080 Y3080_OccNever Y3080_OccYr Y3153 Y3153_cat Y3160 Y3170 Y3155 Y3155_cat

* Keep just the cognitive/psychological variables (and age)
keep aln-mz028b logic_mem-fom_cog_factor1 b916-b921 d842 h151b

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
misstable sum d810-d816, all

* Exposures
misstable sum mz028b-fom_cog_factor1, all


** Also want to get descriptive stats for each exposure split by each RSBB outcome category
foreach var of varlist mz028b-fom_cog_factor1 {
	quietly distinct `var'
	local unique = r(ndistinct)
	
	display ""
	display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	display "Variable " "`var'" " has " `unique' " values."
	
	if `unique' < 10 {
		tab d810 `var', col
		tab d813_grp `var', col
		tab d816 `var', col
	}
	else {
		tab d810, sum(`var')
		tab d813_grp, sum(`var')
		tab d816, sum(`var')
	}
}


**********************************************************************************
*** Correlations between exposures

** Explore correlations between exposures to see how inter-related these factors are
desc mz028b-fom_cog_factor1

* All exposures are continuous, making this easier

* Order variables in order to be analysed
order aln-d816 mz028b ///
logic_mem-fom_cog_factor1 b916-b921 d842 h151b

* Next, rename all these exposures so they are more intuitive and can be read easier on the correlation heatmaps below
rename mz028b ageAtBirth
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

* Associations between cognitive/psychological variables - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
spearman logicMemory-selfEsteem, pw

matrix cor_cog = r(Rho)

heatplot cor_cog, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))

* Save heatmap
graph export ".\Cognitive_Results\G0Mother_corr_heatplot.pdf", replace

* Save matrix as Excel file
putexcel set ".\Cognitive_Results\G0Mother_corrMatrix.xlsx", replace
putexcel A1=matrix(cor_cog), names


************************************************************************************
*** Next, we want to run the actual analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab d810, m

** We want to store both estimates adjusting for age, and also the interaction between age and the exposure, to see whether it's moderated by age. As all cognitive variables are continuous, this makes the script much simpler than for the demographic and SEP script, which includes categorical variables.

** Create a postfile to post results to, then start the loop - Will create three postfiles; one for coefficients and CIs, another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model), and then a third for the (pseduo-)R2 values to assess model fit - NOTE: Have to store pvalues as 'double' format, else really tiny p-values get coded as '0' (as default if float format, which has minimum value of -3.40282346639e+38 [https://blog.stata.com/2012/04/02/the-penultimate-guide-to-precision/]).
capture postclose mother_belief
postfile mother_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\Cognitive_Results\G0Mother_belief_results.dta", replace

capture postclose mother_belief_lr
postfile mother_belief_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G0Mother_belief_results_lr.dta", replace
	
capture postclose mother_belief_r2
postfile mother_belief_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G0Mother_belief_results_r2.dta", replace

foreach var of varlist logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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
	using ".\Cognitive_Results\G0Mother_relig_results.dta", replace

capture postclose mother_relig_lr
postfile mother_relig_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G0Mother_relig_results_lr.dta", replace
	
capture postclose mother_relig_r2
postfile mother_relig_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G0Mother_relig_results_r2.dta", replace

foreach var of varlist logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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

ologit d816_rev ageAtBirth intel_factor, or
brant, detail


** So instead will just run multinomial models with 'not at all' as the baseline/reference category (this also means all outcomes are using the same model, making them easier to compare).


*** Now run the loop to save all the results
capture postclose mother_attend
postfile mother_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\Cognitive_Results\G0Mother_attend_results.dta", replace

capture postclose mother_attend_lr
postfile mother_attend_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G0Mother_attend_results_lr.dta", replace
	
capture postclose mother_attend_r2
postfile mother_attend_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G0Mother_attend_results_r2.dta", replace

foreach var of varlist logicMemory-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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

postclose mother_attend
postclose mother_attend_lr
postclose mother_attend_r2


* Save this file if we want to use it later
save ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", replace



**********************************************************************************
**** Now onto the G0 partner/fathers's analysis
use "G0Partner_PredictorsOfRSBB_B3911.dta", clear

* Drop if mother not enrolled during pregnancy
drop if in_core == 2

*** Descriptive statistics

** Put the RSBB variables at the start of the dataset
order aln pb150 pb153 pb153_grp pb155 FC3000 FC3040 FC3040_grp FC3080 FC3080_OccNever FC3080_OccYr FC3153 FC3153_cat FC3160 FC3170 FC3155 FC3155_cat

* Keep just pregnancy RSBB variables (belief in God, religious affiliation and frequency of church attendance)
drop FC3000 FC3040 FC3040_grp FC3080 FC3080_OccNever FC3080_OccYr FC3153 FC3153_cat FC3160 FC3170 FC3155 FC3155_cat

* Keep just the cognitive/psychological variables (and age)
keep aln-pb155 pb910 pb546-pb551 pa782 esteem_prorated

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
foreach var of varlist pb910 pb546-pb551 pa782 esteem_prorated {
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
misstable sum pb150-pb155, all

* Exposures
misstable sum pb910 pb546-pb551 pa782 esteem_prorated, all


** Also want to get descriptive stats for each exposure split by each RSBB outcome category
foreach var of varlist pb910 pb546-pb551 pa782 esteem_prorated {
	quietly distinct `var'
	local unique = r(ndistinct)
	
	display ""
	display "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	display "Variable " "`var'" " has " `unique' " values."
	
	if `unique' < 10 {
		tab pb150 `var', col
		tab pb153_grp `var', col
		tab pb155 `var', col
	}
	else {
		tab pb150, sum(`var')
		tab pb153_grp, sum(`var')
		tab pb155, sum(`var')
	}
}



**********************************************************************************
*** Correlations between exposures

** Explore correlations between exposures to see how inter-related these factors are
desc pb910 pb546-pb551 pa782 esteem_prorated

* All exposures are continuous, making this easier

* Order variables in order to be analysed
order aln-pb155 pb910 ///
pb546-pb551 pa782 esteem_prorated

* Next, rename all these exposures so they are more intuitive and can be read easier on the correlation heatmaps below
rename pb910 ageInPreg
rename pb546 IPSM_interpAware
rename pb547 IPSM_approval
rename pb548 IPSM_sepAnx
rename pb549 IPSM_timidity
rename pb550 IPSM_fragility
rename pb551 IPSM_total
rename pa782 LoC_external
rename esteem_prorated selfEsteem

* Associations between cognitive/psychological variables - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
spearman IPSM_interpAware-selfEsteem, pw

matrix cor_cog = r(Rho)

heatplot cor_cog, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45)) legend(subtitle(""))

* Save heatmap
graph export ".\Cognitive_Results\G0Partner_corr_heatplot.pdf", replace

* Save matrix as Excel file
putexcel set ".\Cognitive_Results\G0Partner_corrMatrix.xlsx", replace
putexcel A1=matrix(cor_cog), names


************************************************************************************
*** Next, we want to run the actual analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab pb150, m

** We want to store both estimates adjusting for age, and also the interaction between age and the exposure, to see whether it's moderated by age. As all cognitive variables are continuous, this makes the script much simpler than for the demographic and SEP script, which includes categorical variables.

** Create a postfile to post results to, then start the loop - Will create three postfiles; one for coefficients and CIs, another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model), and then a third for the (pseduo-)R2 values to assess model fit - NOTE: Have to store pvalues as 'double' format, else really tiny p-values get coded as '0' (as default if float format, which has minimum value of -3.40282346639e+38 [https://blog.stata.com/2012/04/02/the-penultimate-guide-to-precision/]).
capture postclose partner_belief
postfile partner_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\Cognitive_Results\G0Partner_belief_results.dta", replace

capture postclose partner_belief_lr
postfile partner_belief_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G0Partner_belief_results_lr.dta", replace
	
capture postclose partner_belief_r2
postfile partner_belief_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G0Partner_belief_results_r2.dta", replace

foreach var of varlist IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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

postclose partner_belief
postclose partner_belief_lr
postclose partner_belief_r2


************************************************************************************
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
	using ".\Cognitive_Results\G0Partner_relig_results.dta", replace

capture postclose partner_relig_lr
postfile partner_relig_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G0Partner_relig_results_lr.dta", replace
	
capture postclose partner_relig_r2
postfile partner_relig_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G0Partner_relig_results_r2.dta", replace

foreach var of varlist IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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

ologit pb155_rev ageInPreg LoC_external, or
brant, detail


*** Now run the loop to save all the results
capture postclose partner_attend
postfile partner_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\Cognitive_Results\G0Partner_attend_results.dta", replace

capture postclose partner_attend_lr
postfile partner_attend_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G0Partner_attend_results_lr.dta", replace
	
capture postclose partner_attend_r2
postfile partner_attend_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G0Partner_attend_results_r2.dta", replace

foreach var of varlist IPSM_interpAware-selfEsteem {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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

postclose partner_attend
postclose partner_attend_lr
postclose partner_attend_r2


* Save this file if we want to use it later
save ".\Cognitive_Results\G0Partner_CogPredictorsOfRSBB_B3911_postAnalysis.dta", replace


**********************************************************************************
**** Now onto the G1's analysis
use "G1_PredictorsOfRSBB_B3911.dta", clear


*** Descriptive statistics

** Put the RSBB variables at the start of the dataset
order aln qlet YPG3000 YPG3040 YPG3040_grp YPG3080 YPG3080_OccNever YPG3080_OccYr YPG3153 YPG3153_cat YPG3160 YPG3170 YPG3155 YPG3155_cat

* Keep just pregnancy RSBB variables (belief in God, religious affiliation and frequency of church attendance)
drop YPG3153 YPG3153_cat YPG3160 YPG3170 YPG3155 YPG3155_cat

* Keep just the cognitive/psychological variables (and age and sex)
keep aln-YPG3080_OccYr kz021 YPG8000 f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126


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
foreach var of varlist kz021 YPG8000 f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126 {
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
misstable sum YPG3000-YPG3080_OccYr YPG3080_alt, all

* Exposures
misstable sum male YPG8000 f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126, all


** Also want to get descriptive stats for each exposure split by each RSBB outcome category
foreach var of varlist male YPG8000 f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126 {
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
desc male YPG8000 f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126

* All exposures are continuous or binary, making this easier

* Order variables in order to be analysed
order aln-YPG3080_OccYr YPG8000 male ///
f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126


* Next, rename all these exposures so they are more intuitive and can be read easier on the correlation heatmaps below
rename YPG8000 ageAt28
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
rename fg7364 Openness_age13
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


* Associations between cognitive/psychological variables - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
spearman verbalIQ_age8-globalEsteem_age8, pw

matrix cor_cog = r(Rho)

heatplot cor_cog, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle(""))

* Save heatmap
graph export ".\Cognitive_Results\G1_corr_heatplot.pdf", replace

* And save in .EPS format
heatplot cor_cog, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle("")) ysize(10) xsize(14)
	
graph export ".\Cognitive_Results\G1_corr_heatplot.eps", replace

* Save matrix as Excel file
putexcel set ".\Cognitive_Results\G1_corrMatrix.xlsx", replace
putexcel A1=matrix(cor_cog), names


************************************************************************************
*** Next, we want to run the actual analyses

*** Start with belief in God/divine power - As is a unordered categorical variable, will use multinomial regression (with 'no' as baseline/reference category)
tab YPG3000, m

** We want to store both estimates adjusting for age, and also the interaction between age and the exposure, to see whether it's moderated by age. As all cognitive variables are continuous or binary, this makes the script much simpler than for the demographic and SEP script, which includes categorical variables.

** Create a postfile to post results to, then start the loop - Will create three postfiles; one for coefficients and CIs, another for likelihood ratio tests comparing model fit (first of exposure model to no exposure model, then of interaction model to no interaction model), and then a third for the (pseduo-)R2 values to assess model fit - NOTE: Have to store pvalues as 'double' format, else really tiny p-values get coded as '0' (as default if float format, which has minimum value of -3.40282346639e+38 [https://blog.stata.com/2012/04/02/the-penultimate-guide-to-precision/]).
capture postclose G1_belief
postfile G1_belief str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\Cognitive_Results\G1_belief_results.dta", replace

capture postclose G1_belief_lr
postfile G1_belief_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G1_belief_results_lr.dta", replace
	
capture postclose G1_belief_r2
postfile G1_belief_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G1_belief_results_r2.dta", replace

foreach var of varlist verbalIQ_age8-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\Cognitive_Results\G1_relig_results.dta", replace

capture postclose G1_relig_lr
postfile G1_relig_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G1_relig_results_lr.dta", replace
	
capture postclose G1_relig_r2
postfile G1_relig_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G1_relig_results_r2.dta", replace

foreach var of varlist verbalIQ_age8-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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


** Quick test of whether proportional odds assumption been violated in most basic model (with just age at birth). Ah, it has been violated. 
ologit YPG3080_rev ageAt28, or
brant, detail

ologit YPG3080_rev ageAt28 totalIQ_age8, or
brant, detail


*** Now run the loop to save all the results
capture postclose G1_attend
postfile G1_attend str30 exposure str30 outcome_level str40 exp_level /// 
	n coef lci uci double(p) coef_int lci_int uci_int double(p_int) age_main exp_main ///
	using ".\Cognitive_Results\G1_attend_results.dta", replace

capture postclose G1_attend_lr
postfile G1_attend_lr str30 exposure double(lr_p_main lr_p_int) ///
	using ".\Cognitive_Results\G1_attend_results_lr.dta", replace
	
capture postclose G1_attend_r2
postfile G1_attend_r2 str30 exposure r2_main r2_int ///
	using ".\Cognitive_Results\G1_attend_results_r2.dta", replace

foreach var of varlist verbalIQ_age8-globalEsteem_age8 {
	
	// Save the exposure variable as a macro
	local exp = "`var'"
		
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

postclose G1_attend
postclose G1_attend_lr
postclose G1_attend_r2


* Save this file if we want to use it later
save ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", replace

* And close the log file
log close




***********************************************************************************
***********************************************************************************
***** Next, to start making plots of these results

**** G0 mothers

*** p-value plots

** First outcome - Religious belief
use ".\Cognitive_Results\G0Mother_belief_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Mother_belief_pvalues.dta", replace


** Next outcome - Religious affiliation
use ".\Cognitive_Results\G0Mother_relig_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Mother_relig_pvalues.dta", replace


** Next outcome - Religious attendance
use ".\Cognitive_Results\G0Mother_attend_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'religious attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Mother_attend_pvalues.dta", replace


** Combine all these datasets together
use ".\Cognitive_Results\G0Mother_belief_pvalues.dta", clear
append using ".\Cognitive_Results\G0Mother_relig_pvalues.dta"
append using ".\Cognitive_Results\G0Mother_attend_pvalues.dta"


* Now look at combined results

* Belief/religion/church vars main effects
local bon_thresh = -log10(0.05/15)
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
	ylabel(1(1)15, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(belRelCh_main, replace)

graph export ".\Cognitive_Results\G0Mother_mainEffects_pvalues.pdf", replace

* Belief/religion/church vars interaction effects
local bon_thresh = -log10(0.05/15)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)15, valuelabel labsize(small) angle(0)) ///
	title("Age interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(belRelCh_int, replace)

graph export ".\Cognitive_Results\G0Mother_ageInt_pvalues.pdf", replace


** Combine all these graphs together
graph combine belRelCh_main belRelCh_int, ysize(3) xsize(6)

graph export ".\Cognitive_Results\G0Mother_allData_pvalues.pdf", replace

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

outsheet exp_num-lr_p_intAttend using ".\Cognitive_Results\G0Mother_pvalue_results.csv", comma replace


** And how many exposures were associated with the outcome at both Bonferroni and standard alpha levels?

* Belief in god - main effect
count if lr_p_mainBelief < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainBelief < 0.05
display (r(N) / _N) * 100

* Belief in god - interaction
count if lr_p_intBelief < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intBelief < 0.05
display (r(N) / (_N)) * 100

* Religious affiliation - main effect
count if lr_p_mainReligion < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainReligion < 0.05
display (r(N) / _N) * 100

* Religious affiliation - interaction
count if lr_p_intReligion < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intReligion < 0.05
display (r(N) / (_N)) * 100

* Church attendance - main effect
count if lr_p_mainAttend < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainAttend < 0.05
display (r(N) / _N) * 100

* Church attendance - interaction
count if lr_p_intAttend < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intAttend < 0.05
display (r(N) / (_N)) * 100


*** Pseudo-R2 plots

** First outcome - Religious belief
use ".\Cognitive_Results\G0Mother_belief_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Mother_belief_r2.dta", replace


** Next outcome - Religious affiliation
use ".\Cognitive_Results\G0Mother_relig_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Mother_relig_r2.dta", replace


** Next outcome - Religious attendance
use ".\Cognitive_Results\G0Mother_attend_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'religious attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Mother_attend_r2.dta", replace


** Combine all these datasets together
use ".\Cognitive_Results\G0Mother_belief_r2.dta", clear
append using ".\Cognitive_Results\G0Mother_relig_r2.dta"
append using ".\Cognitive_Results\G0Mother_attend_r2.dta"


* Now look at combined results

* Belief/religion/church vars main effects
twoway (scatter exp_num r2_main if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)15, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(r2_main, replace)

graph export ".\Cognitive_Results\G0Mother_mainEffects_r2.pdf", replace

* Belief/religion/church vars interaction effects
twoway (scatter exp_num r2_int if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)15, valuelabel labsize(small) angle(0)) ///
	title("Age interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(r2_int, replace)

graph export ".\Cognitive_Results\G0Mother_ageInt_r2.pdf", replace


** Combine all these graphs together
graph combine r2_main r2_int, ysize(3) xsize(6)

graph export ".\Cognitive_Results\G0Mother_allData_r2.pdf", replace

graph close _all


** Save these pseudo R2 values as CSV files
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

outsheet exp_num-r2_intAttend using ".\Cognitive_Results\G0Mother_r2_results.csv", comma replace



*** Coefficient plots

** Will read the datasets in then combine together into one single dataset
use ".\Cognitive_Results\G0Mother_belief_results.dta", clear
gen outcome = "Belief"

append using ".\Cognitive_Results\G0Mother_relig_results.dta"
replace outcome = "Relig" if outcome == ""
tab outcome, m

append using ".\Cognitive_Results\G0Mother_attend_results.dta"
replace outcome = "Attend" if outcome == ""
tab outcome, m


** Save these results as CSV files to add to the SI

* Save each result in turn
format coef lci uci coef_int lci_int uci_int %9.3f
format p p_int %9.4f

outsheet exposure-p_int using ".\Cognitive_Results\G0Mother_belief_coefs.csv" if outcome == "Belief", comma replace

outsheet exposure-p_int using ".\Cognitive_Results\G0Mother_relig_coefs.csv" if outcome == "Relig", comma replace

outsheet exposure-p_int using ".\Cognitive_Results\G0Mother_attend_coefs.csv" if outcome == "Attend", comma replace


* Convert format back to default format (so axis on plots display correctly)
format coef lci uci coef_int lci_int uci_int %9.0g
format p p_int %10.0g


** Generate a 'levels' variable which combines all RSBB outcomes together
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


** Plot for cognitive ability factor

* Min and max x-axis values
sum lci uci if exposure == "intel_factor" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "intel_factor", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "intel_factor", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "intel_factor", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "intel_factor", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "intel_factor", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "intel_factor", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Cognitive Ability and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.7 0.8 0.9 1 1.1 1.2 1.3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(cog, replace)
		
graph export ".\Cognitive_Results\G0Mother_CogAbilityResults.pdf", replace


** Plot for IPSM total

* Min and max x-axis values
sum lci uci if exposure == "IPSM_total" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "IPSM_total", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "IPSM_total", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "IPSM_total", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "IPSM_total", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "IPSM_total", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "IPSM_total", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("IPSM total and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(1 1.005 1.01 1.015 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(ipsm, replace)
		
graph export ".\Cognitive_Results\G0Mother_IPSMTotalResults.pdf", replace


** Plot for LoC

* Min and max x-axis values
sum lci uci if exposure == "LoC_external" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "LoC_external", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "LoC_external", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "LoC_external", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("External LoC and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.7 0.8 0.9 1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(loc, replace)
		
graph export ".\Cognitive_Results\G0Mother_LoCResults.pdf", replace


** Plot for self-esteem

* Min and max x-axis values
sum lci uci if exposure == "selfEsteem" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "selfEsteem", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "selfEsteem", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "selfEsteem", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "selfEsteem", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "selfEsteem", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "selfEsteem", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Self-Esteem and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.99 1 1.01 1.02 1.03, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(esteem, replace)
		
graph export ".\Cognitive_Results\G0Mother_selfEsteemResults.pdf", replace


*** Also make a few interaction plots

** Plot for spot word by age interaction

* Min and max x-axis values
sum lci_int uci_int if exposure == "spotWord" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == "spotWord", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == "spotWord", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == "spotWord", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == "spotWord", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == "spotWord", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == "spotWord", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Spot-word task by age interaction", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.99 0.995 1 1.005, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(word_int, replace)
		
graph export ".\Cognitive_Results\G0Mother_wordByAgeInt.pdf", replace


** Plot for cognitive ability factor by age interaction

* Min and max x-axis values
sum lci_int uci_int if exposure == "intel_factor" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == "intel_factor", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == "intel_factor", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == "intel_factor", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == "intel_factor", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == "intel_factor", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == "intel_factor", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Cognitive ability by age interaction", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.98 1 1.02 1.04, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(cog_int, replace)
		
graph export ".\Cognitive_Results\G0Mother_cogAbilityByAgeInt.pdf", replace


** Plot for locus of control factor by age interaction

* Min and max x-axis values
sum lci_int uci_int if exposure == "LoC_external" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == "LoC_external", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == "LoC_external", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == "LoC_external", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("External LoC by age interaction", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.99 1 1.01 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(loc_int, replace)
		
graph export ".\Cognitive_Results\G0Mother_LoCByAgeInt.pdf", replace


*** Predicted probability plots (as multinomial relative risk ratio results not necessarily intuitive to interpret)

** Make one for cognitive ability
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d810 ageAtBirth intel_factor, rrr baseoutcome(3)
margins, at(intel_factor = (-5(1)5))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen intel_factor = fill(-5 -4)
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

twoway (line prob_yes intel_factor, col(black)) ///
	(rarea lci_yes uci_yes intel_factor, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure intel_factor, col(red)) ///
	(rarea lci_notSure uci_notSure intel_factor, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no intel_factor, col(blue)) ///
	(rarea lci_no uci_no intel_factor, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(-5 5)) xlabel(-5(1)5, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Cognitive Ability Factor") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Yes" 3 "Not sure" 5 "No") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(cog_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religious affiliation
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d813_grp ageAtBirth intel_factor, rrr baseoutcome(3)
margins, at(intel_factor = (-5(1)5))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen intel_factor = fill(-5 -4)
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

twoway (line prob_xian intel_factor, col(black)) ///
	(rarea lci_xian uci_xian intel_factor, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other intel_factor, col(red)) ///
	(rarea lci_other uci_other intel_factor, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none intel_factor, col(blue)) ///
	(rarea lci_none uci_none intel_factor, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(-5 5)) xlabel(-5(1)5, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Cognitive Ability Factor") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(cog_relig, replace)
	
	
* Attend church
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d816_rev ageAtBirth intel_factor, rrr baseoutcome(0)
margins, at(intel_factor = (-5(1)5))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen intel_factor = fill(-5 -4)
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

twoway (line prob_no intel_factor, col(black)) ///
	(rarea lci_no uci_no intel_factor, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr intel_factor, col(red)) ///
	(rarea lci_yr uci_yr intel_factor, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth intel_factor, col(blue)) ///
	(rarea lci_mth uci_mth intel_factor, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk intel_factor, col(green)) ///
	(rarea lci_wk uci_wk intel_factor, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(-5 5)) xlabel(-5(1)5, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Cognitive Ability Factor") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Not at all" 3 "1/yr" 5 "1/mth" 7 "1/wk") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(cog_attend, replace)
	
	
* Combine these plots together
graph combine cog_bel cog_relig cog_attend, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G0Mother_CogAbilityPredProbs_combined.pdf", replace

graph close _all


** Also plot interaction between age and cognitive ability - Start with religious belief
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum intel_factor

mlogit d810 c.ageAtBirth##c.intel_factor, rrr baseoutcome(3)
margins, at(intel_factor = (-5(1)5) ageAtBirth = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAtBirth = fill(20 30 40 20 30 40)
egen intel_factor = fill(-5 -5 -5 -4 -4 -4)
gen d810 = 1
predict p1, outcome(1)
sum p1

replace d810 = 2
predict p2, outcome(2)
sum p1 p2

replace d810 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 intel_factor if ageAtBirth == 20, col(black)) ///
	(line p1 intel_factor if ageAtBirth == 30, col(black) lpattern(dash)) ///
	(line p1 intel_factor if ageAtBirth == 40, col(black) lpattern(shortdash)) ///
	(line p2 intel_factor if ageAtBirth == 20, col(red)) ///
	(line p2 intel_factor if ageAtBirth == 30, col(red) lpattern(dash)) ///
	(line p2 intel_factor if ageAtBirth == 40, col(red) lpattern(shortdash)) ///
	(line p3 intel_factor if ageAtBirth == 20, col(blue)) ///
	(line p3 intel_factor if ageAtBirth == 30, col(blue) lpattern(dash)) ///
	(line p3 intel_factor if ageAtBirth == 40, col(blue) lpattern(shortdash)), ///
	xscale(range(-5 5)) xlabel(-5(1)5, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious belief", size(large)) ///
	xtitle("Cognitive Ability Factor") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Yes" 2 "Age 30 - Yes" 3 "Age 40 - Yes" ///
	4 "Age 20 - Not sure" 5 "Age 30 - Not sure" 6 "Age 40 - Not sure"  ///
	7 "Age 20 - No" 8 "Age 30 - No" 9 "Age 40 - No") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(cog_belief_int, replace)

	
* Religious affiliation
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum intel_factor

mlogit d813_grp c.ageAtBirth##c.intel_factor, rrr baseoutcome(3)
margins, at(intel_factor = (-5(1)5) ageAtBirth = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAtBirth = fill(20 30 40 20 30 40)
egen intel_factor = fill(-5 -5 -5 -4 -4 -4)
gen d813_grp = 1
predict p1, outcome(1)
sum p1

replace d813_grp = 2
predict p2, outcome(2)
sum p1 p2

replace d813_grp = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 intel_factor if ageAtBirth == 20, col(black)) ///
	(line p1 intel_factor if ageAtBirth == 30, col(black) lpattern(dash)) ///
	(line p1 intel_factor if ageAtBirth == 40, col(black) lpattern(shortdash)) ///
	(line p2 intel_factor if ageAtBirth == 20, col(red)) ///
	(line p2 intel_factor if ageAtBirth == 30, col(red) lpattern(dash)) ///
	(line p2 intel_factor if ageAtBirth == 40, col(red) lpattern(shortdash)) ///
	(line p3 intel_factor if ageAtBirth == 20, col(blue)) ///
	(line p3 intel_factor if ageAtBirth == 30, col(blue) lpattern(dash)) ///
	(line p3 intel_factor if ageAtBirth == 40, col(blue) lpattern(shortdash)), ///
	xscale(range(-5 5)) xlabel(-5(1)5, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious affiliation", size(large)) ///
	xtitle("Cognitive Ability Factor") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Christian" 2 "Age 30 - Christian" 3 "Age 40 - Christian" ///
	4 "Age 20 - Other" 5 "Age 30 - Other" 6 "Age 40 - Other"  ///
	7 "Age 20 - None" 8 "Age 30 - None" 9 "Age 40 - None") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(cog_relig_int, replace)
	

* Religious attendance
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum intel_factor

mlogit d816_rev c.ageAtBirth##c.intel_factor, rrr baseoutcome(0)
margins, at(intel_factor = (-5(1)5) ageAtBirth = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen ageAtBirth = fill(20 30 40 20 30 40)
egen intel_factor = fill(-5 -5 -5 -4 -4 -4)
gen d816_rev = 0
predict p1, outcome(0)
sum p1

replace d816_rev = 1
predict p2, outcome(1)
sum p1 p2

replace d816_rev = 2
predict p3, outcome(2)
sum p1 p2 p3

replace d816_rev = 3
predict p4, outcome(3)
sum p1 p2 p3 p4

twoway (line p1 intel_factor if ageAtBirth == 20, col(black)) ///
	(line p1 intel_factor if ageAtBirth == 30, col(black) lpattern(dash)) ///
	(line p1 intel_factor if ageAtBirth == 40, col(black) lpattern(shortdash)) ///
	(line p2 intel_factor if ageAtBirth == 20, col(red)) ///
	(line p2 intel_factor if ageAtBirth == 30, col(red) lpattern(dash)) ///
	(line p2 intel_factor if ageAtBirth == 40, col(red) lpattern(shortdash)) ///
	(line p3 intel_factor if ageAtBirth == 20, col(blue)) ///
	(line p3 intel_factor if ageAtBirth == 30, col(blue) lpattern(dash)) ///
	(line p3 intel_factor if ageAtBirth == 40, col(blue) lpattern(shortdash)) ///
	(line p4 intel_factor if ageAtBirth == 20, col(green)) ///
	(line p4 intel_factor if ageAtBirth == 30, col(green) lpattern(dash)) ///
	(line p4 intel_factor if ageAtBirth == 40, col(green) lpattern(shortdash)), ///
	xscale(range(-5 5)) xlabel(-5(1)5, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious attendance", size(large)) ///
	xtitle("Cognitive Ability Factor") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Never" 2 "Age 30 - Never" 3 "Age 40 - Never" ///
	4 "Age 20 - 1/Yr" 5 "Age 30 - 1/Yr" 6 "Age 40 - 1/Yr"  ///
	7 "Age 20 - 1/Mth" 8 "Age 30 - 1/Mth" 9 "Age 40 - 1/Mth" ///
	7 "Age 20 - 1/Wk" 8 "Age 30 - 1/Wk" 9 "Age 40 - 1/Wk") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(cog_attend_int, replace)


* Combine these plots together
graph combine cog_belief_int cog_relig_int cog_attend_int, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G0Mother_CogAbilityPredProbs_ageInt_combined.pdf", replace

graph close _all


** Predicted probability plots for locus of control
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d810 ageAtBirth LoC_external, rrr baseoutcome(3)
margins, at(LoC_external = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen LoC_external = fill(0 1)
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

twoway (line prob_yes LoC_external, col(black)) ///
	(rarea lci_yes uci_yes LoC_external, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure LoC_external, col(red)) ///
	(rarea lci_notSure uci_notSure LoC_external, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no LoC_external, col(blue)) ///
	(rarea lci_no uci_no LoC_external, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("External Locus of Control") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Yes" 3 "Not sure" 5 "No") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(loc_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religious affiliation
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d813_grp ageAtBirth LoC_external, rrr baseoutcome(3)
margins, at(LoC_external = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen LoC_external = fill(0 1)
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

twoway (line prob_xian LoC_external, col(black)) ///
	(rarea lci_xian uci_xian LoC_external, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other LoC_external, col(red)) ///
	(rarea lci_other uci_other LoC_external, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none LoC_external, col(blue)) ///
	(rarea lci_none uci_none LoC_external, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("External Locus of Control") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(loc_relig, replace)
	
	
* Attend church
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit d816_rev ageAtBirth LoC_external, rrr baseoutcome(0)
margins, at(LoC_external = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen LoC_external = fill(0 1)
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

twoway (line prob_no LoC_external, col(black)) ///
	(rarea lci_no uci_no LoC_external, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr LoC_external, col(red)) ///
	(rarea lci_yr uci_yr LoC_external, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth LoC_external, col(blue)) ///
	(rarea lci_mth uci_mth LoC_external, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk LoC_external, col(green)) ///
	(rarea lci_wk uci_wk LoC_external, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("External Locus of Control") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Not at all" 3 "1/yr" 5 "1/mth" 7 "1/wk") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(loc_attend, replace)
	
	
* Combine these plots together
graph combine loc_bel loc_relig loc_attend, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G0Mother_LoCPredProbs_combined.pdf", replace

graph close _all


** Also plot interaction between age and LoC - Start with religious belief
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum LoC_external

mlogit d810 c.ageAtBirth##c.LoC_external, rrr baseoutcome(3)
margins, at(LoC_external = (0(1)12) ageAtBirth = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAtBirth = fill(20 30 40 20 30 40)
egen LoC_external = fill(0 0 0 1 1 1)
gen d810 = 1
predict p1, outcome(1)
sum p1

replace d810 = 2
predict p2, outcome(2)
sum p1 p2

replace d810 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 LoC_external if ageAtBirth == 20, col(black)) ///
	(line p1 LoC_external if ageAtBirth == 30, col(black) lpattern(dash)) ///
	(line p1 LoC_external if ageAtBirth == 40, col(black) lpattern(shortdash)) ///
	(line p2 LoC_external if ageAtBirth == 20, col(red)) ///
	(line p2 LoC_external if ageAtBirth == 30, col(red) lpattern(dash)) ///
	(line p2 LoC_external if ageAtBirth == 40, col(red) lpattern(shortdash)) ///
	(line p3 LoC_external if ageAtBirth == 20, col(blue)) ///
	(line p3 LoC_external if ageAtBirth == 30, col(blue) lpattern(dash)) ///
	(line p3 LoC_external if ageAtBirth == 40, col(blue) lpattern(shortdash)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious belief", size(large)) ///
	xtitle("External LoC") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Yes" 2 "Age 30 - Yes" 3 "Age 40 - Yes" ///
	4 "Age 20 - Not sure" 5 "Age 30 - Not sure" 6 "Age 40 - Not sure"  ///
	7 "Age 20 - No" 8 "Age 30 - No" 9 "Age 40 - No") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(loc_belief_int, replace)

	
* Religious affiliation
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum LoC_external

mlogit d813_grp c.ageAtBirth##c.LoC_external, rrr baseoutcome(3)
margins, at(LoC_external = (0(1)12) ageAtBirth = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAtBirth = fill(20 30 40 20 30 40)
egen LoC_external = fill(0 0 0 1 1 1)
gen d813_grp = 1
predict p1, outcome(1)
sum p1

replace d813_grp = 2
predict p2, outcome(2)
sum p1 p2

replace d813_grp = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 LoC_external if ageAtBirth == 20, col(black)) ///
	(line p1 LoC_external if ageAtBirth == 30, col(black) lpattern(dash)) ///
	(line p1 LoC_external if ageAtBirth == 40, col(black) lpattern(shortdash)) ///
	(line p2 LoC_external if ageAtBirth == 20, col(red)) ///
	(line p2 LoC_external if ageAtBirth == 30, col(red) lpattern(dash)) ///
	(line p2 LoC_external if ageAtBirth == 40, col(red) lpattern(shortdash)) ///
	(line p3 LoC_external if ageAtBirth == 20, col(blue)) ///
	(line p3 LoC_external if ageAtBirth == 30, col(blue) lpattern(dash)) ///
	(line p3 LoC_external if ageAtBirth == 40, col(blue) lpattern(shortdash)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious affiliation", size(large)) ///
	xtitle("External LoC") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Christian" 2 "Age 30 - Christian" 3 "Age 40 - Christian" ///
	4 "Age 20 - Other" 5 "Age 30 - Other" 6 "Age 40 - Other"  ///
	7 "Age 20 - None" 8 "Age 30 - None" 9 "Age 40 - None") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(loc_relig_int, replace)
	

* Religious attendance
use ".\Cognitive_Results\G0Mother_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum LoC_external

mlogit d816_rev c.ageAtBirth##c.LoC_external, rrr baseoutcome(0)
margins, at(LoC_external = (0(1)12) ageAtBirth = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen ageAtBirth = fill(20 30 40 20 30 40)
egen LoC_external = fill(0 0 0 1 1 1)
gen d816_rev = 0
predict p1, outcome(0)
sum p1

replace d816_rev = 1
predict p2, outcome(1)
sum p1 p2

replace d816_rev = 2
predict p3, outcome(2)
sum p1 p2 p3

replace d816_rev = 3
predict p4, outcome(3)
sum p1 p2 p3 p4

twoway (line p1 LoC_external if ageAtBirth == 20, col(black)) ///
	(line p1 LoC_external if ageAtBirth == 30, col(black) lpattern(dash)) ///
	(line p1 LoC_external if ageAtBirth == 40, col(black) lpattern(shortdash)) ///
	(line p2 LoC_external if ageAtBirth == 20, col(red)) ///
	(line p2 LoC_external if ageAtBirth == 30, col(red) lpattern(dash)) ///
	(line p2 LoC_external if ageAtBirth == 40, col(red) lpattern(shortdash)) ///
	(line p3 LoC_external if ageAtBirth == 20, col(blue)) ///
	(line p3 LoC_external if ageAtBirth == 30, col(blue) lpattern(dash)) ///
	(line p3 LoC_external if ageAtBirth == 40, col(blue) lpattern(shortdash)) ///
	(line p4 LoC_external if ageAtBirth == 20, col(green)) ///
	(line p4 LoC_external if ageAtBirth == 30, col(green) lpattern(dash)) ///
	(line p4 LoC_external if ageAtBirth == 40, col(green) lpattern(shortdash)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious attendance", size(large)) ///
	xtitle("External LoC") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Never" 2 "Age 30 - Never" 3 "Age 40 - Never" ///
	4 "Age 20 - 1/Yr" 5 "Age 30 - 1/Yr" 6 "Age 40 - 1/Yr"  ///
	7 "Age 20 - 1/Mth" 8 "Age 30 - 1/Mth" 9 "Age 40 - 1/Mth" ///
	7 "Age 20 - 1/Wk" 8 "Age 30 - 1/Wk" 9 "Age 40 - 1/Wk") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(loc_attend_int, replace)
	
	
* Combine these plots together
graph combine loc_belief_int loc_relig_int loc_attend_int, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G0Mother_LoCPredProbs_ageInt_combined.pdf", replace

graph close _all


	
**** G0 partner results

*** p-value plots

** First outcome - Religious belief
use ".\Cognitive_Results\G0Partner_belief_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Partner_belief_pvalues.dta", replace


** Next outcome - Religious affiliation
use ".\Cognitive_Results\G0Partner_relig_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Partner_relig_pvalues.dta", replace


** Next outcome - Religious attendance
use ".\Cognitive_Results\G0Partner_attend_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'religious attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Partner_attend_pvalues.dta", replace


** Combine all these datasets together
use ".\Cognitive_Results\G0Partner_belief_pvalues.dta", clear
append using ".\Cognitive_Results\G0Partner_relig_pvalues.dta"
append using ".\Cognitive_Results\G0Partner_attend_pvalues.dta"


* Now look at combined results

* Belief/religion/church vars main effects
local bon_thresh = -log10(0.05/8)
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
	ylabel(1(1)8, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(belRelCh_main, replace)

graph export ".\Cognitive_Results\G0Partner_mainEffects_pvalues.pdf", replace

* Belief/religion/church vars interaction effects
local bon_thresh = -log10(0.05/8)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)8, valuelabel labsize(small) angle(0)) ///
	title("Age interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(belRelCh_int, replace)

graph export ".\Cognitive_Results\G0Partner_ageInt_pvalues.pdf", replace


** Combine all these graphs together
graph combine belRelCh_main belRelCh_int, ysize(3) xsize(6)

graph export ".\Cognitive_Results\G0Partner_allData_pvalues.pdf", replace

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

outsheet exp_num-lr_p_intAttend using ".\Cognitive_Results\G0Partner_pvalue_results.csv", comma replace


** And how many exposures were associated with the outcome at both Bonferroni and standard alpha levels?

* Belief in god - main effect
count if lr_p_mainBelief < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainBelief < 0.05
display (r(N) / _N) * 100

* Belief in god - interaction
count if lr_p_intBelief < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intBelief < 0.05
display (r(N) / (_N)) * 100

* Religious affiliation - main effect
count if lr_p_mainReligion < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainReligion < 0.05
display (r(N) / _N) * 100

* Religious affiliation - interaction
count if lr_p_intReligion < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intReligion < 0.05
display (r(N) / (_N)) * 100

* Church attendance - main effect
count if lr_p_mainAttend < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainAttend < 0.05
display (r(N) / _N) * 100

* Church attendance - interaction
count if lr_p_intAttend < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intAttend < 0.05
display (r(N) / (_N)) * 100


*** Pseudo-R2 plots

** First outcome - Religious belief
use ".\Cognitive_Results\G0Partner_belief_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Partner_belief_r2.dta", replace


** Next outcome - Religious affiliation
use ".\Cognitive_Results\G0Partner_relig_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Partner_relig_r2.dta", replace


** Next outcome - Religious attendance
use ".\Cognitive_Results\G0Partner_attend_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'religious attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G0Partner_attend_r2.dta", replace


** Combine all these datasets together
use ".\Cognitive_Results\G0Partner_belief_r2.dta", clear
append using ".\Cognitive_Results\G0Partner_relig_r2.dta"
append using ".\Cognitive_Results\G0Partner_attend_r2.dta"


* Now look at combined results

* Belief/religion/church vars main effects
twoway (scatter exp_num r2_main if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)8, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(r2_main, replace)

graph export ".\Cognitive_Results\G0Partner_mainEffects_r2.pdf", replace

* Belief/religion/church vars interaction effects
twoway (scatter exp_num r2_int if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)8, valuelabel labsize(small) angle(0)) ///
	title("Age interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) size(small)) ///
	name(r2_int, replace)

graph export ".\Cognitive_Results\G0Partner_ageInt_r2.pdf", replace


** Combine all these graphs together
graph combine r2_main r2_int, ysize(3) xsize(6)

graph export ".\Cognitive_Results\G0Partner_allData_r2.pdf", replace

graph close _all


** Save these pseudo R2 values as CSV files
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

outsheet exp_num-r2_intAttend using ".\Cognitive_Results\G0Partner_r2_results.csv", comma replace


*** Coefficient plots

** Will read the datasets in then combine together into one single dataset
use ".\Cognitive_Results\G0Partner_belief_results.dta", clear
gen outcome = "Belief"

append using ".\Cognitive_Results\G0Partner_relig_results.dta"
replace outcome = "Relig" if outcome == ""
tab outcome, m

append using ".\Cognitive_Results\G0Partner_attend_results.dta"
replace outcome = "Attend" if outcome == ""
tab outcome, m


** Save these results as CSV files to add to the SI

* Save each result in turn
format coef lci uci coef_int lci_int uci_int %9.3f
format p p_int %9.4f

outsheet exposure-p_int using ".\Cognitive_Results\G0Partner_belief_coefs.csv" if outcome == "Belief", comma replace

outsheet exposure-p_int using ".\Cognitive_Results\G0Partner_relig_coefs.csv" if outcome == "Relig", comma replace

outsheet exposure-p_int using ".\Cognitive_Results\G0Partner_attend_coefs.csv" if outcome == "Attend", comma replace


* Convert format back to default format (so axis on plots display correctly)
format coef lci uci coef_int lci_int uci_int %9.0g
format p p_int %10.0g


** Generate a 'levels' variable which combines all RSBB outcomes together
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


** Plot for IPSM total

* Min and max x-axis values
sum lci uci if exposure == "IPSM_total" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "IPSM_total", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "IPSM_total", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "IPSM_total", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "IPSM_total", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "IPSM_total", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "IPSM_total", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("IPSM total and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(1 1.005 1.01 1.015 1.02 1.025, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(ipsm, replace)
		
graph export ".\Cognitive_Results\G0Partner_IPSMTotalResults.pdf", replace


** Plot for LoC

* Min and max x-axis values
sum lci uci if exposure == "LoC_external" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "LoC_external", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "LoC_external", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "LoC_external", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("External LoC and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.7 0.8 0.9 1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(loc, replace)
		
graph export ".\Cognitive_Results\G0Partner_LoCResults.pdf", replace


** Plot for self-esteem

* Min and max x-axis values
sum lci uci if exposure == "selfEsteem" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "selfEsteem", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "selfEsteem", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "selfEsteem", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "selfEsteem", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "selfEsteem", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "selfEsteem", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Self-Esteem and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.97 0.98 0.99 1 1.01, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(esteem, replace)
		
graph export ".\Cognitive_Results\G0Partner_selfEsteemResults.pdf", replace


*** Also make a few interaction plots

** Plot for locus of control factor by age interaction

* Min and max x-axis values
sum lci_int uci_int if exposure == "LoC_external" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == "LoC_external", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == "LoC_external", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == "LoC_external", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == "LoC_external", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("External LoC by age interaction", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.99 1 1.01 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(loc_int, replace)
		
graph export ".\Cognitive_Results\G0Partner_LoCByAgeInt.pdf", replace

graph close _all


*** Predicted probability plots (as multinomial relative risk ratio results not necessarily intuitive to interpret)

** Predicted probability plots for locus of control
use ".\Cognitive_Results\G0Partner_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit pb150 ageInPreg LoC_external, rrr baseoutcome(3)
margins, at(LoC_external = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen LoC_external = fill(0 1)
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

twoway (line prob_yes LoC_external, col(black)) ///
	(rarea lci_yes uci_yes LoC_external, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure LoC_external, col(red)) ///
	(rarea lci_notSure uci_notSure LoC_external, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no LoC_external, col(blue)) ///
	(rarea lci_no uci_no LoC_external, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("External Locus of Control") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Yes" 3 "Not sure" 5 "No") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(loc_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religious affiliation
use ".\Cognitive_Results\G0Partner_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit pb153_grp ageInPreg LoC_external, rrr baseoutcome(3)
margins, at(LoC_external = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen LoC_external = fill(0 1)
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

twoway (line prob_xian LoC_external, col(black)) ///
	(rarea lci_xian uci_xian LoC_external, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other LoC_external, col(red)) ///
	(rarea lci_other uci_other LoC_external, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none LoC_external, col(blue)) ///
	(rarea lci_none uci_none LoC_external, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("External Locus of Control") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(loc_relig, replace)
	
	
* Attend church
use ".\Cognitive_Results\G0Partner_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit pb155_rev ageInPreg LoC_external, rrr baseoutcome(0)
margins, at(LoC_external = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen LoC_external = fill(0 1)
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

twoway (line prob_no LoC_external, col(black)) ///
	(rarea lci_no uci_no LoC_external, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr LoC_external, col(red)) ///
	(rarea lci_yr uci_yr LoC_external, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth LoC_external, col(blue)) ///
	(rarea lci_mth uci_mth LoC_external, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk LoC_external, col(green)) ///
	(rarea lci_wk uci_wk LoC_external, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("External Locus of Control") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Not at all" 3 "1/yr" 5 "1/mth" 7 "1/wk") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(loc_attend, replace)
	
	
* Combine these plots together
graph combine loc_bel loc_relig loc_attend, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G0Partner_LoCPredProbs_combined.pdf", replace

graph close _all


** Also plot interaction between age and LoC - Start with religious belief
use ".\Cognitive_Results\G0Partner_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum LoC_external

mlogit pb150 c.ageInPreg##c.LoC_external, rrr baseoutcome(3)
margins, at(LoC_external = (0(1)11) ageInPreg = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageInPreg = fill(20 30 40 20 30 40)
egen LoC_external = fill(0 0 0 1 1 1)
gen pb150 = 1
predict p1, outcome(1)
sum p1

replace pb150 = 2
predict p2, outcome(2)
sum p1 p2

replace pb150 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 LoC_external if ageInPreg == 20, col(black)) ///
	(line p1 LoC_external if ageInPreg == 30, col(black) lpattern(dash)) ///
	(line p1 LoC_external if ageInPreg == 40, col(black) lpattern(shortdash)) ///
	(line p2 LoC_external if ageInPreg == 20, col(red)) ///
	(line p2 LoC_external if ageInPreg == 30, col(red) lpattern(dash)) ///
	(line p2 LoC_external if ageInPreg == 40, col(red) lpattern(shortdash)) ///
	(line p3 LoC_external if ageInPreg == 20, col(blue)) ///
	(line p3 LoC_external if ageInPreg == 30, col(blue) lpattern(dash)) ///
	(line p3 LoC_external if ageInPreg == 40, col(blue) lpattern(shortdash)), ///
	xscale(range(0 11)) xlabel(0(1)11, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious belief", size(large)) ///
	xtitle("External LoC") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Yes" 2 "Age 30 - Yes" 3 "Age 40 - Yes" ///
	4 "Age 20 - Not sure" 5 "Age 30 - Not sure" 6 "Age 40 - Not sure"  ///
	7 "Age 20 - No" 8 "Age 30 - No" 9 "Age 40 - No") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(loc_belief_int, replace)

	
* Religious affiliation
use ".\Cognitive_Results\G0Partner_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum LoC_external

mlogit pb153_grp c.ageInPreg##c.LoC_external, rrr baseoutcome(3)
margins, at(LoC_external = (0(1)11) ageInPreg = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageInPreg = fill(20 30 40 20 30 40)
egen LoC_external = fill(0 0 0 1 1 1)
gen pb153_grp = 1
predict p1, outcome(1)
sum p1

replace pb153_grp = 2
predict p2, outcome(2)
sum p1 p2

replace pb153_grp = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 LoC_external if ageInPreg == 20, col(black)) ///
	(line p1 LoC_external if ageInPreg == 30, col(black) lpattern(dash)) ///
	(line p1 LoC_external if ageInPreg == 40, col(black) lpattern(shortdash)) ///
	(line p2 LoC_external if ageInPreg == 20, col(red)) ///
	(line p2 LoC_external if ageInPreg == 30, col(red) lpattern(dash)) ///
	(line p2 LoC_external if ageInPreg == 40, col(red) lpattern(shortdash)) ///
	(line p3 LoC_external if ageInPreg == 20, col(blue)) ///
	(line p3 LoC_external if ageInPreg == 30, col(blue) lpattern(dash)) ///
	(line p3 LoC_external if ageInPreg == 40, col(blue) lpattern(shortdash)), ///
	xscale(range(0 11)) xlabel(0(1)11, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious affiliation", size(large)) ///
	xtitle("External LoC") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Christian" 2 "Age 30 - Christian" 3 "Age 40 - Christian" ///
	4 "Age 20 - Other" 5 "Age 30 - Other" 6 "Age 40 - Other"  ///
	7 "Age 20 - None" 8 "Age 30 - None" 9 "Age 40 - None") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(loc_relig_int, replace)
	

* Religious attendance
use ".\Cognitive_Results\G0Partner_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum LoC_external

mlogit pb155_rev c.ageInPreg##c.LoC_external, rrr baseoutcome(0)
margins, at(LoC_external = (0(1)11) ageInPreg = (20 30 40))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen ageInPreg = fill(20 30 40 20 30 40)
egen LoC_external = fill(0 0 0 1 1 1)
gen pb155_rev = 0
predict p1, outcome(0)
sum p1

replace pb155_rev = 1
predict p2, outcome(1)
sum p1 p2

replace pb155_rev = 2
predict p3, outcome(2)
sum p1 p2 p3

replace pb155_rev = 3
predict p4, outcome(3)
sum p1 p2 p3 p4

twoway (line p1 LoC_external if ageInPreg == 20, col(black)) ///
	(line p1 LoC_external if ageInPreg == 30, col(black) lpattern(dash)) ///
	(line p1 LoC_external if ageInPreg == 40, col(black) lpattern(shortdash)) ///
	(line p2 LoC_external if ageInPreg == 20, col(red)) ///
	(line p2 LoC_external if ageInPreg == 30, col(red) lpattern(dash)) ///
	(line p2 LoC_external if ageInPreg == 40, col(red) lpattern(shortdash)) ///
	(line p3 LoC_external if ageInPreg == 20, col(blue)) ///
	(line p3 LoC_external if ageInPreg == 30, col(blue) lpattern(dash)) ///
	(line p3 LoC_external if ageInPreg == 40, col(blue) lpattern(shortdash)) ///
	(line p4 LoC_external if ageInPreg == 20, col(green)) ///
	(line p4 LoC_external if ageInPreg == 30, col(green) lpattern(dash)) ///
	(line p4 LoC_external if ageInPreg == 40, col(green) lpattern(shortdash)), ///
	xscale(range(0 11)) xlabel(0(1)11, labsize(small)) ylabel(, labsize(small)) ///
	title("Religious attendance", size(large)) ///
	xtitle("External LoC") ytitle("Predicted probability", margin(small)) ///
	legend(order(1 "Age 20 - Never" 2 "Age 30 - Never" 3 "Age 40 - Never" ///
	4 "Age 20 - 1/Yr" 5 "Age 30 - 1/Yr" 6 "Age 40 - 1/Yr"  ///
	7 "Age 20 - 1/Mth" 8 "Age 30 - 1/Mth" 9 "Age 40 - 1/Mth" ///
	10 "Age 20 - 1/Wk" 11 "Age 30 - 1/Wk" 12 "Age 40 - 1/Wk") cols(3) ///
	colgap(*.5) keygap(*.5) size(vsmall)) ///
	name(loc_attend_int, replace)
	
	
* Combine these plots together
graph combine loc_belief_int loc_relig_int loc_attend_int, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G0Partner_LoCPredProbs_ageInt_combined.pdf", replace

graph close _all



**** G1 offspring

*** p-value plots

** First outcome - Religious belief
use ".\Cognitive_Results\G1_belief_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G1_belief_pvalues.dta", replace


** Next outcome - Religious affiliation
use ".\Cognitive_Results\G1_relig_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\Cognitive_Results\G1_relig_pvalues.dta", replace


** Next outcome - Religious attendance
use ".\Cognitive_Results\G1_attend_results_lr.dta", clear

* Convert string exposure var to numeric
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

* Convert p-values to -log10 p-values
gen logp_main = -log10(lr_p_main)
sum logp_main

gen logp_int = -log10(lr_p_int)
sum logp_int

* Add 'religious attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G1_attend_pvalues.dta", replace


** Combine all these datasets together
use ".\Cognitive_Results\G1_belief_pvalues.dta", clear
append using ".\Cognitive_Results\G1_relig_pvalues.dta"
append using ".\Cognitive_Results\G1_attend_pvalues.dta"


* Now look at combined results

* Belief/religion/church vars main effects
local bon_thresh = -log10(0.05/25)
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
	ylabel(1(1)25, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) colgap(*.5) keygap(*.5) size(small)) ///
	name(belRelCh_main, replace)

graph export ".\Cognitive_Results\G1_mainEffects_pvalues.pdf", replace

* Belief/religion/church vars interaction effects
local bon_thresh = -log10(0.05/25)
local thresh_05 = -log10(0.05)

twoway (scatter exp_num logp_int if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num logp_int if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xline(`bon_thresh', lcol(black) lpattern(dash)) ///
	xline(`thresh_05', lcol(black) lpattern(dot)) ///
	xtitle("-log10 of p-value") ytitle("") ysc(reverse) ///
	ylabel(1(1)25, valuelabel labsize(small) angle(0)) ///
	title("Sex interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) colgap(*.5) keygap(*.5) size(small)) ///
	name(belRelCh_int, replace)

graph export ".\Cognitive_Results\G1_sexInt_pvalues.pdf", replace


** Combine all these graphs together
graph combine belRelCh_main belRelCh_int, ysize(3) xsize(6)

graph export ".\Cognitive_Results\G1_allData_pvalues.pdf", replace

* And convert to EPS format for Wellcome Open Research formatting (also have to enlarge size, else resolution is terrible)
graph combine belRelCh_main belRelCh_int, ysize(10) xsize(20)
graph export ".\Cognitive_Results\G1_allData_pvalues.eps", replace

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

outsheet exp_num-lr_p_intAttend using ".\Cognitive_Results\G1_pvalue_results.csv", comma replace


** And how many exposures were associated with the outcome at both Bonferroni and standard alpha levels?

* Belief in god - main effect
count if lr_p_mainBelief < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainBelief < 0.05
display (r(N) / _N) * 100

* Belief in god - interaction
count if lr_p_intBelief < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intBelief < 0.05
display (r(N) / (_N)) * 100

* Religious affiliation - main effect
count if lr_p_mainReligion < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainReligion < 0.05
display (r(N) / _N) * 100

* Religious affiliation - interaction
count if lr_p_intReligion < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intReligion < 0.05
display (r(N) / (_N)) * 100

* Church attendance - main effect
count if lr_p_mainAttend < 0.05/_N
display (r(N) / _N) * 100

count if lr_p_mainAttend < 0.05
display (r(N) / _N) * 100

* Church attendance - interaction
count if lr_p_intAttend < 0.05/(_N)
display (r(N) / (_N)) * 100

count if lr_p_intAttend < 0.05
display (r(N) / (_N)) * 100


*** Pseudo-R2 plots

** First outcome - Religious belief
use ".\Cognitive_Results\G1_belief_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'belief' as a variable, then save this file (as will merge with other files later on)
gen outcome = "Belief"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G1_belief_r2.dta", replace


** Next outcome - Religious affiliation
use ".\Cognitive_Results\G1_relig_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'religious affiliation' as a variable, then save this file
gen outcome = "Religious affil."
recast str30 outcome
order outcome

save ".\Cognitive_Results\G1_relig_r2.dta", replace


** Next outcome - Religious attendance
use ".\Cognitive_Results\G1_attend_results_r2.dta", clear

* Convert string exposure var to numeric
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

* Add 'religious attendance' as a variable, then save this file
gen outcome = "Church attendance"
recast str30 outcome
order outcome

save ".\Cognitive_Results\G1_attend_r2.dta", replace


** Combine all these datasets together
use ".\Cognitive_Results\G1_belief_r2.dta", clear
append using ".\Cognitive_Results\G1_relig_r2.dta"
append using ".\Cognitive_Results\G1_attend_r2.dta"


* Now look at combined results

* Belief/religion/church vars main effects
twoway (scatter exp_num r2_main if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_main if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)25, valuelabel labsize(small) angle(0)) ///
	title("Main effects") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) colgap(*.5) keygap(*.5) size(small)) ///
	name(r2_main, replace)

graph export ".\Cognitive_Results\G1_mainEffects_r2.pdf", replace

* Belief/religion/church vars interaction effects
twoway (scatter exp_num r2_int if outcome == "Belief", ///
		col(black) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Religious affil.", ///
		col(red) msize(small) msym(D)) ///
	(scatter exp_num r2_int if outcome == "Church attendance", ///
		col(blue) msize(small) msym(D)), ///
	xtitle("Pseudo-R2 value") ytitle("") ysc(reverse) ///
	ylabel(1(1)25, valuelabel labsize(small) angle(0)) ///
	title("Sex interaction") ///
	legend(order(1 "Religious belief" 2 "Religious affiliation" ///
		3 "Religious attendance") rows(1) colgap(*.5) keygap(*.5) size(small)) ///
	name(r2_int, replace)

graph export ".\Cognitive_Results\G1_sexInt_r2.pdf", replace


** Combine all these graphs together
graph combine r2_main r2_int, ysize(3) xsize(6)

graph export ".\Cognitive_Results\G1_allData_r2.pdf", replace

* And convert to EPS format
graph combine r2_main r2_int, ysize(10) xsize(20)
graph export ".\Cognitive_Results\G1_allData_r2.eps", replace

graph close _all


** Save these pseudo R2 values as CSV files
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

outsheet exp_num-r2_intAttend using ".\Cognitive_Results\G1_r2_results.csv", comma replace


*** Coefficient plots

** Will read the datasets in then combine together into one single dataset
use ".\Cognitive_Results\G1_belief_results.dta", clear
gen outcome = "Belief"

append using ".\Cognitive_Results\G1_relig_results.dta"
replace outcome = "Relig" if outcome == ""
tab outcome, m

append using ".\Cognitive_Results\G1_attend_results.dta"
replace outcome = "Attend" if outcome == ""
tab outcome, m


** Save these results as CSV files to add to the SI

* Save each result in turn
format coef lci uci coef_int lci_int uci_int %9.3f
format p p_int %9.4f

outsheet exposure-p_int using ".\Cognitive_Results\G1_belief_coefs.csv" if outcome == "Belief", comma replace

outsheet exposure-p_int using ".\Cognitive_Results\G1_relig_coefs.csv" if outcome == "Relig", comma replace

outsheet exposure-p_int using ".\Cognitive_Results\G1_attend_coefs.csv" if outcome == "Attend", comma replace


* Convert format back to default format (so axis on plots display correctly)
format coef lci uci coef_int lci_int uci_int %9.0g
format p p_int %10.0g


** Generate a 'levels' variable which combines all RSBB outcomes together
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


** Plot for total IQ at age 8

* Min and max x-axis values
sum lci uci if exposure == "totalIQ_age8" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "totalIQ_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "totalIQ_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "totalIQ_age8", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Total IQ age 8 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.98 0.99 1 1.01 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(iq8, replace)
		
graph export ".\Cognitive_Results\G1_totalIQAge8Results.pdf", replace


* And in EPS format
twoway (scatter level_num coef if outcome == "Belief" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "totalIQ_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "totalIQ_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "totalIQ_age8", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Total IQ age 8 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.98 0.99 1 1.01 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(iq8, replace) ysize(10) xsize(14)
		
graph export ".\Cognitive_Results\G1_totalIQAge8Results.eps", replace


** Plot for total IQ at age 15

* Min and max x-axis values
sum lci uci if exposure == "totalIQ_age15" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "totalIQ_age15", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "totalIQ_age15", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "totalIQ_age15", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "totalIQ_age15", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "totalIQ_age15", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "totalIQ_age15", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Total IQ age 15 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.98 0.99 1 1.01 1.02 1.03, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(iq15, replace)
		
graph export ".\Cognitive_Results\G1_totalIQAge15Results.pdf", replace


** Plot for vocabulary test at age 24

* Min and max x-axis values
sum lci uci if exposure == "vocab_age24" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "vocab_age24", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "vocab_age24", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "vocab_age24", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "vocab_age24", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "vocab_age24", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "vocab_age24", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Vocab Test Age 24 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.9 0.95 1 1.05 1.1 1.15, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(vocab24, replace)
		
graph export ".\Cognitive_Results\G1_vocab24Results.pdf", replace


** Plot for personality trait extraversion at age 13

* Min and max x-axis values
sum lci uci if exposure == "extraversion_age13" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "extraversion_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "extraversion_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "extraversion_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "extraversion_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "extraversion_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "extraversion_age13", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Extraversion Age 13 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.94 0.96 0.98 1 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(extra13, replace)
		
graph export ".\Cognitive_Results\G1_extra13Results.pdf", replace


** Plot for personality trait conscientiousness at age 13

* Min and max x-axis values
sum lci uci if exposure == "conscientiousness_age13" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "conscientiousness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "conscientiousness_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "conscientiousness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "conscientiousness_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "conscientiousness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "conscientiousness_age13", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Conscientiousness Age 13 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.94 0.96 0.98 1 1.02 1.04 1.06, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(consc13, replace)
		
graph export ".\Cognitive_Results\G1_consc13Results.pdf", replace


** Plot for personality trait openness to experience at age 13

* Min and max x-axis values
sum lci uci if exposure == "Openness_age13" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "Openness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "Openness_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "Openness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "Openness_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "Openness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "Openness_age13", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Openness Age 13 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.97 0.98 0.99 1 1.01 1.01 1.02 1.03 1.04 1.05, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(open13, replace)
		
graph export ".\Cognitive_Results\G1_open13Results.pdf", replace


** Plot for LoC at age 8

* Min and max x-axis values
sum lci uci if exposure == "loc_age8" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "loc_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "loc_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "loc_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "loc_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "loc_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "loc_age8", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("External LoC Age 8 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.8 0.85 0.9 0.95 1 1.05 1.1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(loc8, replace)
		
graph export ".\Cognitive_Results\G1_LoCage8Results.pdf", replace


** Plot for LoC at age 16

* Min and max x-axis values
sum lci uci if exposure == "loc_age16" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "loc_age16", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "loc_age16", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "loc_age16", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "loc_age16", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "loc_age16", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "loc_age16", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("External LoC Age 16 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.8 0.85 0.9 0.95 1 1.05 1.1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(loc16, replace)
		
graph export ".\Cognitive_Results\G1_LoCage16Results.pdf", replace


** Plot for emotion recogition 'triangles' task at age 13

* Min and max x-axis values
sum lci uci if exposure == "emoRec_triangles_age13" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "emoRec_triangles_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "emoRec_triangles_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "emoRec_triangles_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "emoRec_triangles_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "emoRec_triangles_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "emoRec_triangles_age13", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Emotion Recognition Age 13 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.98 1 1.02 1.04, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(triangles13, replace)
		
graph export ".\Cognitive_Results\G1_triangles13Results.pdf", replace


** Plot for Skuse social cognition at age 8

* Min and max x-axis values
sum lci uci if exposure == "skuseSocCog_age8" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "skuseSocCog_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "skuseSocCog_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "skuseSocCog_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "skuseSocCog_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "skuseSocCog_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "skuseSocCog_age8", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Social Cognition Age 8 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.96 0.98 1 1.02 1.04 1.06 1.08, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(skuse8, replace)
		
graph export ".\Cognitive_Results\G1_skuse8Results.pdf", replace


** Plot for Skuse social cognition at age 16

* Min and max x-axis values
sum lci uci if exposure == "skuseSocCog_age16" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "skuseSocCog_age16", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "skuseSocCog_age16", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "skuseSocCog_age16", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "skuseSocCog_age16", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "skuseSocCog_age16", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "skuseSocCog_age16", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Social Cognition Age 16 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.96 0.98 1 1.02 1.04 1.06 1.08, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(skuse16, replace)
		
graph export ".\Cognitive_Results\G1_skuse16Results.pdf", replace


** Plot for SDQ prosocial sub-scale at age 8

* Min and max x-axis values
sum lci uci if exposure == "SDQ_prosocial_age8" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "SDQ_prosocial_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "SDQ_prosocial_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "SDQ_prosocial_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "SDQ_prosocial_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "SDQ_prosocial_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "SDQ_prosocial_age8", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("SDQ - Prosocial Age 8 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.9 0.95 1 1.05 1.1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(prosocial8, replace)
		
graph export ".\Cognitive_Results\G1_prosocial8Results.pdf", replace


** Plot for SDQ prosocial sub-scale at age 13

* Min and max x-axis values
sum lci uci if exposure == "SDQ_prosocial_age13" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "SDQ_prosocial_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "SDQ_prosocial_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "SDQ_prosocial_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "SDQ_prosocial_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "SDQ_prosocial_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "SDQ_prosocial_age13", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("SDQ - Prosocial Age 13 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.9 1 1.1 1.2 1.3, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(prosocial13, replace)
		
graph export ".\Cognitive_Results\G1_prosocial13Results.pdf", replace


** Plot for SDQ prosocial sub-scale at age 25

* Min and max x-axis values
sum lci uci if exposure == "SDQ_prosocial_age25" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "SDQ_prosocial_age25", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "SDQ_prosocial_age25", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "SDQ_prosocial_age25", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "SDQ_prosocial_age25", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "SDQ_prosocial_age25", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "SDQ_prosocial_age25", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("SDQ - Prosocial Age 25 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.95 1 1.05 1.1 1.15 1.2, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(prosocial25, replace)
		
graph export ".\Cognitive_Results\G1_prosocial25Results.pdf", replace



** Plot for Bachman self-esteem

* Min and max x-axis values
sum lci uci if exposure == "esteem_bachman_age17" & outcome_level != "NA"

twoway (scatter level_num coef if outcome == "Belief" & exposure == "esteem_bachman_age17", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Belief" & exposure == "esteem_bachman_age17", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Relig" & exposure == "esteem_bachman_age17", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Relig" & exposure == "esteem_bachman_age17", ///
			horizontal col(black)) ///
		(scatter level_num coef if outcome == "Attend" & exposure == "esteem_bachman_age17", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci uci level_num if outcome == "Attend" & exposure == "esteem_bachman_age17", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Bachman Self-Esteem Age 17 and RSBB", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.96 0.98 1 1.02 1.04, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(esteem17, replace)
		
graph export ".\Cognitive_Results\G1_selfEsteem17Results.pdf", replace


*** Also make a few interaction plots

** Plot for total IQ at age 15 and sex interaction

* Min and max x-axis values
sum lci_int uci_int if exposure == "totalIQ_age15" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == "totalIQ_age15", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == "totalIQ_age15", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == "totalIQ_age15", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == "totalIQ_age15", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == "totalIQ_age15", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == "totalIQ_age15", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Total IQ Age 15 by Sex interaction", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.94 0.96 0.98 1 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(iq15_int, replace)
		
graph export ".\Cognitive_Results\G1_iq15BySexInt.pdf", replace


** Plot for vocab task at age 24 and sex interaction

* Min and max x-axis values
sum lci_int uci_int if exposure == "vocab_age24" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == "vocab_age24", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == "vocab_age24", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == "vocab_age24", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == "vocab_age24", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == "vocab_age24", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == "vocab_age24", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Vocab Task Age 24 by Sex interaction", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.8 0.9 1 1.1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(vocab24_int, replace)
		
graph export ".\Cognitive_Results\G1_vocab24BySexInt.pdf", replace


** Plot for total IQ at age 8 and sex interaction

* Min and max x-axis values
sum lci_int uci_int if exposure == "totalIQ_age8" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == "totalIQ_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == "totalIQ_age8", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == "totalIQ_age8", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == "totalIQ_age8", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Total IQ Age 8 by Sex interaction", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.97 0.98 0.99 1 1.01 1.02, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(iq8_int, replace)
		
graph export ".\Cognitive_Results\G1_iq8BySexInt.pdf", replace


** Plot for personality agreeableness at age 13 and sex interaction

* Min and max x-axis values
sum lci_int uci_int if exposure == "agreeableness_age13" & outcome_level != "NA"

twoway (scatter level_num coef_int if outcome == "Belief" & exposure == "agreeableness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Belief" & exposure == "agreeableness_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Relig" & exposure == "agreeableness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Relig" & exposure == "agreeableness_age13", ///
			horizontal col(black)) ///
		(scatter level_num coef_int if outcome == "Attend" & exposure == "agreeableness_age13", ///
			col(black) msize(small) msym(D)) ///
		(rspike lci_int uci_int level_num if outcome == "Attend" & exposure == "agreeableness_age13", ///
			horizontal col(black)), ///
		yscale(reverse)	ytitle("") xtitle("Relative risk ratio") ///
		title("Agreeableness Age 13 by Sex interaction", size(medium)) ///
		xline(1, lcol(black) lpattern(shortdash)) xscale(log) ///
		xlabel(0.85 0.9 0.95 1, labsize(small)) ///
		ylabel(0 1 3 4 6 7 8, valuelabel labsize(small) angle(0)) ///
		legend(off) name(vocab13_int, replace)
		
graph export ".\Cognitive_Results\G1_agree13BySexInt.pdf", replace

graph close _all


*** Predicted probability plots (as multinomial relative risk ratio results not necessarily intuitive to interpret)

** Make one for IQ at age 8
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3000 ageAt28 male totalIQ_age8, rrr baseoutcome(3)
margins, at(totalIQ_age8 = (70(1)130))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen totalIQ_age8 = fill(70 71)
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

twoway (line prob_yes totalIQ_age8, col(black)) ///
	(rarea lci_yes uci_yes totalIQ_age8, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure totalIQ_age8, col(red)) ///
	(rarea lci_notSure uci_notSure totalIQ_age8, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no totalIQ_age8, col(blue)) ///
	(rarea lci_no uci_no totalIQ_age8, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Total IQ age 8") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Yes" 3 "Not sure" 5 "No") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(iq8_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religious affiliation
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3040_grp ageAt28 male totalIQ_age8, rrr baseoutcome(3)
margins, at(totalIQ_age8 = (70(1)130))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen totalIQ_age8 = fill(70 71)
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

twoway (line prob_xian totalIQ_age8, col(black)) ///
	(rarea lci_xian uci_xian totalIQ_age8, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other totalIQ_age8, col(red)) ///
	(rarea lci_other uci_other totalIQ_age8, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none totalIQ_age8, col(blue)) ///
	(rarea lci_none uci_none totalIQ_age8, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Total IQ age 8") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(iq8_relig, replace)
	
	
* Attend church
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3080_rev ageAt28 male totalIQ_age8, rrr baseoutcome(0)
margins, at(totalIQ_age8 = (70(1)130))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen totalIQ_age8 = fill(70 71)
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

twoway (line prob_no totalIQ_age8, col(black)) ///
	(rarea lci_no uci_no totalIQ_age8, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr totalIQ_age8, col(red)) ///
	(rarea lci_yr uci_yr totalIQ_age8, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth totalIQ_age8, col(blue)) ///
	(rarea lci_mth uci_mth totalIQ_age8, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk totalIQ_age8, col(green)) ///
	(rarea lci_wk uci_wk totalIQ_age8, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Total IQ age 8") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Not at all" 3 "Occasionally" 5 "1/yr" 7 "1/mth") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(iq8_attend, replace)
	
	
* Combine these plots together
graph combine iq8_bel iq8_relig iq8_attend, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G1_IQ8PredProbs_combined.pdf", replace

graph close _all


*** Repeat for IQ at age 15
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3000 ageAt28 male totalIQ_age15, rrr baseoutcome(3)
margins, at(totalIQ_age15 = (70(1)130))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen totalIQ_age15 = fill(70 71)
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

twoway (line prob_yes totalIQ_age15, col(black)) ///
	(rarea lci_yes uci_yes totalIQ_age15, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure totalIQ_age15, col(red)) ///
	(rarea lci_notSure uci_notSure totalIQ_age15, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no totalIQ_age15, col(blue)) ///
	(rarea lci_no uci_no totalIQ_age15, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Total IQ age 15") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Yes" 3 "Not sure" 5 "No") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(iq15_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religious affiliation
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3040_grp ageAt28 male totalIQ_age15, rrr baseoutcome(3)
margins, at(totalIQ_age15 = (70(1)130))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen totalIQ_age15 = fill(70 71)
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

twoway (line prob_xian totalIQ_age15, col(black)) ///
	(rarea lci_xian uci_xian totalIQ_age15, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other totalIQ_age15, col(red)) ///
	(rarea lci_other uci_other totalIQ_age15, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none totalIQ_age15, col(blue)) ///
	(rarea lci_none uci_none totalIQ_age15, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Total IQ age 15") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(iq15_relig, replace)
	
	
* Attend church
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3080_rev ageAt28 male totalIQ_age15, rrr baseoutcome(0)
margins, at(totalIQ_age15 = (70(1)130))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen totalIQ_age15 = fill(70 71)
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

twoway (line prob_no totalIQ_age15, col(black)) ///
	(rarea lci_no uci_no totalIQ_age15, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr totalIQ_age15, col(red)) ///
	(rarea lci_yr uci_yr totalIQ_age15, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth totalIQ_age15, col(blue)) ///
	(rarea lci_mth uci_mth totalIQ_age15, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk totalIQ_age15, col(green)) ///
	(rarea lci_wk uci_wk totalIQ_age15, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Total IQ age 15") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Not at all" 3 "Occasionally" 5 "1/yr" 7 "1/mth") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(iq15_attend, replace)
	
	
* Combine these plots together
graph combine iq15_bel iq15_relig iq15_attend, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G1_IQ15PredProbs_combined.pdf", replace

graph close _all


** Also plot interaction between sex and total IQ at age 15 - Start with religious belief
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum totalIQ_age15
sum ageAt28
local age = r(mean)

mlogit YPG3000 ageAt28 i.male##c.totalIQ_age15, rrr baseoutcome(3)
margins, at(totalIQ_age15 = (70(1)130) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen totalIQ_age15 = fill(70 70 71 71)
gen YPG3000 = 1
predict p1, outcome(1)
sum p1

replace YPG3000 = 2
predict p2, outcome(2)
sum p1 p2

replace YPG3000 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 totalIQ_age15 if male == 0, col(black)) ///
	(line p1 totalIQ_age15 if male == 1, col(black) lpattern(dash)) ///
	(line p2 totalIQ_age15 if male == 0, col(red)) ///
	(line p2 totalIQ_age15 if male == 1, col(red) lpattern(dash)) ///
	(line p3 totalIQ_age15 if male == 0, col(blue)) ///
	(line p3 totalIQ_age15 if male == 1, col(blue) lpattern(dash)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("IQ at Age 15") ytitle("Predicted probability", margin(small)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Female - Yes" 3 "Female - Not sure" 5 "Female - No" ///
	2 "Male - Yes" 4 "Male - Not sure" 6 "Male - No") cols(3) ///
	colgap(*.5) keygap(*.5)  size(small)) ///
	name(iq15_belief_int, replace)

	
* Religious affiliation
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum totalIQ_age15
sum ageAt28
local age = r(mean)

mlogit YPG3040_grp ageAt28 i.male##c.totalIQ_age15, rrr baseoutcome(3)
margins, at(totalIQ_age15 = (70(1)130) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen totalIQ_age15 = fill(70 70 71 71)
gen YPG3040_grp = 1
predict p1, outcome(1)
sum p1

replace YPG3040_grp = 2
predict p2, outcome(2)
sum p1 p2

replace YPG3040_grp = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 totalIQ_age15 if male == 0, col(black)) ///
	(line p1 totalIQ_age15 if male == 1, col(black) lpattern(dash)) ///
	(line p2 totalIQ_age15 if male == 0, col(red)) ///
	(line p2 totalIQ_age15 if male == 1, col(red) lpattern(dash)) ///
	(line p3 totalIQ_age15 if male == 0, col(blue)) ///
	(line p3 totalIQ_age15 if male == 1, col(blue) lpattern(dash)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("IQ at Age 15") ytitle("Predicted probability", margin(small)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Female - Christian" 3 "Female - Other" 5 "Female - None" ///
	2 "Male - Christian" 4 "Male - Other" 6 "Male - None") cols(3) ///
	colgap(*.5) keygap(*.5) size(small)) ///
	name(iq15_relig_int, replace)
	

* Religious attendance
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum totalIQ_age15
sum ageAt28
local age = r(mean)

mlogit YPG3080_rev ageAt28 i.male##c.totalIQ_age15, rrr baseoutcome(0)
margins, at(totalIQ_age15 = (70(1)130) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen totalIQ_age15 = fill(70 70 71 71)
gen YPG3080_rev = 0
predict p1, outcome(0)
sum p1

replace YPG3080_rev = 1
predict p2, outcome(1)
sum p1 p2

replace YPG3080_rev = 2
predict p3, outcome(2)
sum p1 p2 p3

replace YPG3080_rev = 3
predict p4, outcome(3)
sum p1 p2 p3 p4

twoway (line p1 totalIQ_age15 if male == 0, col(black)) ///
	(line p1 totalIQ_age15 if male == 1, col(black) lpattern(dash)) ///
	(line p2 totalIQ_age15 if male == 0, col(red)) ///
	(line p2 totalIQ_age15 if male == 1, col(red) lpattern(dash)) ///
	(line p3 totalIQ_age15 if male == 0, col(blue)) ///
	(line p3 totalIQ_age15 if male == 1, col(blue) lpattern(dash)) ///
	(line p4 totalIQ_age15 if male == 0, col(green)) ///
	(line p4 totalIQ_age15 if male == 1, col(green) lpattern(dash)), ///
	xscale(range(70 130)) xlabel(70(10)130, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("IQ at Age 15") ytitle("Predicted probability", margin(small)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Female - Never" 2 "Male - Never" 3 "Female - Occasionally" ///
	4 "Male - Occasionally" 5 "Female - 1/Yr" 6 "Male - 1/Yr" ///
	7 "Female - 1/Mth" 8 "Male - 1/Mth") cols(2) colgap(*.5) ///
	keygap(*.5) size(small)) ///
	name(iq15_attend_int, replace)
	
	
* Combine these plots together
graph combine iq15_belief_int iq15_relig_int iq15_attend_int, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G1_IQ15PredProbs_sexInt_combined.pdf", replace

graph close _all


*** Repeat for vocabulary task at age 24
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum vocab_age24 

mlogit YPG3000 ageAt28 male vocab_age24, rrr baseoutcome(3)
margins, at(vocab_age24 = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen vocab_age24 = fill(0 1)
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

twoway (line prob_yes vocab_age24, col(black)) ///
	(rarea lci_yes uci_yes vocab_age24, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure vocab_age24, col(red)) ///
	(rarea lci_notSure uci_notSure vocab_age24, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no vocab_age24, col(blue)) ///
	(rarea lci_no uci_no vocab_age24, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Vocab Task Age 24") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Yes" 3 "Not sure" 5 "No") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(vocab24_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religious affiliation
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3040_grp ageAt28 male vocab_age24, rrr baseoutcome(3)
margins, at(vocab_age24 = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen vocab_age24 = fill(0 1)
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

twoway (line prob_xian vocab_age24, col(black)) ///
	(rarea lci_xian uci_xian vocab_age24, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other vocab_age24, col(red)) ///
	(rarea lci_other uci_other vocab_age24, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none vocab_age24, col(blue)) ///
	(rarea lci_none uci_none vocab_age24, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Vocab Task Age 24") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(vocab24_relig, replace)
	
	
* Attend church
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3080_rev ageAt28 male vocab_age24, rrr baseoutcome(0)
margins, at(vocab_age24 = (0(1)12))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen vocab_age24 = fill(0 1)
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

twoway (line prob_no vocab_age24, col(black)) ///
	(rarea lci_no uci_no vocab_age24, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr vocab_age24, col(red)) ///
	(rarea lci_yr uci_yr vocab_age24, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth vocab_age24, col(blue)) ///
	(rarea lci_mth uci_mth vocab_age24, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk vocab_age24, col(green)) ///
	(rarea lci_wk uci_wk vocab_age24, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Vocab Task Age 24") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Not at all" 3 "Occasionally" 5 "1/yr" 7 "1/mth") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(vocab24_attend, replace)
	
	
* Combine these plots together
graph combine vocab24_bel vocab24_relig vocab24_attend, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G1_vocab24PredProbs_combined.pdf", replace

graph close _all


** Also plot interaction between sex and vocab score at age 24 - Start with religious belief
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum vocab_age24
sum ageAt28
local age = r(mean)

mlogit YPG3000 ageAt28 i.male##c.vocab_age24, rrr baseoutcome(3)
margins, at(vocab_age24 = (0(1)12) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen vocab_age24 = fill(0 0 1 1)
gen YPG3000 = 1
predict p1, outcome(1)
sum p1

replace YPG3000 = 2
predict p2, outcome(2)
sum p1 p2

replace YPG3000 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 vocab_age24 if male == 0, col(black)) ///
	(line p1 vocab_age24 if male == 1, col(black) lpattern(dash)) ///
	(line p2 vocab_age24 if male == 0, col(red)) ///
	(line p2 vocab_age24 if male == 1, col(red) lpattern(dash)) ///
	(line p3 vocab_age24 if male == 0, col(blue)) ///
	(line p3 vocab_age24 if male == 1, col(blue) lpattern(dash)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Vocab Task Age 24") ytitle("Predicted probability", margin(small)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Female - Yes" 3 "Female - Not sure" 5 "Female - No" ///
	2 "Male - Yes" 4 "Male - Not sure" 6 "Male - No") cols(3) ///
	colgap(*.5) keygap(*.5) size(small)) ///
	name(vocab24_belief_int, replace)

	
* Religious affiliation
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum vocab_age24
sum ageAt28
local age = r(mean)

mlogit YPG3040_grp ageAt28 i.male##c.vocab_age24, rrr baseoutcome(3)
margins, at(vocab_age24 = (0(1)12) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen vocab_age24 = fill(0 0 1 1)
gen YPG3040_grp = 1
predict p1, outcome(1)
sum p1

replace YPG3040_grp = 2
predict p2, outcome(2)
sum p1 p2

replace YPG3040_grp = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 vocab_age24 if male == 0, col(black)) ///
	(line p1 vocab_age24 if male == 1, col(black) lpattern(dash)) ///
	(line p2 vocab_age24 if male == 0, col(red)) ///
	(line p2 vocab_age24 if male == 1, col(red) lpattern(dash)) ///
	(line p3 vocab_age24 if male == 0, col(blue)) ///
	(line p3 vocab_age24 if male == 1, col(blue) lpattern(dash)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Vocab Task Age 24") ytitle("Predicted probability", margin(small)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Female - Christian" 3 "Female - Other" 5 "Female - None" ///
	2 "Male - Christian" 4 "Male - Other" 6 "Male - None") cols(3) ///
	colgap(*.5) keygap(*.5) size(small)) ///
	name(vocab24_relig_int, replace)
	

* Religious attendance
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum vocab_age24
sum ageAt28
local age = r(mean)

mlogit YPG3080_rev ageAt28 i.male##c.vocab_age24, rrr baseoutcome(0)
margins, at(vocab_age24 = (0(1)12) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen vocab_age24 = fill(0 0 1 1)
gen YPG3080_rev = 0
predict p1, outcome(0)
sum p1

replace YPG3080_rev = 1
predict p2, outcome(1)
sum p1 p2

replace YPG3080_rev = 2
predict p3, outcome(2)
sum p1 p2 p3

replace YPG3080_rev = 3
predict p4, outcome(3)
sum p1 p2 p3 p4

twoway (line p1 vocab_age24 if male == 0, col(black)) ///
	(line p1 vocab_age24 if male == 1, col(black) lpattern(dash)) ///
	(line p2 vocab_age24 if male == 0, col(red)) ///
	(line p2 vocab_age24 if male == 1, col(red) lpattern(dash)) ///
	(line p3 vocab_age24 if male == 0, col(blue)) ///
	(line p3 vocab_age24 if male == 1, col(blue) lpattern(dash)) ///
	(line p4 vocab_age24 if male == 0, col(green)) ///
	(line p4 vocab_age24 if male == 1, col(green) lpattern(dash)), ///
	xscale(range(0 12)) xlabel(0(1)12, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Vocab Task Age 24") ytitle("Predicted probability", margin(small)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Female - Never" 2 "Male - Never" 3 "Female - Occasionally" ///
	4 "Male - Occasionally" 5 "Female - 1/Yr" 6 "Male - 1/Yr" ///
	7 "Female - 1/Mth" 8 "Male - 1/Mth") cols(2) colgap(*.5) keygap(*.5) ///
	size(small)) ///
	name(vocab24_attend_int, replace)
	
	
* Combine these plots together
graph combine vocab24_belief_int vocab24_relig_int vocab24_attend_int, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G1_vocab24PredProbs_sexInt_combined.pdf", replace

graph close _all



*** Repeat for agreeableness at age 13
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum agreeableness_age13 

mlogit YPG3000 ageAt28 male agreeableness_age13, rrr baseoutcome(3)
margins, at(agreeableness_age13 = (15(1)50))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen agreeableness_age13 = fill(15 16)
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

twoway (line prob_yes agreeableness_age13, col(black)) ///
	(rarea lci_yes uci_yes agreeableness_age13, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_notSure agreeableness_age13, col(red)) ///
	(rarea lci_notSure uci_notSure agreeableness_age13, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_no agreeableness_age13, col(blue)) ///
	(rarea lci_no uci_no agreeableness_age13, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(15 50)) xlabel(15(5)50, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Agreeableness Age 13") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Yes" 3 "Not sure" 5 "No") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(agree13_bel, replace)
	
	
** Repeat for other RSBB outcomes and combine plots together

* Religious affiliation
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3040_grp ageAt28 male agreeableness_age13, rrr baseoutcome(3)
margins, at(agreeableness_age13 = (15(1)50))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen agreeableness_age13 = fill(15 16)
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

twoway (line prob_xian agreeableness_age13, col(black)) ///
	(rarea lci_xian uci_xian agreeableness_age13, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_other agreeableness_age13, col(red)) ///
	(rarea lci_other uci_other agreeableness_age13, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_none agreeableness_age13, col(blue)) ///
	(rarea lci_none uci_none agreeableness_age13, lcol(black) lwidth(vthin) fcol(blue%20)), ///
	xscale(range(15 16)) xlabel(15(5)50, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Agreeableness Age 13") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Christian" 3 "Other" 5 "None") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(agree13_relig, replace)
	
	
* Attend church
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

mlogit YPG3080_rev ageAt28 male agreeableness_age13, rrr baseoutcome(0)
margins, at(agreeableness_age13 = (15(1)50))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen agreeableness_age13 = fill(15 16)
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

twoway (line prob_no agreeableness_age13, col(black)) ///
	(rarea lci_no uci_no agreeableness_age13, lcol(black) lwidth(vthin) fcol(black%20)) ///
	(line prob_yr agreeableness_age13, col(red)) ///
	(rarea lci_yr uci_yr agreeableness_age13, lcol(black) lwidth(vthin) fcol(red%20)) ///
	(line prob_mth agreeableness_age13, col(blue)) ///
	(rarea lci_mth uci_mth agreeableness_age13, lcol(black) lwidth(vthin) fcol(blue%20)) ///
	(line prob_wk agreeableness_age13, col(green)) ///
	(rarea lci_wk uci_wk agreeableness_age13, lcol(black) lwidth(vthin) fcol(green%20)), ///
	xscale(range(15 50)) xlabel(15(5)50, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Agreeableness Age 13") ytitle("Predicted probability") yscale(titlegap(2)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Not at all" 3 "Occasionally" 5 "1/yr" 7 "1/mth") ///
	rows(1) size(small) symxsize(*0.5)) ///
	name(agree13_attend, replace)
	
	
* Combine these plots together
graph combine agree13_bel agree13_relig agree13_attend, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G1_agree13PredProbs_combined.pdf", replace

graph close _all


** Also plot interaction between sex and agreeableness at age 13 - Start with religious belief
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum agreeableness_age13
sum ageAt28
local age = r(mean)

mlogit YPG3000 ageAt28 i.male##c.agreeableness_age13, rrr baseoutcome(3)
margins, at(agreeableness_age13 = (15(1)50) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen agreeableness_age13 = fill(15 15 16 16)
gen YPG3000 = 1
predict p1, outcome(1)
sum p1

replace YPG3000 = 2
predict p2, outcome(2)
sum p1 p2

replace YPG3000 = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 agreeableness_age13 if male == 0, col(black)) ///
	(line p1 agreeableness_age13 if male == 1, col(black) lpattern(dash)) ///
	(line p2 agreeableness_age13 if male == 0, col(red)) ///
	(line p2 agreeableness_age13 if male == 1, col(red) lpattern(dash)) ///
	(line p3 agreeableness_age13 if male == 0, col(blue)) ///
	(line p3 agreeableness_age13 if male == 1, col(blue) lpattern(dash)), ///
	xscale(range(15 50)) xlabel(15(5)50, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Agreeableness Age 13") ytitle("Predicted probability", margin(small)) ///
	title("Religious belief", size(large)) ///
	legend(order(1 "Female - Yes" 3 "Female - Not sure" 5 "Female - No" ///
	2 "Male - Yes" 4 "Male - Not sure" 6 "Male - No") cols(3) ///
	colgap(*.5) keygap(*.5) size(small)) ///
	name(agree13_belief_int, replace)

	
* Religious affiliation
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum agreeableness_age13
sum ageAt28
local age = r(mean)

mlogit YPG3040_grp ageAt28 i.male##c.agreeableness_age13, rrr baseoutcome(3)
margins, at(agreeableness_age13 = (15(1)50) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/3

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen agreeableness_age13 = fill(15 15 16 16)
gen YPG3040_grp = 1
predict p1, outcome(1)
sum p1

replace YPG3040_grp = 2
predict p2, outcome(2)
sum p1 p2

replace YPG3040_grp = 3
predict p3, outcome(3)
sum p1 p2 p3

twoway (line p1 agreeableness_age13 if male == 0, col(black)) ///
	(line p1 agreeableness_age13 if male == 1, col(black) lpattern(dash)) ///
	(line p2 agreeableness_age13 if male == 0, col(red)) ///
	(line p2 agreeableness_age13 if male == 1, col(red) lpattern(dash)) ///
	(line p3 agreeableness_age13 if male == 0, col(blue)) ///
	(line p3 agreeableness_age13 if male == 1, col(blue) lpattern(dash)), ///
	xscale(range(15 50)) xlabel(15(5)50, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Agreeableness Age 13") ytitle("Predicted probability", margin(small)) ///
	title("Religious affiliation", size(large)) ///
	legend(order(1 "Female - Christian" 3 "Female - Other" 5 "Female - None" ///
	2 "Male - Christian" 4 "Male - Other" 6 "Male - None") cols(3) ///
	colgap(*.5) keygap(*.5) size(small)) ///
	name(agree13_relig_int, replace)
	

* Religious attendance
use ".\Cognitive_Results\G1_CogPredictorsOfRSBB_B3911_postAnalysis.dta", clear

sum agreeableness_age13
sum ageAt28
local age = r(mean)

mlogit YPG3080_rev ageAt28 i.male##c.agreeableness_age13, rrr baseoutcome(0)
margins, at(agreeableness_age13 = (15(1)50) male = (0 1))

matrix res = r(table)
matrix list res

local n = colsof(res)/4

clear 
set obs `n'
egen ageAt28 = fill(`age' `age')
egen male = fill(0 1 0 1)
egen agreeableness_age13 = fill(15 15 16 16)
gen YPG3080_rev = 0
predict p1, outcome(0)
sum p1

replace YPG3080_rev = 1
predict p2, outcome(1)
sum p1 p2

replace YPG3080_rev = 2
predict p3, outcome(2)
sum p1 p2 p3

replace YPG3080_rev = 3
predict p4, outcome(3)
sum p1 p2 p3 p4

twoway (line p1 agreeableness_age13 if male == 0, col(black)) ///
	(line p1 agreeableness_age13 if male == 1, col(black) lpattern(dash)) ///
	(line p2 agreeableness_age13 if male == 0, col(red)) ///
	(line p2 agreeableness_age13 if male == 1, col(red) lpattern(dash)) ///
	(line p3 agreeableness_age13 if male == 0, col(blue)) ///
	(line p3 agreeableness_age13 if male == 1, col(blue) lpattern(dash)) ///
	(line p4 agreeableness_age13 if male == 0, col(green)) ///
	(line p4 agreeableness_age13 if male == 1, col(green) lpattern(dash)), ///
	xscale(range(15 50)) xlabel(15(5)50, labsize(small)) ylabel(, labsize(small)) ///
	xtitle("Agreeableness Age 13") ytitle("Predicted probability", margin(small)) ///
	title("Religious attendance", size(large)) ///
	legend(order(1 "Female - Never" 2 "Male - Never" 3 "Female - Occasionally" ///
	4 "Male - Occasionally" 5 "Female - 1/Yr" 6 "Male - 1/Yr" ///
	7 "Female - 1/Mth" 8 "Male - 1/Mth") cols(2) ///
	colgap(*.5) keygap(*.5) size(small)) ///
	name(agree13_attend_int, replace)
	
	
* Combine these plots together
graph combine agree13_belief_int agree13_relig_int agree13_attend_int, iscale(0.5) rows(2)

graph export ".\Cognitive_Results\G1_agree13PredProbs_sexInt_combined.pdf", replace

graph close _all



