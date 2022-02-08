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
corr logicMemory-selfEsteem

matrix cor_cog = r(C)

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
corr IPSM_interpAware-selfEsteem

matrix cor_cog = r(C)

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

* Outcomes
foreach var of varlist YPG3000-YPG3080_OccYr {
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
misstable sum YPG3000-YPG3080_OccYr, all

* Exposures
misstable sum male YPG8000 f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126, all



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


* Associations between cognitive/psychological variables - Then make heat map of correlations (heatplot code adapted from: https://www.stata.com/meeting/germany19/slides/germany19_Jann.pdf)
corr verbalIQ_age8-globalEsteem_age8

matrix cor_cog = r(C)

heatplot cor_cog, color(hcl, diverging intensity(1)) ///
	lower nodiagonal cuts(-1.05(0.1)1.05) xlabel(, angle(45) labsize(vsmall)) ///
	ylabel(, labsize(vsmall)) legend(subtitle(""))

* Save heatmap
graph export ".\Cognitive_Results\G1_corr_heatplot.pdf", replace

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
*** Next, to start making plots of these results

** p-value plots

** Pseudo-R2 value plots

** Coefficient plots

** Predicted probability plots
