*** Predictors of RSBB (B3911) - Data processing script
*** Created 15/11/2021 by Dan Smith
*** Stata v16.0

*** This script reads in the raw ALSPAC data, cleans the variables, and then creates the three datasets used in subsequent analyses.


**********************************************************************************
**** Set working directory, start a log file, read in dataset, and add numlabels

cd "X:\Groups\ARC\DanS\Descriptive_PredictorsOfRSBB_B3911"

capture log close
log using "Desc_RSBB_B3911_DataCleaning_log", replace text

* Make sure no frames in memory
clear frames

use "Desc_RSBB_B3911.dta", clear

numlabel, add


**********************************************************************************
*** Process and tidy all of the variables used in this analysis

*** Start with RSBB outcome variables (note that for G0s belief in God, religious affiliation and church attendance come from pregnancy questionnaires [and only have ~20% missing data for mothers; ~40% for partners/fathers], while the intrinsic, extrinsic and total religiosity variables come from when the study child was aged 28 [so have ~70% missing data for mothers; ~85% for partners/fathers]. As the sources of this data differs, will also analyse G0 data for belief in God, religious affiliation and church attendance at age 28 as well, to compare results. G1 data all comes from the same questionnaire asked at age 28.)

*** G0 mother

** Belief in God/divine power (pregnancy and when study child aged 28) - Code negative values as missing (for age 28 data, recode so consistent with pregnancy data)
tab d810, m

replace d810 = . if d810 < 0
tab d810, m
tab d810

tab Y3000, m

replace Y3000 = . if Y3000 < 0
tab Y3000, m

recode Y3000 (1 = 1) (0 = 3) (2 = 2)
label values Y3000 d810
tab Y3000, m
tab Y3000


** Religious affiliation (pregnancy and when study child aged 28) - Code negative values as missing - Recode to none vs Christian vs other
tab d813, m

replace d813 = . if d813 < 0
tab d813, m

recode d813 (0 = 1) (1/6 = 2) (7/13 = 3), gen(d813_grp)
label define relig_lb 1 "None" 2 "Christian" 3 "Other"
numlabel relig_lb, add
label values d813_grp relig_lb
tab d813_grp, m
tab d813_grp

tab Y3040, m

replace Y3040 = . if Y3040 < 0
tab Y3040, m

recode Y3040 (13 = 1) (1/6 = 2) (7/12 14 = 3), gen(Y3040_grp)
label values Y3040_grp relig_lb
tab Y3040_grp, m
tab Y3040_grp


** Attendance at church/place of worship (pregnancy and when study child aged 28) - Code negative values as missing
tab d816, m

replace d816 = . if d816 < 0
tab d816, m
tab d816

* Age 28 data includes an 'occasionally' response, which pregnancy Q did not have (so data not completely comparable). Will code in three different ways: 1) Leave as is (with additional 'occasionally' category); 2) Code 'occasionally' with 'not at all'; and 3) code 'accasionally' with 'at least once a year'.
tab Y3080, m

replace Y3080 = . if Y3080 < 0
tab Y3080, m

recode Y3080 (0 = 4), gen(Y3080_OccNever)
label define occNever_lb 1 "At least once a week" 2 "At least once a month" 3 "At least once a year" 4 "Occasionally/Not at all"
numlabel occNever_lb, add
label values Y3080_OccNever occNever_lb
tab Y3080_OccNever, m
tab Y3080_OccNever

recode Y3080 (0 = 4) (4 = 3), gen(Y3080_OccYr)
label define occYr_lb 1 "At least once a week" 2 "At least once a month" 3 "At least once a year/Occasionally" 4 "Not at all"
numlabel occYr_lb, add
label values Y3080_OccYr occYr_lb
tab Y3080_OccYr, m
tab Y3080_OccYr

* Also recode 'not at all' so categories are in order
recode Y3080 (0 = 5)
label define attend_lb 1 "At least once a week" 2 "At least once a month" 3 "At least once a year" 4 "Occasionally" 5 "Not at all"
numlabel attend_lb, add
label values Y3080 attend_lb
tab Y3080, m
tab Y3080


** Intrinsic religiosity (age 28 only) - Taken from a three-item subscale of the Duke University Religion Scale (DUREL) asking about experiencing the presence of the divine (e.g., God), whether their religious beliefs lie behind their whole approach to life, and whether they try hard to carry religious over into all other dealings in life. 
tab Y3153, m

* Code negative values as missing
replace Y3153 = . if Y3153 < 0
tab Y3153, m
tab Y3153

sum Y3153

* As distribution very non-normal (spike at 3, then uniform for all higher values), will code into categories and treat as ordinal variable in sensitivity analyses
hist Y3153, freq width(1)

recode Y3153 (3 = 1) (4/7 = 2) (8/11 = 3) (12/15 = 4), gen(Y3153_cat)
label define intrinsic_lb 1 "Lowest IR (3)" 2 "Moderate IR (4-7)" 3 "High IR (8-11)" 4 "Highest IR (12-15)"
numlabel intrinsic_lb, add
label values Y3153_cat intrinsic_lb
tab Y3153_cat, m
tab Y3153_cat


** Extrinsic religiosity (age 28 only) - Based on two extrinsically-loaded questions from the intrinsic/extrinsic religiosity scale developed by Gorsuch and MacPherson (1989) were used: one question was “I attend a place of worship because it helps me to make friends”, and the other asked “I pray mainly to gain relief to protection”. 
tab1 Y3160 Y3170, m

* Code negative values as missing
replace Y3160 = . if Y3160 < 0
tab Y3160, m

replace Y3170 = . if Y3170 < 0
tab Y3170, m

* How to deal with 'not applicable' responses? As the other responses (agree/disagree) assume that participants are religious/engage in these behaviours, these variables are both ordered and unordered categorical variables... Given these complications, I will analyse the two extrinsic outcomes separately
tab Y3160 Y3170, m

tab1 Y3160 Y3170, m
tab1 Y3160 Y3170


** Total DUREL religiosity score, which is a sum of three sub-components: organisational religious activity (attendance at church/place of worship), non-organisational religious activity (time spent in private religious activity, such as prayer) and the intrinsic religiosity score detailed above
tab Y3155, m

* Code negative values as missing
replace Y3155 = . if Y3155 < 0
tab Y3155, m
tab Y3155

sum Y3155

* As distribution very non-normal (spike at 3, then uniform for all higher values), will code into categories and treat as ordinal variable in sensitivity analyses
hist Y3155, freq width(1)

recode Y3155 (5 = 1) (6/10 = 2) (11/15 = 3) (16/20 = 4) (21/26 = 5), gen(Y3155_cat)
label define durel_lb 1 "Lowest DUREL (5)" 2 "6-10 DUREL" 3 "11-15 DUREL" 4 "16-20 DUREL"5 "Highest DUREL (21-26)"
numlabel durel_lb, add
label values Y3155_cat durel_lb
tab Y3155_cat, m


****** Now repeat for the G0 partner/father RSBB data

** Belief in God/divine power (pregnancy and when study child aged 28) - Code negative values as missing (for age 28 data, recode so consistent with pregnancy data)
tab pb150, m

replace pb150 = . if pb150 < 0
tab pb150, m
tab pb150

tab FC3000, m

replace FC3000 = . if FC3000 < 0
tab FC3000, m

recode FC3000 (1 = 1) (0 = 3) (2 = 2)
label values FC3000 pb150
tab FC3000, m
tab FC3000


** Religious affiliation (pregnancy and when study child aged 28) - Code negative values as missing - Recode to none vs Christian vs other
tab pb153, m

replace pb153 = . if pb153 < 0
tab pb153, m

recode pb153 (0 = 1) (1/6 = 2) (7/13 = 3), gen(pb153_grp)
label values pb153_grp relig_lb
tab pb153_grp, m
tab pb153_grp

tab FC3040, m

replace FC3040 = . if FC3040 < 0
tab FC3040, m

recode FC3040 (13 = 1) (1/6 = 2) (7/12 14 = 3), gen(FC3040_grp)
label values FC3040_grp relig_lb
tab FC3040_grp, m
tab FC3040_grp


** Attendance at church/place of worship (pregnancy and when study child aged 28) - Code negative values as missing
tab pb155, m

replace pb155 = . if pb155 < 0
tab pb155, m
tab pb155

* Age 28 data includes an 'occasionally' response, which pregnancy Q did not have (so data not completely comparable). Will code in three different ways: 1) Leave as is (with additional 'occasionally' category); 2) Code 'occasionally' with 'not at all'; and 3) code 'accasionally' with 'at least once a year'.
tab FC3080, m

replace FC3080 = . if FC3080 < 0
tab FC3080, m

recode FC3080 (0 = 4), gen(FC3080_OccNever)
label values FC3080_OccNever occNever_lb
tab FC3080_OccNever, m
tab FC3080_OccNever

recode FC3080 (0 = 4) (4 = 3), gen(FC3080_OccYr)
label values FC3080_OccYr occYr_lb
tab FC3080_OccYr, m
tab FC3080_OccYr

* Also recode 'not at all' so categories are in order
recode FC3080 (0 = 5)
label values FC3080 attend_lb
tab FC3080, m
tab FC3080


** Intrinsic religiosity (age 28 only) - Taken from a three-item subscale of the Duke University Religion Scale (DUREL) asking about experiencing the presence of the divine (e.g., God), whether their religious beliefs lie behind their whole approach to life, and whether they try hard to carry religious over into all other dealings in life. 
tab FC3153, m

* Code negative values as missing
replace FC3153 = . if FC3153 < 0
tab FC3153, m
tab FC3153

sum FC3153

* As distribution very non-normal (spike at 3, then uniform for all higher values), will code into categories and treat as ordinal variable in sensitivity analyses
hist FC3153, freq width(1)

recode FC3153 (3 = 1) (4/7 = 2) (8/11 = 3) (12/15 = 4), gen(FC3153_cat)
label values FC3153_cat intrinsic_lb
tab FC3153_cat, m
tab FC3153_cat


** Extrinsic religiosity (age 28 only) - Based on two extrinsically-loaded questions from the intrinsic/extrinsic religiosity scale developed by Gorsuch and MacPherson (1989) were used: one question was “I attend a place of worship because it helps me to make friends”, and the other asked “I pray mainly to gain relief to protection”. 
tab1 FC3160 FC3170, m

* Code negative values as missing
replace FC3160 = . if FC3160 < 0
tab FC3160, m

replace FC3170 = . if FC3170 < 0
tab FC3170, m

* How to deal with 'not applicable' responses? As the other responses (agree/disagree) assume that participants are religious/engage in these behaviours, these variables are both ordered and unordered categorical variables... Given these complications, I will analyse the two extrinsic outcomes separately
tab FC3160 FC3170, m

tab1 FC3160 FC3170, m
tab1 FC3160 FC3170


** Total DUREL religiosity score, which is a sum of three sub-components: organisational religious activity (attendance at church/place of worship), non-organisational religious activity (time spent in private religious activity, such as prayer) and the intrinsic religiosity score detailed above
tab FC3155, m

* Code negative values as missing
replace FC3155 = . if FC3155 < 0
tab FC3155, m
tab FC3155

sum FC3155

* As distribution very non-normal (spike at 3, then uniform for all higher values), will code into categories and treat as ordinal variable in sensitivity analyses
hist FC3155, freq width(1)

recode FC3155 (5 = 1) (6/10 = 2) (11/15 = 3) (16/20 = 4) (21/26 = 5), gen(FC3155_cat)
label values FC3155_cat durel_lb
tab FC3155_cat, m
tab FC3155_cat


****** Now repeat for the G1 study child RSBB data (all collected at age 28)

** Belief in God/divine power - Code negative values as missing (then recode so consistent with G0 data)
tab YPG3000, m

replace YPG3000 = . if YPG3000 < 0
tab YPG3000, m

recode YPG3000 (1 = 1) (0 = 3) (2 = 2)
label values YPG3000 pb150
tab YPG3000, m
tab YPG3000


** Religious affiliation - Code negative values as missing - Recode to none vs Christian vs other
tab YPG3040, m

replace YPG3040 = . if YPG3040 < 0
tab YPG3040, m

recode YPG3040 (13 = 1) (1/6 = 2) (7/12 14 = 3), gen(YPG3040_grp)
label values YPG3040_grp relig_lb
tab YPG3040_grp, m
tab YPG3040_grp


** Attendance at church/place of worship - Code negative values as missing

* Age 28 data includes an 'occasionally' response, which pregnancy Q did not have (so data not completely comparable). Will code in three different ways: 1) Leave as is (with additional 'occasionally' category); 2) Code 'occasionally' with 'not at all'; and 3) code 'accasionally' with 'at least once a year'.
tab YPG3080, m

replace YPG3080 = . if YPG3080 < 0
tab YPG3080, m

recode YPG3080 (0 = 4), gen(YPG3080_OccNever)
label values YPG3080_OccNever occNever_lb
tab YPG3080_OccNever, m
tab YPG3080_OccNever

recode YPG3080 (0 = 4) (4 = 3), gen(YPG3080_OccYr)
label values YPG3080_OccYr occYr_lb
tab YPG3080_OccYr, m
tab YPG3080_OccYr

* Also recode 'not at all' so categories are in order
recode YPG3080 (0 = 5)
label values YPG3080 attend_lb
tab YPG3080, m
tab YPG3080


** Intrinsic religiosity - Taken from a three-item subscale of the Duke University Religion Scale (DUREL) asking about experiencing the presence of the divine (e.g., God), whether their religious beliefs lie behind their whole approach to life, and whether they try hard to carry religious over into all other dealings in life. 
tab YPG3153, m

* Code negative values as missing
replace YPG3153 = . if YPG3153 < 0
tab YPG3153, m
tab YPG3153

sum YPG3153

* As distribution very non-normal (spike at 3, then uniform for all higher values), will code into categories and treat as ordinal variable in sensitivity analyses
hist YPG3153, freq width(1)

recode YPG3153 (3 = 1) (4/7 = 2) (8/11 = 3) (12/15 = 4), gen(YPG3153_cat)
label values YPG3153_cat intrinsic_lb
tab YPG3153_cat, m
tab YPG3153_cat


** Extrinsic religiosity - Based on two extrinsically-loaded questions from the intrinsic/extrinsic religiosity scale developed by Gorsuch and MacPherson (1989) were used: one question was “I attend a place of worship because it helps me to make friends”, and the other asked “I pray mainly to gain relief to protection”. 
tab1 YPG3160 YPG3170, m

* Code negative values as missing
replace YPG3160 = . if YPG3160 < 0
tab YPG3160, m

replace YPG3170 = . if YPG3170 < 0
tab YPG3170, m

* How to deal with 'not applicable' responses? As the other responses (agree/disagree) assume that participants are religious/engage in these behaviours, these variables are both ordered and unordered categorical variables... Given these complications, I will analyse the two extrinsic outcomes separately
tab YPG3160 YPG3170, m

tab1 YPG3160 YPG3170, m
tab1 YPG3160 YPG3170


** Total DUREL religiosity score, which is a sum of three sub-components: organisational religious activity (attendance at church/place of worship), non-organisational religious activity (time spent in private religious activity, such as prayer) and the intrinsic religiosity score detailed above
tab YPG3155, m

* Code negative values as missing
replace YPG3155 = . if YPG3155 < 0
tab YPG3155, m
tab YPG3155

sum YPG3155

* As distribution very non-normal (spike at 3, then uniform for all higher values), will code into categories and treat as ordinal variable in sensitivity analyses
hist YPG3155, freq width(1)

recode YPG3155 (5 = 1) (6/10 = 2) (11/15 = 3) (16/20 = 4) (21/26 = 5), gen(YPG3155_cat)
label values YPG3155_cat durel_lb
tab YPG3155_cat, m
tab YPG3155_cat


***************************************************************************************
****** Now process all the exposure variables

*** G0 mother exposures

** Demographic variables

* Age at birth
tab mz028b, m

replace mz028b = . if mz028b < 0
tab mz028b, m
sum mz028b

* Age at @28 RSBB questionnaire
tab Y9992, m

replace Y9992 = . if Y9992 < 0
tab Y9992, m
sum Y9992

* Ethnicity - Combine into white vs other than white
tab c800, m

replace c800 = . if c800 < 0
tab c800, m

recode c800 (1 = 0) (2/9 = 1), gen(c800_grp)
label define white_lb 0 "White" 1 "Other than white"
numlabel white_lb, add
label values c800_grp white_lb
tab c800_grp, m
tab c800_grp

* Also want to add in missing ethnicity data collect as part of COVID4 questionnaire
tab c800_grp covid4m_0502, m

replace c800_grp = 0 if covid4m_0502 == 1 & c800_grp == .
replace c800_grp = 1 if covid4m_0502 == 2 & c800_grp == .

tab c800_grp, m
tab c800_grp

* Marital status - Combine widowed/divorced/separated, as well as 1st and 2nd/3rd marriage
tab a525, m

replace a525 = . if a525 < 0
tab a525, m

recode a525 (1 = 1) (5 6 = 2) (2 3 4 = 3), gen(a525_grp)
label define marital_lb 1 "Never married" 2 "Married" 3 "Widowed/Divorced/Separated"
numlabel marital_lb, add
label values a525_grp marital_lb
tab a525_grp, m
tab a525_grp

* Residential mobility - Combine into 0/1/2/3/4/5 or more
tab a005, m

replace a005 = . if a005 < 0
tab a005, m

recode a005 (5/50 = 5), gen(a005_grp)
label define mobility_lb 5 "5 or more"
numlabel mobility_lb, add
label values a005_grp mobility_lb
tab a005_grp, m
tab a005_grp

* Urban/rural status - combine into urban vs rural
tab jan1993ur01ind_M, m

replace jan1993ur01ind_M = . if jan1993ur01ind_M < 0
tab jan1993ur01ind_M, m

recode jan1993ur01ind_M (1 = 1) (2 3 4 = 0), gen(jan1993ur01ind_grp)
label define urban_lb 1 "Urban" 0 "Town/village/Hamlet"
numlabel urban_lb, add
label value jan1993ur01ind_grp urban_lb
tab jan1993ur01ind_grp, m
tab jan1993ur01ind_grp

* Parity - Recode into 0, 1 and 2 or more
tab b032, m

replace b032 = . if b032 < 0
tab b032, m

recode b032 (3/22 = 2), gen(b032_grp)
label define parity_lb 2 "2 or more"
numlabel parity_lb, add
label values b032_grp parity_lb
tab b032_grp, m
tab b032_grp


** Socioeconomic factors/material insecurity

* Highest education
tab c645a, m

replace c645a = . if c645a < 0
tab c645a, m
tab c645a

* Mother's highest education
tab c686a, m

replace c686a = . if c686a < 0
tab c686a, m
tab c686a

* Father's highest education
tab c706a, m

replace c706a = . if c706a < 0
tab c706a, m
tab c706a

* Occupational social class - Will split into low (III manual, IV and V) and high (!, II and III non-manual))
tab c755, m

replace c755 = . if c755 < 0 | c755 == 65
tab c755, m

recode c755 (1 2 3 = 1) (4 5 6 = 0), gen(c755_grp)
label define occ_lb 1 "High (I/II/III non-manual)" 0 "Low (III manual/IV/V)"
numlabel occ_lb, add
label value c755_grp occ_lb
tab c755_grp, m
tab c755_grp

* Maternal occupational social class - Will split into low (III manual, IV and V) and high (!, II and III non-manual))
tab c_sc_mgm, m

replace c_sc_mgm = . if c_sc_mgm < 0 | c_sc_mgm == 65
tab c_sc_mgm, m

recode c_sc_mgm (1 2 3 = 1) (4 5 6 = 0), gen(c_sc_mgm_grp)
label value c_sc_mgm_grp occ_lb
tab c_sc_mgm_grp, m
tab c_sc_mgm_grp

* Paternal occupational social class - Will split into low (III manual, IV and V) and high (!, II and III non-manual))
tab c_sc_mgf, m

replace c_sc_mgf = . if c_sc_mgf < 0 | c_sc_mgf == 65
tab c_sc_mgf, m

recode c_sc_mgf (1 2 3 = 1) (4 5 6 = 0), gen(c_sc_mgf_grp)
label value c_sc_mgf_grp occ_lb
tab c_sc_mgf_grp, m
tab c_sc_mgf_grp

* Household income (weekly - log)
sum logavinceq

hist logavinceq, freq

* IMD
tab jan1993imd2010q5_M, m

replace jan1993imd2010q5_M = . if jan1993imd2010q5_M < 0
tab jan1993imd2010q5_M, m
tab jan1993imd2010q5_M

* Townsend deprivation index
tab jan1993Townsendq5_M, m

replace jan1993Townsendq5_M = . if jan1993Townsendq5_M < 0
tab jan1993Townsendq5_M, m
tab jan1993Townsendq5_M

* Housing status - Recode into owned vs rented vs council/HA vs other
tab a006, m

replace a006 = . if a006 < 0
tab a006, m

recode a006 (0 1 = 1) (3 4  = 2) (2 5 = 3) (6 = 4), gen(a006_grp)
label define housing_lb 1 "Owned/Mortgaged" 2 "Renting" 3 "Council/HA" 4 "Other"
numlabel housing_lb, add
label values a006_grp housing_lb
tab a006_grp, m
tab a006_grp

* Financial difficulties - Combine into yes vs no
tab b594, m

replace b594 = . if b594 < 0
tab b594

recode b594 (5 = 0) (1/4 = 1), gen(b594_grp)
label define fin_lb 0 "No" 1 "Yes"
numlabel fin_lb, add
label values b594_grp fin_lb
tab b594_grp, m
tab b594_grp

* Financial difficulties score
tab c525, m

replace c525 = . if c525 < 0
tab c525

* Traumatic life events in childhood (< 17 years of age) - Both weighted and unweighted scores (weighted = occurred and impact had in participant; unweighted = whether event occurred or not)
tab1 c432 c433, m

replace c432 = . if c432 < 0
replace c433 = . if c433 < 0
tab1 c432 c433, m
tab1 c432 c433

sum c432 c433, d

hist c432, freq width(1)
hist c433, freq width(1)

* Family got poorer in childhood
tab c429a, m

replace c429a = . if c429a < 0
replace c429a = 0 if c429a == 2
label define poor_lb 0 "No" 1 "Yes"
numlabel poor_lb, add
label value c429a poor_lb
tab c429a, m
tab c429a

* Access to car
tab a053, m

replace a053 = . if a053 < 0
replace a053 = 0 if a053 == 2
label define car_lb 0 "No" 1 "Yes"
numlabel car_lb, add
label value a053 car_lb
tab a053, m
tab a053

* Household crowding index
tab a551, m

replace a551 = . if a551 < 0
tab a551, m
tab a551

* Self-reported neighbourhood quality index
tab a636, m

replace a636 = . if a636 < 0
tab a636, m
tab a636

sum a636

* Partner absense in pregnancy (combine two time points - Start with 'C' questionnaire as completed later, then back-fill with 'A' questionnaire if missing)
tab1 a520 c760, m

replace a520 = . if a520 < 0
replace c760 = . if c760 < 0
tab1 a520 c760, m

gen partner_ab = .
replace partner_ab = 0 if c760 == 1
replace partner_ab = 1 if c760 == 2
replace partner_ab = 0 if partner_ab == . & (a520 == 1 | a520 == 2 | a520 == 4)
replace partner_ab = 1 if partner_ab == . & a520 == 3

label define partAb_lb 0 "No - Partner present" 1 "Yes - Partner absent"
numlabel partAb_lb, add
label values partner_ab partAb_lb
tab partner_ab, m


** Cognitive/Psychological factors

* Cognitive tests (from FOM2/3/4 clinics) - includes logic memory, digit backwards, spot-the-word, digit symbol coding, verbal fluency and logic memory delayed tests. Will use PCA to create a single 'general intelligence' factor.
sum fm2cg010-fm4cg015

foreach var of varlist fm2cg010-fm4cg015 {
	replace `var' = . if `var' < 0
}
sum fm2cg010-fm4cg015

* If participants missed FOM2, do they have FOM3 or FOM4 data? Yes, around 700 mothers have FOM3 or FOM4 data, but not FOM 2 - So will use their data if FOM2 is missing
misstable sum fm2cg010 fm3cg010 fm4cg010, all
misstable patterns fm2cg010 fm3cg010 fm4cg010, freq

* Combine this data together, using only complete cases within each FOM clinic
gen fom2_cc = 0
replace fom2_cc = 1 if fm2cg010 < . & fm2cg011 < . & fm2cg012 < . & fm2cg013 < . & fm2cg014 <. & fm2cg015 < .
tab fom2_cc

gen fom3_cc = 0
replace fom3_cc = 1 if fm3cg010 < . & fm3cg011 < . & fm3cg012 < . & fm3cg013 < . & fm3cg014 <. & fm3cg015 < .
tab fom3_cc

gen fom4_cc = 0
replace fom4_cc = 1 if fm4cg010 < . & fm4cg011 < . & fm4cg012 < . & fm4cg013 < . & fm4cg014 <. & fm4cg015 < .
tab fom4_cc

* Logic memory
gen logic_mem = fm2cg010 if fom2_cc == 1
replace logic_mem = fm3cg010 if fom3_cc == 1 & logic_mem == .
replace logic_mem = fm4cg010 if fom4_cc == 1 & logic_mem == .
sum fm2cg010 fm3cg010 fm4cg010 logic_mem

* Digit backwards
gen digit_back = fm2cg011 if fom2_cc == 1
replace digit_back = fm3cg011 if fom3_cc == 1 & digit_back == .
replace digit_back = fm4cg011 if fom4_cc == 1 & digit_back == .
sum fm2cg011 fm3cg011 fm4cg011 digit_back

* Spot the word
gen spot_word = fm2cg012 if fom2_cc == 1
replace spot_word = fm3cg012 if fom3_cc == 1 & spot_word == .
replace spot_word = fm4cg012 if fom4_cc == 1 & spot_word == .
sum fm2cg012 fm3cg012 fm4cg012 spot_word

* Digit symbol
gen digit_symbol = fm2cg013 if fom2_cc == 1
replace digit_symbol = fm3cg013 if fom3_cc == 1 & digit_symbol == .
replace digit_symbol = fm4cg013 if fom4_cc == 1 & digit_symbol == .
sum fm2cg013 fm3cg013 fm4cg013 digit_symbol

* Verbal fluency
gen verbal = fm2cg014 if fom2_cc == 1
replace verbal = fm3cg014 if fom3_cc == 1 & verbal == .
replace verbal = fm4cg014 if fom4_cc == 1 & verbal == .
sum fm2cg014 fm3cg014 fm4cg014 verbal

* Logic memory (delayed)
gen logic_mem_delay = fm2cg015 if fom2_cc == 1
replace logic_mem_delay = fm3cg015 if fom3_cc == 1 & logic_mem_delay == .
replace logic_mem_delay = fm4cg015 if fom4_cc == 1 & logic_mem_delay == .
sum fm2cg015 fm3cg015 fm4cg015 logic_mem_delay

** Run the PCA and extract the first principal component (explains ~43% of the variation in these cognitive tests - Next component explains ~20% of the variation, so first component has strongest association by quite a large amount)
pca logic_mem-logic_mem_delay
predict fom_cog_factor1, score

sum fom_cog_factor1

hist fom_cog_factor1, freq


** Personality (inter-personal sensitivity measure; IPSM) - Includes inter-personal awareness, need for approval, separation anxiety, timidity and fragile inner-self sub-scales; plus total inter-personal sensitivity score
sum b916-b921

* Drop variables indicating number of missing items (ending in 'm')
foreach var of varlist b916-b921 {
	capture drop `var'm
}

* Recode data as missing if negative
foreach var of varlist b916-b921 {
	replace `var' = . if `var' < 0
}

sum b916-b921
sum b916-b921, d

** Locus of control
tab d842, m

replace d842 = . if d842 < 0
tab d842, m
tab d842

sum d842

** Bachman self-esteem score
tab1 h151b, m

replace h151b = . if h151b < 0

sum h151b, d


****** Now on to the G0 partner/father exposures

*** Demographic variables

* Age in pregnancy
tab pb910, m

replace pb910 = . if pb910 < 0
tab pb910, m
sum pb910

* Age at @28 RSBB questionnaire
tab FC9992, m

replace FC9992 = . if FC9992 < 0
tab FC9992, m
sum FC9992

* Ethnicity - Combine into white vs other than white
tab c801, m

replace c801 = . if c801 < 0
tab c801, m

recode c801 (1 = 0) (2/9 = 1), gen(c801_grp)
label values c801_grp white_lb
tab c801_grp, m
tab c801_grp

* Also want to add in missing ethnicity data collect as part of COVID4 questionnaire
tab c801_grp covid4p_0502, m

replace c801_grp = 0 if covid4p_0502 == 1 & c801_grp == .
replace c801_grp = 1 if covid4p_0502 == 2 & c801_grp == .

tab c801_grp, m
tab c801_grp

* Marital status - Combine widowed/divorced/separated, as well as 1st and 2nd/3rd marriage
tab pa065, m

replace pa065 = . if pa065 < 0
tab pa065, m

recode pa065 (1 = 1) (5 6 = 2) (2 3 4 = 3), gen(pa065_grp)
label values pa065_grp marital_lb
tab pa065_grp, m
tab pa065_grp

* Residential mobility - Combine into 0/1/2/3/4/5 or more
tab pa005, m

replace pa005 = . if pa005 < 0
tab pa005, m

recode pa005 (5/50 = 5), gen(pa005_grp)
label values pa005_grp mobility_lb
tab pa005_grp, m
tab pa005_grp


** Socioeconomic/material insecurity factors

* Highest education
tab c666a, m

replace c666a = . if c666a < 0
tab c666a, m
tab c666a

* Highest mother's education
tab pb359a, m

replace pb359a = . if pb359a < 0
tab pb359a, m
tab pb359a

* Highest father's education
tab pb376a, m

replace pb376a = . if pb376a < 0
tab pb376a, m
tab pb376a

* Occupational social class - Will split into low (III manual, IV and V) and high (!, II and III non-manual))
tab c765, m

replace c765 = . if c765 < 0 | c765 == 65
tab c765, m

recode c765 (1 2 3 = 1) (4 5 6 = 0), gen(c765_grp)
label value c765_grp occ_lb
tab c765_grp, m
tab c765_grp

* Mother's occupational social class - Will split into low (III manual, IV and V) and high (!, II and III non-manual))
tab pb_sc_pgm, m

replace pb_sc_pgm = . if pb_sc_pgm < 0 | pb_sc_pgm == 65
tab pb_sc_pgm, m

recode pb_sc_pgm (1 2 3 = 1) (4 5 6 = 0), gen(pb_sc_pgm_grp)
label value pb_sc_pgm_grp occ_lb
tab pb_sc_pgm_grp, m
tab pb_sc_pgm_grp

* Father's occupational social class - Will split into low (III manual, IV and V) and high (!, II and III non-manual))
tab pb_sc_pgf, m

replace pb_sc_pgf = . if pb_sc_pgf < 0 | pb_sc_pgf == 65
tab pb_sc_pgf, m

recode pb_sc_pgf (1 2 3 = 1) (4 5 6 = 0), gen(pb_sc_pgf_grp)
label value pb_sc_pgf_grp occ_lb
tab pb_sc_pgf_grp, m
tab pb_sc_pgf_grp

* Financial difficulties - Combine into yes vs no
tab pb184, m

replace pb184 = . if pb184 < 0
tab pb184

recode pb184 (5 = 0) (1/4 = 1), gen(pb184_grp)
label values pb184_grp fin_lb
tab pb184_grp, m
tab pb184_grp

* Financial difficulties score
tab pd685, m

replace pd685 = . if pd685 < 0
tab pd685

* Traumatic life events in childhood (< 17 years of age) - Both weighted and unweighted scores (weighted = occurred and impact had in participant; unweighted = whether event occurred or not)
tab1 pb481 pb482, m

replace pb481 = . if pb481 < 0
replace pb482 = . if pb482 < 0
tab1 pb481 pb482, m
tab1 pb481 pb482

sum pb481 pb482, d

hist pb481, freq width(1)
hist pb482, freq width(1)

* Family got poorer in childhood
tab pb479a, m

replace pb479a = . if pb479a < 0
replace pb479a = 0 if pb479a == 2
label values pb479a poor_lb
tab pb479a, m
tab pb479a

* Self-reported neighbourhood quality index - Unlike the mother's data on this, for the partners/fathers this scale was changed mid-way through data collection. Most partners answered the same as the mothers (on a scale of 'usually', 'sometimes' or 'never'), while ~1,000 answered on a binary 'yes' vs 'no' response. As not compariable - and more have data - will use the three-item reponse
tab1 pa024-pa029a, m

recode pa024 (1 = 2) (2 = 1) (3 = 0) (-1 = .), gen(lively)
recode pa025 (1 = 2) (2 = 1) (3 = 0) (-1 = .), gen(friendly)
recode pa026 (1 = 0) (2 = 1) (3 = 2) (-1 = .), gen(noisy)
recode pa027 (1 = 2) (2 = 1) (3 = 0) (-1 = .), gen(clean)
recode pa028 (1 = 2) (2 = 1) (3 = 0) (-1 = .), gen(attractive)
recode pa029 (1 = 0) (2 = 1) (3 = 2) (-1 = .), gen(dirty)
tab1 lively-dirty, m

egen neighbour_qual = rowtotal(lively friendly noisy clean attractive dirty)
egen neighbour_miss = rowmiss(lively friendly noisy clean attractive dirty)

tab1 neighbour_qual neighbour_miss, m

replace neighbour_qual = . if neighbour_miss == 6
tab1 neighbour_qual neighbour_miss, m

drop neighbour_miss


** Cognitive/psychological factors

* Personality (inter-personal sensitivity measure; IPSM) - Includes inter-personal awareness, need for approval, separation anxiety, timidity and fragile inner-self sub-scales; plus total inter-personal sensitivity score
sum pb546-pb551

* Drop variables indicating number of missing items (ending in 'm')
foreach var of varlist pb546-pb551 {
	capture drop `var'm
}

* Recode data as missing if negative
foreach var of varlist pb546-pb551 {
	replace `var' = . if `var' < 0
}

sum pb546-pb551
sum pb546-pb551, d

* Locus of control
tab pa782, m

replace pa782 = . if pa782 < 0
tab pa782, m
tab pa782

sum pa782

* Bachman self-esteem score - Need to create this for partners/fathers
tab1 pf3000-pf3010, m

foreach var of varlist pf3000-pf3010 {
	replace `var' = . if `var' < 0
}

tab1 pf3000-pf3010, m

* Create total score and number of missing variables, then prorate, depending on number of missing questions, and include prorated score if <50% missing
recode pf3000 (5 = 0) (4 = 1) (3 = 2) (2 = 3) (1 = 4), gen(esteem1)
recode pf3001 (5 = 0) (4 = 1) (3 = 2) (2 = 3) (1 = 4), gen(esteem2)
recode pf3002 (5 = 0) (4 = 1) (3 = 2) (2 = 3) (1 = 4), gen(esteem3)
recode pf3003 (1 = 0) (2 = 1) (3 = 2) (4 = 3) (5 = 4), gen(esteem4)
recode pf3004 (5 = 0) (4 = 1) (3 = 2) (2 = 3) (1 = 4), gen(esteem5)
recode pf3005 (1 = 0) (2 = 1) (3 = 2) (4 = 3) (5 = 4), gen(esteem6)
recode pf3006 (5 = 0) (4 = 1) (3 = 2) (2 = 3) (1 = 4), gen(esteem7)
recode pf3007 (1 = 0) (2 = 1) (3 = 2) (4 = 3) (5 = 4), gen(esteem8)
recode pf3008 (5 = 0) (4 = 1) (3 = 2) (2 = 3) (1 = 4), gen(esteem9)
recode pf3009 (1 = 0) (2 = 1) (3 = 2) (4 = 3) (5 = 4), gen(esteem10)
sum esteem1-esteem10

egen esteem = rowtotal(esteem1-esteem10)
egen esteem_miss = rowmiss(esteem1-esteem10)

tab1 esteem esteem_miss, m

gen esteem_cc = esteem if esteem_miss == 0
tab1 esteem esteem_cc esteem_miss, m

gen esteem_prorated = round(esteem * (10 / (10 - esteem_miss))) if esteem_miss < 6

tab1 esteem esteem_cc esteem_prorated esteem_miss, m

sum esteem_prorated esteem_cc
sum esteem_prorated esteem_cc, d



****** Now on to the G1 study child exposures

** Demographic variables

* Age of child at @28 RSBB questionnaire
tab YPG8000, m

replace YPG8000 = . if YPG8000 < 0
tab YPG8000, m

* Sex of child
tab kz021, m

replace kz021 = . if kz021 < 0
tab kz021, m

* Ethnicity
tab c804, m

replace c804 = . if c804 < 0
tab c804, m

recode c804 (1 = 0) (2 = 1)

label values c804 white_lb
tab c804, m
tab c804

* Also want to add in missing ethnicity data collect as part of @29 questionnaire
tab c804 YPH2012, m

replace c804 = 0 if YPH2012 == 1 & c804 == .
replace c804 = 1 if YPH2012 == 2 & c804 == .

tab c804, m
tab c804

* Urban/rural status (time-point closest to RSBB questions - Jan 2020) - combine into urban vs rural
tab jan2020ur01ind, m

replace jan2020ur01ind = . if jan2020ur01ind < 0
tab jan2020ur01ind, m

recode jan2020ur01ind (1 = 1) (2 3 4 = 0), gen(jan2020ur01ind_grp)
label value jan2020ur01ind_grp urban_lb
tab jan2020ur01ind_grp, m
tab jan2020ur01ind_grp

* Is a parent (need to combine lots of questionnaires together - Will start with YPG and work backwards to age 20)
tab1 CCU1000 YPA1000 YPB7000 YPC1050 YPD2000 YPE0101 YPF6510 YPG5000, m

foreach var of varlist CCU1000 YPA1000 YPB7000 YPC1050 YPD2000 YPE0101 YPF6510 YPG5000 {
	replace `var' = . if `var' < 0
	tab `var', m
}

recode CCU1000 (1 2 = 1) (3 = 0)
tab CCU1000

recode YPA1000 (1 2 = 1) (3 = 0)
tab YPA1000

recode YPB7000 (1 = 1) (2 = 0)
tab YPB7000

gen parent = .
replace parent = YPG5000 if YPG5000 < .
replace parent = YPF6510 if YPF6510 < . & parent == .
replace parent = YPE0101 if YPE0101 < . & parent == .
replace parent = YPD2000 if YPD2000 < . & parent == .
replace parent = YPC1050 if YPC1050 < . & parent == .
replace parent = YPB7000 if YPB7000 < . & parent == .
replace parent = YPA1000 if YPA1000 < . & parent == .
replace parent = CCU1000 if CCU1000 < . & parent == .
tab parent, m

label define parent_lb 0 "No - Not a parent" 1 "Yes - Is a parent"
numlabel parent_lb, add
label value parent parent_lb
tab parent, m
tab parent

* No measure of marital status for G1, so will use 'living with a partner' at the RSBB time-point as a proxy
tab YPG1052, m

replace YPG1052 = . if YPG1052 < 0
tab YPG1052, m


** Socioeconomic/material insecurity factors

* Highest education (using as-yet-unpublished data)
tab yp_edu, m
tab yp_edu

* IMD (using jan 2020 data)
tab jan2020imd2010q5, m

replace jan2020imd2010q5 = . if jan2020imd2010q5 < 0
tab jan2020imd2010q5, m
tab jan2020imd2010q5

* Townsend deprivation index (again, using Jan 2020 data)
tab jan2020Townsendq5, m

replace jan2020Townsendq5 = . if jan2020Townsendq5 < 0
tab jan2020Townsendq5, m
tab jan2020Townsendq5

* Employed at time of RSBB questionnaire
tab1 YPG1000 YPG1001 YPG1002 YPG1009, m

foreach var of varlist YPG1000 YPG1001 YPG1002 YPG1009 {
	replace `var' = . if `var' < 0
	tab `var', m
}

gen yp_employed = .
replace yp_employed = 0 if YPG1000 < . | YPG1001 < . | YPG1002 < . | YPG1009 < .
replace yp_employed = 1 if YPG1000 == 1 | YPG1001 == 1 | YPG1002 == 1 | YPG1009 == 1
tab yp_employed, m

label define employ_lb 0 "No" 1 "Yes"
numlabel employ_lb, add
label values yp_employed employ_lb
tab yp_employed, m

* Financial difficulties - Combine into yes vs no - Use repeatedly-collected data from 5 recent questionnaires
tab1 YPB6190 YPC2360 YPD1210 YPE6710 YPF6220, m

foreach var of varlist YPB6190 YPC2360 YPD1210 YPE6710 YPF6220 {
	replace `var' = . if `var' < 0
	
	recode `var' (0 5 = 0) (1/4 = 1), gen(`var'_grp)
	tab1 `var' `var'_grp, m
}

gen yp_finDiffs = .
replace yp_finDiffs = YPF6220_grp if YPF6220_grp < .
replace yp_finDiffs = YPE6710_grp if YPE6710_grp < . & yp_finDiffs == .
replace yp_finDiffs = YPD1210_grp if YPD1210_grp < . & yp_finDiffs == .
replace yp_finDiffs = YPC2360_grp if YPC2360_grp < . & yp_finDiffs == .
replace yp_finDiffs = YPB6190_grp if YPB6190_grp < . & yp_finDiffs == .

label values yp_finDiffs fin_lb
tab yp_finDiffs, m
tab yp_finDiffs

* Traumatic life events/adverse childhood experiences in childhood (< 17 years of age; more details on how these variables were derived can be found here: https://wellcomeopenresearch.org/articles/3-106/v1)
desc clon100-clon123

sum clon100-clon123

foreach var of varlist clon100-clon123 {
	replace `var' = . if `var' < 0
}

sum clon100-clon123

tab1 clon100-clon123, m
tab1 clon100-clon123

* YP income (at age 26) - Combine some categories together
tab YPE6020, m

replace YPE6020 = . if YPE6020 < 0
tab YPE6020, m

recode YPE6020 (0 1 = 1) (2 = 2) (3 = 3) (4 = 4) (5 6 7 = 5), gen(yp_income)
tab yp_income, m

label define income_lb 1 "£0-£499" 2 "£500-£999" 3 "£1000-£1499" 4 "£1500-£1999" 5 "£2000 and above"
numlabel income_lb, add
label values yp_income income_lb
tab yp_income, m

* YP housing status (at age 28) - Combine some categories together
tab YPG1060, m

replace YPG1060 = . if YPG1060 < 0
tab YPG1060, m

recode YPG1060 (1 2 = 1) (4 = 2) (3 = 3) (5 = 4), gen(yp_housing)
label values yp_housing housing_lb
tab yp_housing, m

* Biological father absence in childhood (combine time points at ages 3, 4, 8, 10 and 18, working backwards - Will construct both binary 'father absence' variable, and 'age at father absence' categorical variable [never vs 0-4 vs 5 or older]; for the latter, will just use age 8, 10 and 18 data)
tab1 h400 h401 j374 j375 n8050 n8051 q3050 q3051 t1055 t1056, m

foreach var of varlist h400 h401 j374 j375 n8050 n8051 q3050 q3051 t1055 t1056 {
	replace `var' = . if `var' < 0
	tab `var', m
}

* Recode variables so are all consistent (0 = no/lives with bio father vs 1 = yes/not live with bio father)
recode h400 (1 = 0) (2 7 = 1) (9 = .), gen(h400_binFA)
tab1 h400 h400_binFA, m

recode j374 (1 = 1) (2 = 0), gen(j374_binFA)
tab1 j374 j374_binFA, m

recode n8050 (1 = 1) (2 = 0), gen(n8050_binFA)
tab1 n8050 n8050_binFA, m

recode q3050 (1 = 1) (2 = 0), gen(q3050_binFA)
tab1 q3050 q3050_binFA, m

recode t1055 (1 = 0) (2 = 1), gen(t1055_binFA)
tab1 t1055 t1055_binFA, m

gen father_ab = t1055_binFA if t1055_binFA < .
replace father_ab = q3050_binFA if q3050_binFA < . & father_ab == .
replace father_ab = n8050_binFA if n8050_binFA < . & father_ab == .
replace father_ab = j374_binFA if j374_binFA < . & father_ab == .
replace father_ab = h400_binFA if h400_binFA < . & father_ab == .
tab father_ab, m

label define fa_lb 0 "No father absence" 1 "Yes father absence"
numlabel fa_lb, add
label values father_ab fa_lb
tab father_ab, m
tab father_ab

* And age at FA
sum n8051 q3051 t1056

replace n8051 = n8051 / 12
replace q3051 = q3051 / 12
sum n8051 q3051 t1056

gen age_FA = .
replace age_FA = 0 if t1055_binFA == 0
replace age_FA = 1 if t1055_binFA == 1 & t1056 < 5
replace age_FA = 2 if t1055_binFA == 1 & t1056 >=5 & t1056 < . 

replace age_FA = 0 if q3050_binFA == 0 & age_FA == .
replace age_FA = 1 if q3050_binFA == 1 & q3051 < 5 & age_FA == .
replace age_FA = 2 if q3050_binFA == 1 & q3051 >=5 & q3051 < . & age_FA == .

replace age_FA = 0 if n8050_binFA == 0 & age_FA == .
replace age_FA = 1 if n8050_binFA == 1 & n8051 < 5 & age_FA == .
replace age_FA = 2 if n8050_binFA == 1 & n8051 >=5 & n8051 < . & age_FA == .

label define fa_age_lb 0 "No father absence" 1 "FA age 0 to 4" 2 "FA age 5 or older"
numlabel fa_age_lb, add
label values age_FA fa_age_lb
tab age_FA, m
tab age_FA


** Cognitive/psychological factors

* Intelligence/IQ - Based on WISC at age 8 (verbal IQ, performance IQ, and total IQ), WASI at age 15 (total IQ), and WISC digit symbol and vocabulary tasks at age 24
sum f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050

foreach var of varlist f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 {
	replace `var' = . if `var' < 0
}

sum f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050

* Personality (based on big-5) at age 15
sum fg7360-fg7364

foreach var of varlist fg7360-fg7364 {
	replace `var' = . if `var' < 0
}

sum fg7360-fg7364

* Locus of control (at ages 8 and 16; at age 16, need to code into DV)
tab f8lc125, m

replace f8lc125 = . if f8lc125 < 0
tab f8lc125, m
tab f8lc125

tab1 ccs3000-ccs3012, m

* Loop through LOC vars, recode negative values as missing, and recode into 0/1s for the summary score for external LOC (other than for ccs3011, a 'yes' response indicates eLOC)
foreach var of varlist ccs3000-ccs3011 {
	replace `var' = . if `var' < 0
	
	if "`var'" != "ccs3011" {
		recode `var' (1 = 1) (2 = 0), gen(`var'_loc)
		tab1 `var' `var'_loc, m
	}
	
	if "`var'" == "ccs3011" {
		recode `var' (1 = 0) (2 = 1), gen(`var'_loc)
		tab1 `var' `var'_loc, m
	}
	
}

egen loc_age16 = rowtotal(ccs3000_loc-ccs3011_loc), missing
egen loc16_miss = rowmiss(ccs3000_loc-ccs3011_loc)

tab1 loc_age16 loc16_miss, m

replace loc_age16 = . if loc16_miss > 0

tab1 loc_age16 loc16_miss, m

sum loc_age16
sum loc_age16, d


* Negative cognitive style at age 17 (higher values = more negative cognitive style)
sum FJCQ1001

replace FJCQ1001 = . if FJCQ1001 < 0

sum FJCQ1001
sum FJCQ1001, d


* Theory of mind (DANVA faces task) at age 8 - ToM defined using 'at least 7 faces misattributed' - Recode so no = 0
tab f8dv440a, m

replace f8dv440a = . if f8dv440a < 0

replace f8dv440a = 0 if f8dv440a == 2
label define faces_lb 0 "No - < 7 errors" 1 "Yes - 7 or more errors"
numlabel faces_lb, add
label values f8dv440a faces_lb

tab f8dv440a, m
tab f8dv440a


* Theory of mind (triangles task) at age 13 - Coding based on description in methods here: https://www.nature.com/articles/s41598-018-21737-8
desc fg4638-fg4692

* Just keep the scored questions in consecutive order - Drop the response time vars, and move the control Qs to after the scored Qs (are 16 scored questions and 12 control questions)
drop fg4639 fg4641 fg4643 fg4645 fg4647 fg4651 fg4653 fg4655 fg4657 fg4659 fg4661 fg4663 fg4665 fg4667 fg4669 fg4671 fg4673 fg4675 fg4677 fg4679 fg4681 fg4683 fg4685 fg4687 fg4689 fg4691

order aln-fg4638 fg4642 fg4648 fg4650 fg4652 fg4654 fg4660 fg4662 fg4666 fg4668 fg4672 fg4676 fg4680 fg4682 fg4688 fg4692 fg4640 fg4644 fg4646 fg4656 fg4658 fg4664 fg4670 fg4674 fg4678 fg4684 fg4686 fg4690

desc fg4638-fg4692

sum fg4638-fg4692

foreach var of varlist fg4638-fg4692 {
	replace `var' = . if `var' < 0
	tab `var', m
}

* Calculate total score, which is score for positive questions (where mental state matched that of the triangle - e.g., when the triangle was 'happy' the question asked "Was the triangle happy?") substracted by score for negative questions (where mental state did not match that of the triangle - e.g., when the triangle was 'sad' the question asked "Was the triangle scared?"). Were 8 of positve and 8 negative questions, and children answered on a 6 point likerty scale from 0 (not at all) to 5 (extremely). Will only keep data for children who answered all questions.
egen triangles_miss = rowmiss(fg4638-fg4692)
tab triangles_miss

egen triangles_pos = rowtotal(fg4642 fg4648 fg4652 fg4660 fg4668 fg4672 fg4680 fg4688), miss
replace triangles_pos = . if triangles_miss > 0
sum triangles_pos

egen triangles_neg = rowtotal(fg4638 fg4650 fg4654 fg4662 fg4666 fg4676 fg4682 fg4692), miss
replace triangles_neg = . if triangles_miss > 0
sum triangles_neg

gen triangles_total = (triangles_pos - triangles_neg) + 40
sum triangles_total

* Also remove scores where more than 50% of items were the same answer (suggesting they were not attending to the task) - This includes the control questions
foreach var of varlist fg4640-fg4690 {
	replace `var' = . if `var' < 0
	tab `var', m
}

egen all0 = anycount(fg4638-fg4690), values(0)
replace all0 = . if triangles_miss > 0
tab all0

replace triangles_total = . if all0 > 14
drop all0

egen all1 = anycount(fg4638-fg4690), values(1)
replace all1 = . if triangles_miss > 0
tab all1

replace triangles_total = . if all1 > 14
drop all1

egen all2 = anycount(fg4638-fg4690), values(2)
replace all2 = . if triangles_miss > 0
tab all2

replace triangles_total = . if all2 > 14
drop all2

egen all3 = anycount(fg4638-fg4690), values(3)
replace all3 = . if triangles_miss > 0
tab all3

replace triangles_total = . if all3 > 14
drop all3

egen all4 = anycount(fg4638-fg4690), values(4)
replace all4 = . if triangles_miss > 0
tab all4

replace triangles_total = . if all4 > 14
drop all4

egen all5 = anycount(fg4638-fg4690), values(5)
replace all5 = . if triangles_miss > 0
tab all5

replace triangles_total = . if all5 > 14
drop all5

sum triangles_total
hist triangles_total, freq width(1)


* Skuse social cognition scale - At ages 8 and 16 (need to code the age 16 data, though)
tab kr554a, m

replace kr554a = . if kr554a < 0
tab kr554a, m
tab kr554a

sum kr554a
sum kr554a, d

tab1 tc4050-tc4061, m

foreach var of varlist tc4050-tc4061 {
	replace `var' = . if `var' < 0
	recode `var' (1 = 0) (2 = 1) (3 = 2), gen(`var'_edit)
	tab1 `var' `var'_edit, m
}

egen skuse16 = rowtotal(tc4050_edit-tc4061_edit), miss
egen skuse16_miss = rowmiss(tc4050_edit-tc4061_edit)
tab1 skuse16 skuse16_miss, m

replace skuse16 = . if skuse16_miss > 0
sum skuse16


* Autism spectrum quotient scale asked at age 25
desc YPE5510-YPE5559

tab1 YPE5510-YPE5559, m

foreach var of varlist YPE5510-YPE5559 {
	replace `var' = . if `var' < 0
	
	if "`var'" == "YPE5510" | "`var'" == "YPE5512" | "`var'" == "YPE5517" | "`var'" == "YPE5519" | "`var'" == "YPE5520" | "`var'" == "YPE5523" | "`var'" == "YPE5524" | "`var'" == "YPE5526" | "`var'" == "YPE5533" | "`var'" == "YPE5534" | "`var'" == "YPE5536" | "`var'" == "YPE5537" | "`var'" == "YPE5538" | "`var'" == "YPE5539" | "`var'" == "YPE5540" | "`var'" == "YPE5541" | "`var'" == "YPE5543" | "`var'" == "YPE5545" | "`var'" == "YPE5546" | "`var'" == "YPE5547" | "`var'" == "YPE5549" | "`var'" == "YPE5553" | "`var'" == "YPE5556" | "`var'" == "YPE5557" | "`var'" == "YPE5558" | "`var'" == "YPE5559" {
		recode `var' (0 1 = 1) (2 3 = 0), gen(`var'_edit)
		tab1 `var' `var'_edit, m
	} 
	else {
		recode `var' (0 1 = 0) (2 3 = 1), gen(`var'_edit)
		tab1 `var' `var'_edit, m
	}
}

egen autism25 = rowtotal(YPE5510_edit-YPE5559_edit), miss
egen autism25_miss = rowmiss(YPE5510_edit-YPE5559_edit)
tab1 autism25 autism25_miss, m

replace autism25 = . if autism25_miss > 0
sum autism25


* Prosociality subscale of the strengths and difficulties questionnaire - Mother-assessed at ages 8 and 13, and G1-assessed at age 25/26 (need to derive the age 26 data). Am not including the hyperactivity, emotional symptoms, conduct problems and peer problems sub-scales, as these veer into mental health territory, which is not the focus of this paper.
sum kq348a tc4025e

foreach var of varlist kq348a tc4025e {
	replace `var' = . if `var' < 0
}

tab1 kq348a tc4025e, m
tab1 kq348a tc4025e

sum kq348a tc4025e

* For age 25/26 data, just pick out the prosocial items
tab1 YPE1150 YPE1153 YPE1158 YPE1166 YPE1169, m

foreach var of varlist YPE1150 YPE1153 YPE1158 YPE1166 YPE1169 {
	replace `var' = . if `var' < 0
	tab `var', m
}

egen prosocial25 = rowtotal(YPE1150 YPE1153 YPE1158 YPE1166 YPE1169), miss
egen prosocial25_miss = rowmiss(YPE1150 YPE1153 YPE1158 YPE1166 YPE1169)
tab1 prosocial25 prosocial25_miss, m

replace prosocial25 = . if prosocial25_miss > 0
tab prosocial25, m
tab prosocial25

sum prosocial25


* Bachman self-esteem at age 17/18
tab CCXD860a, m

replace CCXD860a = . if CCXD860a < 0
tab CCXD860a, m
tab CCXD860a

sum CCXD860a

* Harter self-esteem at age 8 - Both scholastic competence and global self-worth
tab1 f8se125 f8se126, m

foreach var of varlist f8se125 f8se126 {
	replace `var' = . if `var' < 0
	tab `var', m
}

sum f8se125 f8se126



***********************************************************************************
*** Next, we want to split this data into separate G0 mother, G0 partner/father and G1 study child files ready for analysis

** Will use the 'frames' command for this (new to Stata v.16) so can work with multiple datasets simultaneously

frame rename default main


** G0 mother file (Drop one twin - remove mult_mums - remove data from people who withdrew consent)

* Copy the main dataset/frame
frame copy main mother
frame change mother

* Keep just G0 mother variables we are interested in
keep aln mult_mum_Y in_core ///
d810 d813 d813_grp d816 Y3153 Y3153_cat Y3160 Y3170 Y3155 Y3155_cat ///
Y3000 Y3040 Y3040_grp Y3080_OccNever Y3080_OccYr Y3080 ///
mz028b Y9992 c800_grp a525_grp a005_grp jan1993ur01ind_grp b032_grp ///
c645a c686a c706a c755_grp c_sc_mgf_grp c_sc_mgm_grp logavinceq jan1993imd2010q5_M jan1993Townsendq5_M a006_grp b594_grp c525 c429a c432 c433 a053 a551 a636 partner_ab ///
logic_mem-logic_mem_delay fom_cog_factor1 b916-b921 d842 h151b 

* For mult_mums who have more than one pregnancy enrolled in ALSPAC, just keep one pregnancy (plus drop ALNs where only child enrolled, not the mum)
tab mult_mum_Y, m

drop if mult_mum_Y == 1 | mult_mum_Y == .
drop mult_mum_Y

* Drop one twin
duplicates report aln

duplicates drop

duplicates report aln

* Remove mothers who withdrew consent
tab d810, m

drop if d810 == .a

* Now save this file
save "G0Mother_PredictorsOfRSBB_B3911.dta", replace

* Close this frame and go back to main dataset
frame change main
frame drop mother


** G0 partner/father file (Drop one twin - remove mult_mums - remove data from people who withdrew consent)

* Copy the main dataset/frame
frame copy main partner
frame change partner

* Keep just variables we are interested in
keep aln mult_mum_Y in_core ///
pb150 pb153 pb153_grp pb155 FC3153 FC3153_cat FC3160 FC3170 FC3155 FC3155_cat ///
FC3000 FC3040 FC3040_grp FC3080_OccNever FC3080_OccYr FC3080 ///
pb910 FC9992 c801_grp pa065_grp pa005_grp jan1993ur01ind_grp b032_grp ///
c666a pb359a pb376a c765_grp pb_sc_pgm_grp pb_sc_pgf_grp logavinceq jan1993imd2010q5_M jan1993Townsendq5_M a006_grp pb184_grp pd685 pb479a pb481 pb482 a053 a551 neighbour_qual ///
pb546-pb551 pa782 esteem_prorated 

* For mult_mums who have more than one pregnancy enrolled in ALSPAC, just keep one pregnancy (for this, will assume mult mum pregnancies both have the same partner/father) - Plus drop ALNs where only child enrolled, not the mum)
tab mult_mum_Y, m

drop if mult_mum_Y == 1 | mult_mum_Y == .
drop mult_mum_Y

* Drop one twin
duplicates report aln

duplicates drop

duplicates report aln

* Remove data from partners/fathers and mothers who withdrew consent
tab1 pb150 c666a, m

drop if pb150 == .c
drop if c666a == .a

* Now save this file
save "G0Partner_PredictorsOfRSBB_B3911.dta", replace

* Close this frame and go back to main dataset
frame change main
frame drop partner


** G1 study child file (only keep if alive at 1 year of age - remove data from people who withdrew consent)

* Copy the main dataset/frame
frame copy main g1
frame change g1

* Keep just variables we are interested in
keep aln qlet kz011b ///
YPG3000 YPG3040 YPG3040_grp YPG3080_OccNever YPG3080_OccYr YPG3080 YPG3153 YPG3153_cat YPG3160 YPG3170 YPG3155 YPG3155_cat ///
YPG8000 mz028b kz021 YPG1052 a525_grp c804 a005_grp jan2020ur01ind_grp jan1993ur01ind_grp b032_grp parent ///
yp_edu c645a c666a yp_employed c755_grp c765_grp yp_income logavinceq jan2020imd2010q5 jan1993imd2010q5_M jan2020Townsendq5 jan1993Townsendq5_M yp_housing a006_grp yp_finDiffs c525 clon100-clon123 a053 a551 a636 father_ab age_FA ///
f8ws110 f8ws111 f8ws112 fh6280 FKWI1030 FKWI1050 fg7360-fg7364 f8lc125 loc_age16 FJCQ1001 f8dv440a triangles_total kr554a skuse16 autism25 kq348a tc4025e prosocial25 CCXD860a f8se125 f8se126

* Drop if not alive at 1 year of age
tab kz011b, m

drop if kz011b == 2
drop kz011b

* Remove data from G1 children and mothers who withdrew consent
tab1 YPG3000 kr554a c804, m

drop if YPG3000 == .b
drop if c804 == .a

* Now save this file
save "G1_PredictorsOfRSBB_B3911.dta", replace

* Close this frame and go back to main dataset, then clear all frames
frame change main
frame drop g1

clear frames

* Close the log file
log close
