*** Predictors of RSBB (B3911) - Data processing script
*** Created 15/11/2021 by Dan Smith
*** Stata v16.0

*** This script reads in the raw ALSPAC data, cleans the variables, and then creates the three datasets used in subsequent analyses.


**********************************************************************************
**** Set working directory, start a log file, read in dataset, and add numlabels

cd "X:\Groups\ARC\DanS\Descriptive_PredictorsOfRSBB_B3911"

capture log close
log using "Desc_RSBB_B3911_DataCleaning_log", replace text

use "Desc_RSBB_B3911.dta", clear

numlabel, add


