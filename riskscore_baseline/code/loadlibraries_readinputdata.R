# Activate renv
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))

# load required libraries, functions; Read input data
#library(tidymodels)
library(tidyverse)
library(here)
library(methods)
library(SuperLearner)
library(e1071)
library(glmnet)
library(kyotil)
library(argparse)
library(vimp)
library(nloptr)
library(RhpcBLASctl)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("omp_set_num_threads", "RhpcBLASctl")
library(mice)
library(Hmisc) # wtd.quantile, cut2
library(dplyr)
library(recipes)

if (endsWith(attr(config, "config"), "mock")) {
    if(attr(config, "config")=="moderna_mock") {
      path_to_data <- here("..", paste0("data_raw/moderna/", mapped_data))
    } else {
        # janssen pooled or regions
      path_to_data <- here("..", paste0("data_raw/janssen/", mapped_data))
    } 
} else {
    path_to_data <- mapped_data
}
print(path_to_data)
if (!file.exists(path_to_data)) stop ("while running risk score code (loadlibraries_readinputdata.R): mapped_data not available ===========================================")
inputfileName <- gsub("^.*/", "", mapped_data)
# Define code version to run
# the demo version is simpler and runs faster!
# the production version runs SL with a diverse set of learners
run_prod <- !grepl("Mock", study_name)

# get utility files
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

inputFile <- preprocess(read.csv(path_to_data), study_name)

# Indicator of membership in the cohort included in the analysis that defines the risk score in the placebo arm. It requires:
# 1. baseline SARS-CoV-2 negative, 
# 2. per-protocol, 
# 3. no evidence of SARS-CoV-2 infection or right-censoring up to time point tinterm (2 dose) or tpeak (1 dose)
# 4. lack of missing data on a certain set of baseline input variables (not enfored here because the developer of this script need not have knowledge of risk score requirements)
# no NAs allowed. 
if (study_name == "MockCOVE") {
    # special case, redefined for backward compatibility
    inputFile$Riskscorecohortflag <- with(inputFile, ifelse(Bserostatus==0 & Perprotocol==1, 1, 0))
} else if (study_name == "COVE"){
  inputFile$Riskscorecohortflag <- with(inputFile, ifelse(Perprotocol==1, 1, 0))
} else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")){
    inputFile$Riskscorecohortflag <-
      with(inputFile, ifelse(Bserostatus==0 & Perprotocol==1 & get("EarlyendpointD"%.%timepoints[1]%.%"start1")==0 & get("EventTimePrimaryD"%.%timepoints[1])>=1, 1, 0))

} else if (study_name == "PREVENT19") { # Novavax
  inputFile <- inputFile %>%
    mutate(Riskscorecohortflag = ifelse(Bserostatus==0 & Perprotocol==1, 1, 0),
           RiskscoreAUCflag = ifelse(Trt==1 & Bserostatus==0 & Perprotocol==1 & EarlyendpointD35==0 & EventTimePrimaryD35>=7, 1, 0)
           )
} else if (study_name == "AZD1222") {
    inputFile <- inputFile %>%
      mutate(Riskscorecohortflag = ifelse(Bserostatus==0 & Perprotocol==1, 1, 0),
             RiskscoreAUCflag = ifelse(Trt==1 & Bserostatus==0 & Perprotocol==1 & EarlyendpointD57==0 & EventTimePrimaryD57>=7, 1, 0))
} else if (study_name %in% c("VAT08", "VAT08m", "VAT08b")) { # Sanofi
    inputFile <- inputFile %>%
      mutate(Riskscorecohortflag = ifelse(Perprotocol==1, 1, 0),
             RiskscoreAUCflag = ifelse(Trt==1 & Perprotocol==1 & EarlyinfectionD43==0 & EventTimePrimaryD43>=7, 1, 0))
} else if (study_name %in% c("PROFISCOV")) {
    # Needs Youyi's check; currently do nothing!
} else if (study_name %in% c("COVAIL")) {
  inputFile <- inputFile %>%
    mutate(Riskscorecohortflag = ifelse(Perprotocol == 1, 1, 0),
           RiskscoreAUCflag = ifelse(Riskscorecohortflag == 1, 1, 0))
} else stop("unknown study_name 4")

assertthat::assert_that(
    all(!is.na(inputFile$Riskscorecohortflag)),
    msg = "missing Riskscorecohortflag")
