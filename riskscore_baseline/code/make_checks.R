# Sys.setenv(TRIAL = "janssen_pooled_real")  
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries and functions
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
library(tidymodels)

# Define code version to run
# the demo version is simpler and runs faster!
# the production version runs SL with a diverse set of learners
run_prod <- !grepl("Mock", study_name)

# get utility files
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

############ SETUP INPUT #######################
# Read in data file
inputFile <- read.csv(here::here("..", "data_clean", paste0(attr(config, "config"), "_data_processed.csv"))) 

# Identify the risk demographic variable names that will be used to compute the risk score
# Identify the endpoint variable
if(study_name_code == "COVE"){
  risk_vars <- c(
    "MinorityInd", "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown", 
    "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Other", 
    "Notreported", "Unknown",
    "HighRiskInd", "Sex", "Age", "BMI"
  )
  
  endpoint <- "EventIndPrimaryD57"
  studyName_for_report <- "COVE"
}

if(study_name_code == "ENSEMBLE"){
  input_risk_vars <- c("EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
                       "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
                       "URMforsubcohortsampling", "HighRiskInd", "HIVinfection", 
                       "Sex", "Age", "BMI", 
                       "Country", "Region", "CalendarDateEnrollment")
  
  risk_vars <- c(
    "EthnicityHispanic","EthnicityNotreported", "EthnicityUnknown",
    "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
    "URMforsubcohortsampling", "HighRiskInd", "HIVinfection", 
    "Sex", "Age", "BMI",
    "Country.X1", "Country.X2", "Country.X3", "Country.X4", "Country.X5", "Country.X6", "Country.X7", 
    "Region.X1", "Region.X2", 
    "CalDtEnrollIND.X1"
  )
  
  if(run_prod){
    risk_vars <- append(risk_vars, c("CalDtEnrollIND.X2", "CalDtEnrollIND.X3"))
  }
  
  endpoint <- "EventIndPrimaryIncludeNotMolecConfirmedD29"
  studyName_for_report <- "ENSEMBLE"
  
  # Create binary indicator variables for Country and Region
  inputMod <- inputFile %>%
    drop_na(CalendarDateEnrollment, EventIndPrimaryIncludeNotMolecConfirmedD29) %>%
    mutate(Sex.rand = sample(0:1, n(), replace = TRUE),
           Sex = ifelse(Sex %in% c(2, 3), Sex.rand, Sex), # assign Sex randomly as 0 or 1 if Sex is 2 or 3.
           Country = as.factor(Country),
           Region = as.factor(Region),
           CalDtEnrollIND = case_when(CalendarDateEnrollment < 28 ~ 0,
                                      CalendarDateEnrollment >= 28 & CalendarDateEnrollment < 56 ~ 1,
                                      CalendarDateEnrollment >= 56 & CalendarDateEnrollment < 84 ~ 2,
                                      CalendarDateEnrollment >= 84 & CalendarDateEnrollment < 112 ~ 3,
                                      CalendarDateEnrollment >= 112 & CalendarDateEnrollment < 140 ~ 4,
                                      CalendarDateEnrollment >= 140 & CalendarDateEnrollment < 168 ~ 5),
           CalDtEnrollIND = as.factor(CalDtEnrollIND)) %>%
    select(-Sex.rand)
  
  rec <- recipe(~ Country + Region + CalDtEnrollIND, data = inputMod)
  dummies <- rec %>%
    step_dummy(Country, Region, CalDtEnrollIND) %>%
    prep(training = inputMod)
  inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) %>%
    select(-c(Country, Region, CalDtEnrollIND))
  names(inputMod)<-gsub("\\_",".",names(inputMod))
  
  # # Create interaction variables between Region and CalDtEnrollIND
  # rec <- recipe(EventIndPrimaryIncludeNotMolecConfirmedD29 ~., data = inputMod)
  # int_mod_1 <- rec %>%
  #   step_interact(terms = ~ starts_with("Region"):starts_with("CalDtEnrollIND"))
  # int_mod_1 <- prep(int_mod_1, training = inputMod)
  # inputMod <- bake(int_mod_1, inputMod)
  # names(inputMod)<-gsub("\\_",".",names(inputMod))
  # if(run_prod){
  #   risk_vars <- append(risk_vars, c("Region.X1.x.CalDtEnrollIND.X1", "Region.X1.x.CalDtEnrollIND.X2",
  #                                    "Region.X1.x.CalDtEnrollIND.X3",
  #                                    "Region.X2.x.CalDtEnrollIND.X1", "Region.X2.x.CalDtEnrollIND.X2",
  #                                    "Region.X2.x.CalDtEnrollIND.X3"))
  # }
}


# Check there are no NA values in Riskscorecohortflag!
assertthat::assert_that(
  all(!is.na(inputMod$Riskscorecohortflag)), msg = "NA values present in Riskscorecohortflag!"
)

args <- commandArgs(trailingOnly = TRUE)

if(as.character(args[1]) == "check"){
  source(here("code", "check_if_SL_needs_be_run.R"))
}


if(as.character(args[1]) == "runCVSL"){
  source(here("code", "run_cvsl_riskscore.R"))
}





