###############################################################################
###################  load_libraries.R  ##################################
###############################################################################
###############################################################################

# Activate renv
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))

# load required libraries, functions; Read input data
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
if (!file.exists(path_to_data)) stop ("make dat proc: dataset not available ===========================================")
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
if (study_name %in% c("MockCOVE", "COVE")) {
  # special case, redefined for backward compatibility
  # inputFile$Riskscorecohortflag <- with(inputFile, ifelse(Bserostatus==0 & Perprotocol==1, 1, 0))
  inputFile$Riskscorecohortflag <- with(inputFile, ifelse(Perprotocol==1, 1, 0))  # BSEROPOS CHANGE MADE
  
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
} else if (study_name == "VAT08m") { # Sanofi
  inputFile <- inputFile %>%
    mutate(Riskscorecohortflag = ifelse(Perprotocol==1, 1, 0),
           RiskscoreAUCflag = ifelse(Trt==1 & Perprotocol==1 & EarlyendpointD43==0 & EventTimePrimaryD43>=7, 1, 0))
} else if (study_name %in% c("PROFISCOV")) {
  # Needs Youyi's check; currently do nothing!
  
} else stop("unknown study_name 4")

assertthat::assert_that(
  all(!is.na(inputFile$Riskscorecohortflag)),
  msg = "missing Riskscorecohortflag")


###############################################################################
###################  get_riskscores.R  ##################################
###############################################################################
###############################################################################

if(study_name %in% c("ENSEMBLE", "MockENSEMBLE", "PREVENT19", "AZD1222", "VAT08m", "PROFISCOV")){
  inputFile <- inputFile %>%
    rename(Ptid = Subjectid)
}else if(study_name == "MockCOVE"){
  inputFile <- inputFile %>%
    rename(Ptid = X)
}

# Identify the risk demographic variable names that will be used to compute the risk score
# Identify the endpoint variable
if(study_name %in% c("COVE", "MockCOVE")){
  endpoint <- "EventIndPrimaryD57"
  risk_timepoint <- 57
  studyName_for_report <- "COVE"
  inputMod <- inputFile
  if(study_name %in% c("COVE")){
    risk_vars <- c(
      "MinorityInd", "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl",
      "Multiracial", "Other",
      "Notreported", "Unknown",
      "HighRiskInd", "Sex", "Age", "BMI"
    )
  }
  
  if(study_name %in% c("MockCOVE")){ # as MinorityInd variable is absent in mock!
    risk_vars <- c(
      "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown", 
      "Black", "Asian", "NatAmer", "PacIsl",  
      "Multiracial", "Other", 
      "Notreported", "Unknown",
      "HighRiskInd", "Sex", "Age", "BMI"
    )
  }
  original_risk_vars <- risk_vars
}


if(study_name %in% c("ENSEMBLE", "MockENSEMBLE")){
  risk_vars <- c(
    "EthnicityHispanic","EthnicityNotreported", "EthnicityUnknown",
    "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
    "URMforsubcohortsampling", "HighRiskInd", "HIVinfection", 
    "Sex", "Age", "BMI",
    "Country.X1", "Country.X2", "Country.X3", "Country.X4", "Country.X5", "Country.X6", "Country.X7", 
    "Region.X1", "Region.X2", 
    "CalDtEnrollIND.X1"
  )
  
  # Store original original risk variables as well to check in check_if_SL_needs_be_run.R!
  original_risk_vars <- c(
    "EthnicityHispanic","EthnicityNotreported", "EthnicityUnknown",
    "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
    "URMforsubcohortsampling", "HighRiskInd", "HIVinfection", 
    "Sex", "Age", "BMI",
    "Country", "Region", "CalendarDateEnrollment"
  )
  
  if(run_prod){
    risk_vars <- append(risk_vars, c("CalDtEnrollIND.X2", "CalDtEnrollIND.X3"))
  }
  
  endpoint <- "EventIndPrimaryIncludeNotMolecConfirmedD29"
  risk_timepoint <- 29
  studyName_for_report <- "ENSEMBLE"
  
  # Create binary indicator variables for Country and Region
  inputMod <- inputFile %>%
    drop_na(CalendarDateEnrollment, all_of(endpoint)) %>%
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
  inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) 
  # %>%
  #   select(-c(Country, Region, CalDtEnrollIND))
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

if(study_name == "PREVENT19"){
  inputFile <- inputFile %>%
    mutate(EventIndPrimaryD1rscore = EventIndPrimaryD1,
           EventIndPrimaryD35rauc = ifelse(RiskscoreAUCflag == 1, EventIndPrimaryD35, NA)
    )
  
  risk_vars <- c(
    "Age", "Sex", "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Notreported", "Unknown",
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "Height", "Weight", "BMI", "HighRiskInd"
  )
  original_risk_vars <- risk_vars
  endpoint <- "EventIndPrimaryD1rscore"
  #endpoint <- paste0(endpoint, "rscore")
  riskscore_timepoint <- 1
  vaccAUC_timepoint <- 35
  studyName_for_report <- "PREVENT19"
  inputMod <- inputFile %>%
    filter(Country == 0) # Analysis based off only US subjects 
}

if(study_name == "AZD1222"){
  inputFile <- inputFile %>%
    mutate(EventIndPrimaryD1rscore = EventIndPrimaryD1,
           EventIndPrimaryD57rauc = ifelse(RiskscoreAUCflag == 1, EventIndPrimaryD57, NA))
  risk_vars <- c(
    "Age", "Sex", "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Notreported", "Unknown",
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "BMI", "Country.X1", "Country.X2"
  )
  # Store original original risk variables as well to check in check_if_SL_needs_be_run.R!
  original_risk_vars <- c(
    "Age", "Sex", "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Notreported", "Unknown",
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "BMI", "HighRiskInd", "Country"
  )
  
  endpoint <- "EventIndPrimaryD1rscore"
  riskscore_timepoint <- 1
  vaccAUC_timepoint <- 57
  studyName_for_report <- "AZD1222"
  
  # Create binary indicator variables for Country and Region
  inputMod <- inputFile %>%
    #drop_na(all_of(endpoint)) %>%
    mutate(Country = as.factor(Country))
  
  rec <- recipe(~ Country, data = inputMod)
  dummies <- rec %>%
    step_dummy(Country) %>%
    prep(training = inputMod)
  inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) 
  # %>%
  #   select(-c(Country, Region, CalDtEnrollIND))
  names(inputMod)<-gsub("\\_",".",names(inputMod))
}

if(study_name == "VAT08m"){
  inputFile <- inputFile %>%
    mutate(EventIndPrimaryD1rscore = EventIndPrimaryD1,
           EventIndPrimaryD43rauc = ifelse(RiskscoreAUCflag == 1, EventIndPrimaryD43, NA),
           pooled.age.grp = ifelse(Age >= 60, 1, 0),
           # Pool countries (Japan, Kenya and Nepal) that have sparse endpoints EventIndPrimaryD43)
           Country.pooled = ifelse(Country %in% c(5, 6, 7), 567, Country))
  
  risk_vars <- c(
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
    "URMforsubcohortsampling", "HighRiskInd", "HIVinfection",
    "Sex", "Age", "pooled.age.grp", "BMI", #"BMI.group", "Height", "Weight", 
    "Country.pooled.X2", "Country.pooled.X3", "Country.pooled.X4", "Country.pooled.X8", "Country.pooled.X567", 
    #"USAInd",  
    "CalDtEnrollIND.X1", "CalDtEnrollIND.X2", "CalDtEnrollIND.X3", "CalDtEnrollIND.X4", "CalDtEnrollIND.X5"
  )
  
  # Store original original risk variables as well to check in check_if_SL_needs_be_run.R!
  original_risk_vars <- c(
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
    "URMforsubcohortsampling", "HighRiskInd", "HIVinfection",
    "Sex", "Age", "pooled.age.grp", "BMI", #"BMI.group", "Height", "Weight", 
    "Country", 
    #"USAInd", 
    "CalendarDateEnrollment"
  )
  
  endpoint <- "EventIndPrimaryD1rscore"
  riskscore_timepoint <- 1
  vaccAUC_timepoint <- 43
  studyName_for_report <- "VAT08m"
  
  # Create binary indicator variables for Country and CalendarDateEnrollment
  inputMod <- inputFile %>%
    mutate(Country.pooled = as.factor(Country.pooled),
           CalDtEnrollIND = case_when(CalendarDateEnrollment < 28 ~ 0,
                                      CalendarDateEnrollment >= 28 & CalendarDateEnrollment < 56 ~ 1,
                                      CalendarDateEnrollment >= 56 & CalendarDateEnrollment < 84 ~ 2,
                                      CalendarDateEnrollment >= 84 & CalendarDateEnrollment < 112 ~ 3,
                                      CalendarDateEnrollment >= 112 & CalendarDateEnrollment < 140 ~ 4,
                                      CalendarDateEnrollment >= 140 & CalendarDateEnrollment < 168 ~ 5),
           CalDtEnrollIND = as.factor(CalDtEnrollIND)) 
  
  rec <- recipe(~ Country.pooled + CalDtEnrollIND, data = inputMod)
  dummies <- rec %>%
    step_dummy(Country.pooled, CalDtEnrollIND) %>%
    prep(training = inputMod)
  inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) 
  names(inputMod) <- gsub("\\_", ".", names(inputMod))
}

# Check there are no NA values in Riskscorecohortflag!
if(!study_name %in% c("COVE", "PROFISCOV")){
  assertthat::assert_that(
    all(!is.na(inputMod$Riskscorecohortflag)), msg = "NA values present in Riskscorecohortflag!"
  )
  
  # Save inputFile 
  if(!dir.exists(paste0("output/", Sys.getenv("TRIAL")))){
    dir.create(paste0("output/", Sys.getenv("TRIAL")))
  }
  save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
}

# source(here("code", "check_if_SL_needs_be_run.R"))


###############################################################################
###############################################################################
################### run_cvsl_riskscore.R ###############################
###############################################################################
###############################################################################

inputMod <- inputMod %>%
  drop_na(all_of(endpoint)) 

# Create table of cases in both arms (prior to applying Riskscorecohortflag filter)
  tab <- inputMod %>%
    drop_na(Ptid, Trt, all_of(endpoint)) %>%
    mutate(Trt = ifelse(Trt == 0, "Placebo", "Vaccine"))

  table(tab$Trt, tab %>% pull(endpoint)) %>%
    write.csv(file = here("output", Sys.getenv("TRIAL"), "cases_prior_to_applying_Riskscorecohortflag.csv"))
  rm(tab)


  dat.ph1 <- inputMod %>% filter(Riskscorecohortflag == 1 & Trt == 0)

  dat.ph1 <- dat.ph1 %>%
    # Keep only variables to be included in risk score analyses
    select(Ptid, Trt, Bserostatus, all_of(endpoint), all_of(risk_vars)) %>%
    # Drop any observation with NA values in Ptid, Trt, or endpoint!
    drop_na(Ptid, Trt, all_of(endpoint))

  # Create table of cases in both arms (prior to Risk score analyses)
  tab <- inputMod %>%
    filter(Riskscorecohortflag == 1) %>%
    drop_na(Ptid, Trt, all_of(endpoint)) %>%
    mutate(Trt = ifelse(Trt == 0, "Placebo", "Vaccine"))

  table(tab$Trt, tab %>% pull(endpoint)) %>%
    write.csv(file = here("output", Sys.getenv("TRIAL"), "cases_prior_riskScoreAnalysis.csv"))
  rm(tab)

  # Derive maxVar: the maximum number of variables that will be allowed by SL screens in the models.
  np <- sum(dat.ph1 %>% select(matches(endpoint)))
  maxVar <- max(20, floor(np / 20))
  all_risk_vars <- risk_vars
  
  # Remove a variable if the number of cases in the variable = 1 subgroup is <= 3 or the number of cases in the variable = 0 subgroup is <= 3
  dat.ph1 <- drop_riskVars_with_fewer_0s_or_1s(dat = dat.ph1, 
                                               risk_vars = risk_vars,
                                               np = np)
  
  # Update risk_vars
  risk_vars <- dat.ph1 %>%
    select(-Ptid, -Trt, -Bserostatus, -all_of(endpoint)) %>%  # BSEROPOS CHANGE MADE
    colnames()
  
  # Remove any risk_vars with more than 5% missing values. Impute the missing
  # values for other risk variables using mice package!
  dat.ph1 <- drop_riskVars_with_high_total_missing_values(dat.ph1, risk_vars)
  
  # Update risk_vars
  risk_vars <- dat.ph1 %>%
    select(-Ptid, -Trt, -Bserostatus, -all_of(endpoint)) %>%   # BSEROPOS CHANGE MADE
    colnames()
  
  X_covars2adjust <- dat.ph1 %>%
    select(all_of(risk_vars))
  
  # Save ptids to merge with predictions later
  risk_placebo_ptids <- dat.ph1 %>% filter(Bserostatus == 0) %>% select(Ptid, all_of(endpoint))  # BSEROPOS CHANGE MADE
  
  # Impute missing values in any variable included in risk_vars using the mice package!
  print("Make sure data is clean before conducting imputations!")
  X_covars2adjust <- impute_missing_values(X_covars2adjust, risk_vars)
  
  # # Check for missing values before and after imputation
  # sapply(X_covars2adjust, function(x) sum(is.na(x)))
  
  # Scale X_covars2adjust to have mean 0, sd 1 for all vars
  for (a in colnames(X_covars2adjust)) {
    X_covars2adjust[[a]] <- scale(X_covars2adjust[[a]],
                                  center = mean(X_covars2adjust[[a]], na.rm = T),
                                  scale = sd(X_covars2adjust[[a]], na.rm = T)
    )
  }
  
  X_riskVars <- X_covars2adjust
  
  # BSEROPOS CHANGE MADE
  # Drop Bserostatus == 1 subjects here
  plac_bseropos <- bind_cols(X_riskVars, dat.ph1 %>% select(Ptid, Bserostatus, Trt, all_of(endpoint))) %>% filter(Bserostatus == 1)  # BSEROPOS CHANGE MADE
  plac_bseroneg <- bind_cols(X_riskVars, dat.ph1 %>% select(Ptid, Bserostatus, Trt, all_of(endpoint))) %>% filter(Bserostatus == 0) # BSEROPOS CHANGE MADE
  X_riskVars <- plac_bseroneg %>% select(-c(Ptid, Bserostatus, Trt, all_of(endpoint)))   # BSEROPOS CHANGE MADE
  
  Y <- dat.ph1 %>% filter(Bserostatus == 0) %>% pull(endpoint)  # BSEROPOS CHANGE MADE
  
  # set up outer folds for cv variable importance; do stratified sampling
  V_outer <- 5
  cvControlVar = list(V = V_outer, stratifyCV = TRUE)
  cvControl_quote = quote(list(V = V_outer, stratifyCV = TRUE))
  
  if(study_name != "PREVENT19"){
    #if (np < round(length(Y)*2/100)) {  # Update rule: Do Leave-One-Out CV if number of cases in placebo are less than 2%!
    #V_inner <- 20
    if (np <= 30) {
      V_inner <- length(Y) - 1
    } else {
      V_inner <- 5
    }

    if(V_inner == length(Y) - 1){
      V_inner_quote <- paste0("length(Y) - 1 = ", length(Y) - 1)
    }
    innerCvControlVar = list(list(V = V_inner))
    innerCvControl_quote = quote(list(list(V = V_inner)))
    #}
  }

  if(study_name == "PREVENT19"){
    #if (np < round(length(Y)*2/100)) {  # Update rule: Do Leave-One-Out CV if number of cases in placebo are less than 2%!
      V_inner <- 5
      #V_inner <- length(Y) - 1
      if(V_inner == length(Y) - 1){
        V_inner_quote <- paste0("length(Y) - 1 = ", length(Y) - 1)
        innerCvControlVar = list(list(V = V_inner))
        innerCvControl_quote = quote(list(list(V = V_inner)))
      }else{
        innerCvControlVar = list(list(V = V_inner, stratifyCV = TRUE))
        innerCvControl_quote = quote(list(list(V = V_inner, stratifyCV = TRUE)))
      }
    #}
  }
  
  ## solve cores issue
  #blas_get_num_procs()
  blas_set_num_threads(1)
  #print(blas_get_num_procs())
  stopifnot(blas_get_num_procs() == 1)
  
  # CV.SL inputs
  familyVar = "binomial"
  methodVar = "method.CC_nloglik"
  scaleVar = "identity"
  cvsl_args <- data.frame(matrix(ncol = 2, nrow = 8)) %>%
    rename(Argument = X1,
           Value = X2) %>%
    mutate(Argument = as.character(c("Cases/Total Subjects in placebo group (%)", "family",
                        "method", "scale", "V_outer", "cvControl (outer CV control)",
                        "V_inner", "innerCvControl")),
           Value = as.character(c(paste0(np, "/", length(Y), " (", round(np*100/length(Y), 2), "%)"), familyVar,
                     methodVar, scaleVar, V_outer, cvControl_quote,
                     V_inner, innerCvControl_quote)))

  if(V_inner == length(Y) - 1){
    cvsl_args = cvsl_args %>% mutate(Value = ifelse(Argument == "V_inner", V_inner_quote, Value))
  }

  cvsl_args %>% write.csv(paste0("output/", Sys.getenv("TRIAL"), "/", "cvsl_args.csv"))

  # run super learner ensemble
  fits <- run_cv_sl_once(
    seed = 20210216,
    Y = Y,
    X_mat = X_riskVars,
    family = familyVar,
    method = methodVar,
    scale = scaleVar,
    sl_lib = SL_library,
    cvControl = cvControlVar,
    innerCvControl = innerCvControlVar,
    vimp = FALSE
  )
  
  cvaucs <- list()
  cvaucs[[1]] <- fits$cvaucs
  cvfits <- list()
  cvfits[[1]] <- fits$cvfits
  
  saveRDS(cvaucs, here("output", Sys.getenv("TRIAL"), "cvsl_riskscore_cvaucs.rds"))
  save(cvfits, file = here("output", Sys.getenv("TRIAL"), "cvsl_riskscore_cvfits.rda"))
  save(risk_placebo_ptids, file = here("output", Sys.getenv("TRIAL"), "risk_placebo_ptids.rda"))
  

  if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
    save(run_prod, Y, X_riskVars, weights, inputMod, risk_vars, all_risk_vars, endpoint, maxVar,
         V_outer, V_inner, familyVar, methodVar, scaleVar, studyName_for_report, 
         riskscore_timepoint, vaccAUC_timepoint,
         cvControlVar, inputfileName, mapped_data,
         file = here("output", Sys.getenv("TRIAL"), "objects_for_running_SL.rda"))
  } else {
    save(run_prod, Y, X_riskVars, weights, inputMod, risk_vars, all_risk_vars, endpoint, maxVar,
         V_outer, V_inner, familyVar, methodVar, scaleVar, studyName_for_report, 
         risk_timepoint, 
         cvControlVar, inputfileName, mapped_data, plac_bseropos, plac_bseroneg,
         file = here("output", Sys.getenv("TRIAL"), "objects_for_running_SL.rda"))
  }

  
  #############################################################################
  ########################createRDAfiles_fromSLobjects.R #########################
  ##############################################################################
  ##############################################################################
  # Sys.setenv(TRIAL = "janssen_pooled_realbAb")
  renv::activate(here::here(".."))
  # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
  if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
  source(here::here("..", "_common.R"))
  #-----------------------------------------------
  library(here)
  library(readr)
  require(tidyverse)
  
  # Create fancy/tidy screen names for use in tables and figures
  # @param avgs dataframe containing Screen, Learner, AUCs information as columns
  # @return object containing tidy screen names
  get_fancy_screen_names <- function(avgs) {
    return(avgs %>%
             mutate(fancyScreen = case_when(
               Screen == "screen_highcor_random" ~ "highcor_random",
               Screen == "screen_glmnet" ~ "glmnet",
               Screen == "screen_univariate_logistic_pval" ~ "univar_logistic_pval",
               Screen == "screen_all" ~ "all",
               Screen == "All" ~ "-",
               TRUE ~ as.character(Screen)
             )))
  }
  
  
  # Drop seeds/fits that returned any error
  # @param dat object containing all 10 fits (as lists) from the CV.Superlearner with folds and auc information
  # @return object upon dropping any fit that returned an error
  drop_seeds_with_error <- function(dat) {
    newdat <- vector(mode = "list", length = 1)
    j <- 1
    for (i in seq_along(dat)) {
      if (typeof(dat[[i]][1]) == "list") {
        newdat[[j]] <- dat[[i]]
        j <- j + 1
      }
    }
    newdat
  }
  
  
  # Convert SL object to SL results dataframe
  # @param dat object containing all 10 fits (as lists) from the CV.Superlearner with folds and auc information
  # @return dataframe containing CV-AUCs
  convert_SLobject_to_Slresult_dataframe <- function(dat) {
    # Remove any iteration seeds that returned an error (if any)!
    newdat <- drop_seeds_with_error(dat)
    if (!is.null(newdat[[1]])) {
      as_tibble(do.call(rbind.data.frame, lapply(newdat, function(x) x$aucs))) %>%
        group_by(Learner, Screen) %>%
        summarise(AUC = mean(AUC), se = sqrt(mean(se ^ 2)),
                  .groups = "drop") %>%
        arrange(-AUC) %>%
        mutate(
          ci_ll = AUC - qnorm(0.975) * se, ci_ul = AUC + qnorm(0.975) * se,
          AUCstr = paste0(format(round(AUC, 3), nsmall = 3), " [",
                          format(round(ci_ll, 3), nsmall = 3), ", ",
                          format(round(ci_ul, 3), nsmall = 3), "]"),
          Learner = as.character(Learner),
          Screen = as.character(Screen),
          LearnerScreen = paste(Learner, Screen)
        ) %>%
        get_fancy_screen_names() %>%
        rename(
          Screen_fromRun = Screen,
          Screen = fancyScreen
        )
    }
  }
  
  
  
  # Read in SL objects from folder, get AUCs in a dataframe
  # @param data_file RDS file containing all 10 fits from the CV.Superlearner with folds and auc information
  # @param trt string containing treatment arm (placebo or vaccine)
  # @return dataframe containing CV-AUCs
  readin_SLobjects_fromFolder <- function(data_file, trt) {
    readRDS(data_file) %>%
      convert_SLobject_to_Slresult_dataframe() %>%
      mutate(trt = trt)
  }
  
  
  # Read CV.SL object and save AUCs
  data_file <- here("output", Sys.getenv("TRIAL"), "cvsl_riskscore_cvaucs.rds")
  risk_placebo_cvaucs <- readin_SLobjects_fromFolder(data_file, trt = "placebo")
  save(risk_placebo_cvaucs, file = here("output", Sys.getenv("TRIAL"), "cvsl_risk_placebo_cvaucs.rda"))
  
  #############################################################################
  ##################### tables_figures.R  ##############################
  ##############################################################################
  ##############################################################################
  # Sys.setenv(TRIAL = "janssen_pooled_realbAb")
  # Sys.setenv(TRIAL = "prevent19")
  renv::activate(here::here(".."))
  # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
  if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
  source(here::here("..", "_common.R"))
  #-----------------------------------------------
  ## load libraries and source files #############################################
  library(cvAUC)
  library(conflicted)
  library(tidyverse)
  library(vimp)
  library(kyotil)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library(here)
  conflict_prefer("filter", "dplyr")
  conflict_prefer("summarise", "dplyr")
  conflict_prefer("load", "base")
  source(here("code", "study_specific_functions.R"))
  source(here("code", "utils.R"))
  method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
  ggplot2::theme_set(theme_cowplot())
  
  load(file = here("output", Sys.getenv("TRIAL"), "objects_for_running_SL.rda"))
  rm(Y, X_riskVars, weights, maxVar)
  load(file = here("output", Sys.getenv("TRIAL"), "cvsl_risk_placebo_cvaucs.rda"))
  
  ######## Table of demographic variables used to derive the risk score ##########
  dat <- inputMod %>%
    filter(Riskscorecohortflag == 1 & Trt == 0) %>%
    select(all_of(risk_vars)) 
  
  dat %>%
    map(~ sum(is.na(.))) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Variable Name") %>%
    mutate(V1 = paste0(V1, "/", nrow(dat), " (", format(round((V1 / nrow(dat)) * 100, 1), nsmall = 1), "%)")) %>%
    get_defs_comments_riskVars() %>%
    rename(`Total missing values` = V1) %>%
    select(`Variable Name`, Definition, `Total missing values`, Comments) %>%
    write.csv(here("output", Sys.getenv("TRIAL"), "risk_vars.csv"))
  
  ######## learner-screens #######################################################
  caption <- "All learner-screen combinations (28 in total) used as input to the superlearner."
  
  if (run_prod) {
    tab <- risk_placebo_cvaucs %>%
      filter(!Learner %in% c("SL", "Discrete SL")) %>%
      select(Learner, Screen) %>%
      mutate(
        Screen = fct_relevel(Screen, c(
          "all", "glmnet", "univar_logistic_pval",
          "highcor_random"
        )),
        Learner = as.factor(Learner),
        Learner = fct_relevel(Learner, c(
          "SL.mean", "SL.glm", "SL.glm.interaction",
          "SL.glmnet", "SL.gam", 
          "SL.xgboost", "SL.ranger.imp"
        ))
      ) %>%
      arrange(Learner, Screen) %>%
      distinct(Learner, Screen) %>%
      rename("Screen*" = Screen)
  } else {
    tab <- risk_placebo_cvaucs %>%
      filter(!Learner %in% c("SL", "Discrete SL")) %>%
      select(Learner, Screen) %>%
      mutate(
        Screen = fct_relevel(Screen, c(
          "all", "glmnet", "univar_logistic_pval",
          "highcor_random"
        )),
        Learner = as.factor(Learner),
        Learner = fct_relevel(Learner, c("SL.mean", "SL.glm"))
      ) %>%
      arrange(Learner, Screen) %>%
      distinct(Learner, Screen) %>%
      rename("Screen*" = Screen)
  }
  
  tab %>% write.csv(here("output", Sys.getenv("TRIAL"), "learner-screens.csv"))
  
  ######## SLperformance-plac ####################################################
  if(study_name=="COVE" | study_name=="MockCOVE"){
    caption <- "Performance of Superlearner and all learner-screen combinations (CV-AUCs with 95\\% CIs) for risk score analyses using placebo group and EventIndPrimaryD57 as outcome. Constraint of np/20 is applied to all learners such that no more than 6 input variables were allowed in any model."
  }else if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE"){
    caption <- "Performance of Superlearner and all learner-screen combinations (CV-AUCs with 95\\% CIs) for risk score analyses using placebo group and EventIndPrimaryD29 (including those cases that are not molecularly confirmed) as outcome. Constraint of np/20 is applied to all learners such that no more than 5 input variables were allowed in any model."
  }
  
  sl.perf <- risk_placebo_cvaucs %>%
    mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>%
    select(Learner, Screen, AUCstr)
  
  sl.perf %>% write.csv(here("output", Sys.getenv("TRIAL"), "SLperformance-plac.csv"))
  
  ################################################################################
  # Forest plots for risk_placebo model, yd57 endpoint
  options(bitmapType = "cairo")
  if (run_prod) {
    png(file = here("output", Sys.getenv("TRIAL"), "risk_placebo_cvaucs.png"),
        width = 2000, height = 1100)
    top_learner <- make_forest_plot_prod(risk_placebo_cvaucs)
  } else {
    png(file = here("output", Sys.getenv("TRIAL"), "risk_placebo_cvaucs.png"),
        width = 2000, height = 700)
    top_learner <- make_forest_plot_demo(risk_placebo_cvaucs)
  }
  grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot,
               ncol = 2)
  dev.off()
  
  ################################################################################
  # plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners
  top2_plac <- bind_rows(
    risk_placebo_cvaucs %>% arrange(-AUC) %>%
      filter(!Learner %in% c("SL", "Discrete SL")) %>%
      dplyr::slice(1:2),
    risk_placebo_cvaucs %>%
      filter(Learner == "SL"),
    risk_placebo_cvaucs %>%
      filter(Learner == "Discrete SL")
  ) %>%
    mutate(LearnerScreen = ifelse(Learner == "SL", "Super Learner",
                                  ifelse(Learner == "Discrete SL", Learner,
                                         paste0(Learner, "_", Screen_fromRun))))
  
  # Get cvsl fit and extract cv predictions
  load(file = here("output", Sys.getenv("TRIAL"), "cvsl_riskscore_cvfits.rda"))
  pred <- get_cv_predictions(cv_fit = cvfits[[1]], cvaucDAT = top2_plac)
  
  # plot ROC curve
  options(bitmapType = "cairo")
  png(file = here("output", Sys.getenv("TRIAL"), "ROCcurve_riskscore_plac.png"),
      width = 1000, height = 1000)
  p1 <- plot_roc_curves(pred, cvaucDAT = top2_plac)
  print(p1)
  dev.off()
  
  # plot pred prob plot
  options(bitmapType = "cairo")
  png(file = here("output", Sys.getenv("TRIAL"), "predProb_riskscore_plac.png"),
      width = 1100, height = 1400)
  if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
    p2 <- plot_predicted_probabilities(pred, 1)
  } else {
    p2 <- plot_predicted_probabilities(pred, risk_timepoint)
  }
  print(p2)
  dev.off()
  
  # Use SuperLearner to generate risk scores!
  load(file = here("output", Sys.getenv("TRIAL"), "risk_placebo_ptids.rda"))
  plac <- bind_cols(
    risk_placebo_ptids,
    pred %>% filter(Learner == "SL") %>% select(pred, AUCchar)
  ) %>%
    mutate(risk_score = log(pred / (1 - pred)),
           standardized_risk_score = scale(risk_score,
                                           center = mean(risk_score, na.rm = T),
                                           scale = sd(risk_score, na.rm = T)))
  
  write.csv(plac, here("output", Sys.getenv("TRIAL"), "placebo_ptids_with_riskscores.csv"),
            row.names = FALSE)
  
  save(top2_plac, file = here("output", Sys.getenv("TRIAL"), "plac_top2learners_SL_discreteSL.rda"))
  
  rm(cvfits, pred, p1, p2)
  
  #############################################################################
  #############################################################################
  ##############################################################################
  ##############################################################################
  
  
  
  #############################################################################
  #############################################################################
  ##############################################################################
  ##############################################################################
  
  
  
  #############################################################################
  #############################################################################
  ##############################################################################
  ##############################################################################
  
  
  
  #############################################################################
  #############################################################################
  ##############################################################################
  ##############################################################################