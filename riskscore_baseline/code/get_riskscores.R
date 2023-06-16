# Sys.setenv(TRIAL = "moderna_mock")
# Sys.setenv(TRIAL = "moderna_real")
# Sys.setenv(TRIAL = "janssen_pooled_mock")
# Sys.setenv(TRIAL = "azd1222") # Astra-Zeneca
# Sys.setenv(TRIAL = "prevent19") # Novavax
# Sys.setenv(TRIAL = "vat08m") # Sanofi
# Sys.setenv(TRIAL = "vat08b") # Sanofi
# Sys.setenv(TRIAL = "janssen_pooled_partA") 
# Sys.setenv(TRIAL = "butantan")

print("GET_RISKSCORES.R")

# # Since risk scores are generated for VAT08m and VAT08b using the combined data, 
# # if risk score code is called with study_name = VAT08b, change the study_name to VAT08m
# riskscore_called_using_STUDYNAME <- config$study_name
# if(config$study_name == "VAT08b"){
#   study_name = "VAT08m"
# }

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print("AAA")
source("code/loadlibraries_readinputdata.R")
print("BBB")

if(study_name %in% c("ENSEMBLE", "MockENSEMBLE", "PREVENT19", "AZD1222", "VAT08m", "VAT08b", "PROFISCOV")){
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

if(study_name %in% c("VAT08m", "VAT08b")){
  inputFile <- inputFile %>%
    mutate(EventIndPrimaryD1rscore = EventIndPrimaryD1,
           EventIndPrimaryD43rauc = ifelse(RiskscoreAUCflag == 1, EventIndPrimaryD43, NA),
           pooled.age.grp = ifelse(Age >= 60, 1, 0),
           # # Pool countries (Japan, Kenya and Nepal) that have sparse endpoints EventIndPrimaryD43)
           # Country.pooled = ifelse(Country %in% c(5, 6, 7), 567, Country)
           
           # Pool countries (Japan, Kenya and Nepal) that have sparse endpoints EventIndPrimaryD43)
           # Assign geographic region: Honduras, not Honduras for the Stage 1 trial; India, Mexico, Other/Else for the Stage 2 trial)
           Country.ind = case_when(vacc_trial == "mono" & Country == 3 ~ "Honduras",
                                      vacc_trial == "mono" & Country != 3 ~ "NotHonduras",
                                      vacc_trial == "bi" & Country == 3 ~ "India",
                                      vacc_trial == "bi" & Country == 5 ~ "Mexico",
                                      TRUE ~ "Other")
           )
  
  risk_vars <- c(
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
    "URMforsubcohortsampling", "HighRiskInd", "HIVinfection",
    "Sex", "Age", "pooled.age.grp", "BMI", #"BMI.group", "Height", "Weight", 
    "Country.ind.India", "Country.ind.Mexico", "Country.ind.NotHonduras", "Country.ind.Other", 
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
  studyName_for_report <- "VAT08m_VAT08b_combined"

  # Create binary indicator variables for Country and CalendarDateEnrollment
  inputMod <- inputFile %>%
    mutate(Country.ind = as.factor(Country.ind),
           CalDtEnrollIND = case_when(CalendarDateEnrollment < 28 ~ 0,
                                      CalendarDateEnrollment >= 28 & CalendarDateEnrollment < 56 ~ 1,
                                      CalendarDateEnrollment >= 56 & CalendarDateEnrollment < 84 ~ 2,
                                      CalendarDateEnrollment >= 84 & CalendarDateEnrollment < 112 ~ 3,
                                      CalendarDateEnrollment >= 112 & CalendarDateEnrollment < 140 ~ 4,
                                      CalendarDateEnrollment >= 140 & CalendarDateEnrollment < 168 ~ 5,   
           CalendarDateEnrollment >= 168 & CalendarDateEnrollment < 196 ~ 6),
           CalDtEnrollIND = as.factor(CalDtEnrollIND)) 
  
  rec <- recipe(~ Country.ind + CalDtEnrollIND, data = inputMod)
  dummies <- rec %>%
    step_dummy(Country.ind, CalDtEnrollIND) %>%
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
  if(study_name %in% c("VAT08m", "VAT08b")){
    if(!dir.exists(paste0("output/", Sys.getenv("TRIAL")))){
      dir.create(paste0("output/", Sys.getenv("TRIAL")))
      dir.create(paste0("output/", Sys.getenv("TRIAL"), "/", args[1]))
    }
    if(!dir.exists(paste0("output/", Sys.getenv("TRIAL"), "/", args[1]))){
      dir.create(paste0("output/", Sys.getenv("TRIAL"), "/", args[1]))
    }
    save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
    if(args[1] == "bseroneg"){
      inputFile %>% filter(Bserostatus == 0) %>%
        save(file = paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/inputFile.RData"))
    }else if(args[1] == "bseropos"){
      inputFile %>% filter(Bserostatus == 1) %>%
        save(file = paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/inputFile.RData"))
    }
  }else{
    if(!dir.exists(paste0("output/", Sys.getenv("TRIAL")))){
      dir.create(paste0("output/", Sys.getenv("TRIAL")))
    }
    save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
  }
}

if(study_name == "COVE"){
  save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
}
source(here("code", "check_if_SL_needs_be_run.R"))
