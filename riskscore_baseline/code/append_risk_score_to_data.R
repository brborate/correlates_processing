# Sys.setenv(TRIAL = "janssen_pooled_realbAb")
# Sys.setenv(TRIAL = "prevent19")
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries, cleaned data, and risk score estimates
library(here)
library(tidyverse)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")

print("APPEND_RISK_SCORE_TO_DATA.R")

if(study_name %in% c("VAT08m", "VAT08b")){
  load(paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/objects_for_running_SL.rda"))
  load(paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/inputFile.RData"))
  placebos_risk <- read.csv(here("output", Sys.getenv("TRIAL"), args[1], "placebo_ptids_with_riskscores.csv"))
  vaccinees_risk <- read.csv(here("output", Sys.getenv("TRIAL"), args[1], "vaccine_ptids_with_riskscores.csv"))
}else{
  load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))
  load(paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
  placebos_risk <- read.csv(here("output", Sys.getenv("TRIAL"), "placebo_ptids_with_riskscores.csv"))
  vaccinees_risk <- read.csv(here("output", Sys.getenv("TRIAL"), "vaccine_ptids_with_riskscores.csv"))
}

# merge risk score with cleaned data by IDs, then save updated data file
if(study_name == "COVE"){
  plac_bseropos_risk <- read.csv(here("output", Sys.getenv("TRIAL"), "plac_bseropos_ptids_with_riskscores.csv"))
  
  risk_scores_plac <- bind_rows(placebos_risk, plac_bseropos_risk) %>%
    select(Ptid, risk_score) %>%
    mutate(standardized_risk_score = scale(risk_score,
                                           center = mean(risk_score, na.rm = T),
                                           scale = sd(risk_score, na.rm = T)))
  
  risk_scores_vacc <- bind_rows(vaccinees_risk) %>%
    select(Ptid, risk_score) %>%
    mutate(standardized_risk_score = scale(risk_score,
                                           center = mean(risk_score, na.rm = T),
                                           scale = sd(risk_score, na.rm = T)))
  
  risk_scores <- bind_rows(risk_scores_plac, risk_scores_vacc) 
  
  inputFile <- inputFile %>%
    rename(risk_score_old = risk_score,
           standardized_risk_score_old = standardized_risk_score)
  
  inputFile_with_riskscore <- full_join(inputFile, risk_scores, by = "Ptid") 
  
} else {
  risk_scores <- bind_rows(placebos_risk, vaccinees_risk) %>%
    select(Ptid, risk_score, standardized_risk_score) 
  
  inputFile_with_riskscore <- left_join(inputFile, risk_scores, by = "Ptid") 
}




# For some studies, impute the missing risk scores and standardize them separately for placebo and vaccine groups 
if(study_name == "PREVENT19"){
  if(any(is.na(inputFile_with_riskscore %>% filter(Country == 0 & Riskscorecohortflag == 1) %>% .$risk_score))){
    
    get_imputed_riskscore <- function(dat_with_riskscore, riskVars){
      imp_vars <- dat_with_riskscore %>% 
        select(Ptid, all_of(risk_vars), risk_score) %>%
        tibble::column_to_rownames("Ptid") %>%
        mice(m=5, maxit = 0, method = 'pmm', seed = 500)
      
      complete_vars <- complete(imp_vars, action = 1L) %>%
        rename(risk_score = risk_score) %>%
        mutate(standardized_risk_score = scale(risk_score,
                                               center = mean(risk_score, na.rm = T),
                                               scale = sd(risk_score, na.rm = T))) %>%
        tibble::rownames_to_column("Ptid")
      
      complete_vars %>% select(Ptid, risk_score, standardized_risk_score)
    }
    
    imputed_riskscores <- bind_rows(get_imputed_riskscore(inputFile_with_riskscore %>% filter(Trt == "0"), riskVars = risk_vars), 
                                    get_imputed_riskscore(inputFile_with_riskscore %>% filter(Trt == "1"), riskVars = risk_vars))
    
    inputFile_with_riskscore <- inputFile_with_riskscore %>% 
      select(-c(risk_score, standardized_risk_score)) %>%
      left_join(imputed_riskscores, by = "Ptid")
  }
}
  

# Save inputFile 
if(study_name %in% c("VAT08m", "VAT08b")){
  save(inputFile_with_riskscore, file = paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/inputFile_with_riskscore.RData"))
  save(risk_scores, file = paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/risk_scores.RData"))
}else{
  save(inputFile_with_riskscore, file = paste0("output/", Sys.getenv("TRIAL"), "/inputFile_with_riskscore.RData"))
}

# Create table of cases in both arms (post Risk score analyses)
tab <- inputFile_with_riskscore 
# # For Placebo group, consider all COVID cases occurring after Day 1.
# if(study_name %in% c("PREVENT19")){
#   tab <- tab %>%
#     mutate(EventIndPrimaryD35 = ifelse(Trt == 0 & !is.na(EventIndPrimaryD1) & (EventIndPrimaryD1==1 | EventIndPrimaryD35==1), 1, EventIndPrimaryD35))
# }
# if(study_name %in% c("AZD1222")){
#   tab <- tab %>%
#     mutate(EventIndPrimaryD57 = ifelse(Trt == 0 & !is.na(EventIndPrimaryD1) & (EventIndPrimaryD1==1 | EventIndPrimaryD57==1), 1, EventIndPrimaryD57))
# }
tab <- tab %>%
  filter(Riskscorecohortflag == 1) %>%
  drop_na(Ptid, Trt, all_of(endpoint)) %>%
  mutate(Trt = ifelse(Trt == 0, "Placebo", "Vaccine")) 

if(study_name == "PREVENT19")
  tab <- tab %>% filter(Country == 0)

if(study_name %in% c("VAT08m", "VAT08b")){
  table(tab$Trt, tab %>% pull(endpoint)) %>%
    write.csv(file = here("output", Sys.getenv("TRIAL"), args[1], "cases_post_riskScoreAnalysis.csv"))
}else{
  table(tab$Trt, tab %>% pull(endpoint)) %>%
    write.csv(file = here("output", Sys.getenv("TRIAL"), "cases_post_riskScoreAnalysis.csv"))
}

rm(tab)
