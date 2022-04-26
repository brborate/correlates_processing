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
load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))
load(paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
placebos_risk <- read.csv(here("output", Sys.getenv("TRIAL"), "placebo_ptids_with_riskscores.csv"))
vaccinees_risk <- read.csv(here("output", Sys.getenv("TRIAL"), "vaccine_ptids_with_riskscores.csv"))

# merge risk score with cleaned data by IDs, then save updated data file
risk_scores <- bind_rows(placebos_risk, vaccinees_risk) %>%
  select(Ptid, risk_score, standardized_risk_score)
inputFile_with_riskscore <- left_join(inputFile, risk_scores, by = "Ptid") 

# Save inputFile 
save(inputFile_with_riskscore, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData"))

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
table(tab$Trt, tab %>% pull(endpoint)) %>%
  write.csv(file = here("output", Sys.getenv("TRIAL"), "cases_post_riskScoreAnalysis.csv"))
rm(tab)
