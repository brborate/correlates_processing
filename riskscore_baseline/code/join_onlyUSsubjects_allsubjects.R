# This code is called only for study prevent19 (Novavax), where risk scores are developed for 
# US and Mexican subjects separately.
# For US subjects, the SL used to derive risk scores is trained using only US subjects. 
# For Mexican subjects, the SL used to derive risk scores is trained using both US and Mexican subjects.

args = commandArgs(trailingOnly=TRUE)

source("code/loadlibraries_readinputdata.R")

# load required libraries, cleaned data, and risk score estimates
library(here)
library(tidyverse)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")

print("JOIN_onlyUSsubjects_allsubjects_riskscores.R")

load("output/prevent19/onlyUSsubjects/inputFile_with_riskscore.RData")
onlyUS_riskscores_inputFile = inputFile_with_riskscore
rm(inputFile_with_riskscore)
load("output/prevent19/allsubjects/inputFile_with_riskscore.RData")
all_riskscores_inputFile = inputFile_with_riskscore %>%
  rename(risk_score2 = risk_score,
         standardized_risk_score2 = standardized_risk_score) %>%
  select(Ptid, risk_score2, standardized_risk_score2)
rm(inputFile_with_riskscore)

inputFile_with_riskscore <- onlyUS_riskscores_inputFile %>% 
  left_join(all_riskscores_inputFile, by = "Ptid")

save(inputFile_with_riskscore,
     file = here("output", "prevent19", "inputFile_with_riskscore.RData")) 

source(here("code", "unit_test.R"))
