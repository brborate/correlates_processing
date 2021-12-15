# Sys.setenv(TRIAL = "janssen_pooled_real")
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
load("output/objects_for_running_SL.rda")
load("output/inputFile.RData")
placebos_risk <- read.csv(here("output", "placebo_ptids_with_riskscores.csv"))
vaccinees_risk <- read.csv(here("output", "vaccine_ptids_with_riskscores.csv"))

# merge risk score with cleaned data by IDs, then save updated data file
risk_scores <- rbind(placebos_risk, vaccinees_risk) %>%
  select(Ptid, risk_score, standardized_risk_score)
inputFile_with_riskscore <- left_join(inputFile, risk_scores, by = "Ptid") 
  #%>%
  #filter(Riskscorecohortflag == 1) %>% 
  #drop_na(all_of(endpoint))
# data_name_amended <- paste0(str_remove(paste0(attr(config, "config"), "_data_processed.csv"), ".csv"), "_with_riskscore")

# # Ensure all baseline negative and PP subjects have a risk score!
# if(assertthat::assert_that(
#   all(!is.na(inputFile_with_riskscore %>% .$risk_score)), 
#           msg = "Some baseline negative and PP subjects have NA values in risk score!"
# )){
#   write_csv(inputFile_with_riskscore,
#             here("..", "data_clean", paste0(data_name_amended, ".csv")))
# }

# Save inputFile 
save(inputFile_with_riskscore, file = here("output", "inputFile_with_riskscore.RData"))
# inputFile_with_riskscore %>% 
#   write.csv(here::here("..", "..", "..", "..", "covpn", "p3003", 
#                        "analysis", "correlates", "Part_A_Blinded_Phase_Data", "adata", paste0(attr(config, "config"), "_data_processed_with_riskscoreNEW.csv"))) 

# Create table of cases in both arms (post Risk score analyses)
tab <- inputFile_with_riskscore %>%
  filter(Riskscorecohortflag == 1) %>%
  drop_na(Ptid, Trt, all_of(endpoint)) %>%
  mutate(Trt = ifelse(Trt == 0, "Placebo", "Vaccine")) 
table(tab$Trt, tab %>% pull(endpoint)) %>%
  write.csv(file = here("output", "cases_post_riskScoreAnalysis.csv"))
rm(tab)
