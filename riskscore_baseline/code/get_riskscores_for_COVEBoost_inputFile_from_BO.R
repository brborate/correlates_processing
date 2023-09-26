#Sys.setenv(TRIAL = "moderna_boost")

renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
config <- config::get(config = Sys.getenv("TRIAL"))
for(opt in names(config)){
  eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
}

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

load(paste0("output/moderna_real/objects_for_running_SL.rda"))
load(paste0("output/moderna_real/sl_riskscore_slfits.rda"))
load(paste0("output/moderna_real/X_covars2adjust_scaled_vacc.rda"))

########################################################################
########################################################################
# Read in COVEBoost processed data which has risk scores!
coveboost_processed_data = read.csv("/trials/covpn/p3001/analysis/mapping_immune_correlates/Part_C_Unblinded_Phase_Data/adata/COVID_Moderna_stage2_20230802_withRiskScores.csv")

# Read in Bo's new input file which has additional PTIDS for which risk score has to be derived
inputFile_from_Bo <- read.csv("/trials/covpn/p3001/analysis/mapping_immune_correlates/Part_C_Unblinded_Phase_Data/adata/dat_calendar_time_CoP_0806.csv") %>% 
  rename(Ptid = SUBJID) 

dat_plac2 = inputFile_from_Bo %>%
  filter(!Ptid %in% coveboost_processed_data$Ptid) %>% # Consider only new subjects in Bo's file for which risk scores haven't been derived!
  filter(Trt == 0)

X_covars2adjust_plac2 <- dat_plac2 %>%
  select(all_of(risk_vars))

# Impute missing values in any variable included in risk_vars using the mice package!
print("Make sure data is clean before conducting imputations!")
X_covars2adjust_plac2 <- impute_missing_values(X_covars2adjust_plac2, risk_vars)

# Scale X_covars2adjust_plac according to scale parameters derived from Stage 1 placebo data
X_covars2adjust_plac_scaled2 <- scale(X_covars2adjust_plac2,
                                     center = attr(X_covars2adjust_scaled_plac, "scaled:center"), 
                                     scale = attr(X_covars2adjust_scaled_plac, "scaled:scale"))

attr(X_covars2adjust_plac_scaled2, "scaled:center") <- NULL
attr(X_covars2adjust_plac_scaled2, "scaled:scale") <- NULL

X_covars2adjust_plac_scaled2 <- data.frame(X_covars2adjust_plac_scaled2)

pred_on_plac2 <- predict(sl_riskscore_slfits, newdata = X_covars2adjust_plac_scaled2, onlySL = TRUE)$pred %>%
  as.data.frame()

plac2 <- bind_cols(dat_plac2, pred_on_plac2) %>%
  rename(pred = V1) %>%
  mutate(risk_score = log(pred / (1 - pred))) 

# For ptids used in training SL model in Stage 1, use the risk scores based off CV-predictions in stage 1
load(file = here("output", "moderna_real", "risk_placebo_ptids.rda"))

cove_plac_SL_training <- read.csv(mapped_data) %>%
  select(Ptid, risk_score) %>%
  filter(Ptid %in% risk_placebo_ptids$Ptid) %>%
  rename(risk_score_SLtrain = risk_score)

plac2 <- plac2 %>% left_join(cove_plac_SL_training, by = "Ptid") %>%
  mutate(risk_score = ifelse(!is.na(risk_score_SLtrain), risk_score_SLtrain, risk_score))


########################################################################
########################################################################
# Predict on all stage 2 placebo subjects using SL model from Stage 1
dat_vacc2 = inputFile_from_Bo %>%
  filter(!Ptid %in% coveboost_processed_data$Ptid) %>% # Consider only new subjects in Bo's file for which risk scores haven't been derived!
  filter(Trt == 1)

X_covars2adjust_vacc2 <- dat_vacc2 %>%
  select(all_of(risk_vars))

# Impute missing values in any variable included in risk_vars using the mice package!
print("Make sure data is clean before conducting imputations!")
X_covars2adjust_vacc2 <- impute_missing_values(X_covars2adjust_vacc2, risk_vars)

# Scale X_covars2adjust_vacc according to scale parameters derived from Stage 1 vaccinee data (stored in X_covars2adjust_scaled_vacc)
X_covars2adjust_vacc_scaled2 <- scale(X_covars2adjust_vacc2,
                                            center = attr(X_covars2adjust_scaled_vacc, "scaled:center"), 
                                            scale = attr(X_covars2adjust_scaled_vacc, "scaled:scale"))

attr(X_covars2adjust_vacc_scaled2, "scaled:center") <- NULL
attr(X_covars2adjust_vacc_scaled2, "scaled:scale") <- NULL

X_covars2adjust_vacc_scaled2 <- data.frame(X_covars2adjust_vacc_scaled2)

pred_on_vacc2 <- predict(sl_riskscore_slfits, newdata = X_covars2adjust_vacc_scaled2, onlySL = TRUE)$pred %>%
  as.data.frame()

vacc2 <- bind_cols(dat_vacc2, pred_on_vacc2) %>%
  rename(pred = V1) %>%
  mutate(risk_score = log(pred / (1 - pred))) 

##################################################################
inputFile_with_riskscore <- inputFile_from_Bo %>% 
  left_join(bind_rows(plac2, vacc2) %>% select(Ptid, risk_score), by = "Ptid") %>%
  left_join(coveboost_processed_data %>% select(Ptid, risk_score) %>% rename(risk_score_processeddata = risk_score), by = "Ptid") %>%
  mutate(risk_score = ifelse(is.na(risk_score), risk_score_processeddata, risk_score)) %>%
  select(-c(risk_score_processeddata))

inputFile_with_riskscore %>%
  write.csv("/trials/covpn/p3001/analysis/mapping_immune_correlates/Part_C_Unblinded_Phase_Data/adata/dat_calendar_time_CoP_0806_withRiskScores.csv", row.names = FALSE)


