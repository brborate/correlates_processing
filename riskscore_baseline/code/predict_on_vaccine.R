# Sys.setenv(TRIAL = "janssen_pooled_real")
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
source(here::here("code", "utils.R"))
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
library(aucm)
library(mice)
library(conflicted)
library(gam)
library(xgboost)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

print("PREDICT_ON_VACCINE.R")

# conflict_prefer("omp_set_num_threads", "RhpcBLASctl")
if(study_name %in% c("VAT08m", "VAT08b")){
  load(paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/objects_for_running_SL.rda"))
  load(paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/sl_riskscore_slfits.rda"))
}else{
  load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))
  load(paste0("output/", Sys.getenv("TRIAL"), "/sl_riskscore_slfits.rda"))
}

# load(paste0("output/", Sys.getenv("TRIAL"), "/plac_top2learners_SL_discreteSL.rda"))
# source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
# source(here("code", "utils.R")) # get CV-AUC for all algs

# Predict on vaccine arm
if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
  dat.ph1.vacc <- inputMod %>%
    filter(Riskscorecohortflag == 1 & Trt == 1) %>%
    # Keep only variables to be included in risk score analyses
    select(Ptid, Trt, all_of(endpoint), paste0(sub("1rscore", "", endpoint), paste0(vaccAUC_timepoint, "rauc")), 
           all_of(risk_vars), RiskscoreAUCflag) %>%
    # Drop any observation with NA values in Ptid, Trt, or endpoint!
    drop_na(Ptid, Trt, all_of(endpoint))
} else {
  dat.ph1.vacc <- inputMod %>%
    filter(Riskscorecohortflag == 1 & Trt == 1) %>%
    # Keep only variables to be included in risk score analyses
    select(Ptid, Trt, all_of(endpoint), all_of(risk_vars)) %>%
    # Drop any observation with NA values in Ptid, Trt, or endpoint!
    drop_na(Ptid, Trt, all_of(endpoint))
}

X_covars2adjust_vacc <- dat.ph1.vacc %>%
  select(all_of(risk_vars))

# Impute missing values in any variable included in risk_vars using the mice package!
print("Make sure data is clean before conducting imputations!")
X_covars2adjust_vacc <- impute_missing_values(X_covars2adjust_vacc, risk_vars)

# Scale X_covars2adjust to have mean 0, sd 1 for all vars
X_covars2adjust_scaled_vacc <- get_scaleParams_scaledData(X_covars2adjust_vacc)
X_covars2adjust_scaled_vacc_noattr <- X_covars2adjust_scaled_vacc
attr(X_covars2adjust_scaled_vacc_noattr, "scaled:center") <- NULL
attr(X_covars2adjust_scaled_vacc_noattr, "scaled:scale") <- NULL

# # Scale X_covars2adjust_vacc to have mean 0, sd 1 for all vars
# for (a in colnames(X_covars2adjust_vacc)) {
#   X_covars2adjust_vacc[[a]] <- scale(X_covars2adjust_vacc[[a]],
#     center = mean(X_covars2adjust_vacc[[a]], na.rm = T),
#     scale = sd(X_covars2adjust_vacc[[a]], na.rm = T)
#   )
# }

X_riskVars_vacc <- data.frame(X_covars2adjust_scaled_vacc_noattr)

pred_on_vaccine <- predict(sl_riskscore_slfits, newdata = X_riskVars_vacc, onlySL = TRUE)$pred %>%
  as.data.frame()

if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
  vacc <- dat.ph1.vacc %>% select(Ptid, all_of(endpoint), paste0(sub("1rscore", "", endpoint), paste0(vaccAUC_timepoint, "rauc")), RiskscoreAUCflag)
} else {
  vacc <- dat.ph1.vacc %>% select(Ptid, all_of(endpoint))
}

vacc <- bind_cols(vacc, pred_on_vaccine) %>%
  rename(pred = V1) %>%
  mutate(risk_score = log(pred / (1 - pred)),
         standardized_risk_score = scale(risk_score,
                            center = mean(risk_score, na.rm = T),
                            scale = sd(risk_score, na.rm = T))) 


if(study_name == "COVE"){
  # Predict on baseline seropositives!
  pred_on_plac_bseropos <- predict(sl_riskscore_slfits, 
                                   newdata = plac_bseropos %>% select(names(X_riskVars_vacc)), 
                                   onlySL = TRUE)$pred %>%
    as.data.frame()
  
  plac_bseropos <- bind_cols(plac_bseropos %>% select(Ptid, all_of(endpoint)), pred_on_plac_bseropos) %>%
    rename(pred = V1) %>%
    mutate(risk_score = log(pred / (1 - pred)),
           standardized_risk_score = scale(risk_score,
                                           center = mean(risk_score, na.rm = T),
                                           scale = sd(risk_score, na.rm = T))) 
  
  # Predict on the 6 extra placebo arm + 9 extra vaccine arm subjects in COVEboost (Moderna Stage 2 trial)!
  cove <- read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdata_correlates_originaldatafromModerna_v1.0_Oct28_2021.csv")
  coveboost <- read.csv("/trials/covpn/p3001/analysis/mapping_immune_correlates/Part_C_Unblinded_Phase_Data/adata/COVID_Moderna_stage2_20230710.csv") %>%
    rename(Ptid = Subjectid)
  
  # First, predict on the 6 extra placebo arm subjects in COVEboost trial!
  coveboost_6extra_plac <- coveboost %>% filter(!Ptid %in% cove$Ptid) %>% filter(Trt == 0) 
  coveboost_6extra_plac_riskvars <- coveboost_6extra_plac %>% select(all_of(risk_vars))
  # Scale based off mean and sd for X_covars2adjust_scaled_plac
  scaled_coveboost_6extra_plac_riskvars <- scale(coveboost_6extra_plac_riskvars, 
                                                 center = attr(X_covars2adjust_scaled_plac, "scaled:center"), 
                                                 scale = attr(X_covars2adjust_scaled_plac, "scaled:scale"))
  attr(scaled_coveboost_6extra_plac_riskvars, "scaled:center") <- NULL
  attr(scaled_coveboost_6extra_plac_riskvars, "scaled:scale") <- NULL
  
  # Predict!
  pred_on_scaled_coveboost_6extra_plac_riskvars <- predict(sl_riskscore_slfits, 
                                                           newdata = data.frame(scaled_coveboost_6extra_plac_riskvars),
                                                           onlySL = TRUE)$pred %>%
    as.data.frame()
  
  plac_coveboost_6extra <- bind_cols(coveboost_6extra_plac %>% select(Ptid), pred_on_scaled_coveboost_6extra_plac_riskvars) %>%
    rename(pred = V1) %>%
    mutate(risk_score = log(pred / (1 - pred)),
           standardized_risk_score = scale(risk_score,
                                           center = mean(risk_score, na.rm = T),
                                           scale = sd(risk_score, na.rm = T))) 
  
  # Second, predict on the 9 extra vaccine arm subjects in COVEboost trial!
  coveboost_9extra_vacc <- coveboost %>% filter(!Ptid %in% cove$Ptid) %>% filter(Trt == 1)
  coveboost_9extra_vacc_riskvars <- coveboost_9extra_vacc %>% select(all_of(risk_vars))
  # Scale based off mean and sd for X_covars2adjust_scaled_plac
  scaled_coveboost_9extra_vacc_riskvars <- scale(coveboost_9extra_vacc_riskvars, 
                                                 center = attr(X_covars2adjust_scaled_vacc, "scaled:center"), 
                                                 scale = attr(X_covars2adjust_scaled_vacc, "scaled:scale"))
  
  attr(scaled_coveboost_9extra_vacc_riskvars, "scaled:center") <- NULL
  attr(scaled_coveboost_9extra_vacc_riskvars, "scaled:scale") <- NULL
  
  # Predict!
  pred_on_scaled_coveboost_9extra_vacc_riskvars <- predict(sl_riskscore_slfits, 
                                                           newdata = data.frame(scaled_coveboost_9extra_vacc_riskvars),
                                                           onlySL = TRUE)$pred %>%
    as.data.frame()
  
  vacc_coveboost_9extra <- bind_cols(coveboost_9extra_vacc %>% select(Ptid), pred_on_scaled_coveboost_9extra_vacc_riskvars) %>%
    rename(pred = V1) %>%
    mutate(risk_score = log(pred / (1 - pred)),
           standardized_risk_score = scale(risk_score,
                                           center = mean(risk_score, na.rm = T),
                                           scale = sd(risk_score, na.rm = T))) 
  
}


if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
  # AUC for vaccine arm computed only on cohort with RiskscoreAUCflag==1 and based off endpoint rauc!
  AUCvacc <- vacc %>% 
    filter(RiskscoreAUCflag == 1) %>%
    mutate(AUCchar = format(round(fast.auc(pred, get(paste0(sub("1rscore", "", endpoint), paste0(vaccAUC_timepoint, "rauc")))), 3), nsmall = 3)) %>%
    distinct(AUCchar)
} else {
  AUCvacc <- vacc %>% 
    mutate(AUCchar = format(round(fast.auc(pred, get(endpoint)), 3), nsmall = 3)) %>%
    distinct(AUCchar)
}


vacc <- vacc %>% mutate(AUCchar = AUCvacc$AUCchar) 

if(study_name %in% c("VAT08m", "VAT08b")){
  write.csv(vacc, here("output", Sys.getenv("TRIAL"), args[1], "vaccine_ptids_with_riskscores.csv"), row.names = FALSE)
}else{
  write.csv(vacc, here("output", Sys.getenv("TRIAL"), "vaccine_ptids_with_riskscores.csv"), row.names = FALSE)
}

if(study_name == "COVE"){
  write.csv(plac_bseropos, here("output", Sys.getenv("TRIAL"), "plac_bseropos_ptids_with_riskscores.csv"), row.names = FALSE)
  write.csv(plac_coveboost_6extra, here("output", Sys.getenv("TRIAL"), "plac_coveboost_6extra_ptids_with_riskscores.csv"), row.names = FALSE)
  write.csv(vacc_coveboost_9extra, here("output", Sys.getenv("TRIAL"), "vacc_coveboost_9extra_ptids_with_riskscores.csv"), row.names = FALSE)
  save(X_covars2adjust_scaled_vacc, file = here("output", Sys.getenv("TRIAL"), "X_covars2adjust_scaled_vacc.rda"))
}
