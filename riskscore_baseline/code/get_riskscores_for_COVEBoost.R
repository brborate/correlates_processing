# This code generates risk scores for all subjects in Stage 2 moderna Boost trial 
# and adds them to the mapped_data.
# The risk scores are unchanged for subjects in Stage 2 placebo arm that were involved in training of SL model 
# in Stage 1 data. For all other subjects in Stage 2, the risk scores are predicted using the SL model generated in Stage 1.
# Scaling of Stage 2 risk variables prior to prediction is done separately for placebo and vaccine arms, using the scaling parameters (mean and sd) from Stage 1 data. 

# This code is run from riskscore_baseline directory using the following command: 
# Rscript code/get_riskscores_for_COVEBoost.R

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

print("GET_RISKSCORES_FOR_COVEBoost.R")

print(mapped_data)
if (!file.exists(mapped_data)) stop ("mapped_data dataset not available ===========================================")

# Define code version to run
# the demo version is simpler and runs faster!
# the production version runs SL with a diverse set of learners
run_prod <- !grepl("Mock", study_name)

# get utility files
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

inputFile <- read.csv(mapped_data) %>%
  rename(Ptid = Subjectid)
load(paste0("output/moderna_real/objects_for_running_SL.rda"))
load(paste0("output/moderna_real/sl_riskscore_slfits.rda"))
load(paste0("output/moderna_real/X_covars2adjust_scaled_vacc.rda"))

# Predict on all stage 2 subjects
dat.ph1 <- inputFile %>%
    #filter(Riskscorecohortflag == 1 & Trt == 1) %>%
    # Keep only variables to be included in risk score analyses
    select(Ptid, Trt, all_of(risk_vars)) %>%
    # Drop any observation with NA values in Ptid, Trt, or endpoint!
    drop_na(Ptid, Trt)

########################################################################
########################################################################
# Predict on all stage 2 placebo subjects using SL model from Stage 1
dat.ph1.plac <- dat.ph1 %>% filter(Trt == 0)

X_covars2adjust_plac <- dat.ph1.plac %>%
  select(all_of(risk_vars))

# Impute missing values in any variable included in risk_vars using the mice package!
print("Make sure data is clean before conducting imputations!")
X_covars2adjust_plac <- impute_missing_values(X_covars2adjust_plac, risk_vars)

# Scale X_covars2adjust_plac according to scale parameters derived from Stage 1 placebo data
X_covars2adjust_plac_scaled <- scale(X_covars2adjust_plac,
                                     center = attr(X_covars2adjust_scaled_plac, "scaled:center"), 
                                     scale = attr(X_covars2adjust_scaled_plac, "scaled:scale"))

attr(X_covars2adjust_plac_scaled, "scaled:center") <- NULL
attr(X_covars2adjust_plac_scaled, "scaled:scale") <- NULL

X_covars2adjust_plac_scaled <- data.frame(X_covars2adjust_plac_scaled)

pred_on_plac <- predict(sl_riskscore_slfits, newdata = X_covars2adjust_plac_scaled, onlySL = TRUE)$pred %>%
  as.data.frame()

plac <- bind_cols(dat.ph1.plac, pred_on_plac) %>%
  rename(pred = V1) %>%
  mutate(risk_score = log(pred / (1 - pred))) 

# For ptids used in training SL model in Stage 1, use the risk scores based off CV-predictions in stage 1
load(file = here("output", "moderna_real", "risk_placebo_ptids.rda"))

cove_plac_SL_training <- read.csv(mapped_data) %>%
  select(Ptid, risk_score) %>%
  filter(Ptid %in% risk_placebo_ptids$Ptid) %>%
  rename(risk_score_SLtrain = risk_score)

plac <- plac %>% left_join(cove_plac_SL_training, by = "Ptid") %>%
  mutate(risk_score = ifelse(!is.na(risk_score_SLtrain), risk_score_SLtrain, risk_score))

########################################################################
########################################################################
# Predict on all stage 2 placebo subjects using SL model from Stage 1
dat.ph1.vacc <- dat.ph1 %>% filter(Trt == 1)

X_covars2adjust_vacc <- dat.ph1.vacc %>%
  select(all_of(risk_vars))

# Impute missing values in any variable included in risk_vars using the mice package!
print("Make sure data is clean before conducting imputations!")
X_covars2adjust_vacc <- impute_missing_values(X_covars2adjust_vacc, risk_vars)

# Scale X_covars2adjust_vacc according to scale parameters derived from Stage 1 vaccebo data
X_covars2adjust_vacc_scaled <- scale(X_covars2adjust_vacc,
                                     center = attr(X_covars2adjust_scaled_vacc, "scaled:center"), 
                                     scale = attr(X_covars2adjust_scaled_vacc, "scaled:scale"))

attr(X_covars2adjust_vacc_scaled, "scaled:center") <- NULL
attr(X_covars2adjust_vacc_scaled, "scaled:scale") <- NULL

X_covars2adjust_vacc_scaled <- data.frame(X_covars2adjust_vacc_scaled)

pred_on_vacc <- predict(sl_riskscore_slfits, newdata = X_covars2adjust_vacc_scaled, onlySL = TRUE)$pred %>%
  as.data.frame()

vacc <- bind_cols(dat.ph1.vacc, pred_on_vacc) %>%
  rename(pred = V1) %>%
  mutate(risk_score = log(pred / (1 - pred))) 

##################################################################
inputFile_with_riskscore <- inputFile %>% 
  left_join(bind_rows(plac, vacc) %>% select(Ptid, risk_score), by = "Ptid") 
  
save(inputFile,
     file = here("output", Sys.getenv("TRIAL"), "inputFile.rda"))
save(inputFile_with_riskscore,
     file = here("output", Sys.getenv("TRIAL"), "inputFile_with_riskscore.rda"))

print("Risk scores created for moderna_boost TRIAL ===========================================")

