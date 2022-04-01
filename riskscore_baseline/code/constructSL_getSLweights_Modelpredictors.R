# Sys.setenv(TRIAL = "janssen_pooled_realbAb")
renv::activate(here::here(".."))
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
library(aucm)
library(mice)
library(conflicted)
library(gam)
library(xgboost)
library(ranger)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("omp_set_num_threads", "RhpcBLASctl")
load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))
load(paste0("output/", Sys.getenv("TRIAL"), "/plac_top2learners_SL_discreteSL.rda"))
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs

## solve cores issue
#blas_get_num_procs()
blas_set_num_threads(1)
#print(blas_get_num_procs())
stopifnot(blas_get_num_procs() == 1)

## construct superlearner on placebo arm-----------------------
set.seed(20210216)
sl_riskscore_slfits <- SuperLearner(
  Y = Y, X = X_riskVars, family = familyVar,
  SL.library = SL_library, method = methodVar,
  cvControl = cvControlVar, verbose = FALSE
)

save(sl_riskscore_slfits, file = here("output", Sys.getenv("TRIAL"), "sl_riskscore_slfits.rda"))

# Get Superlearner weights
sl_weights <- sl_riskscore_slfits$coef %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Learner") %>%
  rename(`Weights` = ".") %>%
  arrange(-Weights)

sl_weights %>%
  mutate(Weight = format(round(Weights, 3), nsmall = 3)) %>%
  mutate(
    Screen = paste0("screen_", sapply(strsplit(Learner, "_screen_"), `[`, 2)),
    Learner = sapply(strsplit(Learner, "_screen"), `[`, 1)
  ) %>%
  select(Learner, Screen, Weight) %>%
  write.csv(here("output", Sys.getenv("TRIAL"), "SL_weights.csv"))

top_models <- sl_weights %>%
  filter(Learner != "SL.mean_screen_all") %>%
  .$Learner

# Get predictors selected in the models with highest weights
for (i in seq_along(top_models)) {
  #print(i)
  if(top_models[i] %in% c("SL.glm_screen_univariate_logistic_pval", 
                          "SL.glm.interaction_screen_highcor_random",
                          "SL.glm_screen_all",
                          "SL.glm_screen_glmnet",
                          "SL.glm_screen_highcor_random",
                          "SL.glm.interaction_screen_glmnet",
                          "SL.glm.interaction_screen_univariate_logistic_pval",
                          "SL.gam_screen_glmnet",
                          "SL.gam_screen_univariate_logistic_pval",
                          "SL.gam_screen_highcor_random")) {

    model <- sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object$coefficients %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Predictors") %>%
      rename(`Coefficient` = ".") %>%
      mutate(
        `Odds Ratio` = exp(`Coefficient`),
        Learner = top_models[i])
  }

  if (top_models[i] %in% c("SL.glmnet_screen_all")) {
    model <- coef(sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object, s = "lambda.min") %>%
      as.matrix() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Predictors") %>%
      rename(`Coefficient` = "s1") %>%
      mutate(`Odds Ratio` = exp(`Coefficient`),
             Learner = top_models[i])
  }

  if (top_models[i] %in% c("SL.xgboost_screen_all")) {
    model <- xgboost::xgb.importance(model = sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object) %>%
      as.data.frame() %>%
      mutate(Learner = top_models[i])
  }

  if (top_models[i] %in% c("SL.ranger.imp_screen_all")) {
    model <- sl_riskscore_slfits[["fitLibrary"]][[top_models[i]]]$object$variable.importance %>%
      as.data.frame() %>%
      rename(Importance = ".") %>%
      tibble::rownames_to_column(var = "Predictors") %>%
      mutate(Learner = top_models[i])
  }

  if (top_models[i] == "SL.mean_screen_all")
	next

  if (i == 1) {
    all_models <- model
  } else {
    all_models <- bind_rows(all_models, model)
  }
}

options(scipen=999)

if(run_prod){
  all_models %>%
    left_join(sl_weights, by = "Learner") %>%
    mutate(
      Weight = format(round(Weights, 3), nsmall = 3),
      Coefficient = format(round(Coefficient, 3), nsmall = 3),
      `Odds Ratio` = format(round(`Odds Ratio`, 3), nsmall = 3),
      Importance = format(round(Importance, 3), nsmall = 3),
      Gain = format(round(Gain, 3), nsmall = 3),
      Cover = format(round(Cover, 3), nsmall = 3),
      Frequency = format(round(Frequency, 3), nsmall = 3),
    ) %>%
    mutate(
      Screen = paste0("screen_", sapply(strsplit(Learner, "_screen_"), `[`, 2)),
      Learner = sapply(strsplit(Learner, "_screen"), `[`, 1)
    ) %>%
    select(Learner, Screen, Weight, Predictors, Coefficient, `Odds Ratio`,
           Importance, Feature, Gain, Cover, Frequency) %>%
    write.csv(here("output", Sys.getenv("TRIAL"), "SL_all_models_with_predictors.csv"))
}else{
  all_models %>%
    left_join(sl_weights, by = "Learner") %>%
    mutate(
      Weight = format(round(Weights, 3), nsmall = 3),
      Coefficient = format(round(Coefficient, 3), nsmall = 3),
      `Odds Ratio` = format(round(`Odds Ratio`, 3), nsmall = 3),
    ) %>%
    mutate(
      Screen = paste0("screen_", sapply(strsplit(Learner, "_screen_"), `[`, 2)),
      Learner = sapply(strsplit(Learner, "_screen"), `[`, 1)
    ) %>%
    select(Learner, Screen, Weight, Predictors, Coefficient, `Odds Ratio`) %>%
    write.csv(here("output", Sys.getenv("TRIAL"), "SL_all_models_with_predictors.csv"))
}


rm(sl_riskscore_slfits)

