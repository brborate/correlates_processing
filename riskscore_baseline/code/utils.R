#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
if(!study_name == "COVEBoost"){
  source(here::here("..", "_common.R"))
}

#-----------------------------------------------

# #-----------------------------------------------
# # obligatory to append to the top of each script
# renv::activate(project = here::here())
# # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
# if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
# source(here::here("_common.R"))
# #-----------------------------------------------

# --------------------------------------------------------------------
# Cross-validated predictiveness
# --------------------------------------------------------------------
# get the CV-AUC for a single learner's predicted values
# @param preds the fitted values
# @param Y the outcome
# @param full_y the observed outcome (from the entire dataset, for cross-fitted estimates)
# @param scale what scale should the IPCW correction be applied on?
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale,
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5)
#             for 5-fold cross-validated super learner)
one_auc <- function(preds, Y, full_y = NULL, scale = "identity",
                    weights = rep(1, length(Y)), C = rep(1, length(Y)),
                    Z = NULL, ...) {
  auc_lst <- vimp::measure_auc(
    fitted_values = preds, y = Y, full_y = full_y, C = C, Z = Z,
    ipc_weights = weights,
    ipc_fit_type = "SL", ...
  )
  list(auc = auc_lst$point_est, eif = auc_lst$eif)
}

# get the cross-fitted CV-AUC for a single learner's predicted values
# @param preds the fitted values
# @param Y the outcome
# @param folds the different cv folds that the learner was evaluated on
# @param scale what scale should the IPCW correction be applied on?
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale,
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5)
#             for 5-fold cross-validated super learner)
cv_auc <- function(preds, Y, folds, scale = "identity",
                   weights = rep(1, length(Y)), C = rep(1, length(Y)),
                   Z = NULL, ...) {
  V <- length(folds)
  folds_numeric <- get_cv_sl_folds(folds)
  if (is.null(Z)) {
    folds_z <- folds_numeric
  } else {
    folds_z <- c(folds_numeric, sample(seq_len(V), nrow(Z) - length(folds_numeric),
      replace = TRUE
    ))
  }
  ests_eifs <- lapply(as.list(seq_len(V)), function(v) {
    one_auc(
      preds = preds[folds_numeric == v], Y[folds_numeric == v],
      full_y = Y, scale = scale,
      weights = weights[folds_z == v], C = C[folds_z == v],
      Z = Z[folds_z == v, , drop = FALSE], ...
    )
  })
  est <- mean(unlist(lapply(ests_eifs, function(l) l$auc)))
  var <- mean(unlist(lapply(ests_eifs, function(l) mean(l$eif ^ 2))))
  se <- sqrt(var / length(Y))
  ci <- vimp::vimp_ci(est, se, scale = scale, level = 0.95)
  return(list(auc = est, se = se, ci = ci))
}

# get the folds from a CV.SL object, make them a vector
# @param cv_sl_folds the CV.SL folds (a named list of row numbers)
# @return a vector with the correct folds
get_cv_sl_folds <- function(cv_sl_folds) {
  folds_with_row_nums <- sapply(1:length(cv_sl_folds),
    function(x) {
      list(
        row_nums = cv_sl_folds[[x]],
        fold = rep(x, length(cv_sl_folds[[x]]))
      )
    },
    simplify = FALSE
  )
  folds_df <- data.table::rbindlist(folds_with_row_nums)
  folds_df$fold[order(folds_df$row_nums)]
}
# get the CV-AUC for all learners fit with SL
# @param sl_fit the super learner fit object
# @param scale what scale should the IPCW correction be applied on?
#              Can help with numbers outside of (0, 1)
#             ("identity" denotes the identity scale;
#              "logit" means that AUC is transformed to the logit scale,
#              correction is applied, and then back-transformed)
# @param weights the inverse probability of censoring weights
# @param C the indicator of being observed in phase 2 (1) or not (0)
# @param Z a matrix of predictors observed on all participants in phase 1
#          (can include the outcome)
# @param ... other arguments to measure_auc, but should include at least:
#            a library of learners (using arg "SL.library")
#            and may include control parameters for the super learner
#            (e.g., cvControl = list(V = 5)
#             for 5-fold cross-validated super learner)
get_all_aucs <- function(sl_fit, scale = "identity",
                         weights = rep(1, length(sl_fit$Y)),
                         C = rep(1, length(sl_fit$Y)),
                         Z = NULL, ...) {
  # get the CV-AUC of the SuperLearner predictions
  sl_auc <- cv_auc(
    preds = sl_fit$SL.predict, Y = sl_fit$Y,
    folds = sl_fit$folds,
    scale = scale, weights = weights, C = C, Z = Z, ...
  )
  out <- data.frame(
    Learner = "SL", Screen = "All", AUC = sl_auc$auc,
    se = sl_auc$se, ci_ll = sl_auc$ci[1], ci_ul = sl_auc$ci[2]
  )

  # Get the CV-auc of the Discrete SuperLearner predictions
  discrete_sl_auc <- cv_auc(
    preds = sl_fit$discreteSL.predict, Y = sl_fit$Y,
    folds = sl_fit$folds, scale = scale,
    weights = weights, C = C,
    Z = Z, ...
  )
  out <- rbind(out, data.frame(
    Learner = "Discrete SL", Screen = "All",
    AUC = discrete_sl_auc$auc,
    se = discrete_sl_auc$se,
    ci_ll = discrete_sl_auc$ci[1],
    ci_ul = discrete_sl_auc$ci[2]
  ))

  # Get the cvauc of the individual learners in the library
  get_individual_auc <- function(sl_fit, col, scale = "identity",
                                 weights = rep(1, length(sl_fit$Y)),
                                 C = rep(1, length(sl_fit$Y)), Z = NULL, ...) {
    if (any(is.na(sl_fit$library.predict[, col]))) {
      return(NULL)
    }
    alg_auc <- cv_auc(
      preds = sl_fit$library.predict[, col], Y = sl_fit$Y,
      scale = scale,
      folds = sl_fit$folds, weights = weights,
      C = C, Z = Z, ...
    )
    # get the regexp object
    alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_",
      fixed = TRUE
    )[[1]]
    alg <- tail(alg_screen_string[grepl(".", alg_screen_string,
      fixed = TRUE
    )], n = 1)
    screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string,
      fixed = TRUE
    )],
    collapse = "_"
    )
    data.frame(
      Learner = alg, Screen = screen, AUC = alg_auc$auc,
      se = alg_auc$se,
      ci_ll = alg_auc$ci[1], ci_ul = alg_auc$ci[2]
    )
  }
  other_aucs <- plyr::ldply(
    1:ncol(sl_fit$library.predict),
    function(x) {
      get_individual_auc(
        sl_fit = sl_fit,
        col = x,
        scale = scale,
        weights = weights,
        C = C, Z = Z, ...
      )
    }
  )
  rbind(out, other_aucs)
}


# -------------------------------------------------------------------------
# Run the CV Super Learner
# -------------------------------------------------------------------------
# run CV.SuperLearner for one given random seed
# @param seed the random number seed
# @param Y the outcome
# @param X_mat the covariates
# @param family the family, for super learner (e.g., "binomial")
# @param obsWeights the inverse probability of censoring weights
#                   (for super learner)
# @param sl_lib the super learner library (e.g., "SL.ranger")
# @param method the method for determining the optimal combination
#               of base learners
# @param cvControl a list of control parameters to pass to the
#                  outer super learner
# @param innerCvControl a list of control parameters to pass to the
#                       inner super learners
# @param all_weights the IPC weights for variable importance (full data)
# @param Z the entire (phase 1) dataset (required only if analysis
#                  involves data from phase 2)
# @param C the outcome from the entire (phase 1) dataset
# @param z_lib the learner/s to be used for weighted auc calculations
# @param scale the scale that the IPC correction should be computed on
# @param vimp determines whether or not we save the entire SL fit object
#        (for variable importance, we don't need it so can save some memory
#         by excluding large objects)
run_cv_sl_once <- function(seed = 1, Y = NULL, X_mat = NULL,
                           family = "binomial",
                           obsWeights = rep(1, length(Y)),
                           sl_lib = "SL.ranger",
                           method = "method.CC_nloglik",
                           cvControl = list(V = 5),
                           innerCvControl = list(V = 5),
                           all_weights = rep(1, nrow(Z)),
                           Z = NULL,
                           C = rep(1, length(Y)),
                           z_lib = "SL.glm",
                           scale = "identity",
                           vimp = FALSE) {

  set.seed(seed)
  fit <- SuperLearner::CV.SuperLearner(
    Y = Y, X = X_mat, family = family,
    obsWeights = obsWeights,
    SL.library = sl_lib,
    method = method, cvControl = cvControl,
    innerCvControl = innerCvControl,
    verbose = FALSE
  )

  aucs <- get_all_aucs(sl_fit = fit, scale = scale)

  ret_lst <- list(fit = fit, folds = fit$folds, aucs = aucs)
  if (vimp) {
    ret_lst <- list(fit = fit$SL.predict, folds = fit$folds, aucs = aucs)
  }
  return(list(cvaucs = ret_lst, cvfits = fit))
}

##########################################################################################################
# remove any binary risk variables with fewer than 10 ptids that have a 0 or 1 for that variable
# @param dat the phase 1 dataset
# @param risk_vars the vector of column names of risk variables
# @return a data frame upon removal of any binary risk variables with fewer than 10 ptids that have a 0 or 1 for that variable (COVE analysis)
# @return a data frame upon removal of any binary risk variables with number of cases in the variable = 1 or 0 subgroup is <= 3 (ENSEMBLE analysis)
drop_riskVars_with_fewer_0s_or_1s <- function(dat, risk_vars, np) {

  if(study_name == "COVE"){
    # delete the file drop_riskVars_with_fewer_0s_or_1s.csv
    unlink(here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_fewer_0s_or_1s.csv"))
    # Remove binary risk variables with fewer than 10 ptids that have a 0 or 1 for that variable
    for (i in 1:length(risk_vars)) {
        if ((dat %>% select(starts_with(risk_vars[i])) %>% unique() %>% dim())[1] == 2) {
          if ((dim(dat %>% filter(get(risk_vars[i]) == 1))[1] < 10) | (dim(dat %>% filter(get(risk_vars[i]) == 0))[1] < 10)){
            dat <- dat %>% select(-starts_with(risk_vars[i]))
            print(paste0(risk_vars[i], " dropped from risk score analysis as it had fewer than 10 1's or 0's."))
            # Also print to file
            paste0(risk_vars[i], " dropped from risk score analysis as it had fewer than 10 1's or 0's.") %>%
              write.table(file = here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
          }
        }
      }
  }
  
  if(study_name == "COVE" & !file.exists(here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_fewer_0s_or_1s.csv"))){
    paste0("No binary input variable had fewer than 10 ptids with a 0 or 1 for that variable.") %>%
      write.table(file = here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
  }
  
  #####################################

  if(!study_name %in% c("COVE")){
      # delete the file drop_riskVars_with_fewer_0s_or_1s.csv
      if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
        unlink(here("output", Sys.getenv("TRIAL"), args[1], "drop_riskVars_with_fewer_0s_or_1s.csv"))
      } else {
        unlink(here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_fewer_0s_or_1s.csv"))
      }
      
      # Remove a variable if the number of cases in the variable = 1 subgroup is <= 3 or the number of cases in the variable = 0 subgroup is <= 3
      for (i in 1:length(risk_vars)) {
        print(i)
        if ((dat %>% select(starts_with(risk_vars[i])) %>% unique())[[1]] == c(0,1) | 
            (dat %>% select(starts_with(risk_vars[i])) %>% unique())[[1]] == c(1,0)) {
          if (dat %>% filter(get(risk_vars[i]) == 1) %>% pull(endpoint) %>% sum() <= 3 | dat %>% filter(get(risk_vars[i]) == 0) %>% pull(endpoint) %>% sum() <= 3){
            dat <- dat %>% select(-starts_with(risk_vars[i]))
            print(paste0(risk_vars[i], " dropped from risk score analysis as the number of cases in the variable = 1 or 0 subgroup is <= 3."))
            # Also print to file
            if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
              paste0(risk_vars[i], " dropped from risk score analysis as the number of cases in the variable = 1 or 0 subgroup is <= 3.") %>%
                write.table(file = here("output", Sys.getenv("TRIAL"), args[1], "drop_riskVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
            } else {
              paste0(risk_vars[i], " dropped from risk score analysis as the number of cases in the variable = 1 or 0 subgroup is <= 3.") %>%
                write.table(file = here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
            }
            
          }
        }
      }
  }
  
  if(!study_name %in% c("COVE")){
    if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
      if(!file.exists(here("output", Sys.getenv("TRIAL"), args[1], "drop_riskVars_with_fewer_0s_or_1s.csv"))){
        paste0("No binary input variable had number of cases in the variable = 1 or 0 subgroup <= 3") %>%
          write.table(file = here("output", Sys.getenv("TRIAL"), args[1], "drop_riskVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
      }
    } else {
      if(!file.exists(here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_fewer_0s_or_1s.csv"))){
        paste0("No binary input variable had number of cases in the variable = 1 or 0 subgroup <= 3") %>%
          write.table(file = here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_fewer_0s_or_1s.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
      }
    }
  }

  return(dat)
}



# remove any binary risk variables that have more than 5% values missing
# @param X dataframe containing all risk variables
# @param riskVars the vector of column names of risk variables
# @return a data frame upon removal of any binary risk variables that have more than 5% values missing
drop_riskVars_with_high_total_missing_values <- function(X, riskVars) {
  
  # delete the file drop_riskVars_with_high_total_missing_values.csv
  if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
    unlink(here("output", Sys.getenv("TRIAL"), args[1], "drop_riskVars_with_high_total_missing_values.csv"))
  }else{
    unlink(here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_high_total_missing_values.csv"))
  }
  
  covars_highNAvalues <- vector()
  for (i in 1:length(riskVars)) {
    total_NAs <- sum(is.na(X %>% pull(riskVars[i])))
    percent_NAs <- total_NAs / length(X %>% pull(riskVars[i]))

    if (percent_NAs > 0.05) {
      print(paste0("WARNING: ", riskVars[i], " variable has more than 5% values missing! This variable will be dropped from SuperLearner analysis."))
      # Also print to file
      if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
        paste0(riskVars[i], " variable has more than 5% values missing and was dropped from risk score analysis.") %>%
          write.table(file = here("output", Sys.getenv("TRIAL"), args[1], "drop_riskVars_with_high_total_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
      }else{
        paste0(riskVars[i], " variable has more than 5% values missing and was dropped from risk score analysis.") %>%
          write.table(file = here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_high_total_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
      }
      covars_highNAvalues <- riskVars[i]
    }
  }
  
  if(!file.exists(here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_high_total_missing_values.csv"))){
    if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
      paste0("No variables had more than 5% values missing.") %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), args[1], "drop_riskVars_with_high_total_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
    }else{
      paste0("No variables had more than 5% values missing.") %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), "drop_riskVars_with_high_total_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F) 
    }
  }
  
  if(length(covars_highNAvalues) == 0)
    return(X)
  else
    return(X %>% select(-all_of(covars_highNAvalues)))
}


# impute missing values in the risk variables using mice package
# @param X dataframe containing all risk variables
# @param riskVars the vector of column names of risk variables
# @return a data frame updated with imputations for missing values of risk variables
impute_missing_values <- function(X, riskVars) {
  # delete the file impute_missing_values.csv
  if(study_name %in% c("VAT08m", "VAT08")){
    unlink(here("output", Sys.getenv("TRIAL"), "bseroneg", "impute_missing_values.csv"))
  }else if(study_name == "PREVENT19"){
    unlink(here("output", Sys.getenv("TRIAL"), args[1], "impute_missing_values.csv"))
  }else{
    unlink(here("output", Sys.getenv("TRIAL"), "impute_missing_values.csv"))
  }
  
  covars <- vector()
  # First identify risk demographic variables having missing values
  for (i in 1:length(riskVars)) {
    total_NAs <- sum(is.na(X %>% pull(riskVars[i])))
    percent_NAs <- total_NAs / length(X %>% pull(riskVars[i]))

    if (percent_NAs > 0.05) {
      print(paste0("WARNING: ", riskVars[i], " variable has more than 5% values missing; all missing values will be imputed!"))
    }

    if (total_NAs > 0) {
      if (i == 1) {
        covars <- riskVars[i]
      } else {
        covars <- c(covars, riskVars[i])
      }
    }
  }

  if (length(covars) == 0) {
    print("No missing values to impute in any risk variables!")
    # Also print to file
    if(study_name %in% c("VAT08m", "VAT08")){
      paste("Imputing missing values in following variables: None") %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), "bseroneg", "impute_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
    }else if(study_name == "PREVENT19"){
      paste("Imputing missing values in following variables: None") %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), args[1], "impute_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
    }else{
      paste("Imputing missing values in following variables: None") %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), "impute_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
    }
  } else {
    print(paste("Imputing missing values in following variables: ", paste(as.character(covars), collapse = ", ")))
    # Also print to file
    if(study_name %in% c("VAT08m", "VAT08")){
      paste("Missing values were imputed for the following variables:", paste(as.character(covars), collapse = ", ")) %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), "bseroneg", "impute_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
    }else if(study_name == "PREVENT19"){
      paste("Missing values were imputed for the following variables:", paste(as.character(covars), collapse = ", ")) %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), args[1], "impute_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
    }else{
      if(!dir.exists(paste0("output/", Sys.getenv("TRIAL")))){
        dir.create(paste0("output/", Sys.getenv("TRIAL")))
      }
      paste("Missing values were imputed for the following variables:", paste(as.character(covars), collapse = ", ")) %>%
        write.table(file = here("output", Sys.getenv("TRIAL"), "impute_missing_values.csv"), sep=",", append = TRUE, row.names = F, col.names = F)
    }
    
    n.imp <- 1
    impVars <- X %>% select(all_of(covars))
    
    # deal with constant variables
    for (a in names(impVars)) {
      print(a)
      if (all(impVars[[a]]==min(impVars[[a]], na.rm=TRUE), na.rm=TRUE)) {
        X[[a]] = min(impVars[[a]], na.rm=TRUE)
        covars = covars[!covars %in% a]
      }
    }
    
    noimpVars <- X %>% select(-all_of(covars)) %>% colnames()

    init = mice(X, maxit=0) 
    meth = init$method
    predM = init$predictorMatrix 
    meth[c(noimpVars)] = ""
    
    # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
    X <- mice(X, method=meth, predictorMatrix=predM, m=n.imp,
                   printFlag = FALSE, seed=1, diagnostics = FALSE, remove_collinear = FALSE)
    X <- mice::complete(X, action = 1L)
    
  }
  return(X)
}


# choose learners for getting predictions
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @return dataframe sorted according to descending AUC for individual learners followed by SL and Discrete SL
choose_learners <- function(cvaucDAT) {
  cvaucDAT %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    arrange(-AUC) %>%
    .[1:2, ] %>%
    bind_rows(cvaucDAT %>%
      filter(Learner %in% c("SL", "Discrete SL")))
}

# get CV predictions
# @param cv_fit fit from running CV.Superlearner
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @return a dataframe with predictions for each subject for each Learner/Learner-Screen combination present in cvaucDAT
get_cv_predictions <- function(cv_fit, cvaucDAT) {
  top3 <- choose_learners(cvaucDAT)

  predict <- cv_fit[["library.predict"]] %>%
    as.data.frame() %>%
    bind_cols(cv_fit[["discreteSL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Discrete SL"))) %>%
    bind_cols(cv_fit[["SL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Super Learner"))) %>%
    bind_cols(cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y"))) %>%
    gather("algo", "pred", -Y) %>%
    filter(algo %in% c(top3$LearnerScreen))

  predict %>%
    left_join(top3 %>% select(Screen_fromRun, Learner, Screen, AUC, LearnerScreen), by = c("algo" = "LearnerScreen")) %>%
    mutate(
      learnerScreen = paste0(Learner, "_", Screen),
      learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
      AUCchar = format(round(AUC, 3), nsmall = 3),
      learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
      learnerScreen = reorder(learnerScreen, -AUC)
    )
}



# Plot ROC curves for SL, discrete.SL and topRanking learner-screen combinations
# @param predict dataframe returned by get_cv_predictions function
# @param cvaucDAT a dataframe containing Learner, Screen, and AUC
# @return ggplot object containing the ROC curves
plot_roc_curves <- function(predict, cvaucDAT) {
  top3 <- choose_learners(cvaucDAT)

  roc.obj <- predict %>%
    group_by(algo) %>%
    nest() %>%
    mutate(
      pred.obj = purrr::map(data, ~ ROCR::prediction(.x$pred, .x$Y)),
      perf.obj = purrr::map(pred.obj, ~ ROCR::performance(.x, "tpr", "fpr")),
      roc.dat = purrr::map(perf.obj, ~ tibble(
        xval = .x@x.values[[1]],
        yval = .x@y.values[[1]]
      ))
    )

  roc.obj %>%
    unnest(roc.dat) %>%
    select(algo, xval, yval) %>%
    ungroup() %>%
    left_join(top3 %>% select(Screen_fromRun, Learner, Screen, AUC, LearnerScreen), by = c("algo" = "LearnerScreen")) %>%
    mutate(
      learnerScreen = paste0(Learner, "_", Screen),
      learnerScreen = ifelse(Learner %in% c("SL", "Discrete SL"), algo, learnerScreen),
      AUCchar = format(round(AUC, 3), nsmall = 3),
      learnerScreen = paste0(learnerScreen, " (", AUCchar, ")"),
      learnerScreen = reorder(learnerScreen, -AUC)
    ) %>%
    ggplot(aes(x = xval, y = yval, col = learnerScreen)) +
    geom_step(lwd = 2) +
    theme(
      legend.position = "top",
      legend.direction = "vertical",
      legend.box = "horizontal",
      legend.title=element_text(size=20),
      legend.text=element_text(size=20),
      axis.ticks.length = unit(.35, "cm"),
      axis.text = element_text(size = 23),
      axis.title = element_text(size = 30)
    ) +
    labs(x = "Cross-Validated False Positive Rate", y = "Cross-Validated True Positive Rate", col = "Model (CV-AUC)") +
    geom_abline(intercept = 0, slope = 1)
}


# Plot predicted probability plots for SL, discrete.SL and topRanking learner-screen combinations
# @param pred dataframe returned by get_cv_predictions function
# @return ggplot object containing the predicted probability plots
plot_predicted_probabilities <- function(pred, day) {
  if(study_name == "COVE" | study_name == "MockCOVE"){
    cases = "Post Day 57 Cases"
  }
  if(study_name == "ENSEMBLE" | study_name == "MockENSEMBLE"){
    cases = "Post Day 29 Cases"
  }
  
  pred %>%
    mutate(Ychar = ifelse(Y == 0, "Non-Cases", paste0("Post Day ", day, " Cases"))) %>%
    ggplot(aes(x = Ychar, y = pred, color = Ychar)) +
    geom_jitter(width = 0.06, size = 3, shape = 21, fill = "white") +
    geom_violin(alpha = 0.05, color = "black", lwd=1.5) +
    geom_boxplot(alpha = 0.05, width = 0.15, color = "black", outlier.size = NA, outlier.shape = NA, lwd=1.5) +
    theme_bw() +
    #scale_color_manual(values = c("#56B4E9", "#E69F00")) +
    scale_color_manual(values = c("#00468B", "#8B0000")) +
    facet_wrap(vars(learnerScreen), ncol = 1) +
    labs(y = "CV estimated predicted probability of COVID-19 disease", x = "") +
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 25),
      axis.text = element_text(size = 23),
      axis.ticks.length = unit(.35, "cm"),
      axis.title.y = element_text(size = 30)
    )
}



# Create forest plot for demo run
# @param avgs dataframe containing Screen, Learner, AUC estimates and CIs as columns
# @return list of 2 ggplot objects: one containing forest plot and the other containing labels (Screen, Learner and CV-AUCs)
make_forest_plot_demo <- function(avgs) {
  lowestXTick <- floor(min(avgs$ci_ll) * 10) / 10
  top_learner_plot <- ggplot() +
    geom_pointrange(avgs %>% mutate(LearnerScreen = fct_reorder(LearnerScreen, AUC, .desc = F)), mapping = aes(x = LearnerScreen, y = AUC, ymin = ci_ll, ymax = ci_ul), size = 1, color = "blue", fill = "blue", shape = 20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, 1, 0.1), labels = seq(lowestXTick, 1, 0.1), limits = c(lowestXTick, 1)) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_blank(),
      plot.margin = unit(c(1, -0.15, 1, -0.15), "cm"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    )

  total_learnerScreen_combos <- length(avgs$LearnerScreen)

  avgs_withCoord <- avgs %>%
    select(Learner, Screen, AUCstr) %>%
    gather("columnVal", "strDisplay") %>%
    mutate(
      xcoord = case_when(
        columnVal == "Learner" ~ 1,
        columnVal == "Screen" ~ 1.5,
        columnVal == "AUCstr" ~ 2
      ),
      ycoord = rep(total_learnerScreen_combos:1, 3)
    )

  top_learner_nms_plot <- ggplot(avgs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust = 1, vjust = 0, size = 5) +
    xlim(0.7, 2) +
    theme(
      plot.margin = unit(c(2.2, -0.15, 3, -0.15), "cm"),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 2, color = "white"),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}



# Create forest plot for prod run
# @param avgs dataframe containing Screen, Learner, AUC estimates and CIs as columns
# @return list of 2 ggplot objects: one containing forest plot and the other containing labels (Screen, Learner and CV-AUCs)
make_forest_plot_prod <- function(avgs) {
  lowestXTick <- floor(min(avgs$ci_ll) * 10) / 10
  top_learner_plot <- ggplot() +
    geom_pointrange(avgs %>% mutate(LearnerScreen = fct_reorder(LearnerScreen, AUC, .desc = F)), mapping = aes(x = LearnerScreen, y = AUC, ymin = ci_ll, ymax = ci_ul), size = 3, color = "blue", fill = "blue", shape = 20) +
    coord_flip() +
    scale_y_continuous(breaks = seq(lowestXTick, 1, 0.1), labels = seq(lowestXTick, 1, 0.1), limits = c(lowestXTick, 1)) +
    theme_bw() +
    labs(y = "CV-AUC [95% CI]", x = "") +
    theme(
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 30),
      axis.title = element_text(size = 30),
      axis.ticks.length = unit(.35, "cm"),
      axis.text.y = element_blank(),
      plot.margin = unit(c(0.8, -0.15, 0.2, -0.15), "cm"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    )

  total_learnerScreen_combos <- length(avgs$LearnerScreen)

  avgs_withCoord <- avgs %>%
    select(Learner, Screen, AUCstr) %>%
    gather("columnVal", "strDisplay") %>%
    mutate(
      xcoord = case_when(
        columnVal == "Learner" ~ 1,
        columnVal == "Screen" ~ 1.5,
        columnVal == "AUCstr" ~ 2
      ),
      ycoord = rep(total_learnerScreen_combos:1, 3)
    )

  top_learner_nms_plot <- ggplot(avgs_withCoord, aes(x = xcoord, y = ycoord, label = strDisplay)) +
    geom_text(hjust = 1, vjust = 0, size = 10) +
    xlim(0.7, 2) +
    theme(
      plot.margin = unit(c(1.0, -0.15, 1.7, -0.15), "cm"),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 2, color = "white"),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}



# This function scales each column in a numeric dataset and also stores the 
# scaling parameters as attributes to be accessed later
get_scaleParams_scaledData <- function(inputdat, groupStr){
  # Calculate mean and standard deviation for each column
  means <- colMeans(inputdat, na.rm=T)
  stds <- apply(inputdat, 2, sd, na.rm=T)
  
  # Scale each column using mean and standard deviation
  scaled_inputdat <- scale(inputdat, center = means, scale = stds)
  
  # Scaling params can be accessed using:
  # attr(scaled_inputdat ,"scaled:center")
  # attr(scaled_inputdat, "scaled:scale")
  # These params are accessed for scaling the extra 15 subjects from COVEBoost!
  scaled_inputdat
}






#' Check if a DataFrame is Standardized
#'
#' This function checks whether all columns in a dataframe are standardized,
#' meaning they have a mean close to 0 and a standard deviation close to 1.
#'
#' @param df A numeric dataframe to check.
#' @param tol A numeric tolerance level for checking if mean is approximately 0 and standard deviation is approximately 1. Default is `1e-6`.
#'
#' @return A logical value: `TRUE` if all columns are standardized, `FALSE` otherwise.
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   A = rnorm(1000, mean = 0, sd = 1),
#'   B = rnorm(1000, mean = 0, sd = 1),
#'   C = rnorm(1000, mean = 0, sd = 1)
#' )
#' check_standardized(df)  # Should return TRUE
#'
#' df_non_standardized <- data.frame(
#'   A = rnorm(1000, mean = 5, sd = 2),
#'   B = rnorm(1000, mean = 0, sd = 1)
#' )
#' check_standardized(df_non_standardized)  # Should return FALSE
#'
#' @export
check_standardized <- function(df, tol = 1e-6) {
  means <- colMeans(df)
  stds <- apply(df, 2, sd)
  
  mean_check <- all(abs(means) < tol)
  std_check <- all(abs(stds - 1) < tol)
  
  return(mean_check && std_check)
}
