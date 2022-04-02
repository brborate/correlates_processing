inputMod <- inputMod %>%
  drop_na(all_of(endpoint)) 

# Create table of cases in both arms (prior to applying Riskscorecohortflag filter)
  tab <- inputMod %>%
    drop_na(Ptid, Trt, all_of(endpoint)) %>%
    mutate(Trt = ifelse(Trt == 0, "Placebo", "Vaccine")) 
  
  table(tab$Trt, tab %>% pull(endpoint)) %>%
    write.csv(file = here("output", Sys.getenv("TRIAL"), "cases_prior_to_applying_Riskscorecohortflag.csv"))
  rm(tab)
  
  
  dat.ph1 <- inputMod %>% filter(Riskscorecohortflag == 1 & Trt == 0)
  
  dat.ph1 <- dat.ph1 %>%
    # Keep only variables to be included in risk score analyses
    select(Ptid, Trt, all_of(endpoint), all_of(risk_vars)) %>%
    # Drop any observation with NA values in Ptid, Trt, or endpoint!
    drop_na(Ptid, Trt, all_of(endpoint))
  
  # Create table of cases in both arms (prior to Risk score analyses)
  tab <- inputMod %>%
    filter(Riskscorecohortflag == 1) %>%
    drop_na(Ptid, Trt, all_of(endpoint)) %>%
    mutate(Trt = ifelse(Trt == 0, "Placebo", "Vaccine")) 
  
  table(tab$Trt, tab %>% pull(endpoint)) %>%
    write.csv(file = here("output", Sys.getenv("TRIAL"), "cases_prior_riskScoreAnalysis.csv"))
  rm(tab)
  
  # Derive maxVar: the maximum number of variables that will be allowed by SL screens in the models.
  np <- sum(dat.ph1 %>% select(matches(endpoint)))
  maxVar <- max(20, floor(np / 20))
  all_risk_vars <- risk_vars
  
  # Remove a variable if the number of cases in the variable = 1 subgroup is <= 3 or the number of cases in the variable = 0 subgroup is <= 3
  dat.ph1 <- drop_riskVars_with_fewer_0s_or_1s(dat = dat.ph1, 
                                               risk_vars = risk_vars,
                                               np = np)
  
  # Update risk_vars
  risk_vars <- dat.ph1 %>%
    select(-Ptid, -Trt, -all_of(endpoint)) %>%
    colnames()
  
  # Remove any risk_vars with more than 5% missing values. Impute the missing
  # values for other risk variables using mice package!
  dat.ph1 <- drop_riskVars_with_high_total_missing_values(dat.ph1, risk_vars)
  
  # Update risk_vars
  risk_vars <- dat.ph1 %>%
    select(-Ptid, -Trt, -all_of(endpoint)) %>%
    colnames()
  
  X_covars2adjust <- dat.ph1 %>%
    select(all_of(risk_vars))
  
  # Save ptids to merge with predictions later
  risk_placebo_ptids <- dat.ph1 %>% select(Ptid, all_of(endpoint))
  
  # Impute missing values in any variable included in risk_vars using the mice package!
  print("Make sure data is clean before conducting imputations!")
  X_covars2adjust <- impute_missing_values(X_covars2adjust, risk_vars)
  
  # # Check for missing values before and after imputation
  # sapply(X_covars2adjust, function(x) sum(is.na(x)))
  
  # Scale X_covars2adjust to have mean 0, sd 1 for all vars
  for (a in colnames(X_covars2adjust)) {
    X_covars2adjust[[a]] <- scale(X_covars2adjust[[a]],
                                  center = mean(X_covars2adjust[[a]], na.rm = T),
                                  scale = sd(X_covars2adjust[[a]], na.rm = T)
    )
  }
  
  X_riskVars <- X_covars2adjust
  Y <- dat.ph1 %>% pull(endpoint)
  
  # set up outer folds for cv variable importance; do stratified sampling
  V_outer <- 5
  cvControlVar = list(V = V_outer, stratifyCV = TRUE)
  cvControl_quote = quote(list(V = V_outer, stratifyCV = TRUE))
  
  if(study_name != "PREVENT19"){
    #if (np < round(length(Y)*2/100)) {  # Update rule: Do Leave-One-Out CV if number of cases in placebo are less than 2%!
    #V_inner <- 20
    if (np <= 30) {
      V_inner <- length(Y) - 1
    } else {
      V_inner <- 5
    }

    if(V_inner == length(Y) - 1){
      V_inner_quote <- paste0("length(Y) - 1 = ", length(Y) - 1)
    }
    innerCvControlVar = list(list(V = V_inner))
    innerCvControl_quote = quote(list(list(V = V_inner)))
    #}
  }

  if(study_name == "PREVENT19"){
    #if (np < round(length(Y)*2/100)) {  # Update rule: Do Leave-One-Out CV if number of cases in placebo are less than 2%!
      #V_inner <- 5
      V_inner <- length(Y) - 1
      if(V_inner == length(Y) - 1){
        V_inner_quote <- paste0("length(Y) - 1 = ", length(Y) - 1)
        innerCvControlVar = list(list(V = V_inner))
        innerCvControl_quote = quote(list(list(V = V_inner)))
      }else{
        innerCvControlVar = list(list(V = V_inner, stratifyCV = TRUE))
        innerCvControl_quote = quote(list(list(V = V_inner, stratifyCV = TRUE)))
      }
    #}
  }
  
  ## solve cores issue
  #blas_get_num_procs()
  blas_set_num_threads(1)
  #print(blas_get_num_procs())
  stopifnot(blas_get_num_procs() == 1)
  
  # CV.SL inputs
  familyVar = "binomial"
  methodVar = "method.CC_nloglik"
  scaleVar = "identity"
  cvsl_args <- data.frame(matrix(ncol = 2, nrow = 8)) %>%
    rename(Argument = X1,
           Value = X2) %>%
    mutate(Argument = as.character(c("Cases/Total Subjects in placebo group (%)", "family",
                        "method", "scale", "V_outer", "cvControl (outer CV control)",
                        "V_inner", "innerCvControl")),
           Value = as.character(c(paste0(np, "/", length(Y), " (", round(np*100/length(Y), 2), "%)"), familyVar,
                     methodVar, scaleVar, V_outer, cvControl_quote,
                     V_inner, innerCvControl_quote)))

  if(V_inner == length(Y) - 1){
    cvsl_args = cvsl_args %>% mutate(Value = ifelse(Argument == "V_inner", V_inner_quote, Value))
  }

  cvsl_args %>% write.csv(paste0("output/", Sys.getenv("TRIAL"), "/", "cvsl_args.csv"))

  # run super learner ensemble
  fits <- run_cv_sl_once(
    seed = 20210216,
    Y = Y,
    X_mat = X_riskVars,
    family = familyVar,
    method = methodVar,
    scale = scaleVar,
    sl_lib = SL_library,
    cvControl = cvControlVar,
    innerCvControl = innerCvControlVar,
    vimp = FALSE
  )
  
  cvaucs <- list()
  cvaucs[[1]] <- fits$cvaucs
  cvfits <- list()
  cvfits[[1]] <- fits$cvfits
  
  saveRDS(cvaucs, here("output", Sys.getenv("TRIAL"), "cvsl_riskscore_cvaucs.rds"))
  save(cvfits, file = here("output", Sys.getenv("TRIAL"), "cvsl_riskscore_cvfits.rda"))
  save(risk_placebo_ptids, file = here("output", Sys.getenv("TRIAL"), "risk_placebo_ptids.rda"))
  save(run_prod, Y, X_riskVars, weights, inputMod, risk_vars, all_risk_vars, endpoint, maxVar,
       V_outer, V_inner, familyVar, methodVar, scaleVar, studyName_for_report, risk_timepoint, 
       cvControlVar, inputfileName, 
       file = here("output", Sys.getenv("TRIAL"), "objects_for_running_SL.rda"))

