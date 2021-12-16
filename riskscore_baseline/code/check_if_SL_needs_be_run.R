# Check if only change in input dataset is marker data. 
# If so, simply pull the earlier risk scores and add to the new processed dataset!
generate_new_riskscores <- function(){
  source(here("code", "clean_output_and_figs_dirs.R"))
  source(here("code", "run_cvsl_riskscore.R"))
  source(here("code", "createRDAfiles_fromSLobjects.R"))
  source(here("code", "tables_figures.R"))
  source(here("code", "constructSL_predict_on_vaccine.R"))
  source(here("code", "get_SLweights_Modelpredictors.R"))
  source(here("code", "append_risk_score_to_data.R"))
}


if(study_name %in% c("ENSEMBLE", "MockENSEMBLE")){
  if(file.exists(paste0("output/", attr(config, "config"), "_inputFile_with_riskscore.RData"))){
    load(paste0("output/", attr(config, "config"), "_inputFile_with_riskscore.RData"))
    old_processed <- inputFile_with_riskscore %>% 
      select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), EthnicityHispanic,EthnicityNotreported, EthnicityUnknown,
             Black, Asian, NatAmer, PacIsl, Multiracial, Notreported, Unknown,
             URMforsubcohortsampling, HighRiskInd, HIVinfection, 
             Sex, Age, BMI, Country, Region, CalendarDateEnrollment) 
    
    new_processed <- inputFile %>% 
      select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), EthnicityHispanic,EthnicityNotreported, EthnicityUnknown,
             Black, Asian, NatAmer, PacIsl, Multiracial, Notreported, Unknown,
             URMforsubcohortsampling, HighRiskInd, HIVinfection, 
             Sex, Age, BMI, Country, Region, CalendarDateEnrollment) 
    
    if(all.equal(old_processed, new_processed) == TRUE){
       message("No change in input data. Superlearner will not be run. Risk scores from earlier run appended to raw data!")
    }else{
      message("There is change in input data. Superlearner needs to be run and new risk scores generated!")
      generate_new_riskscores()
    }
  }
  if(!file.exists(paste0("output/", attr(config, "config"), "_inputFile_with_riskscore.RData"))){
    message(paste0("riskscore_baseline/", paste0("output/", attr(config, "config"), "_inputFile_with_riskscore.RData"), " does not exist. Superlearner needs to be run and new risk scores generated!"))
    generate_new_riskscores()
    }
}



if(study_name %in% c("COVE", "MockCOVE")){
  if(file.exists(paste0("output/", attr(config, "config"), "_inputFile_with_riskscore.RData"))){
    load(paste0("output/", attr(config, "config"), "_inputFile_with_riskscore.RData"))
    old_processed <- inputFile_with_riskscore %>% 
      select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), all_of(risk_vars)) 
    
    new_processed <- inputFile %>% 
      select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), all_of(risk_vars)) 
    
    if(all.equal(old_processed, new_processed) == TRUE){
      message("No change in input data. Superlearner will not be run. Risk scores from earlier run appended to raw data!")
    }else{
      message("There is change in input data. Superlearner needs to be run and new risk scores generated!")
      generate_new_riskscores()
    }
  }
  if(!file.exists(paste0("output/", attr(config, "config"), "_inputFile_with_riskscore.RData"))){
    message(paste0("riskscore_baseline/", paste0("output/", attr(config, "config"), "_inputFile_with_riskscore.RData"), " does not exist. Superlearner needs to be run and new risk scores generated!"))
    generate_new_riskscores()
  }
}


# 
# # Append risk scores from janssen_pooled_real to janssen_pooled_realADCP
# if(study_name == "ENSEMBLE" & attr(config, "config") == "janssen_pooled_realADCP"){
#   real_processed <- read.csv(here::here("..", "data_clean", "janssen_pooled_real_data_processed_with_riskscore.csv")) %>%
#     select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), EthnicityHispanic,EthnicityNotreported, EthnicityUnknown,
#            Black, Asian, NatAmer, PacIsl, Multiracial, Notreported, Unknown,
#            URMforsubcohortsampling, HighRiskInd, HIVinfection, 
#            Sex, Age, BMI, Country, Region, CalendarDateEnrollment) 
#   
#   realADCP_processed <- inputFile %>% 
#     filter(Riskscorecohortflag == 1) %>%
#     select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), EthnicityHispanic,EthnicityNotreported, EthnicityUnknown,
#            Black, Asian, NatAmer, PacIsl, Multiracial, Notreported, Unknown,
#            URMforsubcohortsampling, HighRiskInd, HIVinfection, 
#            Sex, Age, BMI, Country, Region, CalendarDateEnrollment) 
#   
#   if(all.equal(real_processed, realADCP_processed)){
#     dat_with_riskscore <- inputFile %>% left_join(read.csv(here::here("..", "data_clean", "janssen_pooled_real_data_processed_with_riskscore.csv")) %>%
#                                                     select(Ptid, risk_score, standardized_risk_score), 
#                                                   by = "Ptid") %>%
#       filter(Riskscorecohortflag == 1) %>% 
#       drop_na(all_of(endpoint))
#     data_name_amended <- paste0(str_remove(paste0(attr(config, "config"), "_data_processed.csv"), ".csv"), "_with_riskscore")
#     
#     # Ensure all baseline negative and PP subjects have a risk score!
#     if(assertthat::assert_that(
#       all(!is.na(dat_with_riskscore %>% .$risk_score)), 
#       msg = "Some baseline negative and PP subjects have NA values in risk score!"
#     )){
#       write_csv(dat_with_riskscore,
#                 here("..", "data_clean", paste0(data_name_amended, ".csv")))
#       message("No change in input data. Superlearner will not be run. Risk scores from earlier run appended to data_processed.csv!")
#     }
#   }else{
#     message("There is change in input data. Superlearner needs to be run and new risk scores generated!")
#   }
# }
# 
# 
