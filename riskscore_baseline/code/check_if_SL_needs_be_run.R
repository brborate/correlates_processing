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


if(file.exists(paste0("output/", Sys.getenv("TRIAL"), "/", attr(config, "config"), "_inputFile_with_riskscore.RData")) |
   (study_name == "ENSEMBLE" & file.exists(paste0("output/janssen_pooled_realbAb/janssen_pooled_realbAb_inputFile_with_riskscore.RData")))){
  
  if(study_name == "ENSEMBLE"){
    load("output/janssen_pooled_realbAb/janssen_pooled_realbAb_inputFile_with_riskscore.RData")
  }else{
    load(paste0("output/", Sys.getenv("TRIAL"), "/", attr(config, "config"), "_inputFile_with_riskscore.RData"))
  } 
    
  old_processed <- inputFile_with_riskscore %>%
    select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), all_of(original_risk_vars), risk_score, standardized_risk_score)

  new_processed <- inputFile %>%
    select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), all_of(original_risk_vars))

  if(all.equal(old_processed %>% select(-c(risk_score, standardized_risk_score)), new_processed) == TRUE){
    message("No change in input data. Superlearner will not be run. Risk scores from earlier run will be appended to raw data!")
    inputFile_with_riskscore <- left_join(inputFile,
                                          old_processed %>%
                                            select(Ptid, risk_score, standardized_risk_score), by = "Ptid")

    save(inputFile_with_riskscore, file = paste0("output/", Sys.getenv("TRIAL"), "/", attr(config, "config"), "_inputFile_with_riskscore.RData"))
    }else{
      message("There is change in input data. Superlearner needs to be run and new risk scores generated!")
      generate_new_riskscores()
    }
}else{
  message(paste0("riskscore_baseline/", paste0("output/", Sys.getenv("TRIAL"), "/", attr(config, "config"), "_inputFile_with_riskscore.RData"), " does not exist. Superlearner needs to be run and new risk scores generated!"))
  generate_new_riskscores()
}

# Perform sanity check if mock data!
source(here("code", "unit_test.R"))
