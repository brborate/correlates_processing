# Check if only change in input dataset is marker data. 
# If so, simply pull the earlier risk scores and add to the new processed dataset!
if(study_name == "ENSEMBLE"){
    old_processed <- read.csv(here::here("..", "..", "previous_data_clean", paste0(attr(config, "config"), "_data_processed.csv"))) %>%
      select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), EthnicityHispanic,EthnicityNotreported, EthnicityUnknown,
             Black, Asian, NatAmer, PacIsl, Multiracial, Notreported, Unknown,
             URMforsubcohortsampling, HighRiskInd, HIVinfection, 
             Sex, Age, BMI, Country, Region, CalendarDateEnrollment) 
    
    new_processed <- inputFile %>% 
      select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), EthnicityHispanic,EthnicityNotreported, EthnicityUnknown,
             Black, Asian, NatAmer, PacIsl, Multiracial, Notreported, Unknown,
             URMforsubcohortsampling, HighRiskInd, HIVinfection, 
             Sex, Age, BMI, Country, Region, CalendarDateEnrollment) 
    
    if(all.equal(old_processed, new_processed)){
      dat_with_riskscore <- inputFile %>% left_join(read.csv(here::here("..", "..", "previous_data_clean", paste0(str_remove(paste0(attr(config, "config"), "_data_processed.csv"), ".csv"), "_with_riskscore.csv"))) %>%
                                                      select(Ptid, risk_score, standardized_risk_score), 
                                                    by = "Ptid") %>%
        filter(Riskscorecohortflag == 1) %>% 
        drop_na(all_of(endpoint))
      data_name_amended <- paste0(str_remove(paste0(attr(config, "config"), "_data_processed.csv"), ".csv"), "_with_riskscore")
      
      # Ensure all baseline negative and PP subjects have a risk score!
      if(assertthat::assert_that(
        all(!is.na(dat_with_riskscore %>% .$risk_score)), 
        msg = "Some baseline negative and PP subjects have NA values in risk score!"
      )){
        write_csv(dat_with_riskscore,
                  here("..", "data_clean", paste0(data_name_amended, ".csv")))
        message("No change in input data. Superlearner will not be run. Risk scores from earlier run appended to data_processed.csv!")
      }
    }else{
      message("There is change in input data. Superlearner needs to be run and new risk scores generated!")
    }
}


