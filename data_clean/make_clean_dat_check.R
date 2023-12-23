#Sys.setenv(TRIAL = "vat08_combined")

source(here::here("_common.R"))


# same as deploy script
if (attr(config, "config") %in% c("prevent19", "moderna_real", "moderna_boost", "janssen_partA_VL", "vat08_combined", "covail")) {
  data_name_amended <- c(paste0(attr(config, "config"), "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv"))
  
} else if(attr(config, "config") %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA")) {
  data_name_amended <- c( paste0(attr(config, "config"), "_data_processed_with_riskscore"), 
                          paste0(attr(config, "config"), "senior_data_processed_with_riskscore"),
                          paste0(attr(config, "config"), "nonsenior_data_processed_with_riskscore"))
  
} else {
  data_name_amended <- c(paste0(attr(config, "config"), "_data_processed_with_riskscore"))
}

# load data and rename first column (ID)
dat_clean <- read.csv(here::here("data_clean/csv", data_name_amended)) 


# leave comments below for checks implemented in make_dat_proc.R
## missing markers imputed properly in each stratum
## imputed values of missing markers merged properly for all individuals in the two phase sample

## presence of values lower than LLOD / 2
assays_plusN = c(assays, "bindN")

failed_llod_check <- NULL
for (a in assays_plusN) {
  for (t in c("B", paste0("Day", config$timepoints))) {
    pass <- all(dat_clean[[paste0(t,a)]] >= ifelse(llods[a]/2>=1, .999, 1.001) * log10(llods[a] / 2), na.rm = TRUE) #.999 is a necessary fudge factor
    if(!pass){
        failed_llod_check <- c(failed_llod_check, paste0(t,a))
    }
  }
}

if(length(failed_llod_check) > 1){
    stop(paste0("Values of assays less than LLOD / 2 for: ", 
                paste(failed_llod_check, sep = ", "), "\n"))
}



## missing values in variables that should have no missing values
variables_with_no_missing <- 
    c(
      c(outer(c("ph1.D", "ph2.D", "EarlyendpointD", "TwophasesampIndD", "EarlyinfectionD"), config$timepoints, paste0)),
      "age.geq.65", "MinorityInd",
      "ph1.immuno",
      "ph2.immuno"
      )

failed_variables_missing <- failed_variables_01 <- NULL
for(variable in variables_with_no_missing){
    pass <- all(!is.na(dat_clean[[variable]]))
    if(!pass){
        failed_variables_missing <- c(failed_variables_missing, variable)
    }
}

if(length(failed_variables_missing) > 0){
    stop(paste0("Unexpected missingness in: ", paste(failed_variables_missing, collapse = ", ")))   
}


if (study_name=="VAT08") {
  # verify event indicator and time variables for Omicron using first and second event variables
  
  # For the 7 participants with first event a known-lineage Omicron case, this first event is counted as an endpoint and the second event is ignored;
  # For the 2 participants with first event a missing-lineage case in 2021 and second event known-lineage Omicron, the first event is ignored and the second event is counted as an endpoint; and
  # For the 9 participants with first event a missing-lineage case in 2022 and second event known-lineage Omicron, the first event is counted as an endpoint (with hotdeck imputation employed as is generally done) and the second event is ignored.
  
  with(dat_clean, table(EventTypeFirstInfection, EventTypeSecondInfection, useNA = "ifany"))

  with(subset(dat_clean, !is.na(EventIndFirstInfectionD22)), table(EventTypeFirstInfection, EventTypeSecondInfection, useNA = "ifany"))
  with(subset(dat_clean, !is.na(EventIndFirstInfectionD43)), table(EventTypeFirstInfection, EventTypeSecondInfection, useNA = "ifany"))
  
  # EventTypeFirstInfection and seq1.variant have the same info
  with(dat_clean, table(EventTypeFirstInfection, seq1.variant, useNA = "ifany"))
  
  with(dat_clean, table(seq1.variant.hotdeck1, seq1.variant, useNA = "ifany"))
  
  with(dat_clean, table(EventIndKnownLineageOmicronD22, EventIndOmicronD22M12hotdeck1, useNA = "ifany"))
  
  with(dat_clean, table(EventIndKnownLineageOmicronD22, EventTypeFirstInfection, useNA = "ifany"))
  
  with(dat_clean, table(EventIndKnownLineageOmicronD1, EventTypeFirstInfection, useNA = "ifany"))
  
  with(dat_clean, table(EventIndFirstInfectionD22, EventIndFirstInfectionD1, useNA = "ifany"))
  with(dat_clean, table(EventIndFirstInfectionD22, EventIndFirstInfectionD43, useNA = "ifany"))
}