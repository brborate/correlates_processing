#Sys.setenv(TRIAL = "moderna_mock")
#-----------------------------------------------
renv::activate(here::here())
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("_common.R"))
#-----------------------------------------------
library(here)


if (endsWith(attr(config, "config"), "mock")) {
    if(attr(config, "config")=="moderna_mock") {
      path_to_data <- here(".", paste0("data_raw/moderna/", mapped_data))
    } else {
        # janssen pooled or regions
      path_to_data <- here(".", paste0("data_raw/janssen/", mapped_data))
    } 
} else {
    path_to_data <- mapped_data
}
print(path_to_data)
if (!file.exists(path_to_data)) stop ("make raw dat check: dataset not available ===========================================")
dat_proc <- read.csv(path_to_data)

# load data and rename first column (ID)
colnames(dat_proc)[1] <- "Ptid"


## missing values in variables that should have no missing values
## binary variables only take values 0/1
bin_variables_with_no_missing <-
    c("Trt", "Bserostatus", #"Age", # age is not binary
      "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl", "Multiracial",
      "Other", "Notreported", "Unknown",
      "SubcohortInd",
      "EventIndPrimaryD1", 
      if(study_name_code=="ENSEMBLE") c("HIVinfection"))
failed_variables_missing <- failed_variables_01 <- NULL
for(variable in bin_variables_with_no_missing){
    pass <- all(!is.na(dat_proc[[variable]]))
    if(!pass){
        failed_variables_missing <- c(failed_variables_missing, variable)
    }
    pass <- all(dat_proc[[variable]] %in% c(0,1))
    if(!pass){
        failed_variables_01 <- c(failed_variables_01, variable)
    }
}
if(length(failed_variables_missing) > 0){
    stop(paste0("Unexpected missingness in: ", paste(failed_variables_missing,
                                                     collapse = ", ")))
}
if(length(failed_variables_01) > 0){
    stop(paste0("Unexpected values in: ", paste(failed_variables_01,
                                                collapse = ", ")))
}



## missing values in variables that should have no missing values
## binary variables only take values 0/1
variables_with_no_missing <-
    c("EventTimePrimaryD1", 
      if(study_name_code=="ENSEMBLE") c("HIVinfection"))
failed_variables_missing <- NULL
for(variable in variables_with_no_missing){
    pass <- all(!is.na(dat_proc[[variable]]))
    if(!pass){
        failed_variables_missing <- c(failed_variables_missing, variable)
    }
}
if(length(failed_variables_missing) > 0){
    stop(paste0("Unexpected missingness in: ", paste(failed_variables_missing,
                                                     collapse = ", ")))
}





# check failure times for sanity
## EventIndPrimaryDtp2==1 implies EventIndPrimaryDtp1==1
if(two_marker_timepoints) {
    pass <- with(dat_proc, {
        idx = get("EventIndPrimaryD"%.%timepoints[2]) == 1
        idx[is.na(idx)]=FALSE
        all(get("EventIndPrimaryD"%.%timepoints[1])[idx] == 1)
    })
    if(is.na(pass) | !pass){
        stop(paste0("Some individuals with qualifying events for Day tp1 analysis are labeled ",
                    "as having no event for the Day tp2 analysis."))
    }
    
    ## cases that qualify for both events have shorter follow up for Day tp2 analysis
    pass <- with(dat_proc, {
        idx <- get("EventIndPrimaryD"%.%timepoints[2]) == 1 & get("EventIndPrimaryD"%.%timepoints[1]) == 1
        idx[is.na(idx)]=FALSE
        all(get("EventTimePrimaryD"%.%timepoints[2])[idx] < get("EventTimePrimaryD"%.%timepoints[1])[idx])
    })
    if(is.na(pass) | !pass){
        stop(paste0("Amongst individuals who have events that qualify for both Day tp1 and Day tp2 ",
                    "some follow up times are *longer* for Day tp2 than for Day tp1"))
    }

    ## consistency between event time variables for the cases
    if (study_name != "MockCOVE") {
        pass <- with(dat_proc, {
            tmp = get("NumberdaysD1toD"%.%timepoints[2]) - get("NumberdaysD1toD"%.%timepoints[1]) == get("EventTimePrimaryD"%.%timepoints[1]) - get("EventTimePrimaryD"%.%timepoints[2])
            all(tmp | is.na(tmp))
        })
        if(!pass){
            stop(paste0("NumberdaysD1toDtp2 - NumberdaysD1toDtp1 == EventTimePrimaryDtp1 - EventTimePrimaryDtp2 fails for some rows"))
        }
    }



}
