# Sys.setenv(TRIAL = "janssen_pooled_real")  
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries, cleaned data, and risk score estimates
library(here)
library(tidyverse)
data_name_amended <- paste0(str_remove(paste0(attr(config, "config"), "_data_processed.csv"), ".csv"), "_with_riskscore")
# Move current copy of risk score dataset to archive

# Copy current copy of risk score dataset to archive 
# Remove current copy of risk score dataset from adata
# Copy new copy of risk score dataset to adata
if(study_name_code == "ENSEMBLE"){
  file.copy(from =  here("..", "..", "p3003", "analysis", "correlates", "Part_A_Blinded_Phase_Data", "adata", paste0(data_name_amended, ".csv")),
            to =  here("..", "..", "p3003", "analysis", "correlates", "Part_A_Blinded_Phase_Data", "adata", "archive", 
                       paste0(data_name_amended, 
                              file.info(here("..", "..", "previous_data_clean", paste0(data_name_amended, ".csv")))$ctime, 
                              ".csv")),
            copy.date = TRUE)
  
  file.remove(from =  here("..", "..", "p3003", "analysis", "correlates", "Part_A_Blinded_Phase_Data", "adata", paste0(data_name_amended, ".csv")))
  
  file.copy(from =  here("..", "data_clean", paste0(data_name_amended, ".csv")),
            to = here("..", "..", "p3003", "analysis", "correlates", "Part_A_Blinded_Phase_Data", "adata", paste0(data_name_amended, ".csv")),
            copy.date = TRUE)
}



