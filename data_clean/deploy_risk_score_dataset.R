# Sys.setenv(TRIAL = "janssen_pooled_realbAb")  
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here())
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("_common.R"))
#-----------------------------------------------

# load required libraries, cleaned data, and risk score estimates
library(here)
library(tidyverse)
data_name_amended <- paste0(str_remove(paste0(attr(config, "config"), "_data_processed.csv"), ".csv"), "_with_riskscore")

# Copy current deployed copy of risk score dataset in adata to archive 
# Remove current deployed copy of risk score dataset from adata
# Copy new copy of risk score dataset from data_clean to adata
if(study_name_code == "ENSEMBLE"){
  if(file.exists(paste0("/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended, ".csv"))){
    file.copy(from = paste0("/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended, ".csv"),
              to =  paste0("/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/archive/", data_name_amended, "_",
                           str_replace(str_replace_all(file.info(paste0("/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended, ".csv"))$mtime,
                                                       ":",
                                                       "."), " ", " time "), 
                           ".csv"),
              copy.date = TRUE)
    
    file.remove(from = paste0("/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended, ".csv"))
  }
    
  file.copy(from = paste0("data_clean/", data_name_amended, ".csv"),
            to = paste0("/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended, ".csv"),
            copy.date = TRUE)
}



