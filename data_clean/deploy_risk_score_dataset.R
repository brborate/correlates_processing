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
data_name_amended <- c(paste0(attr(config, "config"), "_data_processed_with_riskscore"), 
                       paste0(attr(config, "config"), "_data_processed_for_immunogenicity"))

# Request reason for adata update from deployer!
update_reason <- paste0(Sys.Date(), " ", readline(prompt = "Enter reason for updating adata (this text will be added to adata/README change log): ")) 

# Copy current deployed copy of risk score dataset in adata to archive 
# Remove current deployed copy of risk score dataset from adata
# Copy new copy of risk score dataset from data_clean to adata
for(j in 1:length(data_name_amended)){
  if(file.exists(paste0("data_clean/", data_name_amended[j], ".csv"))){
    if(file.exists(paste0(strsplit(mapped_data, "mapping")[[1]][1],
                          "correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended[j], ".csv"))){
      
      if(!file.exists(paste0(strsplit(mapped_data, "mapping")[[1]][1],
                             "correlates/Part_A_Blinded_Phase_Data/adata/archive"))){
        
        dir.create(paste0(strsplit(mapped_data, "mapping")[[1]][1],
                          "correlates/Part_A_Blinded_Phase_Data/adata/archive"))
      }
      
      file.copy(from = paste0(strsplit(mapped_data, "mapping")[[1]][1], "correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended[j], ".csv"),
                to =  paste0(strsplit(mapped_data, "mapping")[[1]][1], "correlates/Part_A_Blinded_Phase_Data/adata/archive/", data_name_amended[j], "_",
                             str_replace(str_replace_all(file.info(paste0(strsplit(mapped_data, "mapping")[[1]][1], "correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended[j], ".csv"))$mtime,
                                                         ":",
                                                         "."), " ", " time "), 
                             ".csv"),
                copy.date = TRUE)
      
      file.remove(from = paste0(strsplit(mapped_data, "mapping")[[1]][1], "correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended[j], ".csv"))
    }
    
    file.copy(from = paste0("data_clean/", data_name_amended[j], ".csv"),
              to = paste0(strsplit(mapped_data, "mapping")[[1]][1], "correlates/Part_A_Blinded_Phase_Data/adata/", data_name_amended[j], ".csv"),
              copy.date = TRUE)
    
    # Add reason for adata update to README file!
    write(update_reason, file = paste0(strsplit(mapped_data, "mapping")[[1]][1], "correlates/Part_A_Blinded_Phase_Data/adata/readme.txt"), append=TRUE)
    
  }else{
    print(paste0("data_clean/", data_name_amended[j], ".csv not found!"))
  }
}
