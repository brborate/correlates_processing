# Sys.setenv(TRIAL = "janssen_pooled_realbAb")  
# Sys.setenv(TRIAL = "moderna_real")  
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here())
#-----------------------------------------------

config <- config::get(config = Sys.getenv("TRIAL"))

library(stringr)

if (attr(config, "config")=="moderna_boost") {
  data_name_amended <- c(paste0(attr(config, "config"), "_data_processed_", format(Sys.Date(), "%Y%m%d")))
  # Request reason for adata update from deployer!
  cat("Enter reason for updating adata without quotes (this text will be added to adata/README change log): ")
  args <- readLines(con = "stdin", n = 1)
  update_reason <- paste0(Sys.Date(), " ", args[[1]]) 
  deploy_path =      "/trials/covpn/p3001/analysis/correlates/Part_C_Unblinded_Phase_Data/adata/"
  
} else if (attr(config, "config")=="janssen_partA_VL") {
  data_name_amended <- c(paste0(attr(config, "config"), "_data_processed_", format(Sys.Date(), "%Y%m%d")))
  # Request reason for adata update from deployer!
  cat("Enter reason for updating adata without quotes (this text will be added to adata/README change log): ")
  args <- readLines(con = "stdin", n = 1)
  update_reason <- paste0(Sys.Date(), " ", args[[1]]) 
  deploy_path =      "/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
  
} else {
  # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
  if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
  source(here::here("_common.R"))
  
  # load required libraries, cleaned data, and risk score estimates
  library(here)
  library(tidyverse)
  if(attr(config, "config") %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA")) {
    data_name_amended <- c( paste0(attr(config, "config"), "_data_processed_with_riskscore"), 
                            paste0(attr(config, "config"), "senior_data_processed_with_riskscore"),
                            paste0(attr(config, "config"), "nonsenior_data_processed_with_riskscore"))
  } else {
    data_name_amended <- c(paste0(attr(config, "config"), "_data_processed_with_riskscore"))
  }
  
  cat("Enter reason for updating adata without quotes (this text will be added to adata/README change log): ")
  args <- readLines(con = "stdin", n = 1)
  
  # Request reason for adata update from deployer!
  update_reason <- paste0(Sys.Date(), " ", args[[1]]) 
  
  # set deployment path for each study
  deploy_path <- switch(study_name,
                        COVE =      "/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                        AZD1222 =   "/trials/covpn/p3002/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                        ENSEMBLE =  "/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                        PREVENT19 = "/trials/covpn/p3004/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                        VAT08m =    "/trials/covpn/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                        PROFISCOV = "/networks/cavd/Objective 4/GH-VAP/ID127-Gast/correlates/adata/",
                        stop("study_name not supported 1"))  
  
}


  
# Copy current deployed copy of risk score dataset in adata to archive 
# Remove current deployed copy of risk score dataset from adata
# Copy new copy of risk score dataset from data_clean to adata
for(j in 1:length(data_name_amended)){
  if(file.exists(paste0("data_clean/", data_name_amended[j], ".csv"))){
        if(file.exists(paste0(deploy_path, data_name_amended[j], ".csv"))){
        
              if(!file.exists(paste0(deploy_path, "archive"))){
                dir.create(paste0(deploy_path, "archive"))
              }
              
              file.copy(from = paste0(deploy_path, data_name_amended[j], ".csv"),
                        to =  paste0(deploy_path, "archive/", data_name_amended[j], "_", 
                                     str_replace(str_replace_all(file.info(paste0(deploy_path, data_name_amended[j], ".csv"))$mtime,
                                                                 ":",
                                                                 "."), " ", " time "), 
                                     ".csv"),
                        copy.date = TRUE)
              
              file.remove(from = paste0(deploy_path, data_name_amended[j], ".csv"))
        } 
    
    file.copy(from = paste0("data_clean/", data_name_amended[j], ".csv"),
              to = paste0(deploy_path, data_name_amended[j], ".csv"),
              copy.date = TRUE)
    
    # Add reason for adata update to README file!
    write(update_reason, file = paste0(deploy_path, "readme.txt"), append=TRUE)
    
  }else{
    print(paste0("data_clean/", data_name_amended[j], ".csv not found!"))
  }
}
