# Sys.setenv(TRIAL = "janssen_pooled_realbAb")  
# Sys.setenv(TRIAL = "moderna_real")  
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here())
#-----------------------------------------------

config <- config::get(config = Sys.getenv("TRIAL"))

if(Sys.getenv("TRIAL")=="") stop("Environmental variable TRIAL not defined!!!!!!!!!!!!!!")
TRIAL=Sys.getenv("TRIAL")

config <- config::get(config = TRIAL)
for(opt in names(config))eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))


library(stringr)
library(here)
library(tidyverse)


# set deployment path for each study
deploy_path <- switch(study_name,
                      COVE =      "/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                      COVEBoost =      "/trials/covpn/p3001/analysis/correlates/Part_C_Unblinded_Phase_Data/adata/",
                      AZD1222 =   "/trials/covpn/p3002/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                      ENSEMBLE =  "/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                      PREVENT19 = "/trials/covpn/p3004/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                      VAT08 =    "/trials/covpn/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/adata/",
                      PROFISCOV = "/networks/cavd/Objective 4/GH-VAP/ID127-Gast/correlates/adata/",
                      IARC_HPV = "/networks/cavd/Objective 4/GH-VAP/ID27-Sankaranarayanan/analysis/correlates/adata/",
                      COVAIL =    "/trials/covpn/COVAILcorrelates/analysis/correlates/adata/",
                      NVX_UK302 =    "/trials/covpn/p3004/analysis/correlates/UK302/adata/",
                      stop("study_name not supported 1"))  



if (TRIAL %in% c()) {
  data_name_amended <- c(paste0(attr(config, "config"), "_data_processed_with_riskscore"))
  
} else if(attr(config, "config") %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA")) {
  data_name_amended <- c( paste0(attr(config, "config"), "_data_processed_with_riskscore"), 
                          paste0(attr(config, "config"), "senior_data_processed_with_riskscore"),
                          paste0(attr(config, "config"), "nonsenior_data_processed_with_riskscore"))
  
} else {
  # c("prevent19", "moderna_real", "moderna_boost", "janssen_partA_VL", "vat08_combined", "vat08_nAb", "id27hpv", "covail", "nvx_uk302")) {
  data_name_amended <- c(paste0(attr(config, "config"), "_data_processed_", format(Sys.Date(), "%Y%m%d")))
}

cat("Enter reason for updating adata without quotes (this text will be added to adata/changelog): ")
args <- readLines(con = "stdin", n = 1)


# Request reason for adata update from deployer!
update_reason <- paste0(Sys.Date(), " ", args[[1]]) 


  
# Copy current deployed copy of risk score dataset in adata to archive 
# Remove current deployed copy of risk score dataset from adata
# Copy new copy of risk score dataset from data_clean/csv to adata
for(j in 1:length(data_name_amended)){
  if(file.exists(paste0("data_clean/csv/", data_name_amended[j], ".csv"))){
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
    
    file.copy(from = paste0("data_clean/csv/", data_name_amended[j], ".csv"),
              to = paste0(deploy_path, data_name_amended[j], ".csv"),
              copy.date = TRUE)
    cat(paste0("Deployed to: ", deploy_path, data_name_amended[j], ".csv\n"))
    
    # Add reason for adata update to changelog!
    write(update_reason, file = paste0(deploy_path, "changelog.txt"), append=TRUE)
    
  }else{
    print(paste0("data_clean/csv/", data_name_amended[j], ".csv not found!"))
  }
}
