print("CLEAN_OUTPUT_DIR.R")

if(study_name %in% c("VAT08m", "VAT08b", "PREVENT19")){
  files_to_remove <- grep(list.files(path=paste0("output/", Sys.getenv("TRIAL"), "/", args[1])), 
                          pattern='inputFile.RData', invert=TRUE, value=TRUE)
  
  if(length(files_to_remove) > 0){
    file.remove(paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/", files_to_remove))
  }
}else{
  files_to_remove <- grep(list.files(path=paste0("output/", Sys.getenv("TRIAL"))), 
                          pattern='inputFile.RData', invert=TRUE, value=TRUE)
  
  if(length(files_to_remove) > 0){
    file.remove(paste0("output/", Sys.getenv("TRIAL"), "/", files_to_remove))
  }
}


#file.remove(list.files(paste0("output/", Sys.getenv("TRIAL")), full.names = TRUE))
#file.remove((list.files(here("figs"), full.names = TRUE)))







