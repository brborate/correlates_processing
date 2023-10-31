# This code is called only for study vat08m (Sanofi), where risk scores are developed for 
# baseline seronegative and baseline seropositive subjects separately. 

# load required libraries, cleaned data, and risk score estimates
library(here)
library(tidyverse)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")

print("STACK_bseroneg_bseropos_riskscores.R")

load("output/vat08_combined/bseroneg/inputFile.RData")
dneg_inputFile = inputFile
rm(inputFile)
load("output/vat08_combined/bseropos/inputFile.RData")
dpos_inputFile = inputFile
rm(inputFile)

# load("output/vat08m/bseroneg/inputFile_with_riskscore.RData")
# bseroneg_with_riskscore <- inputFile_with_riskscore
# rm(inputFile_with_riskscore)
# load("output/vat08m/bseropos/inputFile_with_riskscore.RData")
# bseropos_with_riskscore <- inputFile_with_riskscore
# rm(inputFile_with_riskscore)

load("output/vat08_combined/bseroneg/risk_scores.RData")
bseroneg_riskscores <- risk_scores
rm(risk_scores)
load("output/vat08_combined/bseropos/risk_scores.RData")
bseropos_riskscores <- risk_scores
rm(risk_scores)

risk_scores <- bind_rows(bseroneg_riskscores, bseropos_riskscores)

if(!identical(dneg_inputFile, dpos_inputFile)){
  stop("Execution stopped. For study VAT08, dneg_inputFile is not same as dpos_inputFile!")
}else{
  inputFile_with_riskscore <- left_join(dneg_inputFile, risk_scores, by = "Ptid") 
} 
  
save(inputFile_with_riskscore,
     file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData")) 

source("code/get_riskscores.R")
