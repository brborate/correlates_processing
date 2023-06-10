load("output/vat08m/bseroneg/inputFile_with_riskscore.RData")
bseroneg_with_riskscore <- inputFile_with_riskscore
rm(inputFile_with_riskscore)
load("output/vat08m/bseropos/inputFile_with_riskscore.RData")
bseropos_with_riskscore <- inputFile_with_riskscore
rm(inputFile_with_riskscore)

inputFile_bseroneg_bseropos_with_riskscore <- dplyr::bind_rows(bseroneg_with_riskscore, 
                                                        bseropos_with_riskscore)

save(inputFile_bseroneg_bseropos_with_riskscore,
     file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_bseroneg_bseropos_with_riskscore.RData")) 
