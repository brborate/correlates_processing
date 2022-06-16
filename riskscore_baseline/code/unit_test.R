# library(digest)
# if(attr(config, "config") %in% c("janssen_pooled_mock", "moderna_mock") & Sys.getenv ("NOCHECK")=="") {
#   assertthat::assert_that(
#     digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == ifelse(attr(config, "config")=="janssen_pooled_mock", 
#                                                                                        "f7a8225eb5fa8cc9a5426211988b9d95", 
#                                                                                        "95368009ca10fc4b2e075885442e6e31"),
#     msg = "Failed risk score digest check. New digest for inputFile_with_riskscore is "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))
#   print("======================= Passed risk score digest check. Digest for inputFile_with_riskscore is same: "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))])%.%" =======================")  
# }


library(digest)
if (Sys.getenv("NOCHECK") == "") {
  if (attr(config, "config") == "moderna_mock") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "95368009ca10fc4b2e075885442e6e31", 
                            msg = "failed risk_score digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))]))    
  } else if (attr(config, "config") == "janssen_pooled_mock") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "f7a8225eb5fa8cc9a5426211988b9d95", 
                            msg = "failed risk_score digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))]))    
  } else if (attr(config, "config") == "prevent19") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "ae6dd4697c65a4b7511f890c219d17be", 
                            msg = "failed risk_score digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))]))    
  } else if (attr(config, "config") == "vat08m") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "232efb271e32526c042350921fa74041", 
                            msg = "failed risk_score digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))]))    
  } 
  print("======================= Passed risk_score digest check =======================")    
}