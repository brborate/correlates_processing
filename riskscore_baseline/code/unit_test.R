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
  if (attr(config, "config") == "moderna_real") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "4d1e96b334113ce3099ca35f5999637a", 
                            msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))    
  } else if (attr(config, "config") == "moderna_mock") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "95368009ca10fc4b2e075885442e6e31", 
                            msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))    
  } else if (attr(config, "config") == "janssen_pooled_mock") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "f7a8225eb5fa8cc9a5426211988b9d95", 
                            msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))    
  } else if (attr(config, "config") == "prevent19") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "53f991ec3b75c8643f2c90fe7252c25e", 
                            msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))    
  } #else if (attr(config, "config") == "vat08m") {
    #assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "3cc3c3c9d1536be807bdde83788bddd7", 
    #                        msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))    
  #} 
  else if (attr(config, "config") == "azd1222") {
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == "90bc4f44f626a160d866f1cf7369c645", 
                            msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))    
  } 
  print("======================= Passed risk_score digest check =======================")    
}