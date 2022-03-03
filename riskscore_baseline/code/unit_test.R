library(digest)
if(attr(config, "config") %in% c("janssen_pooled_mock", "moderna_mock") & Sys.getenv ("NOCHECK")=="") {
  assertthat::assert_that(
    digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == ifelse(attr(config, "config")=="janssen_pooled_mock", 
                                                                                       "f7a8225eb5fa8cc9a5426211988b9d95", 
                                                                                       "95368009ca10fc4b2e075885442e6e31"),
    msg = "Failed risk score digest check. New digest for inputFile_with_riskscore is "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))
  print("======================= Passed risk score digest check. Digest for inputFile_with_riskscore is same: "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))])%.%" =======================")  
}

