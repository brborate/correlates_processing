library(digest)

tmp <- switch(attr(config, "config"),
             moderna_real = "e0eed0edd95569334059e09ac374ca50",
             moderna_mock = "9df4cd6639381811e763c2dddc0a12fd",
             janssen_pooled_mock = "f7a8225eb5fa8cc9a5426211988b9d95",
             janssen_pooled_real = "c38fb43e2c87cf2d392757840af68bba",
             azd1222 = "23cee1e1ba96ef85326dadceb34b3c6f",
             #azd1222_bAb = "",
             prevent19 = "16ec3b6ee4e5b2f9e03755e0d2b033d7",
             vat08m = "477d501bbcbc1d002fe2e75cd5e46146",
             NA)    


if (Sys.getenv("NOCHECK") == "" &
    all.equal(names(inputFile_with_riskscore %>% select(Ptid, risk_score, standardized_risk_score)), c("Ptid", "risk_score", "standardized_risk_score"))) {
  
  if (!is.na(tmp)) assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == tmp, 
                                           msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))])) 
  
  print("======================= Passed risk_score digest check =======================")    
}
