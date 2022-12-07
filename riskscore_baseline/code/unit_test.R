library(digest)

tmp <- switch(attr(config, "config"),
             moderna_real = "e0eed0edd95569334059e09ac374ca50",
             moderna_mock = "9df4cd6639381811e763c2dddc0a12fd",
             janssen_pooled_mock = "f7a8225eb5fa8cc9a5426211988b9d95",
             janssen_pooled_real = "c38fb43e2c87cf2d392757840af68bba",
             azd1222 = "8dd08a2f6738091ebdd4748d5b631682",
             #azd1222_bAb = "",
             prevent19 = "53f991ec3b75c8643f2c90fe7252c25e",
             NA)    


if (Sys.getenv("NOCHECK") == "" &
    all.equal(names(inputFile_with_riskscore %>% select(Ptid, risk_score, standardized_risk_score)), c("Ptid", "risk_score", "standardized_risk_score"))) {
  
  if (!is.na(tmp)) assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == tmp, 
                                           msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))])) 
  
  print("======================= Passed risk_score digest check =======================")    
}
