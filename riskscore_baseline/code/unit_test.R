library(digest)

if(attr(config, "config") == "vat08_combined"){
  tmp <- switch(args[1],
                bseroneg = "95d865388d3e7dfa35900d2c170b1671", #"b6044e6b8ca8ddfadf0708e89550dfb9",
                bseropos = "99c40506a225b1144f0218493947a44f", #"66b6a2e3c7dc48dcfbfa8f2990a39b5f",
                SLnotrun = "c4482f985c4ad6d50b2c07ad5c03999a",
                NA) 
} else if(attr(config, "config") == "prevent19"){
  tmp <- switch(args[1],
                onlyUSsubjects = "432e24718a92a9ee3c6fe4a71af177f8", 
                allsubjects = "41588e90ae5325c87b5304b2b9087a3b", 
                SLnotrun = "795a3b72a0dde6aff71b7487f738839c",
                NA) 
} else {
  tmp <- switch(attr(config, "config"),
                moderna_real = "9f238311b6f252204336eecdc85efdc1",
                moderna_mock = "9df4cd6639381811e763c2dddc0a12fd",
                janssen_pooled_mock = "f7a8225eb5fa8cc9a5426211988b9d95",
                janssen_pooled_real = "c38fb43e2c87cf2d392757840af68bba",
                azd1222 = "23cee1e1ba96ef85326dadceb34b3c6f",
                NA)   
}


if (Sys.getenv("NOCHECK") == "" &
    all.equal(names(inputFile_with_riskscore %>% select(Ptid, risk_score, standardized_risk_score)), c("Ptid", "risk_score", "standardized_risk_score"))) {
  
  if (!is.na(tmp)) assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == tmp, 
                                           msg = "failed risk_score digest check. new digest "%.%digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))])) 
  
  print("======================= Passed risk_score digest check =======================")    
}
