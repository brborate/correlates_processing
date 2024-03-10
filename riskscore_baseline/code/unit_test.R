library(digest)

if(attr(config, "config") == "prevent19"){
  tmp <- switch(args[1],
                onlyUSsubjects = "64709ad36f1e982b3a1c9fa3c0d20532", 
                allsubjects = "eed29375dfb27bdd9c25abfe769d9a1f", 
                SLnotrun = "0cffe738dca698637ec40c2c3cc2a476",
                NA) 
} else {
  tmp <- switch(attr(config, "config"),
                moderna_real = "9f238311b6f252204336eecdc85efdc1",
                moderna_mock = "9df4cd6639381811e763c2dddc0a12fd",
                janssen_pooled_mock = "f7a8225eb5fa8cc9a5426211988b9d95",
                janssen_pooled_real = "c38fb43e2c87cf2d392757840af68bba",
                azd1222 = "23cee1e1ba96ef85326dadceb34b3c6f",
                vat08_combined = "e6cb5d601d6bd45e26401d991affba15",
                NA)   
}


if (Sys.getenv("NOCHECK") == "" &
    all.equal(names(inputFile_with_riskscore %>% select(Ptid, risk_score, standardized_risk_score)), c("Ptid", "risk_score", "standardized_risk_score"))) {
  
  if (!is.na(tmp)) 
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == tmp, 
                            msg = paste0("failed risk_score digest check. Old digest ", tmp, ". New digest ", digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))) 
  
  print("======================= Passed risk_score digest check =======================")    
}
