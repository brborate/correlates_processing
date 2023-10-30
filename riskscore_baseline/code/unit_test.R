library(digest)

if(attr(config, "config") == "vat08_combined"){
  tmp <- switch(args[1],
                bseroneg = "f8a532f9d53a6815b15d38f8e906f3e7", 
                bseropos = "1640fc96414a3af75bc6cf098d49113c", 
                SLnotrun = "8a6896b15935b9ee207f369416dd3fc0",
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
  
  if (!is.na(tmp)) 
    assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == tmp, 
                            msg = paste0("failed risk_score digest check. Old digest ", tmp, ". New digest ", digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]))) 
  
  print("======================= Passed risk_score digest check =======================")    
}
