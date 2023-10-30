library(digest)

if(attr(config, "config") == "vat08_combined"){
  tmp <- switch(args[1],
                bseroneg = "aa894cc045426014202293c5ac7da4e1",
                bseropos = "c3dba271f211a97f90005ab8d6ff1422", 
                SLnotrun = "8baa0b7080c1dac718e58d31327e47ab",
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
