library(digest)

if(attr(config, "config") == "prevent19"){
  tmp <- switch(args[1],
                onlyUSsubjects = "76405bdd76df7565dc6fa697a0e34d0e", 
                allsubjects = "07e22b6acf911b1aa7993722332971d2", 
                SLnotrun = "76405bdd76df7565dc6fa697a0e34d0e",
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



if(Sys.getenv("NOCHECK") == "" &
   all.equal(names(inputFile_with_riskscore %>% select(Ptid, risk_score, standardized_risk_score)), c("Ptid", "risk_score", "standardized_risk_score"))){
  
  if(!is.na(tmp)){
    if(attr(config, "config") %in% c("moderna_real", "moderna_mock", "janssen_pooled_mock", "janssen_pooled_real", "azd1222", "vat08_combined")){
      assertthat::assert_that(digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))]) == tmp, 
                              msg = paste0("failed risk_score digest check. Old digest ", tmp, ". New digest ", digest(inputFile_with_riskscore[order(names(inputFile_with_riskscore))])))
      
      print("======================= Passed risk_score digest check =======================") 
    } else {
      only_riskscores <- inputFile_with_riskscore %>% select(Ptid, risk_score, standardized_risk_score) 
      assertthat::assert_that(digest(only_riskscores[order(names(only_riskscores))]) == tmp, 
                              msg = paste0("failed risk_score digest check. Old digest ", tmp, ". New digest ", digest(only_riskscores[order(names(only_riskscores))]))) 
      
      print("======================= Passed risk_score digest check =======================")
    }
  }
}
  



