library(digest)
if(attr(config, "config") %in% c("janssen_pooled_mock", "moderna_mock") & Sys.getenv ("NOCHECK")=="") {
  assertthat::assert_that(
    digest(inputFile_with_riskscore)==ifelse(attr(config, "config")=="janssen_pooled_mock", "2e68d1710090ba47dbf22420afb4dac7", "6e61dc1bedfae9969f4a377b13f1109a"),
    msg = "failed sanity check: inputFile_with_riskscore_"%.%digest(inputFile_with_riskscore))
  print("Passed sanity check: inputFile_with_riskscore_"%.%digest(inputFile_with_riskscore))
}

