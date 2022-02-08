library(digest)
if(attr(config, "config") %in% c("janssen_pooled_mock", "moderna_mock")) {
  assertthat::assert_that(
    # digest(inputFile_with_riskscore)==ifelse(attr(config, "config")=="janssen_pooled_mock", "5cb73bf9810a11ae0b3d6638753bb690", "76bae5ccbfce93842efcd8b93fa3f38f"),
    # msg = "failed sanity check")
    digest(inputFile_with_riskscore)==ifelse(attr(config, "config")=="janssen_pooled_mock", "2e68d1710090ba47dbf22420afb4dac7", "76bae5ccbfce93842efcd8b93fa3f38f"),
    msg = "failed sanity check")
  print("Passed sanity check")
}
