library(digest)
if(attr(config, "config") %in% c("janssen_pooled_mock", "moderna_mock")) {
  assertthat::assert_that(
    digest(inputFile_with_riskscore)==ifelse(attr(config, "config")=="janssen_pooled_mock", "fa4b46ae6c39fdda27dccc944938aaad", "43895d21d723439f96d183c8898be370"),
    msg = "failed sanity check")
  print("Passed sanity check")
}
