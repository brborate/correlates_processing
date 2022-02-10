# Activate renv
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))

# load required libraries, functions; Read input data
library(tidyverse)
library(here)
library(methods)
library(SuperLearner)
library(e1071)
library(glmnet)
library(kyotil)
library(argparse)
library(vimp)
library(nloptr)
library(RhpcBLASctl)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("omp_set_num_threads", "RhpcBLASctl")
library(mice)
library(tidymodels)
library(Hmisc) # wtd.quantile, cut2
library(dplyr)
library(recipes)

if(startsWith(tolower(study_name), "mock")) {
  path_to_data <- here("..", "data_raw", data_raw_dir, data_in_file)
} else {
  path_to_data <- data_in_file
}
print(path_to_data)
if (!file.exists(path_to_data)) stop ("make dat proc: dataset not available ===========================================")
inputfileName <- gsub("^.*/", "", data_in_file)
# Define code version to run
# the demo version is simpler and runs faster!
# the production version runs SL with a diverse set of learners
run_prod <- !grepl("Mock", study_name)

# get utility files
source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
source(here("code", "utils.R")) # get CV-AUC for all algs
inputFile <- preprocess.for.risk.score(read.csv(path_to_data), study_name)