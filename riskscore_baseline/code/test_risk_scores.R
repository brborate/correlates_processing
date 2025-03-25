# Sys.setenv(TRIAL = "janssen_pooled_realbAb")
# Sys.setenv(TRIAL = "prevent19")
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries, cleaned data, and risk score estimates
library(here)
library(tidyverse)
library(conflicted)
library(gridExtra)
conflicted::conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")

print("TEST_RISK_SCORES.R")

if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  output_path = paste0("output/", Sys.getenv("TRIAL"), "/", args[1])
} else {
  output_path = paste0("output/", Sys.getenv("TRIAL"), "/")
}

if(study_name %in% c("VAT08m", "VAT08")){
  load(paste0(output_path, "/inputFile_with_riskscore.RData"))
  if(args[1] == "bseroneg"){
    dat = inputFile_with_riskscore %>% 
      filter(Bserostatus == 0) 
  } else if(args[1] == "bseropos"){
    dat = inputFile_with_riskscore %>% 
      filter(Bserostatus == 1) 
  }
} else {
  load(paste0(output_path, "/inputFile_with_riskscore.RData"))
  dat = inputFile_with_riskscore 
}

# Density plot: risk scores & standardized risk scores grouped by Trt
options(bitmapType = "cairo")
if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  png(file = here(output_path, paste0("density_plot_riskscores_", args[1], ".png")),
      width = 1000, height = 1000)
}else{
  png(file = here(output_path, "density_plot_riskscores.png"),
      width = 1000, height = 500)
}

p1 <- grid.arrange(dat %>% 
                     ggplot(aes(x = risk_score, fill = factor(Trt))) +
                     geom_density(alpha = 0.3) +  # Add transparency for overlapping areas
                     labs(
                       title = paste0(args[1], ": Density Plot of risk scores grouped by Trt"),
                       x = "Risk scores",
                       y = "Density",
                       fill = "Trt"
                     ), 
                   dat %>%
                     ggplot(aes(x = standardized_risk_score, fill = factor(Trt))) +
                     geom_density(alpha = 0.3) +  # Add transparency for overlapping areas
                     labs(
                       title = paste0(args[1], ": Density Plot of standardized risk scores grouped by Trt"),
                       x = "Standardized risk scores",
                       y = "Density",
                       fill = "Trt"
                     ), 
                   ncol = 1)
 
print(p1)
dev.off()


# Ask user to check the created density plots! 
# Convert input to numeric (since readline() returns a string)

cat(paste0("CHECK THE DENSITY PLOTS CREATED IN riskscore_baseline/", output_path, "\n"))



