# Sys.setenv(TRIAL = "janssen_pooled_realbAb")
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries and functions
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
library(aucm)
library(mice)
library(conflicted)
library(gam)
library(xgboost)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))

if(!any(sapply(c("COVE", "janssen"), grepl, study_name)))
  endpoint <- paste0(sub("rscore", "", endpoint), "rauc")

vacc <- read.csv(here("output", Sys.getenv("TRIAL"), "vaccine_ptids_with_riskscores.csv"))

# plot ROC curve on vaccinees
pred.obj <- ROCR::prediction(vacc$pred, vacc %>% pull(endpoint))
perf.obj <- ROCR::performance(pred.obj, "tpr", "fpr")

options(bitmapType = "cairo")
png(file = here("output", Sys.getenv("TRIAL"), "ROCcurve_riskscore_vacc_onlySL.png"),
    width = 1000, height = 1000)

print(data.frame(xval = perf.obj@x.values[[1]],
           yval = perf.obj@y.values[[1]],
           learner = paste0("Superlearner (", unique(vacc$AUCchar), ")")) %>% 
  ggplot(aes(x = xval, y = yval, col = learner)) +
  geom_step(lwd = 2) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.title=element_text(size=20),
    legend.text=element_text(size=20),
    axis.ticks.length = unit(.35, "cm"),
    axis.text = element_text(size = 23),
    axis.title = element_text(size = 30)
  ) +
  labs(x = "False Positive Rate", y = "True Positive Rate", col = "Model (AUC)") +
  geom_abline(intercept = 0, slope = 1) + 
  scale_color_manual(values = "purple"))

dev.off()

# plot pred prob plot on vaccinees
options(bitmapType = "cairo")
png(file = here("output", Sys.getenv("TRIAL"), "predProb_riskscore_vacc_onlySL.png"),
    width = 1100, height = 700)
# if(study_name_code == "COVE"){
#   cases = "Post Day 57 Cases"
# }
# if(study_name_code == "ENSEMBLE"){
#   cases = "Post Day 29 Cases"
# }
print(vacc %>%
  mutate(Ychar = ifelse(get(endpoint) == 0, "Non-Cases", paste0("Post Day ", risk_timepoint, " Cases"))) %>%
  ggplot(aes(x = Ychar, y = pred, color = Ychar)) +
  geom_jitter(width = 0.06, size = 3, shape = 21, fill = "white") +
  geom_violin(alpha = 0.05, color = "black", lwd=1.5) +
  geom_boxplot(alpha = 0.05, width = 0.15, color = "black", outlier.size = NA, outlier.shape = NA, lwd=1.5) +
  theme_bw() +
  #scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#00468B", "#8B0000")) +
  labs(y = "Predicted probability of COVID-19 disease", x = "") +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 25),
    axis.text = element_text(size = 23),
    axis.ticks.length = unit(.35, "cm"),
    axis.title.y = element_text(size = 30)
  ))
dev.off()

