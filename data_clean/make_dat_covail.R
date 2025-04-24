Sys.setenv(TRIAL = "covail")
source(here::here("_common.R"))
# no need to run renv::activate(here::here()) b/c .Rprofile exists

{
library(tidyverse)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)
library(here)
library(mdw)

begin=Sys.time()
}


########################################################################################################
# 1. read mapped data with risk score added

{
  # # load risk score from running risk analysis
  # load(file = paste0('riskscore_baseline/output/',TRIAL,'/inputFile_with_riskscore.RData'))
  # dat_proc <- inputFile_with_riskscore    
  
  # load risk score from a file
  dat.risk = read.csv("/trials/covpn/COVAILcorrelates/analysis/correlates/adata/risk_score.csv")
  # read mapped data
  dat_raw = read.csv(mapped_data)
  dat_proc = preprocess(dat_raw, study_name)   
  names(dat_proc)[[1]]="Ptid"
  dat_proc$risk_score = dat.risk$risk_score[match(dat_proc$Pti, dat.risk$Ptid)]
  dat_proc$standardized_risk_score = dat.risk$standardized_risk_score[match(dat_proc$Ptid, dat.risk$Ptid)]
  
  # bring in imputed variant column
  dat.lineage = read.csv('/trials/covpn/COVAILcorrelates/analysis/correlates/adata/lineages/covail_lineages_export_v1.csv')
  dat_proc$COVIDlineage = dat.lineage$inf1.lineage[match(dat_proc$Ptid, dat.lineage$ptid)]
  dat_proc$COVIDlineageObserved = !dat.lineage$inf1.imputed[match(dat_proc$Ptid, dat.lineage$ptid)]
  # check NA
  stopifnot(!any(is.na(dat_proc$COVIDlineage[dat_proc$ph1.D15==1 & dat_proc$COVIDIndD22toD181==1])))
  stopifnot(!any(is.na(dat_proc$COVIDlineage[dat_proc$ph1.D29==1 & dat_proc$COVIDIndD36toD181==1])))
  # this is not true: !any(is.na(dat_proc$COVIDlineage[dat_proc$ph1.D15==1 & dat_proc$AsympInfectIndD15to181==1]))
  
  # bring in FOI
  dat.foi = read.csv('/trials/covpn/COVAILcorrelates/analysis/correlates/adata/covail_foi_v2.csv')
  dat_proc$FOIoriginal = dat.foi$foi[match(dat_proc$Ptid, dat.foi$ptid)]
  dat_proc$FOIstandardized = scale(dat_proc$FOIoriginal)
  # check NA
  stopifnot(!any(is.na(dat_proc$FOI[dat_proc$ph1.D15==1])))
  stopifnot(!any(is.na(dat_proc$FOI[dat_proc$ph1.D29==1])))
  
}


########################################################################################################
# define Senior and race/ethnicity
{
  colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
  dat_proc <- dat_proc %>% mutate(age.geq.65 = as.integer(Age >= 65))
  dat_proc$Senior = as.integer(dat_proc$Age>=switch(study_name, COVE=65, MockCOVE=65, ENSEMBLE=60, MockENSEMBLE=60, PREVENT19=65, AZD1222=65, VAT08=60, PROFISCOV=NA, COVAIL=65, NVX_UK302=65, stop("unknown study_name 1")))
  
  # ethnicity labeling
  dat_proc$ethnicity <- ifelse(dat_proc$EthnicityHispanic == 1, labels.ethnicity[1], labels.ethnicity[2])
  dat_proc$ethnicity[dat_proc$EthnicityNotreported == 1 | dat_proc$EthnicityUnknown == 1] <- labels.ethnicity[3]
  dat_proc$ethnicity <- factor(dat_proc$ethnicity, levels = labels.ethnicity)
  
  dat_proc$race = 0 # not applicable, but has to define a value so that the next chunk of code can run

  dat_proc$WhiteNonHispanic <- NA
  # WhiteNonHispanic=1 IF race is White AND ethnicity is not Hispanic
  dat_proc$WhiteNonHispanic <-
    ifelse(dat_proc$race == "White" &
             dat_proc$ethnicity == "Not Hispanic or Latino", 1,
           dat_proc$WhiteNonHispanic
    )
  # WhiteNonHispanic=0 IF race is not "white or unknown" OR ethnicity is Hispanic
  dat_proc$WhiteNonHispanic <-
    ifelse(!dat_proc$race %in% c(labels.race[1], last(labels.race)) |
             dat_proc$ethnicity == "Hispanic or Latino", 0,
           dat_proc$WhiteNonHispanic
    )
  dat_proc$MinorityInd = 1-dat_proc$WhiteNonHispanic
  # set NA to 0 in both WhiteNonHispanic and MinorityInd. This means the opposite things for how NA's are interpreted and that gave us the option to use one or the other
  dat_proc$WhiteNonHispanic[is.na(dat_proc$WhiteNonHispanic)] = 0
  dat_proc$MinorityInd[is.na(dat_proc$MinorityInd)] = 0

}




###############################################################################
# 3. stratum variables

# there are no demographics stratum for subcohort sampling
dat_proc$Bstratum =  1 

names(Bstratum.labels) <- Bstratum.labels
  

# demo.stratum: correlates sampling strata
# Moderna: 1 ~ 6 defines the 6 baseline strata within trt/serostatus
# may have NA b/c URMforsubcohortsampling may be NA
dat_proc$demo.stratum = 1 # # there are no demographics stratum for subcohort sampling
dat_proc <- dat_proc %>% mutate(tps.stratum = arm)
if (!is.null(dat_proc$tps.stratum)) table(dat_proc$tps.stratum)


# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+4. 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# A case will have a Wstratum even if its tps.stratum is NA
# The case is defined using EventIndPrimaryD29
dat_proc$Wstratum = dat_proc$tps.stratum

if (!is.null(dat_proc$Wstratum)) table(dat_proc$Wstratum) # variables may be named other than Wstratum



################################################################################
# 4. Define ph1, ph2, and weights
# Note that Wstratum may have NA if any variables to form strata has NA

# define must_have_assays for ph2 definition
must_have_assays <- NULL
# the whole cohort is treated as ph1 and ph2
dat_proc$TwophasesampIndD15 = dat_proc$ph1.D15 
dat_proc$TwophasesampIndD29 = dat_proc$ph1.D29
  

# PP = no violation + marker available at d1 and d15
# Immunemarkerset = PP & no infection between enrollment and D15+6
# ph1.D15 = Immunemarkerset & arm!=3
dat_proc[["ph2.D15"]]=dat_proc$ph1.D15
dat_proc[["wt.D15"]] = 1
dat_proc[["ph2.D92"]]=dat_proc$ph1.D92
dat_proc[["wt.D92"]] = 1
dat_proc[["ph2.D29"]]=dat_proc$ph1.D29
dat_proc[["wt.D29"]] = 1




###############################################################################
# 5. impute missing biomarkers in ph2 (assay imputation)
#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, each time point, but not Delta

# loop through the time points
# first impute (B, D29, D57) among TwophasesampIndD57==1
# next impute (B, D29) among TwophasesampIndD29==1


#### impute D15 BA.4/BA.5 nAb among ph1.D15

dat.tmp.impute <- subset(dat_proc, ph1.D15 == 1)

# first, do it based on a linear regression with D29 BA.4/BA.5 nAb and number of days between the D15 and D29 visit
fit = lm(Day15pseudoneutid50_BA.4.BA.5~Day29pseudoneutid50_BA.4.BA.5, dat.tmp.impute)
predicted = predict(fit, subset(dat.tmp.impute, select=Day29pseudoneutid50_BA.4.BA.5))
dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5 = ifelse(is.na(dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5), predicted, dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5)

# there are still some ptids with missing D15 BA4.5 b/c not all have D29 BA4.5. impute with BA1: mypairs(dat.tmp.impute["Day15"%.%assays[1:5]])
fit = lm(Day15pseudoneutid50_BA.4.BA.5~Day15pseudoneutid50_BA.1, dat.tmp.impute)
predicted = predict(fit, subset(dat.tmp.impute, select=Day15pseudoneutid50_BA.1))
dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5 = ifelse(is.na(dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5), predicted, dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5)

# save a copy for Lauren before assigning imputed values
dat_proc$Day15pseudoneutid50_BA.4.BA.5_unimputed = dat_proc$Day15pseudoneutid50_BA.4.BA.5
# populate dat_proc with the imputed values
imp.markers='Day15pseudoneutid50_BA.4.BA.5'
dat_proc[dat_proc[["ph1.D15"]]==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["ph1.D15"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]

assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc[["ph1.D15"]] == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)


#### none missing to impute D29 among those at risk at D29
dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD29 & AsympInfectIndD15to29==0)
with(dat.tmp.impute, print(table(!is.na(get("Day29"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
# thus, no missingness actually

#### none missing to impute D91 among those at risk at D91
dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD91 & AsympInfectIndD15to91==0)
with(dat.tmp.impute, print(table(!is.na(get("Day91"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
# thus, no missingness actually

#### none missing to impute D181 among those at risk at D181
dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD181 & AsympInfectIndD15to181==0)
with(dat.tmp.impute, print(table(!is.na(get("Day181"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
# thus, no missingness actually



###############################################################################
# 6. transformation of the markers
# e.g. converting binding variables from AU to IU for binding assays

assays.includeN = assays


###############################################################################
# 7. add mdw scores

kp = dat_proc$ph1.D15==1

# myboxplot(dat_proc[kp, c("B"%.%assays[1:5], "Day15"%.%assays[1:5])], names=sub("pseudoneutid50_", "", rep(assays[1:5],2)))
# mypairs(dat_proc["Day15"%.%assays[1:5]])
# mypairs(dat_proc["B"%.%assays[1:5]])
# corplot(Day15pseudoneutid50_D614G~Bpseudoneutid50_D614G, dat_proc)
# sapply(dat_proc[kp, c("B"%.%assays[1:5], "Day15"%.%assays[1:5])], sd)

nAb = setdiff(assays[startsWith(assays, "pseudoneutid50_")], c('pseudoneutid50_MDW'))

# use Day15 to derive weights
mdw.wt.nAb=tryCatch({
  tree.weight(cor(dat_proc[kp, "Day15"%.%nAb], use='complete.obs'))
}, error = function(err) {
  print(err$message)
  rep(1/length(nAb), length(nAb))
})
# using Day15 and using B lead to similar weights: BA1 and BA4 BA5 together takes about half of the weight
print(mdw.wt.nAb)
write.csv(mdw.wt.nAb, file = here("data_clean", "csv", TRIAL%.%"_nAb_mdw_weights.csv"))

# apply to all time points
for (t in c("B", "Day15", "Day29", "Day91", "Day181")) {
  dat_proc[, t%.%'pseudoneutid50_MDW'] = as.matrix(dat_proc[, t%.%nAb]) %*% mdw.wt.nAb
}


###############################################################################
# 8. add fold change markers
# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta

assays1=assays
# mdw scores delta are computed as weighted average of delta, not as delta of mdw
# assays1 = assays.includeN[!endsWith(assays.includeN, "_mdw")]

tmp=list()
for (a in assays1) {
  for (t in c("B", paste0("Day", config$timepoints)) ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat_proc[[t %.% a]])
  }
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

for (tp in rev(timepoints)) {
  dat_proc["Delta"%.%tp%.%"overB" %.% assays1] <- tmp["Day"%.%tp %.% assays1] - tmp["B" %.% assays1]
}   

if(two_marker_timepoints) {
  dat_proc["Delta"%.%timepoints[2]%.%"over"%.%timepoints[1] %.% assays1] <- tmp["Day"%.% timepoints[2]%.% assays1] - tmp["Day"%.%timepoints[1] %.% assays1]
}


# also need D29 delta for sanofi arms
assays1=setdiff(assays, "pseudoneutid50Duke_BA.2.12.1")

tmp=list()
for (a in assays1) {
  for (t in c("B", "Day29") ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat_proc[[t %.% a]])
  }
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

dat_proc["Delta29overB" %.% assays1] <- tmp["Day29" %.% assays1] - tmp["B" %.% assays1]



###############################################################################
# 9. add discrete/trichotomized markers

# mRNA arms
dat_proc$tmp = with(dat_proc, ph1.D15 & TrtonedosemRNA==1) 
assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
all.markers = c("B"%.%assays, "Day15"%.%assays, "Delta15overB"%.%assays)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D15")

# # Sanofi arms
dat_proc$tmp = with(dat_proc, ph1.D29 & TrtSanofi==1)
assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
all.markers = c("Day29"%.%assays, "Delta29overB"%.%assays)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D29")

# remove the temp ph2 column
dat_proc$tmp = NULL


# subset on subset_variable
if(!is.null(config$subset_variable) & !is.null(config$subset_value)){
  if(subset_value != "All") {
    include_in_subset <- dat_proc[[subset_variable]] == subset_value
    dat_proc <- dat_proc[include_in_subset, , drop = FALSE]
  }
}


###############################################################################
# 10. impute covariates if necessary



###############################################################################
# special handling 




###############################################################################
# digest check

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         covail = "bc72cec5055c250ba9bfca8513ab3344",
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, 
      msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))



print("run time: "%.%format(Sys.time()-begin, digits=1))
