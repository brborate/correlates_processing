Sys.setenv(TRIAL = "azd1222_stage2")
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

dat_raw=read.csv(mapped_data)
dat_proc = preprocess(dat_raw, study_name)   
colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 

# borrow risk score from azd1222
load(file = 'riskscore_baseline/output/azd1222/inputFile_with_riskscore.RData')
# stage 2 dataset has fewer rows than stage 1
dat_proc$risk_score = inputFile_with_riskscore$risk_score[match(dat_proc$Ptid, inputFile_with_riskscore$Ptid)]
dat_proc$standardized_risk_score = inputFile_with_riskscore$standardized_risk_score[match(dat_proc$Ptid, inputFile_with_riskscore$Ptid)]

# subset to vaccine and seroneg, necessary for getting weights correctly
dat_proc = subset(dat_proc, Trt==1 & Bserostatus==0)

# there are 11 ptids without risk_score. these ptids are absent from stage 1 risk score dataset and from stage 1 mapped dataset
# will impute them later in this script
# subset(dat_proc, is.na(risk_score), c(Ptid, Country))
# subset(inputFile_with_riskscore, Ptid=="D8110C00001/E20052580178")
# subset(dat_mapped_stage1, Subjectid=="D8110C00001/E2005329019")


########################################################################################################
# 2. define Senior and race/ethnicity
{
  colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
  dat_proc <- dat_proc %>% mutate(age.geq.65 = as.integer(Age >= 65))
  dat_proc$Senior = as.integer(dat_proc$Age>=switch(study_name, COVE=65, MockCOVE=65, ENSEMBLE=60, MockENSEMBLE=60, PREVENT19=65, AZD1222=65, VAT08=60, PROFISCOV=NA, COVAIL=65, NVX_UK302=65, stop("unknown study_name 1")))
  
  # ethnicity labeling
  dat_proc$ethnicity <- ifelse(dat_proc$EthnicityHispanic == 1, labels.ethnicity[1], labels.ethnicity[2])
  dat_proc$ethnicity[dat_proc$EthnicityNotreported == 1 | dat_proc$EthnicityUnknown == 1] <- labels.ethnicity[3]
  dat_proc$ethnicity <- factor(dat_proc$ethnicity, levels = labels.ethnicity)
  
  # race labeling
  if (study_name %in% c("COVE", "MockCOVE")) {
    dat_proc <- dat_proc %>%
      mutate(
        race = labels.race[1],
        race = case_when(
          Black == 1 ~ labels.race[2],
          Asian == 1 ~ labels.race[3],
          NatAmer == 1 ~ labels.race[4],
          PacIsl == 1 ~ labels.race[5],
          Multiracial == 1 ~ labels.race[6],
          Other == 1 ~ labels.race[7],
          Notreported == 1 | Unknown == 1 ~ labels.race[8], # labels.race has 8 levels for COVE and MockCOVE
          TRUE ~ labels.race[1]
        ),
        race = factor(race, levels = labels.race)
      )
    
  } else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE", "PREVENT19", "AZD1222")) {
    # remove the Other category
    dat_proc <- dat_proc %>%
      mutate(
        race = labels.race[1],
        race = case_when(
          Black == 1 ~ labels.race[2],
          Asian == 1 ~ labels.race[3],
          NatAmer == 1 ~ labels.race[4],
          PacIsl == 1 ~ labels.race[5],
          Multiracial == 1 ~ labels.race[6],
          Notreported == 1 | Unknown == 1 ~ labels.race[7],
          TRUE ~ labels.race[1]
        ),
        race = factor(race, levels = labels.race)
      )
    
  } else if (study_name %in% c("PROFISCOV", "VAT08", "COVAIL", "NVX_UK302")) {
    dat_proc$race = 0 # not applicable, but has to define a value so that the next chunk of code can run
    
  } else stop("unknown study_name 2")
  
  
  
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
  
  # set MinorityInd to 0 for latin america and south africa
  if (study_name %in% c("COVE", "MockCOVE", "PROFISCOV", "VAT08", "COVAIL", "NVX_UK302")) {
    # nothing to do
    # COVE only has US data
    
  } else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    dat_proc$MinorityInd[dat_proc$Region!=0] = 0
    
  } else if (study_name=="PREVENT19") {
    dat_proc$MinorityInd[dat_proc$Country!=0] = 0 # 0 is US
    
  } else if (study_name=="AZD1222") {
    dat_proc$MinorityInd[dat_proc$Country!=2] = 0 # 2 is US
    
  } else stop("unknown study_name 3")  
  
  # for the mock datasets, hardcode AnyinfectionD1 
  if (study_name %in% c("MockENSEMBLE", "MockCOVE")) dat_proc$AnyinfectionD1=0

}



###############################################################################
# 3. stratum variables

# Bstratum: randomization strata
dat_proc$Bstratum =  with(dat_proc, Senior + 1)
names(Bstratum.labels) <- Bstratum.labels


# demo.stratum: correlates sampling strata
# may have NA b/c URMforsubcohortsampling may be NA
# make the order same as Table 3 in the SAP
#    US, >65, non-Minority
#    US, <65, non-Minority
#    US, >65, Minority
#    US, <65, Minority
#    non-US, >65
#    non-US, <65
dat_proc$demo.stratum = with(dat_proc, strtoi(paste0(URMforsubcohortsampling, 1-Senior), base = 2)) + 1
dat_proc$demo.stratum = with(dat_proc, ifelse(Country==2, demo.stratum, ifelse(Senior, 5, 6))) # 2 is US


# tps stratum, used in tps regression and to define Wstratum
# since we only care about vacc and seroneg, we set it to demo.stratum
dat_proc <- dat_proc %>% mutate(tps.stratum = demo.stratum)
mytable(dat_proc$tps.stratum)

  
# Wstratum is used to compute sampling weights. 
# Note that a case should still have a Wstratum even if its tps.stratum is NA

dat_proc$Wstratum = dat_proc$tps.stratum

# Delta/ancestral/minor variants are case-sampling strata

cond=!is.na(dat_proc$KnownOrImputedDeltaCOVIDIndD57_7toD360) & dat_proc$KnownOrImputedDeltaCOVIDIndD57_7toD360==1
dat_proc$Wstratum[cond] = 100+dat_proc$tps.stratum[cond]

cond=!is.na(dat_proc$MinorVariantsCOVIDIndD57_7toD360) & dat_proc$MinorVariantsCOVIDIndD57_7toD360==1
dat_proc$Wstratum[cond] = 300+dat_proc$tps.stratum[cond]

cond=!is.na(dat_proc$AncestralCOVIDIndD57_7toD360) & dat_proc$AncestralCOVIDIndD57_7toD360==1
dat_proc$Wstratum[cond] = 200+dat_proc$tps.stratum[cond]

# Severe cases form one stratum and are not stratified by demographics groups because we sample all severe cases
# severe has to come last, 
dat_proc$Wstratum[!is.na(dat_proc$SevereCOVIDIndD57_7toD360) & dat_proc$SevereCOVIDIndD57_7toD360==1] = 401

mytable(dat_proc$Wstratum)

# merge 105 (KnownOrImputedDelta, Non-US Age ≥ 65) with 101 (KnownOrImputedDelta, US Age ≥ 65 URM)
# this ptid has unknown lineage, which we assume to be Delta, which is why the ptid was not sampled at design time
dat_proc$Wstratum[dat_proc$Wstratum==105]=101


###############################################################################
# 4. Define ph1, ph2, and weights
# Note that Wstratum may have NA if any variables to form strata has NA

tp='57_120'

dat_proc[["ph1.D"%.%tp]] = with(dat_proc, 
  Perprotocol==1 & 
    get("AnyInfectionD1toD"%.%tp)==0 & # no evidence of any infection by D57_120
    COVIDTimeD57to21Dec10 >= 120 & # COVID time or censor time is after D57_120
  # either Delta COVID, severe COVID, or no evidence of infection till 22Mar26
  (KnownOrImputedDeltaCOVIDIndD57_7toD360==1 | SevereCOVIDIndD57_7toD360 | AnyInfectionD1toD57_7==0)
)

# two separate indicators, one for bAb and one for nAb

# nAb, we only use Delta to impute ancestral
dat_proc[["TwophasesampIndD57nAb"]] = complete.cases(dat_proc[,c("Day57pseudoneutid50_Delta")])      
dat_proc[["ph2.D"%.%tp%.%"nAb"]] = dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD57nAb"]]
dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp, ph2="ph2.D"%.%tp%.%"nAb", Wstratum="Wstratum", wt="wt.D"%.%tp%.%"nAb", verbose=T) 

# bAb
dat_proc[["TwophasesampIndD57bAb"]] = complete.cases(dat_proc[,c("Day57bindSpike_Alpha")])      
dat_proc[["ph2.D"%.%tp%.%"bAb"]] = dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD57bAb"]]
dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp, ph2="ph2.D"%.%tp%.%"bAb", Wstratum="Wstratum", wt="wt.D"%.%tp%.%"bAb", verbose=T) 

mytable(dat_proc$TwophasesampIndD57nAb, dat_proc$TwophasesampIndD57bAb)
mytable(dat_proc$TwophasesampIndD57nAb, dat_proc$TwophasesampIndD57bAb, dat_proc$AnyInfectionD1)


###############################################################################
# 5. impute missing biomarkers in ph2 (assay imputation)
#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, each time point, but not Delta

n.imp <- 1

# impute nAb
dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD57nAb") == 1)
imp.markers=paste0("Day57", c("pseudoneutid50_D614G", "pseudoneutid50_Delta"))

trt=1; sero=0

imp <- dat.tmp.impute %>% dplyr::filter(Trt == trt & Bserostatus==sero) %>% select(all_of(imp.markers))         
if(any(is.na(imp))) {
  # if there is no variability, fill in NA with constant values
  for (a in names(imp)) {
    if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
  }            
  # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  dat.tmp.impute[dat.tmp.impute$Trt == trt & dat.tmp.impute$Bserostatus == sero , imp.markers] <- mice::complete(imp, action = 1)
}                

# missing markers imputed properly?
assertthat::assert_that(
  all(complete.cases(dat.tmp.impute[, imp.markers])),
  msg = "missing markers imputed properly?"
)    

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc[["TwophasesampIndD57nAb"]]==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["TwophasesampIndD57nAb"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]

assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD57nAb"]] == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)


  
  
###############################################################################
# 6. transformation of the markers if necessary
# e.g., convert binding variables from AU to IU for binding assays



###############################################################################
# 7. add mdw scores



###############################################################################
# 8. add fold change markers
# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta



###############################################################################
# 9. add discrete/trichotomized markers

bAb.assays = subset(assay_metadata, panel=="bindSpike", assay, drop=T)
dat_proc$tmp = with(dat_proc, Trt==1 & Bserostatus==0 & get("ph2.D57_120bAb")) 
dat_proc = add.trichotomized.markers (dat_proc, c("Day57"%.%bAb.assays), ph2.col.name="tmp", wt.col.name="wt.D57_120bAb")
dat_proc$tmp = NULL
  
nAb.assays = subset(assay_metadata, panel=="id50", assay, drop=T)
dat_proc$tmp = with(dat_proc, Trt==1 & Bserostatus==0 & get("ph2.D57_120nAb")) 
dat_proc = add.trichotomized.markers (dat_proc, c("Day57"%.%nAb.assays), ph2.col.name="tmp", wt.col.name="wt.D57_120nAb")
dat_proc$tmp = NULL




###############################################################################
# 10. impute covariates if necessary

# some ptids have missing risk score
summary (subset(dat_proc, Trt==1 & Bserostatus==0, risk_score))


n.imp <- 1
dat.tmp.impute <- dat_proc

imp.markers=c("risk_score", "Age", "Sex", "BMI", "URMforsubcohortsampling", "HighRiskInd", 
              "USAInd", "HIVinfection", "CalendarGrp", "Bserostatus")
  


imp <- dat.tmp.impute %>%  select(all_of(imp.markers))         
if(any(is.na(imp))) {
    # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
    imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
    dat.tmp.impute[, imp.markers] <- mice::complete(imp, action = 1)
}                

# missing markers imputed properly?
assertthat::assert_that(
    all(complete.cases(dat.tmp.impute[, imp.markers])),
    msg = "missing covariates imputed properly?"
)    

# populate dat_proc imp.markers with the imputed values
dat_proc[, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[, "Ptid"], dat.tmp.impute$Ptid), ]

# imputed values of missing markers merged properly for all individuals in the two phase sample?
assertthat::assert_that(
  all(complete.cases(dat_proc[, imp.markers])),
  msg = "imputed values of missing covariates merged properly for all individuals?"
)




###############################################################################
# special handling
# e.g., merging datasets




###############################################################################
# digest check
###############################################################################

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         azd1222_stage2 = "ac88cfcaebec995068d5278f69a747d0",
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, 
      msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))



print("run time: "%.%format(Sys.time()-begin, digits=1))
