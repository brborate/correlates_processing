#Sys.setenv(TRIAL = "moderna_real")
#Sys.setenv(TRIAL = "janssen_partA_VL")
#Sys.setenv(TRIAL = "vat08_combined")
#Sys.setenv(TRIAL = "id27hpv")

# no need to run renv::activate(here::here()) b/c .Rprofile exists


# vat08_combined: Stage is incorporated into demo.stratum

source(here::here("_common.R"))

if (TRIAL=="moderna_boost") stop("For moderna_boost, run make_dat_moderna_boost.R") 

library(tidyverse)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)
library(here)
library(mdw)

begin=Sys.time()


########################################################################################################
# read mapped data with risk score added

if (make_riskscore) {
  if (TRIAL=="janssen_partA_VL") {
    # read hot deck data
    tmp = sub(".csv","_hotdeck.csv",mapped_data)
    if (file.exists(tmp)) dat_raw=read.csv(tmp) else stop("hotdeck file not exists, run hotdeck R script first")
    dat_proc = preprocess(dat_raw, study_name)   
    colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
    
    # borrow risk score from janssen_pooled_partA
    load(file = 'riskscore_baseline/output/janssen_pooled_partA/inputFile_with_riskscore.RData')
    stopifnot(all(dat_proc$Ptid==inputFile_with_riskscore$Ptid))
    dat_proc$risk_score = inputFile_with_riskscore$risk_score
    dat_proc$standardized_risk_score = inputFile_with_riskscore$standardized_risk_score
    
    # create a new event indicator variable that censors cases without VL, which also include all non-molec confirmed cases
    # required for defining weights
    dat_proc$EventIndPrimaryHasVLD29 = dat_proc$EventIndPrimaryIncludeNotMolecConfirmedD29
    dat_proc$EventIndPrimaryHasVLD29[is.na(dat_proc$seq1.log10vl) & !is.na(dat_proc$EventIndPrimaryHasVLD29)] = 0 # if EventIndPrimaryHasVLD29 is NA, this will remain NA
    dat_proc$EventIndPrimaryD29 = dat_proc$EventIndPrimaryHasVLD29
    # EventTimePrimaryD29 is set to EventIndPrimaryIncludeNotMolecConfirmedD29 in preprocess()
    
   } else if (TRIAL=="vat08_combined") {
     # read hot deck data
     tmp = sub(".csv","_hotdeck.csv",mapped_data)
     if (file.exists(tmp)) dat_raw=read.csv(tmp) else stop("hotdeck file not exists, run hotdeck R script first")
     
     dat_proc = preprocess(dat_raw, study_name)   
     colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
     
     # add risk score
     load(file = paste0('riskscore_baseline/output/',TRIAL,'/inputFile_with_riskscore.RData'))
     stopifnot(all(inputFile_with_riskscore$Ptid==dat_proc$Ptid))
     dat_proc$risk_score = inputFile_with_riskscore$risk_score
     dat_proc$standardized_risk_score = inputFile_with_riskscore$standardized_risk_score
     
     # ptids with missing Bserostatus already filtered out in preprocess
     
     # define event indicator and event time variables based on seq1.variant.hotdeck1 etc
     for (t in c(1,22,43)) {
       for (i in 1:10) {
         dat_proc[[paste0("EventIndOmicronD",t,"hotdeck",i)]] = 
           ifelse(!is.na(dat_proc[["seq1.variant.hotdeck"%.%i]]) & dat_proc[["seq1.variant.hotdeck"%.%i]]=="Omicron" & !is.na(dat_proc[["EventIndFirstInfectionD"%.%t]]),
                  1,
                  0)
         dat_proc[[paste0("EventTimeOmicronD",t,"hotdeck",i)]] = ifelse(dat_proc[[paste0("EventIndOmicronD",t,"hotdeck",i)]] ==1, 
                                                                        min(dat_proc[["EventTimeKnownLineageOmicronD"%.%t]],    dat_proc[["EventTimeMissingLineageD"%.%t]]),
                                                                        max(dat_proc[["EventTimeKnownLineageNonOmicronD"%.%t]], dat_proc[["EventTimeMissingLineageD"%.%t]]))
       }
     }
     
     

    } else {
      # load inputFile_with_riskscore.Rdata, a product of make riskscore_analysis, which calls preprocess and makes risk scores
      load(file = paste0('riskscore_baseline/output/',TRIAL,'/inputFile_with_riskscore.RData'))
      dat_proc <- inputFile_with_riskscore    
    }
  
} else {
    dat_raw=read.csv(mapped_data)
    dat_proc = preprocess(dat_raw, study_name)   
    colnames(dat_proc)[1] <- "Ptid" 
}

# define new variables
dat_proc <- dat_proc %>% mutate(age.geq.65 = as.integer(Age >= 65))
dat_proc$Senior = as.integer(dat_proc$Age>=switch(study_name, COVE=65, MockCOVE=65, ENSEMBLE=60, MockENSEMBLE=60, PREVENT19=65, AZD1222=65, VAT08=60, PROFISCOV=NA, IARC_HPV=15, stop("unknown study_name 1")))
  
# for the mock datasets, hardcode AnyinfectionD1 
if (study_name %in% c("MockENSEMBLE", "MockCOVE")) dat_proc$AnyinfectionD1=0

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
      
} else if (study_name %in% c("PROFISCOV", "VAT08")) {
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
if (study_name %in% c("COVE", "MockCOVE", "PROFISCOV", "VAT08")) {
    # nothing to do
    # COVE only has US data
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    dat_proc$MinorityInd[dat_proc$Region!=0] = 0
    
} else if (study_name=="PREVENT19") {
    dat_proc$MinorityInd[dat_proc$Country!=0] = 0 # 0 is US
    
} else if (study_name=="AZD1222") {
    dat_proc$MinorityInd[dat_proc$Country!=2] = 0 # 2 is US
    
} else stop("unknown study_name 3")



# check coding via tables
#table(dat_proc$race, useNA = "ifany")
#table(dat_proc$WhiteNonHispanic, useNA = "ifany")
#table(dat_proc$race, dat_proc$WhiteNonHispanic, useNA = "ifany")


## URMforsubcohortsampling is defined in data_raw, the following is a check
#dat_proc$URMforsubcohortsampling <- NA
#dat_proc$URMforsubcohortsampling[dat_proc$Black==1] <- 1
#dat_proc$URMforsubcohortsampling[dat_proc$ethnicity=="Hispanic or Latino"] <- 1
#dat_proc$URMforsubcohortsampling[dat_proc$NatAmer==1] <- 1
#dat_proc$URMforsubcohortsampling[dat_proc$PacIsl==1] <- 1
#dat_proc$URMforsubcohortsampling[dat_proc$Asian==1 & dat_proc$EthnicityHispanic==0 & dat_proc$EthnicityUnknown==0 & dat_proc$EthnicityNotreported==0] <- 0
#dat_proc$URMforsubcohortsampling[dat_proc$Multiracial==1 & dat_proc$EthnicityHispanic==0 & dat_proc$EthnicityUnknown==0 & dat_proc$EthnicityNotreported==0] <- 0
#dat_proc$URMforsubcohortsampling[dat_proc$Other==1 & dat_proc$EthnicityHispanic==0 & dat_proc$EthnicityUnknown==0 & dat_proc$EthnicityNotreported==0] <- 0
## Add observed White Non Hispanic:
#dat_proc$URMforsubcohortsampling[dat_proc$EthnicityHispanic==0 & dat_proc$EthnicityUnknown==0 & dat_proc$EthnicityNotreported==0 & dat_proc$Black==0 & dat_proc$Asian==0 & dat_proc$NatAmer==0 & dat_proc$PacIsl==0 &
#dat_proc$Multiracial==0 & dat_proc$Other==0 & dat_proc$Notreported==0 & dat_proc$Unknown==0] <- 0
#
#
## nonURM=1 IF race is White AND ethnicity is not Hispanic
#dat_proc$nonURM <- NA
#dat_proc$nonURM <-
#  ifelse(dat_proc$race %in% c("White","Asian","Other","Multiracial") &
#    dat_proc$ethnicity == "Not Hispanic or Latino", 1,
#  dat_proc$nonURM
#  )
## nonURM=0 IF race is not "white or unknown" OR ethnicity is Hispanic
#dat_proc$nonURM <-
#  ifelse(!dat_proc$race %in% c("White","Asian","Other","Multiracial","Not reported and unknown") |
#    dat_proc$ethnicity == "Hispanic or Latino", 0,
#    dat_proc$nonURM
#  )
#dat_proc$URM.2 = 1-dat_proc$nonURM
#
#with(dat_proc, table(URMforsubcohortsampling, URM.2, useNA="ifany"))




###############################################################################
# stratum variables
# The code for Bstratum is trial specifc
# The code for tps.stratum and Wstratum are not trial specific since they are constructed on top of Bstratum
###############################################################################

# Bstratum: randomization strata
# Moderna: 1 ~ 3, defines the 3 baseline strata within trt/serostatus
if (study_name=="COVE" | study_name=="MockCOVE" ) {
    dat_proc$Bstratum = with(dat_proc, ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3)))
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
    dat_proc$Bstratum =  with(dat_proc, strtoi(paste0(Senior, HighRiskInd), base = 2)) + 1
    
} else if (study_name %in% c("PREVENT19", "AZD1222", "VAT08")) {
    dat_proc$Bstratum = with(dat_proc, ifelse(Senior, 2, 1))
    
} else if (study_name %in% c("PROFISCOV")) {
    dat_proc$Bstratum = 1 # there are no demographics stratum for subcohort sampling
    
} else stop("unknown study_name 4")

names(Bstratum.labels) <- Bstratum.labels

#with(dat_proc, table(Bstratum, Senior, HighRiskInd))


# demo.stratum: correlates sampling strata
# Moderna: 1 ~ 6 defines the 6 baseline strata within trt/serostatus
# may have NA b/c URMforsubcohortsampling may be NA
if (study_name=="COVE" | study_name=="MockCOVE" ) {
    dat_proc$demo.stratum = with(dat_proc, ifelse (URMforsubcohortsampling==1, ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3)), 3+ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3))))
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
    # first step, stratify by age and high risk
    dat_proc$demo.stratum =  with(dat_proc, strtoi(paste0(Senior, HighRiskInd), base = 2)) + 1
    # second step, stratify by country
    dat_proc$demo.stratum=with(dat_proc, ifelse(Region==0 & URMforsubcohortsampling==0, demo.stratum + 4, demo.stratum)) # US, non-URM
    dat_proc$demo.stratum[dat_proc$Region==1] = dat_proc$demo.stratum[dat_proc$Region==1] + 8 # Latin America
    dat_proc$demo.stratum[dat_proc$Region==2] = dat_proc$demo.stratum[dat_proc$Region==2] + 12 # Southern Africa
    # the above sequence ends up setting US URM=NA to NA
    
    assertthat::assert_that(
        all(!with(dat_proc, xor(is.na(demo.stratum),  Region==0 & is.na(URMforsubcohortsampling) ))),
        msg = "demo.stratum is na if and only if URM is NA and north america")
    
} else if (study_name=="PREVENT19" ) {
    dat_proc$demo.stratum = with(dat_proc, strtoi(paste0(URMforsubcohortsampling, Senior, HighRiskInd), base = 2)) + 1
    dat_proc$demo.stratum = with(dat_proc, ifelse(Country==0, demo.stratum, ifelse(!Senior, 9, 10))) # 0 is US
        
} else if (study_name=="AZD1222" ) {
#    US, <65, non-Minority
#    US, >65, non-Minority
#    US, <65, Minority
#    US, >65, Minority
#    non-US, <65
#    non-US, >65
    dat_proc$demo.stratum = with(dat_proc, strtoi(paste0(URMforsubcohortsampling, Senior), base = 2)) + 1
    dat_proc$demo.stratum = with(dat_proc, ifelse(Country==2, demo.stratum, ifelse(!Senior, 5, 6))) # 2 is US
        
} else if (study_name=="VAT08" ) {
#    Stage 1, Not HND, Not senior
#    Stage 1, Not HND, senior
#    Stage 1, HND,     Not senior
#    Stage 1, HND,     senior
#    Stage 2, Not senior
#    Stage 2, senior
    dat_proc$demo.stratum = dat_proc$Bstratum
    dat_proc$demo.stratum = with(dat_proc, ifelse(Trialstage==1 & Country==3, demo.stratum+2, demo.stratum)) # Stage 1, HND
    dat_proc$demo.stratum = with(dat_proc, ifelse(Trialstage==2, demo.stratum+4, demo.stratum)) # Stage 2
    
} else if (study_name=="PROFISCOV" ) {
    dat_proc$demo.stratum = 1 # # there are no demographics stratum for subcohort sampling
    
} else stop("unknown study_name 5")  
  
# names(demo.stratum.labels) <- demo.stratum.labels

with(dat_proc, table(demo.stratum))

# tps stratum, 1 ~ 4*max(demo.stratum), used in tps regression
dat_proc <- dat_proc %>%
  mutate(
    tps.stratum = demo.stratum + strtoi(paste0(Trt, Bserostatus), base = 2) * max(demo.stratum,na.rm=T)
  )

# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+4. 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# A case will have a Wstratum even if its tps.stratum is NA
# The case is defined using EventIndPrimaryD29

max.tps=max(dat_proc$tps.stratum,na.rm=T)
dat_proc$Wstratum = dat_proc$tps.stratum
tps.cnt=max.tps+1
if(TRIAL %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA",
                                 "janssen_partA_VL")) {
  # cases sampling weights are also conditional on region and age group
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==0 & Region==0 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==0 & Region==0 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==0 & Region==1 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==0 & Region==1 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==0 & Region==2 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==0 & Region==2 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==1 & Region==0 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==1 & Region==0 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==1 & Region==1 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==1 & Region==1 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==1 & Region==2 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==1 & Region==2 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==0 & Region==0 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==0 & Region==0 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==0 & Region==1 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==0 & Region==1 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==0 & Region==2 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==0 & Region==2 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==1 & Region==0 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==1 & Region==0 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==1 & Region==1 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==1 & Region==1 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==1 & Region==2 & Senior==0)]=tps.cnt; tps.cnt=tps.cnt+1
  dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==1 & Region==2 & Senior==1)]=tps.cnt; tps.cnt=tps.cnt+1
  
} else if (study_name %in% c("COVE", "MockCOVE", "ENSEMBLE", "MockENSEMBLE", "AZD1222")) {
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==0)]=max.tps+1
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==1)]=max.tps+2
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==0)]=max.tps+3
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==1)]=max.tps+4
    
} else if (study_name == "PREVENT19") {
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD21==1 & Trt==0 & Bserostatus==0)]=max.tps+1
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD21==1 & Trt==0 & Bserostatus==1)]=max.tps+2
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD21==1 & Trt==1 & Bserostatus==0)]=max.tps+3
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD21==1 & Trt==1 & Bserostatus==1)]=max.tps+4
    
} else if (study_name == "VAT08") {
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==0)]=max.tps+1
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==1)]=max.tps+2
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==0)]=max.tps+3
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==1)]=max.tps+4
    
    # collapse Senior and non-Senior in the naive for both trial stage 1 and 2 due to empty cells in senior in 
    # Stage 2, naive
    # Stage 1, nnaive, Honduras
    dat_proc$Wstratum[dat_proc$Wstratum %in% c(17,18)]=17
    dat_proc$Wstratum[dat_proc$Wstratum %in% c(21,22)]=21
    
} else if (study_name == "PROFISCOV") {
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD43==1 & Trt==0 & Bserostatus==0)]=max.tps+1
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD43==1 & Trt==0 & Bserostatus==1)]=max.tps+2
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD43==1 & Trt==1 & Bserostatus==0)]=max.tps+3
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD43==1 & Trt==1 & Bserostatus==1)]=max.tps+4
    
} else stop("unknown study_name 6")


#subset(dat_proc, Trt==1 & Bserostatus==1 & EventIndPrimaryD29 == 1)[1:3,]

with(dat_proc, table(tps.stratum))
with(dat_proc, table(Wstratum))




###############################################################################
# observation-level weights
# Note that Wstratum may have NA if any variables to form strata has NA
###############################################################################


# define must_have_assays for ph2 definition
if (study_name %in% c("COVE", "MockCOVE")) {
    must_have_assays <- c("bindSpike", "bindRBD")
    
} else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")) {
    if (endsWith(TRIAL,"ADCP")) {
        must_have_assays <- c("ADCP")    
    } else {
        must_have_assays <- c("bindSpike", "bindRBD")
    }    
    
} else if (study_name %in% c("PREVENT19")) {
    must_have_assays <- c("bindSpike")
    
} else if (study_name %in% c("AZD1222")) {
    if (TRIAL=="azd1222") {
        must_have_assays <- c("pseudoneutid50")
        
    } else if (TRIAL=="azd1222_bAb") {
        must_have_assays <- c("bindSpike")
        
    } else stop("need to define must_have_assays")
    
} else if (study_name %in% c("PROFISCOV")) {
    if (TRIAL=="profiscov") {
        must_have_assays <- c("bindSpike")
        
    } else if (TRIAL=="profiscov_lvmn") {
        must_have_assays <- c("bindSpike")
        
    } else stop("need to define must_have_assays")
    
} else if (study_name %in% c("VAT08")) {
    must_have_assays <- c("bindSpike")
        
} else stop("unknown study_name 7")


# examine missingness
# # sanofi 
# with(subset(dat_proc, Trt==1), table(!is.na(BbindSpike)))
# with(subset(dat_proc, Trt==1), table(!is.na(Day22bindSpike)))
# with(subset(dat_proc, Trt==1), table(!is.na(Day43bindSpike)))
# with(subset(dat_proc, Trt==1), table(!is.na(BbindSpike) & !is.na(Day22bindSpike) & !is.na(Day43bindSpike)))
# with(subset(dat_proc, Trt==1), table(!is.na(Day22bindSpike), !is.na(Day43bindSpike), !is.na(BbindSpike) ))
# with(subset(dat_proc, Trt==1 & EventIndPrimaryD22), table(!is.na(Day22bindSpike), !is.na(Day43bindSpike), !is.na(BbindSpike) ))
# with(subset(dat_proc, Trt==1 & EventIndPrimaryD22), table(!is.na(Day22bindSpike), !is.na(Day43bindSpike), Stage ))



# TwophasesampInd: be in the subcohort or a case after time point 1  &  have the necessary markers
if (study_name %in% c("COVE", "MockCOVE", "MockENSEMBLE", "PREVENT19")) {
    if (two_marker_timepoints) {
    # require baseline and timpoint 1
        dat_proc[["TwophasesampIndD"%.%timepoints[2]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("B"%.%must_have_assays, "Day"%.%timepoints[1]%.%must_have_assays, "Day"%.%timepoints[2]%.%must_have_assays)])      
    }
    # require baseline
    dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("B"%.%must_have_assays, "Day"%.%timepoints[1]%.%must_have_assays)])      
        
} else if (study_name %in% c("VAT08")) {
  # require baseline but not time point 1
  dat_proc[["TwophasesampIndD"%.%timepoints[2]]] = 
    with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
    complete.cases(dat_proc[,c("B"%.%must_have_assays, "Day"%.%timepoints[2]%.%must_have_assays)])      

  # same as tp2
  dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = dat_proc[["TwophasesampIndD"%.%timepoints[2]]]

  
} else if (study_name %in% c("AZD1222", "PROFISCOV")) {
    if (two_marker_timepoints) {
        # does not require baseline or time point 1
        dat_proc[["TwophasesampIndD"%.%timepoints[2]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("Day"%.%timepoints[2]%.%must_have_assays)])        
        
        # does not require baseline
        dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            # adding | is because if D57 is present, D29 will be imputed if missing
            (complete.cases(dat_proc[,c("Day"%.%timepoints[1]%.%must_have_assays)]) | complete.cases(dat_proc[,c("Day"%.%timepoints[2]%.%must_have_assays)])) 
    } else {
        # does not require baseline
        dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            (complete.cases(dat_proc[,c("Day"%.%timepoints[1]%.%must_have_assays)]) ) 
    }
        
} else if (study_name=="ENSEMBLE") {
    if (contain(TRIAL, "EUA")) {
        # require baseline
        dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("B"%.%must_have_assays, "Day"%.%timepoints[1]%.%must_have_assays)])      
        
    } else if (contain(TRIAL, "partA")) {
        # does not require baseline
        dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("Day"%.%timepoints[1]%.%must_have_assays)])
        
        # define TwophasesampIndD29variant, which is same as TwophasesampIndD29 except that for the non-cases, it is limited to a subset with variant ID50
        # need this anymore b/c even though we impute, for example, Beta ID50 for non-Beta cases in RSA, not all non-cases have variant ID50s
        if (TRIAL=="janssen_partA_VL"){
          # an indicator for non-cases in the cohort
          select = with(dat_proc, SubcohortInd & (is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0))
          # an indicator for being ancestral in US, beta in RSA, and G/L/M/Z in LatAm
          tmp = complete.cases(dat_proc[,c("Day"%.%timepoints[1]%.%must_have_assays)]) &
                with(dat_proc, Region==0 & !is.na(Day29pseudoneutid50) |
                               Region==2 & !is.na(Day29pseudoneutid50_Beta) |
                               Region==1 & (!is.na(Day29pseudoneutid50_Gamma) | !is.na(Day29pseudoneutid50_Lambda) | !is.na(Day29pseudoneutid50_Mu) | !is.na(Day29pseudoneutid50_Zeta)) )
          dat_proc$TwophasesampIndD29variant = dat_proc$TwophasesampIndD29 
          dat_proc$TwophasesampIndD29variant[select] = tmp [select]
        }
        
    }
    
} else stop("unknown study_name 8")



# weights 
for (tp in rev(timepoints)) { # rev is just so that digest passes
    tmp = with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp) >= 7)
    wts_table <- with(dat_proc[tmp,], table(Wstratum, get("TwophasesampIndD"%.%tp)))
    wts_norm <- rowSums(wts_table) / wts_table[, 2]
    dat_proc[["wt.D"%.%tp]] <- wts_norm[dat_proc$Wstratum %.% ""]
    # the step above assigns weights for some subjects outside ph1. the next step makes them NA
    dat_proc[["wt.D"%.%tp]] = ifelse(with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp)>=7), dat_proc[["wt.D"%.%tp]], NA) 
    dat_proc[["ph1.D"%.%tp]]=!is.na(dat_proc[["wt.D"%.%tp]])
    dat_proc[["ph2.D"%.%tp]]=dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD"%.%tp]]
    
    assertthat::assert_that(
        all(!is.na(subset(dat_proc, tmp & !is.na(Wstratum))[["wt.D"%.%tp]])),
        msg = "missing wt.D for D analyses ph1 subjects")
}

## additional weights

if (TRIAL=="janssen_partA_VL"){
  # weights that depend on TwophasesampIndD29variant and Wstratum
  
  tmp = with(dat_proc, get("EarlyendpointD29")==0 & Perprotocol==1 & get("EventTimePrimaryD29") >= 7)
  wts_table <- with(dat_proc[tmp,], table(Wstratum, TwophasesampIndD29variant))
  wts_norm <- rowSums(wts_table) / wts_table[, 2]
  dat_proc[["wt.D29variant"]] <- wts_norm[dat_proc$Wstratum %.% ""]
  # the step above assigns weights for some subjects outside ph1. the next step makes them NA
  dat_proc[["wt.D29variant"]] = ifelse(with(dat_proc, get("EarlyendpointD29")==0 & Perprotocol==1 & get("EventTimePrimaryD29")>=7), dat_proc[["wt.D29variant"]], NA) 
  #no need to define ph1.D29variant since it is identical to ph1.D29
  dat_proc[["ph2.D29variant"]] = dat_proc[["ph1.D29"]] & dat_proc$TwophasesampIndD29variant
  
  assertthat::assert_that(
    all(!is.na(subset(dat_proc, tmp & !is.na(Wstratum))[["wt.D29variant"]])),
    msg = "missing wt.D for D analyses ph1 subjects")
  
  # RSA
  # non-cases
  # everyone who has Beta also has ancestral
  with(subset(dat_proc,Trt==1 & ph1.D29 & Region==2 & !EventIndPrimaryHasVLD29), table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Beta)))
  
  # LatAm
  # non-cases
  # everyone who has Gamma also has ancestral
  with(subset(dat_proc,Trt==1 & ph1.D29 & Region==1 & !EventIndPrimaryHasVLD29), table(!is.na(Day29pseudoneutid50), !is.na(Day29pseudoneutid50_Gamma)))
  # same set of ptids have Gamma etc
  with(subset(dat_proc,Trt==1 & ph1.D29 & Region==1 & !EventIndPrimaryHasVLD29), table(!is.na(Day29pseudoneutid50_Gamma), !is.na(Day29pseudoneutid50_Mu)))
  with(subset(dat_proc,Trt==1 & ph1.D29 & Region==1 & !EventIndPrimaryHasVLD29), table(!is.na(Day29pseudoneutid50_Gamma), !is.na(Day29pseudoneutid50_Zeta)))
  with(subset(dat_proc,Trt==1 & ph1.D29 & Region==1 & !EventIndPrimaryHasVLD29), table(!is.na(Day29pseudoneutid50_Gamma), !is.na(Day29pseudoneutid50_Lambda)))
  # cases
  with(subset(dat_proc,Trt==1 & ph1.D29 & Region==1), table(!is.na(Day29pseudoneutid50_Gamma), !is.na(Day29pseudoneutid50_Mu), EventIndPrimaryHasVLD29, useNA="ifany"))
  
}

# Starting at 1 day post D29 visit
if(study_name %in% c("ENSEMBLE", "MockENSEMBLE")) {
    wts_table2 <-  dat_proc %>% dplyr::filter(EarlyendpointD29start1==0 & Perprotocol==1 & EventTimePrimaryD29>=1) %>%
      with(table(Wstratum, TwophasesampIndD29))
    wts_norm2 <- rowSums(wts_table2) / wts_table2[, 2]
    dat_proc$wt.D29start1 <- wts_norm2[dat_proc$Wstratum %.% ""]
    dat_proc$wt.D29start1 = ifelse(with(dat_proc,  EarlyendpointD29start1==0 & Perprotocol==1 & EventTimePrimaryD29>=1), dat_proc$wt.D29start1, NA)
    dat_proc$ph1.D29start1=!is.na(dat_proc$wt.D29start1)
    dat_proc$ph2.D29start1=with(dat_proc, ph1.D29start1 & TwophasesampIndD29)
    
    assertthat::assert_that(
        all(!is.na(subset(dat_proc,           EarlyendpointD29start1==0 & Perprotocol==1 & EventTimePrimaryD29>=1 & !is.na(Wstratum), select=wt.D29start1, drop=T))),
        msg = "missing wt.D29start1 for D29start1 analyses ph1 subjects")
} 

# weights for intercurrent cases
if(two_marker_timepoints) {
    tp=timepoints[1]
    tmp = with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp) >= 7 & get("EventIndPrimaryD"%.%tp)==1 
        & get("EventTimePrimaryD"%.%tp) <= 6 + get("NumberdaysD1toD"%.%timepoints[2]) - get("NumberdaysD1toD"%.%tp))
    wts_table2 <- with(dat_proc[tmp,], table(Wstratum, get("TwophasesampIndD"%.%tp)))
    wts_norm2 <- rowSums(wts_table2) / wts_table2[, 2]
    dat_proc$wt.intercurrent.cases <- wts_norm2[dat_proc$Wstratum %.% ""]
    dat_proc$wt.intercurrent.cases = ifelse(tmp, 
                                            dat_proc$wt.intercurrent.cases, 
                                            NA)
    dat_proc$ph1.intercurrent.cases=!is.na(dat_proc$wt.intercurrent.cases)
    dat_proc$ph2.intercurrent.cases=with(dat_proc, ph1.intercurrent.cases & get("TwophasesampIndD"%.%tp))    
    
    assertthat::assert_that(
        all(!is.na(subset(dat_proc, tmp & !is.na(Wstratum), select=wt.intercurrent.cases, drop=T))),
        msg = "missing wt.intercurrent.cases for intercurrent analyses ph1 subjects")
}

# weights for immunogenicity analyses that use subcohort only and are not enriched by cases outside subcohort
tp=timepoints[ifelse(two_marker_timepoints, 2, 1)]
tmp = with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1)
wts_table <- with(dat_proc[tmp,], table(tps.stratum, get("TwophasesampIndD"%.%tp) & SubcohortInd))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc$wt.subcohort <- wts_norm[dat_proc$tps.stratum %.% ""]
dat_proc$wt.subcohort = ifelse(tmp, dat_proc$wt.subcohort, NA)
dat_proc$ph1.immuno=!is.na(dat_proc$wt.subcohort)
dat_proc$ph2.immuno=with(dat_proc, ph1.immuno & SubcohortInd & get("TwophasesampIndD"%.%tp))

assertthat::assert_that(
    all(!is.na(subset(dat_proc, tmp & !is.na(tps.stratum), select=wt.subcohort, drop=T))), 
    msg = "missing wt.subcohort for immuno analyses ph1 subjects")
        


## check missing risk score 
#with(subset(dat_proc, ph1.D35 & Trt==1 & Bserostatus==0), table(is.na(risk_score), Riskscorecohortflag, ph2.D35, EventIndPrimaryD35))
#with(subset(dat_proc, Country==0 & Bserostatus==0 & Trt==1 & ph1.D35 & !is.na(BbindSpike) & !is.na(Day35bindSpike)), table(is.na(Day35pseudoneutid50), EventIndPrimaryD35))
#with(subset(dat_proc, Country==0 & Bserostatus==0 & Trt==1 & ph2.D35), table(is.na(Day35pseudoneutid50), EventIndPrimaryD35))
#with(subset(dat_proc, Country==0 & Bserostatus==0 & Trt==1 & ph2.immuno & !is.na(BbindSpike) & !is.na(Day35bindSpike)), table(is.na(Day35pseudoneutid50), EventIndPrimaryD35))
#with(subset(dat_proc, Country==0 & Bserostatus==0 & Trt==1 & ph2.immuno), table(EventIndPrimaryD35, useNA="ifany"))


###############################################################################
# impute missing biomarkers in ph2 (assay imputation)
#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, each time point, but not Delta
###############################################################################

# loop through the time points
# first impute (B, D29, D57) among TwophasesampIndD57==1
# next impute (B, D29) among TwophasesampIndD29==1
for (tp in rev(timepoints)) {    
    n.imp <- 1
    dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD"%.%tp) == 1)
    
    if (TRIAL=="janssen_partA_VL") {
      # assays include many more assays than what are needed here
      imp.markers=c(outer(c("B", "Day"%.%tp), c("bindSpike", "bindRBD", "pseudoneutid50"), "%.%"))
      
    } else {
      if(two_marker_timepoints) {
        imp.markers=c(outer(c("B", if(tp==timepoints[2]) "Day"%.%timepoints else "Day"%.%tp), assays, "%.%"))
      } else {
        imp.markers=c(outer(c("B", "Day"%.%tp), assays, "%.%"))
      }
      # mdw markers are not imputed
      imp.markers=imp.markers[!endsWith(imp.markers, "_mdw")]
    }
        
    for (trt in unique(dat_proc$Trt)) {
    for (sero in unique(dat_proc$Bserostatus)) {    
        #summary(subset(dat.tmp.impute, Trt == 1 & Bserostatus==0)[imp.markers])      
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
    }
    }
    
    # missing markers imputed properly?
    assertthat::assert_that(
        all(complete.cases(dat.tmp.impute[, imp.markers])),
        msg = "missing markers imputed properly?"
    )    
    
    # populate dat_proc imp.markers with the imputed values
    dat_proc[dat_proc[["TwophasesampIndD"%.%tp]]==1, imp.markers] <-
      dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["TwophasesampIndD"%.%tp]]==1, "Ptid"], dat.tmp.impute$Ptid), ]
    
    assertthat::assert_that(
      all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD"%.%tp]] == 1, imp.markers])),
      msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
    )
}


# for janssen_partA_VL, impute 10 copies of ID50 variant markers
if (TRIAL=="janssen_partA_VL") {
  
  n.imp <- 10
  
  # RSA, impute Beta for non-Beta cases
  select = with(dat_proc, Trt == 1 & Bserostatus==0 & Region==2 & TwophasesampIndD29variant==1, drop=T)
  imp.markers = c("Day29bindSpike", "Day29bindRBD", "Day29pseudoneutid50", "Day29pseudoneutid50_Beta")
  imp <- dat_proc[select,] %>% select(all_of(imp.markers)) %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  # add 10 new columns for each of the variants to the dataset
  for (i in 1:10) {
    dat_proc[["Day29pseudoneutid50_Beta"%.%i]] = NA
    dat_proc[select, "Day29pseudoneutid50_Beta"%.%i] = mice::complete(imp, action=i)[,4]
  } 
  assertthat::assert_that(
    all(complete.cases(dat_proc[select, "Day29pseudoneutid50_Beta"%.%1:10])),
    msg = "janssen_partA_VL: imputed values of missing markers merged properly for all individuals in the two phase sample?"
  )
  
  # LatAm, impute four variants
  select = with(dat_proc, Trt == 1 & Bserostatus==0 & Region==1 & TwophasesampIndD29variant==1, drop=T)
  imp.markers = c("Day29bindSpike", "Day29bindRBD", "Day29pseudoneutid50", "Day29pseudoneutid50_Gamma", "Day29pseudoneutid50_Lambda", "Day29pseudoneutid50_Mu", "Day29pseudoneutid50_Zeta")
  imp <- dat_proc[select, ] %>% select(all_of(imp.markers)) %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  # add 10 new columns for each of the variants to the dataset
  for (i in 1:10) {
    dat_proc[["Day29pseudoneutid50_Gamma"%.%i]]  = NA
    dat_proc[["Day29pseudoneutid50_Lambda"%.%i]] = NA
    dat_proc[["Day29pseudoneutid50_Mu"%.%i]]     = NA
    dat_proc[["Day29pseudoneutid50_Zeta"%.%i]]   = NA
    dat_proc[select, "Day29pseudoneutid50_Gamma"%.%i]  = mice::complete(imp, action=i)[,4]
    dat_proc[select, "Day29pseudoneutid50_Lambda"%.%i] = mice::complete(imp, action=i)[,5]
    dat_proc[select, "Day29pseudoneutid50_Mu"%.%i]     = mice::complete(imp, action=i)[,6]
    dat_proc[select, "Day29pseudoneutid50_Zeta"%.%i]   = mice::complete(imp, action=i)[,7]
  } 
  assertthat::assert_that(
    all(complete.cases(dat_proc[select, c("Day29pseudoneutid50_Gamma"%.%1:10, "Day29pseudoneutid50_Lambda"%.%1:10, "Day29pseudoneutid50_Mu"%.%1:10, "Day29pseudoneutid50_Zeta"%.%1:10)])),
    msg = "janssen_partA_VL: imputed values of missing markers merged properly for all individuals in the two phase sample?"
  )
  
}  
  
  
###############################################################################
# transformation of the markers
###############################################################################


# converting binding variables from AU to IU for binding assays
# COVE only 
# moderna_real immune.csv file is not on international scale and other mapped data files are

# some trials have N some don't
if (is.null(config$assay_metadata)) {
  includeN = switch(study_name, COVE=1, MockCOVE=1, ENSEMBLE=1, MockENSEMBLE=1, PREVENT19=0, AZD1222=0, VAT08=0, PROFISCOV=1, stop("unknown study_name 9"))
  assays.includeN=c(assays, if(includeN==1) "bindN")
} else {
  if (TRIAL=="janssen_partA_VL") {
    assays.includeN = c("bindSpike", "bindRBD", "pseudoneutid50", "bindN")
  } else {
    assays.includeN = assays
  }
}


if(study_name=="COVE"){
    # conversion is only done for COVE for backward compatibility
    convf=c(bindSpike=0.0090, bindRBD=0.0272, bindN=0.0024, pseudoneutid50=0.242, pseudoneutid80=1.502)    
    for (a in assays.includeN) {
      for (t in c("B", paste0("Day", config$timepoints)) ) {
          dat_proc[[t %.% a]] <- dat_proc[[t %.% a]] + log10(convf[a])
      }
    }
}


# censoring values below LLOD
# COVE and mock only
# after COVE, the mapped data comes censored

if(study_name %in% c("COVE", "MockCOVE")){
    for (a in assays.includeN) {
      for (t in c("B", paste0("Day", config$timepoints)) ) {
        dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(llods[a]), log10(llods[a] / 2), dat_proc[[t %.% a]])
      }
    }
} else if(study_name %in% c("MockENSEMBLE")){
    for (a in assays.includeN) {
      for (t in c("B", paste0("Day", config$timepoints)) ) {
        dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(lloqs[a]), log10(lloqs[a] / 2), dat_proc[[t %.% a]])
      }
    }
}



# define delta for dat_proc
# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta

# mdw scores delta are computed as weighted aveage of delta, not as delta of mdw
assays1 = assays.includeN[!endsWith(assays.includeN, "_mdw")]

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




# define mdw scores 
if(TRIAL == "vat08_combined") {

  for (iAna in 1:4) {
    if (iAna==1) kp = dat_proc$Trialstage==1 & dat_proc$Bserostatus==1
    if (iAna==2) kp = dat_proc$Trialstage==2 & dat_proc$Bserostatus==0
    if (iAna==3) kp = dat_proc$Trialstage==2 & dat_proc$Bserostatus==1
    if (iAna==4) kp = dat_proc$Trialstage==1 & dat_proc$Bserostatus==0 # not used for correlates, but do it for completeness
    
    bAb = setdiff(assays[startsWith(assays, "bindSpike")], c('bindSpike_mdw'))
    nAb = setdiff(assays[startsWith(assays, "pseudoneutid50")], c('pseudoneutid50_mdw'))
    
    for (t in c("B", "Day"%.%timepoints, "Delta"%.%timepoints%.%"overB", "Delta43over22")) {
      # bAb
      # mdw weights are computed from vaccine arm
      mdw.weights=tryCatch({
        tree.weight(cor(dat_proc[kp & dat_proc$Trt==1, t%.%bAb], use='complete.obs'))
      }, error = function(err) {
        print(err$message)
        rep(1/length(bAb), length(bAb))
      })
      print(mdw.weights)  
      # center and scale are computed from vaccine arm
      centers=attr(scale(dat_proc[kp & dat_proc$Trt==1, t%.%bAb]), "scaled:center")
      scales=attr(scale(dat_proc[kp & dat_proc$Trt==1, t%.%bAb]), "scaled:scale")
      # transform both vaccine and placebo
      dat_proc[kp, t%.%'bindSpike_mdw'] = scale(dat_proc[kp, t%.%bAb], center=centers, scale=scales) %*% mdw.weights

      # nAb
      # mdw weights are computed from vaccine arm
      mdw.weights=tryCatch({
        tree.weight(cor(dat_proc[kp & dat_proc$Trt==1, t%.%nAb], use='complete.obs'))
      }, error = function(err) {
        print(err$message)
        rep(1/length(nAb), length(nAb))
      })
      print(mdw.weights)  
      # center and scale are computed from vaccine arm
      centers=attr(scale(dat_proc[kp & dat_proc$Trt==1, t%.%nAb]), "scaled:center")
      scales =attr(scale(dat_proc[kp & dat_proc$Trt==1, t%.%nAb]), "scaled:scale")
      # transform both vaccine and placebo
      dat_proc[kp, t%.%'pseudoneutid50_mdw'] = scale(dat_proc[kp, t%.%nAb], center=centers, scale=scales) %*% mdw.weights
    }
  }
}


# add two synthetic ID50 markers by region for ensemble

if(TRIAL %in% c("janssen_pooled_EUA", "janssen_na_EUA", "janssen_la_EUA", "janssen_sa_EUA")) {
    dat_proc$Day29pseudoneutid50la = ifelse(dat_proc$Day29pseudoneutid50-0.124 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50-0.124 )
    dat_proc$Day29pseudoneutid50sa = ifelse(dat_proc$Day29pseudoneutid50-0.556 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50-0.556 )
} else if(TRIAL %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA")) {
    dat_proc$Day29pseudoneutid50la = ifelse(dat_proc$Day29pseudoneutid50-0.255 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50-0.255 )
    dat_proc$Day29pseudoneutid50sa = ifelse(dat_proc$Day29pseudoneutid50-0.458 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50-0.458 )
}



###############################################################################
# subset on subset_variable
###############################################################################

if(!is.null(config$subset_variable) & !is.null(config$subset_value)){
    if(subset_value != "All") {
        include_in_subset <- dat_proc[[subset_variable]] == subset_value
        dat_proc <- dat_proc[include_in_subset, , drop = FALSE]
    }
}


###############################################################################
# impute covariates if necessary
# do this last so as not to change earlier values
###############################################################################

if (TRIAL %in% c("profiscov", "profiscov_lvmn", "vat08_combined")) {
    # no risk score for profiscov, but some have missing BMI
    n.imp <- 1
    dat.tmp.impute <- dat_proc
    
    if (TRIAL %in% c("profiscov", "profiscov_lvmn")) {
      imp.markers=c("HighRiskInd", "Sex", "Age", "BMI")
    } else if (TRIAL %in% c("vat08_combined")) {     
      imp.markers=c("FOI", "risk_score")
    }
        
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
} 


###############################################################################
# special handling 
###############################################################################

if(TRIAL == "moderna_real") {
    # modernal is a special case because how the code and manuscripts co-evolve 
    # special handling is required to preserve the imputed values used for manuscripts
    # the following steps treat seroneg and seropos populations separately and combine them to form one dataset
    
    # For the baseline seronegative population, use P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.csv, which was used for the lvmn manuscript
    dat_proc.tmp=read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.csv")
    # EarlyendpointD29start1 is deprecated at some point after P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.csv was made
    dat_proc.tmp=subset(dat_proc.tmp, select=-EarlyinfectionD29start1)
    dat_proc.tmp=subset(dat_proc.tmp, select=-EarlyendpointD29start1)
    # sort columns to combine with dat_proc
    dat_proc.tmp=dat_proc.tmp[,sort(names(dat_proc.tmp))]
    
    # For the baseline seropos population, impute lvmn data by calling the following script
    dat_proc=subset(dat_proc, Bserostatus==1, select=-c(risk_score_old, standardized_risk_score_old))
    source(here::here("data_clean", "add_lvmn_to_cove_analysisreadydataset.R"))
    # sort columns to combine with dat_proc.tmp
    dat_proc=dat_proc[,sort(names(dat_proc))]
    
    # combine baseline neg and pos
    stopifnot(all(names(dat_proc) == names(dat_proc.tmp)))
    stopifnot(all(dat_proc$Bserostatus.tmp == 0))
    stopifnot(all(dat_proc$Bserostatus == 1))
    dat_proc=rbind(dat_proc, dat_proc.tmp)
    
    # calendar time variable
    adsl_partc <- read.csv("/trials/covpn/p3001/analysis/mapping_immune_correlates/Part_C_Unblinded_Phase_Data/qdata/20230515/adsl.csv", stringsAsFactors = F)
    adsl_partc$CalendarDateEnrollment <- as.Date(adsl_partc$TRTSDT) - min(as.Date(adsl_partc$TRTSDT), na.rm=T)
    dat_proc=dat_proc %>% left_join(adsl_partc %>% 
                  dplyr::mutate(Ptid = gsub("mRNA-1273-P301-|-", "", USUBJID)) %>% 
                  dplyr::select(Ptid, CalendarDateEnrollment),
                by = c("Ptid"))
    
} else if(TRIAL %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA",
                                        "janssen_partA_VL")) {

    # add bin numbers associated with the biweekly calendar period of each endpoint 
    dat.eventtime.bin = read.csv("/trials/covpn/p3003/analysis/post_covid/sieve/Part_A_Blinded_Phase_Data/adata/omnibus/cpn3003_time_to_event_v6a.csv")
    if (TRIAL %in% c("janssen_pooled_partA", "janssen_partA_VL")) {
        v.name="endpointDate.Bin.Pooled"
    } else if (TRIAL %in% c("janssen_na_partA")) {
        v.name="endpointDate.Bin.usa"
    } else if (TRIAL %in% c("janssen_la_partA")) {
        v.name="endpointDate.Bin.latin.america"
    } else if (TRIAL %in% c("janssen_sa_partA")) {
        v.name="endpointDate.Bin.rsa"
    } 
    dat_proc$endpointDate.Bin = dat.eventtime.bin[[v.name]][match(dat_proc$Ptid, dat.eventtime.bin$USUBJID)]
    
    # add bin numbers associated with the biweekly calendar period of enrollment
    binTime = 13 # since bin period is biweekly!
    binVec = 1:8
    breaksVec = c(0, binTime*binVec+(binVec-1))    
    if(max(breaksVec) > max(dat_proc$CalendarDateEnrollment)){
      dat_proc <- dat_proc %>% 
        mutate(BinnedEnrollmentBiweekly = cut(CalendarDateEnrollment, breaks = breaksVec, 
                                              include.lowest = T, 
                                              labels = c("0", "1", "2", "3", "4", "5", "6", "7")))
    } else {
      "The maximum value in CalendarDateEnrollment is higher than max(breaksVec)! Update binVec!"
    }
    rm(binTime, binVec, breaksVec)
    
    # # add event time from sieve analysis
    # dat_proc$sieve.time = dat.eventtime.bin$time[match(dat_proc$Ptid, dat.eventtime.bin$USUBJID)]
    # dat_proc$sieve.status = dat.eventtime.bin$status[match(dat_proc$Ptid, dat.eventtime.bin$USUBJID)]
    

} else if(TRIAL == "prevent19") {
    # first round submission lacks RBD
    # for revision, RBD is added. To reproduce results from the first revision, we add RBD to the analysis-ready dataset, instead of reprocessing all three markers together, which will lead to changes in imputed values in the first two markers
    source(here::here("data_clean", "add_rbd_to_prevent19_analysisreadydataset.R"))
    
} else if(TRIAL == "azd1222") {
  # add bindSpike data for multivariable modles
  source(here::here("data_clean", "add_bindSpike_to_azd1222ID50_analysisreadydataset.R"))
  
} else if(TRIAL == "vat08_combined") {
  # add region variable for regression
  # Honduras (3), not Honduras for the Stage 1 trial nnaive
  # Mexico (9), Other/Else country for the Stage 2 trial naive
  # India (4), Mexico (9), Other/Else country for the Stage 2 trial nnaive
  dat_proc$Region = NA
  
  dat_proc$Region[dat_proc$Trialstage==1 & dat_proc$Bserostatus==0] = "Stage1naive"
  
  dat_proc$Region[dat_proc$Trialstage==1 & dat_proc$Bserostatus==1 & dat_proc$Country==3] = "HON_Stage1Nnaive"
  dat_proc$Region[dat_proc$Trialstage==1 & dat_proc$Bserostatus==1 & dat_proc$Country!=3] = "NotHON_Stage1Nnaive"
  
  dat_proc$Region[dat_proc$Trialstage==2 & dat_proc$Bserostatus==0 & dat_proc$Country==9] = "MEX_Stage2naive"
  dat_proc$Region[dat_proc$Trialstage==2 & dat_proc$Bserostatus==0 & dat_proc$Country!=9] = "NotMEX_Stage2naive"
  
  dat_proc$Region[dat_proc$Trialstage==2 & dat_proc$Bserostatus==1 & dat_proc$Country==9] = "MEX_Stage2Nnaive"
  dat_proc$Region[dat_proc$Trialstage==2 & dat_proc$Bserostatus==1 & dat_proc$Country==4] = "IND_Stage2Nnaive"
  dat_proc$Region[dat_proc$Trialstage==2 & dat_proc$Bserostatus==1 & dat_proc$Country!=9 & dat_proc$Country!=4] = "NotMEXIND_Stage2Nnaive"


}




###############################################################################
# digest check
###############################################################################

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         moderna_mock = "34e297fd1a736f9320573ff1d2944904",
         moderna_real = "684004027985464d02c4e157a86667e8",
         janssen_pooled_mock = "f3e286effecf1581eec34707fc4d468f",
         janssen_pooled_EUA = "c38fb43e2c87cf2d392757840af68bba",
         azd1222 = "f573e684800003485094c18120361663",
         azd1222_bAb = "fc3851aff1482901f079fb311878c172",
         prevent19 = "a4c1de3283155afb103261ce6ff8cec2",
         janssen_pooled_partA = "335d2628adb180d3d07745304d7bf603",
         janssen_partA_VL = "e7925542e4a1ccc1cc94c0e7a118da95", 
         vat08_combined = "00b5e96ee8fcdb89b87b619a1c62da48", 
         NA)    
    if (!is.na(tmp)) assertthat::assert_that(digest(dat_proc[order(names(dat_proc))])==tmp, msg = "failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))]))    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))


# split into Senior and non-Senior for ENSEMBLE partA
if(TRIAL %in% c("janssen_pooled_partA", "janssen_na_partA", "janssen_la_partA", "janssen_sa_partA")) {
    write_csv(subset(dat_proc, Age>=60), file = here("data_clean", paste0(TRIAL, "senior_data_processed_with_riskscore.csv")))
    write_csv(subset(dat_proc, Age< 60),  file = here("data_clean", paste0(TRIAL, "nonsenior_data_processed_with_riskscore.csv")))
}


print("run time: "%.%format(Sys.time()-begin, digits=1))
