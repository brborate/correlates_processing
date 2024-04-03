Sys.setenv(TRIAL = "prevent19_stage2")
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

# borrow risk score from prevent19, even though the clinical database changed slightly
load(file = 'riskscore_baseline/output/prevent19/inputFile_with_riskscore.RData')
nrow(inputFile_with_riskscore)
nrow(dat_proc)

# ptids are formatted differently in the two datasets, need to be transformed
inputFile_with_riskscore$Ptid = sub("2019nCoV301","2019nCoV-301",inputFile_with_riskscore$Ptid)

# stage 2 dataset has fewer rows than stage 1 b/c it is vaccine only, baseline sero-negative only
dat_proc$risk_score = inputFile_with_riskscore$risk_score[match(dat_proc$Ptid, inputFile_with_riskscore$Ptid)]
dat_proc$standardized_risk_score = inputFile_with_riskscore$standardized_risk_score[match(dat_proc$Ptid, inputFile_with_riskscore$Ptid)]

# subset to vaccine and seroneg, necessary for getting weights correctly
dat_proc = subset(dat_proc, Trt==1 & Bserostatus==0)



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
dat_proc$demo.stratum = with(dat_proc, strtoi(paste0(URMforsubcohortsampling, Senior, HighRiskInd), base = 2)) + 1
dat_proc$demo.stratum = with(dat_proc, ifelse(Country==0, demo.stratum, ifelse(!Senior, 9, 10))) # 0 is US
with(dat_proc, table(demo.stratum))


# tps stratum, used in tps regression and to define Wstratum
dat_proc <- dat_proc %>% mutate(tps.stratum = demo.stratum + strtoi(paste0(Trt, Bserostatus), base = 2) * max(demo.stratum,na.rm=T))
if (!is.null(dat_proc$tps.stratum)) table(dat_proc$tps.stratum)


# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# A case will have a Wstratum even if its tps.stratum is NA

dat_proc$Wstratum = dat_proc$tps.stratum
max.tps=max(dat_proc$tps.stratum,na.rm=T)
# cannot use KnownOrImputedDeltaCOVIDInd21Apr19to22Mar26 b/c cases are not defined using that
dat_proc$Wstratum[with(dat_proc, KnownOrImputedDeltaCOVIDIndD35_108to21Dec10==1 & Trt==1 & Bserostatus==0)]=max.tps+2
# severe has to come second to overwrite delta
dat_proc$Wstratum[with(dat_proc, SevereCOVIDIndD35_108to21Dec10 ==1 & Trt==1 & Bserostatus==0)]=max.tps+2

table(dat_proc$Wstratum) 



################################################################################
# 4. Define ph1, ph2, and weights
# Note that Wstratum may have NA if any variables to form strata has NA

tp='35_108'

# ph1
dat_proc[["ph1.D"%.%tp]] = with(dat_proc, 
                                Perprotocol==1 & 
                                  get("AnyInfectionD1toD"%.%tp)==0 & # no evidence of any infection by D35_108
                                  COVIDTimeD35to21Dec10 >= 108 & # COVID time or censor time is after D35_108
                                  # either Delta COVID, severe COVID, or no evidence of infection till 22Mar26
                                  (KnownOrImputedDeltaCOVIDInd21Apr19to22Mar26==1 | SevereCOVIDInd21Apr19to22Mar26==1 | AnyInfectionD1to22Mar26==0)
)


# TwophasesampInd: be in the subcohort or a case after time point 1  &  have the necessary markers
# requires both bAb and nAb
dat_proc$TwophasesampIndD35 = with(dat_proc, 
                                     # bAb is all or none
                                    !is.na(Day35bindSpike_D614) & 
                                     # nAb will be used to impute each other
                                    (!is.na(Day35pseudoneutid50_D614G) | !is.na(Day35pseudoneutid50_Delta) )
)

# remove three ptids from TwophasesampIndD35 because NVX programmer used the investigator’s 
# assessment of severity (variable SEV) rather than the final severity assessment ASEV. 
dat_proc$TwophasesampIndD35[dat_proc$Ptid %in% c(
  "2019nCoV-301-US228-0105", "2019nCoV-301-US232-0013", # non-severe Mu and Gamma cases in the final determination
  "2019nCoV-301-US179-0082") # a less than mild Delta case in the final determination
  ] = F 

  
# ph2
dat_proc[["ph2.D"%.%tp]] = dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD35"]] & 
  # remove cases outside D35_108to21Dec10
  !(dat_proc$KnownOrImputedDeltaCOVIDInd21Apr19to22Mar26==1 & dat_proc$KnownOrImputedDeltaCOVIDIndD35_108to21Dec10==0)

# weights
dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp, ph2="ph2.D"%.%tp, Wstratum="Wstratum", wt="wt.D"%.%tp, verbose=F) 
  


###############################################################################
# 5. impute missing biomarkers in ph2 (assay imputation)
#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, each time point, but not Delta

# loop through the time points
# first impute (B, D29, D57) among TwophasesampIndD57==1
# next impute (B, D29) among TwophasesampIndD29==1

for (tp in rev(timepoints)) {    
    n.imp <- 1
    
    dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD"%.%tp) == 1)

    # markers can be TRIAL-specific
    # no need to impute bAb
    imp.markers=paste0("Day"%.%tp, c("pseudoneutid50_D614G", "pseudoneutid50_Delta"))
    # mdw markers are not imputed
    imp.markers=imp.markers[!endsWith(imp.markers, "_mdw")]

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


  
###############################################################################
# 6. transformation of the markers

# converting binding variables from AU to IU for binding assays
# COVE only 
# moderna_real immune.csv file is not on international scale and other mapped data files are

# some trials have N some don't
if (is.null(config$assay_metadata)) {
  includeN = switch(study_name, COVE=1, MockCOVE=1, ENSEMBLE=1, MockENSEMBLE=1, PREVENT19=0, AZD1222=0, VAT08=0, PROFISCOV=1, COVAIL=1, stop("unknown study_name 9"))
  assays.includeN=c(assays, if(includeN==1) "bindN")
  
} else if (TRIAL=="janssen_partA_VL") {
    assays.includeN = c("bindSpike", "bindRBD", "pseudoneutid50", "bindN")
    
} else {
  assays.includeN = assays
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


###############################################################################
# 7. add mdw scores



###############################################################################
# 8. add fold change markers
# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta


# skipping b/c there is no baseline data


###############################################################################
# 9. add discrete/trichotomized markers

dat_proc$tmp = with(dat_proc, Trt==1 & Bserostatus==0 & get("ph2.D35_108")) 
dat_proc = add.trichotomized.markers (dat_proc, c("Day35"%.%assays), ph2.col.name="tmp", wt.col.name="wt.D35_108")
dat_proc$tmp = NULL
  


###############################################################################
# 10. impute covariates if necessary



###############################################################################
# special handling 




###############################################################################
# digest check

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         prevent19_stage2 = "0197ec25780673886bd7b985a4ffdba4",
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))



print("run time: "%.%format(Sys.time()-begin, digits=1))
