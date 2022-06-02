#Sys.setenv(TRIAL = "azd1222")
renv::activate(here::here())
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("_common.R"))

library(tidyverse)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)
library(here)


#mapped.dat=read.csv(mapped_data)


########################################################################################################
# read mapped data with risk score added

# inputFile_with_riskscore has the same number of rows as dat_raw, but adds columns related to earlyendpoint, earlyinfection (preprocess.for.risk.score) and risk scores
load(file = paste0("riskscore_baseline/output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData"))
dat_proc <- inputFile_with_riskscore

#with(dat_proc[dat_proc$Trt==1,], table(!is.na(Day29bindSpike), !is.na(Day29bindRBD), EventIndPrimaryIncludeNotMolecConfirmedD29)) # same missingness

# hardcode AnyinfectionD1 for the mock datasets
if (study_name %in% c("MockENSEMBLE", "MockCOVE")) {
    dat_proc$AnyinfectionD1=0
}


colnames(dat_proc)[1] <- "Ptid" 
dat_proc <- dat_proc %>% mutate(age.geq.65 = as.integer(Age >= 65))
dat_proc$Senior = as.integer(dat_proc$Age>=switch(study_name, COVE=65, MockCOVE=65, ENSEMBLE=60, MockENSEMBLE=60, PREVENT19=65, AZD1222=65, stop("unknown study_name")))
  

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
      
} else stop("unknown study_name")

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
if (study_name %in% c("COVE", "MockCOVE")) {
    # nothing to do, US data only
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    dat_proc$MinorityInd[dat_proc$Region!=0] = 0

} else if (study_name=="PREVENT19") {
    dat_proc$MinorityInd[dat_proc$Country!=0] = 0 # 0 is US

} else if (study_name=="AZD1222") {
    dat_proc$MinorityInd[dat_proc$Country!=2] = 0 # 2 is US

} else stop("unknown study_name")



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
    
} else if (study_name %in% c("PREVENT19", "AZD1222")) {
    dat_proc$Bstratum = with(dat_proc, ifelse(Senior, 1, 0))
    
} else stop("unknown study_name")

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
        
} else stop("unknown study_name")  
  
names(demo.stratum.labels) <- demo.stratum.labels


# tps stratum, 1 ~ 4*max(demo.stratum), used in tps regression
dat_proc <- dat_proc %>%
  mutate(
    tps.stratum = demo.stratum + strtoi(paste0(Trt, Bserostatus), base = 2) * length(demo.stratum.labels)
  )

# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+4. 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# A case will have a Wstratum even if its tps.stratum is NA
# The case is defined using EventIndPrimaryD29

max.tps=max(dat_proc$tps.stratum,na.rm=T)
dat_proc$Wstratum = dat_proc$tps.stratum
if (study_name %in% c("COVE", "MockCOVE", "ENSEMBLE", "MockENSEMBLE", "AZD1222")) {
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==0)]=max.tps+1
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==0 & Bserostatus==1)]=max.tps+2
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==0)]=max.tps+3
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD29==1 & Trt==1 & Bserostatus==1)]=max.tps+4
    
} else if (study_name == "PREVENT19") {
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD21==1 & Trt==0 & Bserostatus==0)]=max.tps+1
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD21==1 & Trt==0 & Bserostatus==1)]=max.tps+2
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD21==1 & Trt==1 & Bserostatus==0)]=max.tps+3
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD21==1 & Trt==1 & Bserostatus==1)]=max.tps+4
    
} else stop("unknown study_name")


#subset(dat_proc, Trt==1 & Bserostatus==1 & EventIndPrimaryD29 == 1)[1:3,]

with(dat_proc, table(tps.stratum))



###############################################################################
# observation-level weights
# Note that Wstratum may have NA if any variables to form strata has NA
###############################################################################


# define must have assays for ph2 definition
if (study_name %in% c("COVE", "MockCOVE")) {
    must_have_assays <- c("bindSpike", "bindRBD")
    
} else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")) {
    if (endsWith(attr(config, "config"),"ADCP")) {
        must_have_assays <- c("ADCP")    
    } else {
        must_have_assays <- c("bindSpike", "bindRBD")
    }    
} else if (study_name %in% c("PREVENT19")) {
    must_have_assays <- c("bindSpike")
    
} else if (study_name %in% c("AZD1222")) {
    if (attr(config, "config")=="azd1222") {
        must_have_assays <- c("pseudoneutid50")
    } else if (attr(config, "config")=="azd1222_bAb") {
        must_have_assays <- c("bindSpike")
    } else stop("need to define must_have_assays")
    
    
} else stop("unknown study_name")


# TwophasesampInd: be in the case or subcohort  &  have the necessary markers
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

} else if (study_name=="AZD1222") {
    # does not require baseline or time point 1
    dat_proc[["TwophasesampIndD"%.%timepoints[2]]] = 
        with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
        complete.cases(dat_proc[,c("Day"%.%timepoints[2]%.%must_have_assays)])        
    # does not require baseline
    dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
        with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
        # adding | is because if D57 is present, D29 will be imputed if missing
        (complete.cases(dat_proc[,c("Day"%.%timepoints[1]%.%must_have_assays)]) | complete.cases(dat_proc[,c("Day"%.%timepoints[2]%.%must_have_assays)])) 
        
} else if (study_name=="ENSEMBLE") {
    if (contain(attr(config, "config"), "real")) {
        # require baseline
        dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("B"%.%must_have_assays, "Day"%.%timepoints[1]%.%must_have_assays)])      
    } else if (contain(attr(config, "config"), "partA")) {
        # does not require baseline
        dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("Day"%.%timepoints[1]%.%must_have_assays)])
    }
    
} else stop("unknown study_name")



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
        


###############################################################################
# impute missing neut biomarkers in ph2
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
    
    if(two_marker_timepoints) {
        imp.markers=c(outer(c("B", if(tp==timepoints[2]) "Day"%.%timepoints else "Day"%.%tp), assays, "%.%"))
    } else {
        imp.markers=c(outer(c("B", "Day"%.%tp), assays, "%.%"))
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
    
    # imputed values of missing markers merged properly for all individuals in the two phase sample?
    assertthat::assert_that(
      all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD"%.%tp]] == 1, imp.markers])),
      msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
    )
}


assays.includeN=c(assays, if(!study_name %in% c("PREVENT19","AZD1222")) "bindN")


###############################################################################
# converting binding variables from AU to IU for binding assays
# moderna_real immune.csv file is not on international scale and other mapped data files are
###############################################################################

# COVE only 
if(study_name=="COVE"){
    # conversion is only done for COVE for backward compatibility
    convf=c(bindSpike=0.0090, bindRBD=0.0272, bindN=0.0024, pseudoneutid50=0.242, pseudoneutid80=1.502)    
    for (a in assays.includeN) {
      for (t in c("B", paste0("Day", config$timepoints)) ) {
          dat_proc[[t %.% a]] <- dat_proc[[t %.% a]] + log10(convf[a])
      }
    }
}


###############################################################################
# censoring values below LLOD
# after COVE, the mapped data comes censored
###############################################################################

# COVE and mock only
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
# define delta for dat_proc
###############################################################################

tmp=list()
# delta is computed after lloq censoring
for (a in assays.includeN) {
  for (t in c("B", paste0("Day", config$timepoints)) ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(lloqs[a]), log10(lloqs[a] / 2), dat_proc[[t %.% a]])
  }
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

# some trials have N some don't
for (tp in rev(timepoints)) {
    dat_proc["Delta"%.%tp%.%"overB" %.% assays.includeN] <- tmp["Day"%.%tp %.% assays.includeN] - tmp["B" %.% assays.includeN]
}   
if(two_marker_timepoints) {
    dat_proc["Delta"%.%timepoints[2]%.%"over"%.%timepoints[1] %.% assays.includeN] <- tmp["Day"%.% timepoints[2]%.% assays.includeN] - tmp["Day"%.%timepoints[1] %.% assays.includeN]
}



###############################################################################
# add two synthetic ID50 markers by region for ensemble
###############################################################################

if(attr(config, "config") %in% c("janssen_pooled_real", "janssen_na_real", "janssen_la_real", "janssen_sa_real")) {
    dat_proc$Day29pseudoneutid50sa = ifelse(dat_proc$Day29pseudoneutid50sa – 0.470 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50sa – 0.470)
    dat_proc$Day29pseudoneutid50la = ifelse(dat_proc$Day29pseudoneutid50la – 0.271 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50la – 0.271)
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
# bundle data sets and save as CSV
###############################################################################

library(digest)
if(attr(config, "config") %in% c("janssen_pooled_mock", "moderna_mock") & Sys.getenv ("NOCHECK")=="") {
    assertthat::assert_that(
        digest(dat_proc[order(names(dat_proc))])==
        ifelse(attr(config, "config")=="janssen_pooled_mock", 
            "028549acb994980fbb6832198308ec4a", 
            "1902296f23fca88c4757ed7a84fbe7d7"),
        msg = "failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))]))    
    print("======================= Passed make_dat_proc digest check =======================")    
}



write_csv(dat_proc, file = here("data_clean", paste0(attr(config, "config"), "_data_processed_with_riskscore.csv")))
