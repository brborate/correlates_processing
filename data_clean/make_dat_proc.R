#Sys.setenv(TRIAL = "moderna_real")
#Sys.setenv(TRIAL = "vat08_combined")
#Sys.setenv(TRIAL = "covail")
#Sys.setenv(TRIAL = "janssen_partA_VL")
#Sys.setenv(TRIAL = "azd1222_stage2")

# no need to run renv::activate(here::here()) b/c .Rprofile exists

source(here::here("_common.R"))

if (TRIAL %in% c("moderna_boost","id27hpv")) stop("For moderna_boost, run make_dat_moderna_boost.R") 

library(tidyverse)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)
library(here)
library(mdw)

begin=Sys.time()



########################################################################################################
# read mapped data with risk score added

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
  # define here b/c it is needed to define weights
  dat_proc$EventIndPrimaryHasVLD29 = dat_proc$EventIndPrimaryIncludeNotMolecConfirmedD29
  # if EventIndPrimaryHasVLD29 is NA, this will remain NA
  dat_proc$EventIndPrimaryHasVLD29[is.na(dat_proc$seq1.log10vl) & !is.na(dat_proc$EventIndPrimaryHasVLD29)] = 0 
  dat_proc$EventIndPrimaryD29 = dat_proc$EventIndPrimaryHasVLD29
  # EventTimePrimaryD29 is already set to EventIndPrimaryIncludeNotMolecConfirmedD29 in preprocess()
  
  country.codes=c("USA", "ARG", "BRA", "CHL", "COL", "MEX", "PER", "ZAF")
  dat_proc$cc=country.codes[dat_proc$Country+1]
  
  # update SubcohortInd to 1) include the Colombians sampled for the variants study 
  dat_proc$COL_variants_study = with(dat_proc, cc=="COL" & EventIndPrimaryIncludeNotMolecConfirmedD1==0 & SubcohortInd!=1 & !is.na(Day29bindSpike_D614))
  mytable(dat_proc$COL_variants_study)
  dat_proc$SubcohortInd = ifelse(dat_proc$SubcohortInd | dat_proc$COL_variants_study, 1, 0)
  # 2) remove ARV users
  dat_proc$SubcohortInd = ifelse(dat_proc$SubcohortInd & dat_proc$ARVuseDay29==0, 1, 0)
  
  # create a mdw score
  # the reason to do it in the beginning is that this score replaces the five delta markers completely
  delta_markers = c("bindSpike_AY.2", "bindSpike_B.1.617.2_AY.4", "bindSpike_AY.12", "bindSpike_AY.1", "bindSpike_B.1.617.2" )
  mdw.wt=tryCatch({
    tree.weight(cor(dat_proc["Day29"%.%delta_markers], use='complete.obs'))
  }, error = function(err) {
    print(err$message)
    rep(1/length(delta_markers), length(delta_markers))
  })
  write.csv(mdw.wt, file = here("data_clean", "csv", TRIAL%.%"_delta_score_mdw_weights.csv"))
  t="Day29"
  dat_proc[, t%.%'bindSpike_DeltaMDW'] = c(as.matrix(dat_proc[, t%.%delta_markers]) %*% mdw.wt)
  # remove delta_markers
  for (a in delta_markers) dat_proc[, t%.%a] = NULL

  
  
 } else if (study_name=="VAT08") {

   # read hot deck data
   tmp = sub(".csv","_hotdeck.csv",mapped_data)
   if (file.exists(tmp)) dat_raw=read.csv(tmp) else stop("hotdeck file not exists, run hotdeck R script first")
   
   dat_proc = preprocess(dat_raw, study_name)   
   colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
   
   # scale FOI
   dat_proc$FOI = scale(log10(dat_proc$FOI+1))[,1]
   
   # add country code
   country.codes=c("Colombia", "Ghana", "Honduras", "India", "Japan", "Kenya", "Nepal", "United States", "Mexico", "Uganda", "Ukraine")
   continents=c("Colombia"=1, "Ghana"=2, "Honduras"=1, "India"=3, "Japan"=3, "Kenya"=2, "Nepal"=3, "United States"=5, "Mexico"=1, "Uganda"=2, "Ukraine"=4)
   dat_proc$cc = country.codes[dat_proc$Country]
   dat_proc$continent = continents[dat_proc$cc]
   table(dat_proc$continent, dat_proc$cc)
   
   # add risk score
   load(file = paste0('riskscore_baseline/output/vat08_combined/inputFile_with_riskscore.RData'))
   stopifnot(all(inputFile_with_riskscore$Ptid==dat_proc$Ptid))
   dat_proc$risk_score = inputFile_with_riskscore$risk_score
   dat_proc$standardized_risk_score = inputFile_with_riskscore$standardized_risk_score
   
   # ptids with missing Bserostatus already filtered out in preprocess
   
   # define event indicator and event time variables based on seq1.variant.hotdeck1 etc
   for (t in c(1,22,43)) {
     for (i in 1:10) {
       dat_proc[[paste0("EventIndOmicronD",t,"M12hotdeck",i)]]  = 
         ifelse(!is.na(dat_proc[["seq1.variant.hotdeck"%.%i]]) & dat_proc[["seq1.variant.hotdeck"%.%i]]=="Omicron" & !is.na(dat_proc[["EventIndFirstInfectionD"%.%t]]),
                1,
                0)
       
       dat_proc[[paste0("EventTimeOmicronD",t,"M12hotdeck",i)]] = ifelse(dat_proc[[paste0("EventIndOmicronD",t,"M12hotdeck",i)]] ==1, 
                                                                      pmin(dat_proc[["EventTimeKnownLineageOmicronD"%.%t]],    dat_proc[["EventTimeMissingLineageD"%.%t]]),
                                                                      pmax(dat_proc[["EventTimeKnownLineageNonOmicronD"%.%t]], dat_proc[["EventTimeMissingLineageD"%.%t]]))
     }
   }
   
   # create event time and indicator variables censored on 180 days post dose 2
   for (t in c(1,22,43)) {
     for (i in 1:10) {
       dat_proc[[paste0("EventIndOmicronD",t,"M6hotdeck",i)]]  = ifelse (dat_proc[[paste0("EventTimeOmicronD43M12hotdeck",i)]]>180-21, 0,   dat_proc[[paste0("EventIndOmicronD",t,"M12hotdeck",i)]])
       dat_proc[[paste0("EventTimeOmicronD",t,"M6hotdeck",i)]] = ifelse (dat_proc[[paste0("EventTimeOmicronD43M12hotdeck",i)]]>180-21, 180-21, dat_proc[[paste0("EventTimeOmicronD",t,"M12hotdeck",i)]])
     }
   }
   
   
} else if (TRIAL == "covail") {
  # load risk score
  load(file = paste0('riskscore_baseline/output/',TRIAL,'/inputFile_with_riskscore.RData'))
  dat_proc <- inputFile_with_riskscore    

  # bring in imputed variant column
  dat.lineage = read.csv('/trials/covpn/COVAILcorrelates/analysis/correlates/adata/lineages/covail_lineages_export_v1.csv')
  dat_proc$COVIDlineage = dat.lineage$inf1.lineage[match(dat_proc$Ptid, dat.lineage$ptid)]
  dat_proc$COVIDlineageObserved = dat.lineage$inf1.observed[match(dat_proc$Ptid, dat.lineage$ptid)]
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
  
  
  
} else if (TRIAL == "azd1222_stage2") {
  dat_raw=read.csv(mapped_data)
  dat_proc = preprocess(dat_raw, study_name)   
  colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
  
  # borrow risk score from azd1222_stage2
  load(file = 'riskscore_baseline/output/azd1222/inputFile_with_riskscore.RData')
  # stage 2 dataset has fewer rows than stage 1
  dat_proc$risk_score = inputFile_with_riskscore$risk_score[match(dat_proc$Ptid, inputFile_with_riskscore$Ptid)]
  dat_proc$standardized_risk_score = inputFile_with_riskscore$standardized_risk_score[match(dat_proc$Ptid, inputFile_with_riskscore$Ptid)]
  
  
} else {
  if (make_riskscore) {
    # load inputFile_with_riskscore.Rdata, a product of make riskscore_analysis, which calls preprocess and makes risk scores
    load(file = paste0('riskscore_baseline/output/',TRIAL,'/inputFile_with_riskscore.RData'))
    dat_proc <- inputFile_with_riskscore    
  } else {
    dat_raw=read.csv(mapped_data)
    dat_proc = preprocess(dat_raw, study_name)   
  }
  
}


# define new variables
{ # use this to navigate faster in Rstudio
  colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
  dat_proc <- dat_proc %>% mutate(age.geq.65 = as.integer(Age >= 65))
  dat_proc$Senior = as.integer(dat_proc$Age>=switch(study_name, COVE=65, MockCOVE=65, ENSEMBLE=60, MockENSEMBLE=60, PREVENT19=65, AZD1222=65, VAT08=60, PROFISCOV=NA, COVAIL=65, stop("unknown study_name 1")))
  
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
    
  } else if (study_name %in% c("PROFISCOV", "VAT08", "COVAIL")) {
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
  if (study_name %in% c("COVE", "MockCOVE", "PROFISCOV", "VAT08", "COVAIL")) {
    # nothing to do
    # COVE only has US data
    
  } else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    dat_proc$MinorityInd[dat_proc$Region!=0] = 0
    
  } else if (study_name=="PREVENT19") {
    dat_proc$MinorityInd[dat_proc$Country!=0] = 0 # 0 is US
    
  } else if (study_name=="AZD1222") {
    dat_proc$MinorityInd[dat_proc$Country!=2] = 0 # 2 is US
    
  } else stop("unknown study_name 3")  
}



###############################################################################
# stratum variables
# The code for Bstratum is trial specifc
# The code for tps.stratum and Wstratum are not trial specific since they are constructed on top of Bstratum
###############################################################################

# Bstratum: randomization strata
# e.g., Moderna: 1 ~ 3, defines the 3 baseline strata within trt/serostatus
if (study_name=="COVE" | study_name=="MockCOVE" ) {
    dat_proc$Bstratum = with(dat_proc, ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3)))
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
    dat_proc$Bstratum =  with(dat_proc, strtoi(paste0(Senior, HighRiskInd), base = 2)) + 1
    
} else if (study_name %in% c("PREVENT19", "AZD1222")) {
  dat_proc$Bstratum =  with(dat_proc, Senior + 1)
  
} else if (study_name %in% c("VAT08")) {
  dat_proc$Bstratum =  1 
  
} else if (study_name %in% c("PROFISCOV", "COVAIL")) {
  dat_proc$Bstratum = 1 # there are no demographics stratum for subcohort sampling
  
} else stop("unknown study_name 4")

names(Bstratum.labels) <- Bstratum.labels


# demo.stratum: correlates sampling strata
# Moderna: 1 ~ 6 defines the 6 baseline strata within trt/serostatus
# may have NA b/c URMforsubcohortsampling may be NA
if (study_name=="COVE" | study_name=="MockCOVE" ) {
    dat_proc$demo.stratum = with(dat_proc, ifelse (URMforsubcohortsampling==1, ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3)), 3+ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3))))
    
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
    # first step, stratify by age and high risk
    dat_proc$demo.stratum =  with(dat_proc, strtoi(paste0(Senior, HighRiskInd), base = 2)) + 1
    # second step, stratify by region
    dat_proc$demo.stratum=with(dat_proc, ifelse(Region==0 & URMforsubcohortsampling==0, demo.stratum + 4, demo.stratum)) # US, non-URM
    dat_proc$demo.stratum[dat_proc$Region==1] = dat_proc$demo.stratum[dat_proc$Region==1] + 8 # Latin America
    dat_proc$demo.stratum[dat_proc$Region==2] = dat_proc$demo.stratum[dat_proc$Region==2] + 12 # Southern Africa
    # the above sequence ends up setting US URM=NA to NA
    
    # for the variants study, in LatAm, COL becomes part of the demographics sampling strata
    # we add 8 for COL b/c +4 is for RSA
    if (TRIAL=="janssen_partA_VL") {
      dat_proc$demo.stratum[dat_proc$Region==1] = with(subset(dat_proc,Region==1), 
        ifelse(cc=="COL", demo.stratum+8, demo.stratum)
      )
    }
    
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
    nCountries = 10 # 10 is easier to work with than length(unique(dat_proc$Country))
    
    # for nAb
    # use continent for naive and country for nnaive
    dat_proc$demo.stratum = ifelse(dat_proc$Bserostatus==0, dat_proc$continent, dat_proc$Country)
    unique(subset(dat_proc, cc=='Japan', c(demo.stratum, Bserostatus)))
    # move JPN to be in the US group for Naive
    dat_proc$demo.stratum[dat_proc$Bserostatus==0 & dat_proc$cc=='Japan']=5
    unique(subset(dat_proc, cc=='Japan', c(demo.stratum, Bserostatus)))
    
    # for bAb
    # use continent for all
    dat_proc$demo.stratum2 = dat_proc$continent
    # move JPN to be in the US group 
    dat_proc$demo.stratum2[dat_proc$cc=='Japan']=5

    
} else if (study_name %in% c("PROFISCOV", "COVAIL") ) {
    dat_proc$demo.stratum = 1 # # there are no demographics stratum for subcohort sampling

        
} else stop("unknown study_name 5")  
# names(demo.stratum.labels) <- demo.stratum.labels

with(dat_proc, table(demo.stratum))



###########################################
# tps stratum, used in tps regression and to define Wstratum
if (study_name=="VAT08" ) {
  # for nAb markers. include Senior in the list of stratification variables
  dat_proc <- dat_proc %>% mutate(tps.stratum.nAb = demo.stratum + 
                            strtoi(paste0(Trialstage-1, Bserostatus, Trt, Senior), base = 2) * nCountries)
  
  # not include Senior in the list of stratification variables, used for bAb markers
  dat_proc <- dat_proc %>% mutate(tps.stratum.bAb = demo.stratum2 + 
                                strtoi(paste0(Trialstage-1, Bserostatus, Trt, Senior), base = 2) * nCountries)
  
  # # original, Oct 30 version
  # dat_proc <- dat_proc %>% mutate(tps.stratum.original = strtoi(paste0(Trt, Bserostatus, Trialstage-1, Senior), base = 2))
  

} else if (study_name=="COVAIL" ) {
  dat_proc <- dat_proc %>% mutate(tps.stratum = arm)
  
} else {
  dat_proc <- dat_proc %>% mutate(tps.stratum = demo.stratum + strtoi(paste0(Trt, Bserostatus), base = 2) * max(demo.stratum,na.rm=T))
}

if (!is.null(dat_proc$tps.stratum)) table(dat_proc$tps.stratum)



###########################################
# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+4. 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# A case will have a Wstratum even if its tps.stratum is NA
# The case is defined using EventIndPrimaryD29

if (study_name == "VAT08") {
  max.tps=max(dat_proc$tps.stratum.nAb,na.rm=T) # 160
  myprint(max.tps)

  # nAb
  dat_proc$Wstratum.nAb = dat_proc$tps.stratum.nAb
  dat_proc$Wstratum.nAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==0 & Trialstage==1)]=max.tps+1
  dat_proc$Wstratum.nAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==1 & Trialstage==1)]=max.tps+2
  dat_proc$Wstratum.nAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==0 & Trialstage==1)]=max.tps+3
  dat_proc$Wstratum.nAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==1 & Trialstage==1)]=max.tps+4
  dat_proc$Wstratum.nAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==0 & Trialstage==2)]=max.tps+5
  dat_proc$Wstratum.nAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==1 & Trialstage==2)]=max.tps+6
  dat_proc$Wstratum.nAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==0 & Trialstage==2)]=max.tps+7
  dat_proc$Wstratum.nAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==1 & Trialstage==2)]=max.tps+8
  
  
  # bAb
  dat_proc$Wstratum.bAb = dat_proc$tps.stratum.bAb
  dat_proc$Wstratum.bAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==0 & Trialstage==1)]=max.tps+1
  dat_proc$Wstratum.bAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==1 & Trialstage==1)]=max.tps+2
  dat_proc$Wstratum.bAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==0 & Trialstage==1)]=max.tps+3
  dat_proc$Wstratum.bAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==1 & Trialstage==1)]=max.tps+4
  dat_proc$Wstratum.bAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==0 & Trialstage==2)]=max.tps+5
  dat_proc$Wstratum.bAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==1 & Trialstage==2)]=max.tps+6
  dat_proc$Wstratum.bAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==0 & Trialstage==2)]=max.tps+7
  dat_proc$Wstratum.bAb[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==1 & Trialstage==2)]=max.tps+8
  
  # original
  # dat_proc$Wstratum.original = dat_proc$tps.stratum.original
  # dat_proc$Wstratum.original[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==0 & Trialstage==1)]=max.tps+1
  # dat_proc$Wstratum.original[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==1 & Trialstage==1)]=max.tps+2
  # dat_proc$Wstratum.original[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==0 & Trialstage==1)]=max.tps+3
  # dat_proc$Wstratum.original[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==1 & Trialstage==1)]=max.tps+4
  # dat_proc$Wstratum.original[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==0 & Trialstage==2)]=max.tps+5
  # dat_proc$Wstratum.original[with(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Bserostatus==1 & Trialstage==2)]=max.tps+6
  # dat_proc$Wstratum.original[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==0 & Trialstage==2)]=max.tps+7
  # dat_proc$Wstratum.original[with(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Bserostatus==1 & Trialstage==2)]=max.tps+8
  
} else if (study_name=="COVAIL" ) {
  dat_proc$Wstratum = dat_proc$tps.stratum
  
} else {
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
    
  } else if (study_name == "PROFISCOV") {
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD43==1 & Trt==0 & Bserostatus==0)]=max.tps+1
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD43==1 & Trt==0 & Bserostatus==1)]=max.tps+2
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD43==1 & Trt==1 & Bserostatus==0)]=max.tps+3
    dat_proc$Wstratum[with(dat_proc, EventIndPrimaryD43==1 & Trt==1 & Bserostatus==1)]=max.tps+4
    
  } else stop("unknown study_name 6")  
  
}  

if (!is.null(dat_proc$Wstratum)) table(dat_proc$Wstratum) # variables may be named other than Wstratum



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
      
  } else if (TRIAL=="janssen_partA_VL") {
    # we change from bindSpike to pseudoneutid50+bindSpike
    # 4 ptids have Day29bindSpike but no Day29pseudoneutid50, they are all vaccinees, baseline neg, non-cases, US
    # 1 case in US has Day29pseudoneutid50 but no Day29bindSpike
    must_have_assays <- c("pseudoneutid50", "bindSpike")
    
  } else {
    # bindSpike and bindRBD are all or none
    must_have_assays <- c("bindSpike","bindRBD")
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
  
    
} else if (study_name %in% c("VAT08", "COVAIL")) {
  must_have_assays <- NULL

  
} else stop("unknown study_name 7")



# TwophasesampInd: be in the subcohort or a case after time point 1  &  have the necessary markers
if (study_name %in% c("COVE", "MockCOVE", "MockENSEMBLE", "PREVENT19")) {
    if (two_marker_timepoints) {
    # require baseline and timepoint 1
        dat_proc[["TwophasesampIndD"%.%timepoints[2]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("B"%.%must_have_assays, "Day"%.%timepoints[1]%.%must_have_assays, "Day"%.%timepoints[2]%.%must_have_assays)])      
    }
    # require baseline
    dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
            with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
            complete.cases(dat_proc[,c("B"%.%must_have_assays, "Day"%.%timepoints[1]%.%must_have_assays)])      
        
    
} else if (study_name=="ENSEMBLE") {
  if (endsWith(TRIAL, "EUA")) {
    # require baseline
    dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
      with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
      complete.cases(dat_proc[,c("B"%.%must_have_assays, "Day"%.%timepoints[1]%.%must_have_assays)])      
    
  } else if (endsWith(TRIAL, "partA")) {
    # does not require baseline
    dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
      with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
      complete.cases(dat_proc[,c("Day"%.%timepoints[1]%.%must_have_assays)])
    
  } else if (TRIAL=="janssen_partA_VL") {
    # does not require baseline
    dat_proc[["TwophasesampIndD"%.%timepoints[1]]] = 
      with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
      (complete.cases(dat_proc[,"Day29"%.%must_have_assays]) | dat_proc$COL_variants_study) # add 98 COL to ph2
    
    # define TwophasesampIndD29variant, which is same as TwophasesampIndD29 except that 
    # for the non-cases, it is limited to a subset with variant ID50
    dat_proc$TwophasesampIndD29variant = dat_proc$TwophasesampIndD29 
    noncases = with(dat_proc, SubcohortInd & (is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0))
    # tmp is an indicator for having ancestral marker in US, beta in RSA, and Gamma/Lambda/Mu/Zeta in LatAm (G/L/M/Z is all or none, so only G is actually sufficient)
    tmp = (complete.cases(dat_proc[,"Day29"%.%must_have_assays]) | dat_proc$COL_variants_study) &
      with(dat_proc, 
           Region==0 & !is.na(Day29pseudoneutid50) |
           Region==1 & !is.na(Day29pseudoneutid50_Gamma) | # G/L/M/Z
           Region==2 & !is.na(Day29pseudoneutid50_Beta) 
      ) 
    dat_proc$TwophasesampIndD29variant[noncases] = tmp [noncases]
    # check if they all have variants bAb
    # 6 from RSA and 1 from LatAm needs bAb imputed
    # all those have variants bAb have variants nAb
    with(subset(dat_proc, SubcohortInd & EventIndPrimaryIncludeNotMolecConfirmedD1==0), 
         table(!is.na(Day29bindSpike_D614), TwophasesampIndD29variant, Region))
    
  }
  

} else if (study_name %in% c("VAT08")) {

  # baseline any nAb
  dat_proc$baseline.nAb = with(dat_proc, 
             !is.na(Bpseudoneutid50) | 
               !is.na(Bpseudoneutid50_B.1.351) | 
               !is.na(Bpseudoneutid50_BA.1) | 
               !is.na(Bpseudoneutid50_BA.2) | 
               !is.na(Bpseudoneutid50_BA.4.5)) 
  # D43 any nAb
  dat_proc$D43.nAb = with(dat_proc, 
             !is.na(Day43pseudoneutid50) | 
               !is.na(Day43pseudoneutid50_B.1.351) | 
               !is.na(Day43pseudoneutid50_BA.1) | 
               !is.na(Day43pseudoneutid50_BA.2) | 
               !is.na(Day43pseudoneutid50_BA.4.5)) 
  # D22 any nAb
  dat_proc$D22.nAb = with(dat_proc, 
             !is.na(Day22pseudoneutid50) | 
               !is.na(Day22pseudoneutid50_B.1.351) | 
               !is.na(Day22pseudoneutid50_BA.1) | 
               !is.na(Day22pseudoneutid50_BA.2) | 
               !is.na(Day22pseudoneutid50_BA.4.5)) 
  
  dat_proc[["TwophasesampIndD43nAb"]] = dat_proc$baseline.nAb & dat_proc$D43.nAb
  dat_proc[["TwophasesampIndD22nAb"]] = dat_proc$baseline.nAb & dat_proc$D22.nAb
  

  # baseline any bAb
  dat_proc$baseline.bAb = with(dat_proc, 
                                 !is.na(BbindSpike) | 
                                 !is.na(BbindSpike_beta) | 
                                 !is.na(BbindSpike_alpha) | 
                                 !is.na(BbindSpike_gamma) | 
                                 !is.na(BbindSpike_delta1) | 
                                 !is.na(BbindSpike_delta2) | 
                                 !is.na(BbindSpike_delta3) | 
                                 !is.na(BbindSpike_omicron)) 
  # D43 any bAb
  dat_proc$D43.bAb = with(dat_proc, 
                                 !is.na(Day43bindSpike) | 
                                 !is.na(Day43bindSpike_beta) | 
                                 !is.na(Day43bindSpike_alpha) | 
                                 !is.na(Day43bindSpike_gamma) | 
                                 !is.na(Day43bindSpike_delta1) | 
                                 !is.na(Day43bindSpike_delta2) | 
                                 !is.na(Day43bindSpike_delta3) | 
                                 !is.na(Day43bindSpike_omicron)) 
  # D22 any bAb
  dat_proc$D22.bAb = with(dat_proc, 
                            !is.na(Day22bindSpike) | 
                            !is.na(Day22bindSpike_beta) | 
                            !is.na(Day22bindSpike_alpha) | 
                            !is.na(Day22bindSpike_gamma) | 
                            !is.na(Day22bindSpike_delta1) | 
                            !is.na(Day22bindSpike_delta2) | 
                            !is.na(Day22bindSpike_delta3) | 
                            !is.na(Day22bindSpike_omicron)) 
  
  dat_proc[["TwophasesampIndD43bAb"]] = dat_proc$baseline.bAb & dat_proc$D43.bAb
  dat_proc[["TwophasesampIndD22bAb"]] = dat_proc$baseline.bAb & dat_proc$D22.bAb

    
  # # the original
  # # require baseline but not time point 1
  # dat_proc[["TwophasesampIndD"%.%timepoints[2]%.%"original"]] =
  #   with(dat_proc, SubcohortInd | !(is.na(get("EventIndPrimaryD"%.%timepoints[1])) | get("EventIndPrimaryD"%.%timepoints[1]) == 0)) &
  #   (dat_proc$EventIndFirstInfectionD1==0 & complete.cases(dat_proc[,c("B", "Day"%.%timepoints[2])%.%'bindSpike']) |
  #      dat_proc$EventIndFirstInfectionD1==1 & (complete.cases(dat_proc[,c("B", "Day"%.%timepoints[2])%.%'pseudoneutid50']) |
  #                                                complete.cases(dat_proc[,c("B", "Day"%.%timepoints[2])%.%'bindSpike'])))
  # 
  # # same as tp2
  # dat_proc[["TwophasesampIndD"%.%timepoints[1]%.%"original"]] = dat_proc[["TwophasesampIndD"%.%timepoints[2]%.%"original"]]
  
  
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
        
  
} else if (study_name=="COVAIL" ) {
  # the whole cohort is treated as ph1 and ph2
  dat_proc$TwophasesampIndD15 = dat_proc$ph1.D15 
  dat_proc$TwophasesampIndD29 = dat_proc$ph1.D29
  
  
} else stop("unknown study_name 8")



# weights 
if (TRIAL=='vat08_combined') {
  # two sets of weights, one for bAb and one for nAb, the one for bAb is also used for multivariate analyses
  for (tp in rev(timepoints)) { # rev is just so that digest passes
    
    # # for correlates studies, filter out 
    # # stage 1 naive and 
    # # countries based on lack of cases in both vacc and plac arms
    # # country 4 for stage 2 nnaive
    # # country 5, 6 for stage 1 nnaive
    # # country 2, 4 for stage 2 naive
    # kp = kp & with(dat_proc, 
    #   !(
    #         Trialstage==2 & Bserostatus==1 & Country %in% c(4)   | 
    #         Trialstage==1 & Bserostatus==1 & Country %in% c(5,6) |
    #         Trialstage==2 & Bserostatus==0 & Country %in% c(2,4) |
    #         Trialstage==1 & Bserostatus==0 
    #   ) 
    # )
    
    kp = with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp) >= 7)

    # nAb
    
    wts_table <- with(dat_proc[kp,], table(Wstratum.nAb, get("TwophasesampIndD"%.%tp%.%"nAb")))
    print(t(wts_table))
    
    unique(subset(dat_proc, Wstratum.nAb%in%c(66,76,130,140,150,160), c(tps.stratum.nAb, Wstratum.nAb, Trialstage, Bserostatus, Trt, cc, Senior))) 
    unique(subset(dat_proc, Wstratum.nAb%in%c(103,113), c(tps.stratum.nAb, Wstratum.nAb, Trialstage, Bserostatus, Trt, continent, Senior))) 
    
    # collapse nonseninor (i) and senior strata (i+10) to remove empty ph2 cells
    # note that since we merge Wstratum but not tps.stratum. tsp.stratum is not in sync with Wstratum
    i=66;  dat_proc$Wstratum.nAb[dat_proc$Wstratum.nAb %in% c(i,i+10)]=i   # Stage 1, NN, Vacc, Kenya
    i=130; dat_proc$Wstratum.nAb[dat_proc$Wstratum.nAb %in% c(i,i+10)]=i   # Stage 2, NN, Plac, Uganda
    i=150; dat_proc$Wstratum.nAb[dat_proc$Wstratum.nAb %in% c(i,i+10)]=i   # Stage 2, NN, Vacc, Uganda
    i=103; dat_proc$Wstratum.nAb[dat_proc$Wstratum.nAb %in% c(i,i+10)]=i   # Stage 2, N,  Vacc, Asia
    
    wts_table <- with(dat_proc[kp,], table(Wstratum.nAb, get("TwophasesampIndD"%.%tp%.%"nAb")))
    print(t(wts_table))
    # check again for empty cells
    stopifnot(all(wts_table[,2]!=0))
    
    wts_norm <- rowSums(wts_table) / wts_table[, 2]
    dat_proc[["wt.D"%.%tp%.%".nAb"]] <- wts_norm[dat_proc$Wstratum.nAb %.% ""]
    # the step above assigns weights for some subjects outside ph1. the next step makes them NA
    dat_proc[["wt.D"%.%tp%.%".nAb"]] = ifelse(kp, dat_proc[["wt.D"%.%tp%.%".nAb"]], NA) 
    dat_proc[["ph1.D"%.%tp]]=!is.na(dat_proc[["wt.D"%.%tp%.%".nAb"]])
    dat_proc[["ph2.D"%.%tp%.%".nAb"]]=dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD"%.%tp%.%"nAb"]]
    
    assertthat::assert_that(
      all(!is.na(subset(dat_proc, kp & !is.na(Wstratum.nAb))[["wt.D"%.%tp%.%'.nAb']])),
      msg = "missing wt.D for D analyses ph1 subjects")
    
    
    # bAb
    wts_table <- with(dat_proc[kp,], table(Wstratum.bAb, get("TwophasesampIndD"%.%tp%.%"bAb")))
    print(t(wts_table))
    
    unique(subset(dat_proc, Wstratum.bAb%in%c(2,13,71,93,113), c(tps.stratum.bAb, Wstratum.bAb, Trialstage, Bserostatus, Trt, continent, Senior))) 
    unique(subset(dat_proc, Wstratum.bAb%in%c(12,3,61,83,103), c(tps.stratum.bAb, Wstratum.bAb, Trialstage, Bserostatus, Trt, continent, Senior))) 
    
    # collapse nonseninor (i) and senior strata (i+10) to remove empty ph2 cells
    # note that since we merge Wstratum but not tps.stratum. tsp.stratum is not in sync with Wstratum
    i=2;   dat_proc$Wstratum.bAb[dat_proc$Wstratum.bAb %in% c(i,i+10)]=i   # Stage 1, N, Plac, Africa
    i=3;   dat_proc$Wstratum.bAb[dat_proc$Wstratum.bAb %in% c(i,i+10)]=i   # Stage 1, N, Plac, Asia
    i=61;  dat_proc$Wstratum.bAb[dat_proc$Wstratum.bAb %in% c(i,i+10)]=i   # Stage 1, NN, Vacc, LatAm
    i=83;  dat_proc$Wstratum.bAb[dat_proc$Wstratum.bAb %in% c(i,i+10)]=i   # Stage 2, N, Plac, Asia
    i=103; dat_proc$Wstratum.bAb[dat_proc$Wstratum.bAb %in% c(i,i+10)]=i   # Stage 2, N, Vacc, Asia
    
    wts_table <- with(dat_proc[kp,], table(Wstratum.bAb, get("TwophasesampIndD"%.%tp%.%"bAb")))
    print(t(wts_table))
    # check again for empty cells
    stopifnot(all(wts_table[,2]!=0))
    
    wts_norm <- rowSums(wts_table) / wts_table[, 2]
    dat_proc[["wt.D"%.%tp%.%".bAb"]] <- wts_norm[dat_proc$Wstratum.bAb %.% ""]
    # the step above assigns weights for some subjects outside ph1. the next step makes them NA
    dat_proc[["wt.D"%.%tp%.%".bAb"]] = ifelse(kp, dat_proc[["wt.D"%.%tp%.%".bAb"]], NA) 
    dat_proc[["ph1.D"%.%tp]]=!is.na(dat_proc[["wt.D"%.%tp%.%".bAb"]])
    dat_proc[["ph2.D"%.%tp%.%".bAb"]]=dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD"%.%tp%.%"bAb"]]
    
    assertthat::assert_that(
      all(!is.na(subset(dat_proc, kp & !is.na(Wstratum.bAb))[["wt.D"%.%tp%.%'.bAb']])),
      msg = "missing wt.D for D analyses ph1 subjects")
    
    
    
    # # original
    # tmp = with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp) >= 7)
    # wts_table <- with(dat_proc[tmp,], table(Wstratum.original, get("TwophasesampIndD"%.%tp%.%"original")))
    # print(wts_table)
    # stopifnot(all(wts_table[,2]!=0))
    # wts_norm <- rowSums(wts_table) / wts_table[, 2]
    # dat_proc[["wt.D"%.%tp%.%".original"]] <- wts_norm[dat_proc$Wstratum.original %.% ""]
    # # the step above assigns weights for some subjects outside ph1. the next step makes them NA
    # dat_proc[["wt.D"%.%tp%.%".original"]] = ifelse(with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp)>=7), 
    #                                                dat_proc[["wt.D"%.%tp%.%".original"]], NA) 
    # dat_proc[["ph1.D"%.%tp]]=!is.na(dat_proc[["wt.D"%.%tp%.%".original"]])
    # dat_proc[["ph2.D"%.%tp%.%".original"]]=dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD"%.%tp%.%"original"]]
  }
  
  
} else if (TRIAL=="covail" ) {
  # PP = no violation + marker available at d1 and d15
  # Immunemarkerset = PP & no infection between enrollment and D15+6
  # ph1.D15 = Immunemarkerset & arm!=3
  dat_proc[["ph1.D15"]]=dat_proc$ph1.D15 
  dat_proc[["ph2.D15"]]=dat_proc$ph2.D15
  dat_proc[["wt.D15"]] = 1
  dat_proc[["ph1.D92"]]=dat_proc$ph1.D92
  dat_proc[["ph2.D92"]]=dat_proc$ph2.D92
  dat_proc[["wt.D92"]] = 1
  dat_proc[["ph1.D29"]]=dat_proc$ph1.D29
  dat_proc[["ph2.D29"]]=dat_proc$ph2.D29
  dat_proc[["wt.D29"]] = 1
  

} else if (TRIAL=="janssen_partA_VL") {
  tp=29
  tmp = with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp) >= 7)
  wts_table <- with(dat_proc[tmp,], table(Wstratum, get("TwophasesampIndD"%.%tp)))
  print(wts_table)
  wts_norm <- rowSums(wts_table) / wts_table[, 2]
  dat_proc[["wt.D"%.%tp]] <- wts_norm[dat_proc$Wstratum %.% ""]
  # the step above assigns weights for some subjects outside ph1. the next step makes them NA
  dat_proc[["wt.D"%.%tp]] = ifelse(with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp)>=7), dat_proc[["wt.D"%.%tp]], NA) 
  dat_proc[["ph1.D"%.%tp]]=!is.na(dat_proc[["wt.D"%.%tp]])
  dat_proc[["ph2.D"%.%tp]]=dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD"%.%tp]]
  
  assertthat::assert_that(
    all(!is.na(subset(dat_proc, tmp & !is.na(Wstratum))[["wt.D"%.%tp]])),
    msg = "missing wt.D for D analyses ph1 subjects")
  

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
  

} else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")) {
  # EUA and PartA
  tp=29
  tmp = with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp) >= 7)
  wts_table <- with(dat_proc[tmp,], table(Wstratum, get("TwophasesampIndD"%.%tp)))
  print(wts_table)
  wts_norm <- rowSums(wts_table) / wts_table[, 2]
  dat_proc[["wt.D"%.%tp]] <- wts_norm[dat_proc$Wstratum %.% ""]
  # the step above assigns weights for some subjects outside ph1. the next step makes them NA
  dat_proc[["wt.D"%.%tp]] = ifelse(with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp)>=7), dat_proc[["wt.D"%.%tp]], NA) 
  dat_proc[["ph1.D"%.%tp]]=!is.na(dat_proc[["wt.D"%.%tp]])
  dat_proc[["ph2.D"%.%tp]]=dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampIndD"%.%tp]]
  
  assertthat::assert_that(
    all(!is.na(subset(dat_proc, tmp & !is.na(Wstratum))[["wt.D"%.%tp]])),
    msg = "missing wt.D for D analyses ph1 subjects")

  
  # Starting at 1 day post D29 visit
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

    
} else {
  # the default
  for (tp in rev(timepoints)) { # rev is just so that digest passes
    tmp = with(dat_proc, get("EarlyendpointD"%.%tp)==0 & Perprotocol==1 & get("EventTimePrimaryD"%.%tp) >= 7)
    wts_table <- with(dat_proc[tmp,], table(Wstratum, get("TwophasesampIndD"%.%tp)))
    print(wts_table)
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
}  
  

# immunogenicity weights and intercurrent weights
if (!TRIAL %in% c('vat08_combined','covail')) {
  
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
  
}


###############################################################################
# impute missing biomarkers in ph2 (assay imputation)
#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, each time point, but not Delta
###############################################################################

if (study_name%in%c("COVAIL")) {
  
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
  
  # populate dat_proc with the imputed values
  imp.markers='Day15pseudoneutid50_BA.4.BA.5'
  dat_proc[dat_proc[["ph1.D15"]]==1, imp.markers] <-
    dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["ph1.D15"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]
  
  assertthat::assert_that(
    all(complete.cases(dat_proc[dat_proc[["ph1.D15"]] == 1, imp.markers])),
    msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
  )
  
  
  #### impute D29 among those at risk at D29
  dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD29 & AsympInfectIndD15to29==0)
  with(dat.tmp.impute, print(table(!is.na(get("Day29"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
  # thus, no missingness actually
  
  #### impute D91 among those at risk at D91
  dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD91 & AsympInfectIndD15to91==0)
  with(dat.tmp.impute, print(table(!is.na(get("Day91"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
  # thus, no missingness actually
  
  #### impute D181 among those at risk at D181
  dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD181 & AsympInfectIndD15to181==0)
  with(dat.tmp.impute, print(table(!is.na(get("Day181"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
  # thus, no missingness actually
  
  
  
} else if(study_name%in%c("VAT08")) {

  n.imp <- 1

  # impute each time point separately
  for (i in 1:2) {
    for (tp in c("B","Day22","Day43")) {
      
      if (i==1) {
        # nAb
        if (tp=="B") kp = dat_proc$baseline.nAb
        if (tp=="Day22") kp = dat_proc$D22.nAb
        if (tp=="Day43") kp = dat_proc$D43.nAb
        imp.markers=tp%.%assays[startsWith(assays,"pseudoneut") & !endsWith(assays,"mdw")]
      } else {
        # bAb
        if (tp=="B") kp = dat_proc$baseline.bAb
        if (tp=="Day22") kp = dat_proc$D22.bAb
        if (tp=="Day43") kp = dat_proc$D43.bAb
        imp.markers=tp%.%assays[startsWith(assays,"bindSpike") & !endsWith(assays,"mdw")]
      }
      
      dat.tmp.impute <- dat_proc[kp,]
      
      
      for (trt in unique(dat_proc$Trt)) {
        for (sero in unique(dat_proc$Bserostatus)) {    
          # for (stage in unique(dat_proc$Trialstage)) {
          imp <- dat.tmp.impute %>% dplyr::filter(Trt == trt & Bserostatus==sero) %>% select(all_of(imp.markers))         
          if(any(is.na(imp))) {
            # if there is no variability, fill in NA with constant values
            for (a in names(imp)) {
              if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
            }            
            # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
            imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
            dat.tmp.impute[dat.tmp.impute$Trt == trt & dat.tmp.impute$Bserostatus == sero, imp.markers] <- mice::complete(imp, action = 1)
          } 
          # }
        }
      }
      
      # missing markers imputed properly?
      assertthat::assert_that(
        all(complete.cases(dat.tmp.impute[, imp.markers])),
        msg = "missing markers imputed properly?"
      )    
      
      # populate dat_proc imp.markers with the imputed values
      dat_proc[kp==1, imp.markers] <-
        dat.tmp.impute[imp.markers][match(dat_proc[kp==1, "Ptid"], dat.tmp.impute$Ptid), ]
      
      assertthat::assert_that(
        all(complete.cases(dat_proc[kp == 1, imp.markers])),
        msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
      )
      
    }    
  }

  
} else if (TRIAL=="janssen_partA_VL") {
  # in the data before the variants study, all those have ancestral nAb have ancestral bAb
  
  assays_panel19 = assays[startsWith(assays, "bindSpike_")]
  
  # three imputation steps are done

  # impute 98 non-case bindSpike values (COL) by 521 bindSpike_D614 values for the variants study samples
  n.imp = 1
  select = with(dat_proc, TwophasesampIndD29variant==1 & !is.na(Day29bindSpike_D614))
  imp.markers = c("Day29bindSpike", "Day29bindSpike_D614")
  imp <- dat_proc[select,] %>% select(all_of(imp.markers)) 
  summary(imp)
  nrow(tmp)
  imp = imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  dat_proc[select, "Day29bindSpike"] = mice::complete(imp, action=1)[,"Day29bindSpike"]
  assertthat::assert_that(
    all(complete.cases(dat_proc[select, "Day29bindSpike"])),
    msg = "janssen_partA_VL: imputed values of missing markers merged properly for all individuals in the two phase sample?"
  )
  
  
  # impute occasional non-case ptids with variants nAb measured but missing variants bAb
  
  # in LatAm, build model with 243 non-cases ancestral/G/L/M/Z nAb + ancestral bAb to predict 1 ptid's missing variants bAb
  n.imp = 1
  select = with(dat_proc, Region==1 & TwophasesampIndD29variant==1 & EventIndPrimaryIncludeNotMolecConfirmedD1==0)
  imp.markers = "Day29" %.% setdiff(assays, c("pseudoneutid50_Delta", "pseudoneutid50_Beta", "bindRBD", "bindN"))
  imp <- dat_proc[select,] %>% select(all_of(imp.markers))
  summary(imp)
  nrow(imp)
  subset(dat_proc, select & is.na(Day29bindSpike_D614), c(Region, Trt, Bserostatus, ph2.D29))
  
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)   
  for (a in assays_panel19) {
    dat_proc[select, "Day29"%.%a] = mice::complete(imp, action=1)[,"Day29"%.%a]
    assertthat::assert_that(
      all(complete.cases(dat_proc[select, "Day29"%.%a])),
      msg = "janssen_partA_VL: imputed values of missing markers merged properly for all individuals in the two phase sample?"
    )
  }
  
  # in RSA, we will use 94 non-cases ancestral/Beta nAb + ancestral bAb to predict 6 ptids' variants bAb 
  n.imp = 1
  select = with(dat_proc, Region==2 & TwophasesampIndD29variant==1 & EventIndPrimaryIncludeNotMolecConfirmedD1==0)
  imp.markers = "Day29" %.% setdiff(assays, c("pseudoneutid50_Zeta",      "pseudoneutid50_Mu",
                                              "pseudoneutid50_Gamma",     "pseudoneutid50_Lambda", "bindRBD", "bindN"))
  imp <- dat_proc[select,] %>% select(all_of(imp.markers)) 
  nrow(imp)
  summary(imp)
  
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)   
  for (a in assays_panel19) {
    dat_proc[select, "Day29"%.%a] = mice::complete(imp, action=1)[,"Day29"%.%a]
    assertthat::assert_that(
      all(complete.cases(dat_proc[select, "Day29"%.%a])),
      msg = "janssen_partA_VL: imputed values of missing markers merged properly for all individuals in the two phase sample?"
    )
  }
  
  
  # impute 10 copies of variants nAb and bAb for all cases with ancestral nAb
  
  n.imp <- 10
  
  # create placeholders for imputed ab markers
  for (i in 1:10) {
    dat_proc[["Day29pseudoneutid50_Beta"%.%"_"%.%i]] = NA
    dat_proc[["Day29pseudoneutid50_Gamma"%.%"_"%.%i]]  = NA
    dat_proc[["Day29pseudoneutid50_Lambda"%.%"_"%.%i]] = NA
    dat_proc[["Day29pseudoneutid50_Mu"%.%"_"%.%i]]     = NA
    dat_proc[["Day29pseudoneutid50_Zeta"%.%"_"%.%i]]   = NA
    for (a in assays_panel19) dat_proc[["Day29"%.%a%.%"_"%.%i]] = NA
  }
  
  # LatAm, impute four variants
  # Day29bindRBD not in the list b/c there are 98 missing values and there is no need to impute it 
  # TwophasesampIndD29variant includes all cases
  select = with(dat_proc, Trt == 1 & Bserostatus==0 & Region==1 & TwophasesampIndD29variant==1)
  imp.markers = c("Day29bindSpike", "Day29pseudoneutid50", "Day29pseudoneutid50_Gamma", "Day29pseudoneutid50_Lambda", "Day29pseudoneutid50_Mu", "Day29pseudoneutid50_Zeta", "Day29"%.%assays_panel19)
  imp <- dat_proc[select, ] %>% select(all_of(imp.markers)) 
  summary(imp)
  nrow(imp)
  
  imp = imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  # add 10 new columns for each of the variants to the dataset
  for (i in 1:n.imp) {
    dat_proc[select, "Day29pseudoneutid50_Gamma"%.%"_"%.%i]  = mice::complete(imp, action=i)[,"Day29pseudoneutid50_Gamma"]
    dat_proc[select, "Day29pseudoneutid50_Lambda"%.%"_"%.%i] = mice::complete(imp, action=i)[,"Day29pseudoneutid50_Lambda"]
    dat_proc[select, "Day29pseudoneutid50_Mu"%.%"_"%.%i]     = mice::complete(imp, action=i)[,"Day29pseudoneutid50_Mu"]
    dat_proc[select, "Day29pseudoneutid50_Zeta"%.%"_"%.%i]   = mice::complete(imp, action=i)[,"Day29pseudoneutid50_Zeta"]
    for (a in assays_panel19) dat_proc[select, "Day29"%.%a%.%"_"%.%i] = mice::complete(imp, action=i)[,"Day29"%.%a]
  } 
  assertthat::assert_that(
    all(complete.cases(dat_proc[select, c("Day29pseudoneutid50_Gamma"%.%"_"%.%1:10, "Day29pseudoneutid50_Lambda"%.%"_"%.%1:10, 
                                          "Day29pseudoneutid50_Mu"%.%"_"%.%1:10, "Day29pseudoneutid50_Zeta"%.%"_"%.%1:10)])),
    msg = "janssen_partA_VL: imputed values of missing markers merged properly for all individuals in the two phase sample?"
  )
  

  # RSA, impute Beta nAb for non-Beta cases, impute Delta nAb for non-Delta cases (the latter can be skipped if desired b/c not used for correlates)
  select = with(dat_proc, Trt == 1 & Bserostatus==0 & Region==2 & TwophasesampIndD29variant==1)
  imp.markers = c("Day29bindSpike", "Day29bindRBD", "Day29pseudoneutid50", "Day29pseudoneutid50_Beta", "Day29pseudoneutid50_Delta", "Day29"%.%assays_panel19)
  imp <- dat_proc[select,] %>% select(all_of(imp.markers)) 
  summary(imp)
  nrow(imp)
  
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  # add 10 new columns for each of the variants to the dataset
  for (i in 1:10) {
    dat_proc[select, "Day29pseudoneutid50_Beta"%.%"_"%.%i] = mice::complete(imp, action=i)[,"Day29pseudoneutid50_Beta"]
    dat_proc[select, "Day29pseudoneutid50_Delta"%.%"_"%.%i] = mice::complete(imp, action=i)[,"Day29pseudoneutid50_Delta"]
    for (a in assays_panel19) dat_proc[select, "Day29"%.%a%.%"_"%.%i] = mice::complete(imp, action=i)[,"Day29"%.%a]
  } 
  assertthat::assert_that(
    all(complete.cases(dat_proc[select, c("Day29pseudoneutid50_Beta"%.%"_"%.%1:10)])),
    all(complete.cases(dat_proc[select, c("Day29pseudoneutid50_Delta"%.%"_"%.%1:10)])),
    msg = "janssen_partA_VL: imputed values of missing markers merged properly for all individuals in the two phase sample?"
  )

    
  # remove Day29bindSpike_D614 from dataset so that there is only one ancestral bAb variable in the dataset
  dat_proc = subset(dat_proc, select=-Day29bindSpike_D614)
  
  
  # add 10 identical copies of ancestral markers to make programming easier (every marker has 10 copies)
  for (i in 1:10) {
    dat_proc[, "Day29pseudoneutid50"%.%"_"%.%i] = dat_proc[, "Day29pseudoneutid50"]
    dat_proc[, "Day29bindSpike"%.%"_"%.%i]      = dat_proc[, "Day29bindSpike"]
  } 
  
  
} else {
  # loop through the time points
  # first impute (B, D29, D57) among TwophasesampIndD57==1
  # next impute (B, D29) among TwophasesampIndD29==1
  
  for (tp in rev(timepoints)) {    
      n.imp <- 1
      
      # .nAb does not need to be further imputed
      if (TRIAL=='vat08_combined') dat_proc[["TwophasesampIndD"%.%tp]] <- dat_proc[["TwophasesampIndD"%.%tp%.%'bAb']]
  
      dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD"%.%tp) == 1)
      
      if(two_marker_timepoints) {
        imp.markers=c(outer(c("B", if(tp==timepoints[2]) "Day"%.%timepoints else "Day"%.%tp), assays, "%.%"))
      } else {
        imp.markers=c(outer(c("B", "Day"%.%tp), assays, "%.%"))
      }
      # mdw markers are not imputed
      imp.markers=imp.markers[!endsWith(imp.markers, "_mdw")]

      for (trt in unique(dat_proc$Trt)) {
      for (sero in unique(dat_proc$Bserostatus)) {    
        
        if (study_name=="VAT08") {
          # further separate by trial stage
          for (stage in unique(dat_proc$Trialstage)) {
            imp <- dat.tmp.impute %>% dplyr::filter(Trt == trt & Bserostatus==sero & Trialstage==stage) %>% select(all_of(imp.markers))         
            if(any(is.na(imp))) {
              # if there is no variability, fill in NA with constant values
              for (a in names(imp)) {
                if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
              }            
              # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
              imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
              dat.tmp.impute[dat.tmp.impute$Trt == trt & dat.tmp.impute$Bserostatus == sero & dat.tmp.impute$Trialstage == stage, imp.markers] <- mice::complete(imp, action = 1)
            } 
          }
          
        } else {
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
}


  
  
###############################################################################
# transformation of the markers
###############################################################################


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
# add mdw scores

if(study_name == "COVAIL") {
  
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
  
  
} else if(study_name == "VAT08") {
  # should be the same for vat08_combined and vat08_nAb
  bAb = setdiff(assays[startsWith(assays, "bindSpike")], c('bindSpike_mdw'))
  nAb = setdiff(assays[startsWith(assays, "pseudoneutid50")], c('pseudoneutid50_mdw'))
  
  # use D43, stage 2, nnaive, vaccine to derive weights
  kp = dat_proc$Trialstage==2 & dat_proc$Bserostatus==1 & dat_proc$Trt==1
  
  # bAb
  mdw.wt.bAb=tryCatch({
    tree.weight(cor(dat_proc[kp, "Day43"%.%bAb], use='complete.obs'))
  }, error = function(err) {
    print(err$message)
    rep(1/length(bAb), length(bAb))
  })
  write.csv(mdw.wt.bAb, file = here("data_clean", "csv", TRIAL%.%"_bAb_mdw_weights.csv"))
  # apply to all time points, to both stages, naive/nnaive, vaccine and placebo
  for (t in c("B", "Day"%.%timepoints)) { # , "Delta"%.%timepoints%.%"overB", "Delta43over22"
    dat_proc[, t%.%'bindSpike_mdw'] = as.matrix(dat_proc[, t%.%bAb]) %*% mdw.wt.bAb
  }
  
  # nAb
  mdw.wt.nAb=tryCatch({
    tree.weight(cor(dat_proc[kp, "Day43"%.%nAb], use='complete.obs'))
  }, error = function(err) {
    print(err$message)
    rep(1/length(nAb), length(nAb))
  })
  write.csv(mdw.wt.nAb, file = here("data_clean", "csv", TRIAL%.%"_nAb_mdw_weights.csv"))
  # apply to all time points, to both stages, naive/nnaive, vaccine and placebo
  for (t in c("B", "Day"%.%timepoints)) { # , "Delta"%.%timepoints%.%"overB", "Delta43over22"
    dat_proc[, t%.%'pseudoneutid50_mdw'] = as.matrix(dat_proc[, t%.%nAb]) %*% mdw.wt.nAb
  }
  
}




###############################################################################
# add delta for dat_proc
# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta


if (TRIAL %in% c("janssen_partA_VL")) {
  # skipping b/c there is no baseline data
  
} else {

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
  
  
  if (TRIAL=="covail") {
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
  }
  
}



###############################################################################
# add discrete markers

if (TRIAL=="covail") {
  # mRNA arms
  dat_proc$tmp = with(dat_proc, ph1.D15 & TrtonedosemRNA==1) 
  assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
  all.markers = c("B"%.%assays, "Day15"%.%assays, "Delta15overB"%.%assays)
  dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D15")
  
  # Sanofi arms
  dat_proc$tmp = with(dat_proc, ph1.D29 & TrtSanofi==1)
  assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
  all.markers = c("Day29"%.%assays, "Delta29overB"%.%assays)
  dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D29")
  
  # remove the temp ph2 column
  dat_proc$tmp = NULL

  
} else if (TRIAL=="janssen_partA_VL") {
  
  # use Latin America for ancestral
  # define a temp ph2 column
  dat_proc$tmp = with(dat_proc, Trt==1 & Bserostatus==0 & ph1.D29 & Region==1) 
  assays = c("pseudoneutid50", "bindSpike")
  all.markers = c("Day29"%.%assays)
  dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D29")
  
  # use Latin America for Lambda, Mu, Zeta, Gamma
  # define a temp ph2 column
  dat_proc$tmp = with(dat_proc, Trt==1 & Bserostatus==0 & ph2.D29variant & Region==1) 
  assays = c("pseudoneutid50_Zeta", "pseudoneutid50_Mu", "pseudoneutid50_Gamma", "pseudoneutid50_Lambda",
                                    "bindSpike_B.1.621", "bindSpike_P.1", "bindSpike_C.37")
  all.markers = c("Day29"%.%assays)
  dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="ph2.D29variant", wt.col.name="wt.D29variant")
  # cut the imputed markers
  cutpoints=attr(dat_proc, "marker.cutpoints")
  for (a in all.markers) {
    for (i in 1:10) dat_proc[[a%.%"_"%.%i%.%"cat"]]=factor(cut(dat_proc[[a%.%"_"%.%i]], breaks = c(-Inf, cutpoints[[a]], Inf)))
  }
  
  # use South Africa for Beta, Delta
  # define a temp ph2 column
  dat_proc$tmp = with(dat_proc, Trt==1 & Bserostatus==0 & ph2.D29variant & Region==2) 
  assays = c("pseudoneutid50_Beta", "pseudoneutid50_Delta",
             "bindSpike_B.1.351", "bindSpike_DeltaMDW")
  all.markers = c("Day29"%.%assays)
  dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D29variant")
  # cut the imputed markers
  cutpoints=attr(dat_proc, "marker.cutpoints")
  for (a in all.markers) {
    for (i in 1:10) dat_proc[[a%.%"_"%.%i%.%"cat"]]=factor(cut(dat_proc[[a%.%"_"%.%i]], breaks = c(-Inf, cutpoints[[a]], Inf)))
  }
  
  # remove the temp ph2 column
  dat_proc$tmp = NULL
  
  # add 10 identical copies of ancestral markers to make programming easier (every marker has 10 copies)
  for (i in 1:10) {
    dat_proc[, "Day29pseudoneutid50"%.%"_"%.%i%.%"cat"] = dat_proc[, "Day29pseudoneutid50cat"]
    dat_proc[, "Day29bindSpike"%.%"_"%.%i%.%"cat"]      = dat_proc[, "Day29bindSpikecat"]
  } 
}


###############################################################################
# add two synthetic ID50 markers by region for ensemble

if (study_name=="ENSEMBLE") {
  if(endsWith(TRIAL, "EUA")) {
    dat_proc$Day29pseudoneutid50la = ifelse(dat_proc$Day29pseudoneutid50-0.124 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50-0.124 )
    dat_proc$Day29pseudoneutid50sa = ifelse(dat_proc$Day29pseudoneutid50-0.556 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50-0.556 )
    
  } else if(endsWith(TRIAL, "partA")) {
    dat_proc$Day29pseudoneutid50la = ifelse(dat_proc$Day29pseudoneutid50-0.255 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50-0.255 )
    dat_proc$Day29pseudoneutid50sa = ifelse(dat_proc$Day29pseudoneutid50-0.458 < log10(lloqs["pseudoneutid50"]), log10(lloqs["pseudoneutid50"]/2), dat_proc$Day29pseudoneutid50-0.458 )
  }
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

if (TRIAL %in% c("profiscov", "profiscov_lvmn", "vat08_combined", "vat08_nAb")) {
    # no risk score for profiscov, but some have missing BMI
    n.imp <- 1
    dat.tmp.impute <- dat_proc
    
    if (TRIAL %in% c("profiscov", "profiscov_lvmn")) {
      imp.markers=c("HighRiskInd", "Sex", "Age", "BMI")
    } else if (TRIAL %in% c("vat08_combined", "vat08_nAb")) {     
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
  
} else if(study_name == "VAT08") {
  # add region variable for regression
  # Honduras (3), not Honduras for the Stage 1 trial nnaive
  # Mexico (9), Other/Else country for the Stage 2 trial naive
  # India (4), Mexico (9), Other/Else country for the Stage 2 trial nnaive
  dat_proc$Region = NA
  
  dat_proc$Region[dat_proc$Trialstage==1 & dat_proc$Bserostatus==0] = "Stage1naive"
  
  dat_proc$Region[dat_proc$Trialstage==1 & dat_proc$Bserostatus==1 & dat_proc$Country==3] = "HON_Stage1Nnaive"
  dat_proc$Region[dat_proc$Trialstage==1 & dat_proc$Bserostatus==1 & dat_proc$Country!=3] = "NotHON_Stage1Nnaive"
  
  dat_proc$Region[dat_proc$Trialstage==2 & dat_proc$Country==9] = "MEX_Stage2"
  dat_proc$Region[dat_proc$Trialstage==2 & dat_proc$Country==4] = "IND_Stage2"
  dat_proc$Region[dat_proc$Trialstage==2 & !dat_proc$Country%in%c(4,9)] = "NotMEXIND_Stage2"

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
         janssen_pooled_partA = "335d2628adb180d3d07745304d7bf603",
         janssen_partA_VL = "be70e58897d461c242f930d09bbbcd0a", 
         azd1222 = "f573e684800003485094c18120361663",
         azd1222_bAb = "fc3851aff1482901f079fb311878c172",
         prevent19 = "a4c1de3283155afb103261ce6ff8cec2",
         vat08_combined = "d82e4d1b597215c464002962d9bd01f7", 
         covail = "6f75d31eff6089a930784373e56ed8ae", 
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
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
