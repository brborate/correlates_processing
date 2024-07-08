Sys.setenv(TRIAL = "vat08_combined")
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

bAb.1 = assays[startsWith(assays, "bindSpike")]
nAb.1 = assays[startsWith(assays, "pseudoneutid50")]
bAb = setdiff(bAb.1, c('bindSpike_mdw'))
nAb = setdiff(nAb.1, c('pseudoneutid50_mdw'))

}


########################################################################################################
# 1. read mapped data with risk score added

{
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
region.1 = c( # stage 1
 "United States" = 1, "Japan" = 1, 
 "Colombia" = 2, "Honduras" = 2, 
 "Ghana" = 3, "Kenya" = 3, 
 "Nepal" = 4, "India" = 4)
region.2 = c( # stage 2
 "Colombia" = 1, "Mexico" = 1, 
 "Ghana" = 2, "Kenya" = 2, "Uganda" = 2,
 "Nepal" = 3, "India" = 3)
# first set it to stage 1 region, then change the region for stage 2 countries
dat_proc$region = region.1[dat_proc$cc] 
dat_proc$region = ifelse(dat_proc$Trialstage==2, region.2[dat_proc$cc], dat_proc$region)


# add risk score
load(file = paste0('riskscore_baseline/output/vat08_combined/inputFile_with_riskscore.RData'))
stopifnot(all(inputFile_with_riskscore$Ptid==dat_proc$Ptid))
dat_proc$risk_score = inputFile_with_riskscore$risk_score
dat_proc$standardized_risk_score = inputFile_with_riskscore$standardized_risk_score

# ptids with missing Bserostatus already filtered out in preprocess

# define ten copies of imputed Omicron indicator and event time variables based on seq1.variant.hotdeck1 etc
for (t in c(1,22,43)) {
 for (i in 1:10) {
   dat_proc[[paste0("EventIndOmicronD",t,"M12hotdeck",i)]]  = 
     ifelse(!is.na(dat_proc[["seq1.variant.hotdeck"%.%i]]) & dat_proc[["seq1.variant.hotdeck"%.%i]]=="Omicron" & !is.na(dat_proc[["EventIndFirstInfectionD"%.%t]]),
            1,
            0)
   
   dat_proc[[paste0("EventTimeOmicronD",t,"M12hotdeck",i)]] = 
     ifelse(dat_proc[[paste0("EventIndOmicronD",t,"M12hotdeck",i)]] ==1, 
        pmin(dat_proc[["EventTimeKnownLineageOmicronD"%.%t]],    dat_proc[["EventTimeMissingLineageD"%.%t]]),
        pmax(dat_proc[["EventTimeKnownLineageNonOmicronD"%.%t]], dat_proc[["EventTimeMissingLineageD"%.%t]]))
 }
}

# create event time and indicator variables censored after M6 (instead of M12) post dose 2
for (t in c(1,22,43)) {
 for (i in 1:10) {
   dat_proc[[paste0("EventIndOmicronD",t,"M6hotdeck",i)]]  = ifelse (dat_proc[[paste0("EventTimeOmicronD43M12hotdeck",i)]]>180-21, 0,   dat_proc[[paste0("EventIndOmicronD",t,"M12hotdeck",i)]])
   dat_proc[[paste0("EventTimeOmicronD",t,"M6hotdeck",i)]] = ifelse (dat_proc[[paste0("EventTimeOmicronD43M12hotdeck",i)]]>180-21, 180-21, dat_proc[[paste0("EventTimeOmicronD",t,"M12hotdeck",i)]])
 }
}
   
}



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

{
# Bstratum: randomization strata
dat_proc$Bstratum =  with(dat_proc, Senior + 1)
names(Bstratum.labels) <- Bstratum.labels

# demographics
dat_proc$demo.stratum = dat_proc$region


# tps stratum, used in tps regression and to define Wstratum

# Wstratum
# Used to compute sampling weights. 
# Differs from tps stratum in that cases are separated out
# A case will have a Wstratum even if its tps.stratum is NA

dat_proc$tps.stratum = dat_proc$region
# 1-4: Non-senior; 5-8: Senior
# for stage 2, strata 4, 8 etc are empty because there are only three regions
cond=dat_proc$Senior==1
dat_proc$tps.stratum[cond] = dat_proc$tps.stratum[cond] + 4

# make cases 10 in Wstratum
dat_proc$Wstratum =  dat_proc$tps.stratum
cond=dat_proc$EventIndPrimaryD22==1
dat_proc$Wstratum[cond] = 10

# cross with trt and stage

# 1-10: vaccine; 11-20: placebo
cond = dat_proc$Trt==0
dat_proc$tps.stratum[cond] = dat_proc$tps.stratum[cond] + 10
dat_proc$Wstratum[cond] = dat_proc$Wstratum[cond] + 10
# 1-20: stage 1; 51-70: stage 2
cond = dat_proc$Trialstage==2
dat_proc$tps.stratum[cond] = dat_proc$tps.stratum[cond] + 50
dat_proc$Wstratum[cond] = dat_proc$Wstratum[cond] + 50

# cross with naive status

# Bserostatus==0: <100

# Bserostatus==1 & RAPDIAG==0: 101-170
cond = dat_proc$Bserostatus==1 & dat_proc$RAPDIAG=="NEGATIVE"
dat_proc$tps.stratum[cond] = dat_proc$tps.stratum[cond] + 100
dat_proc$Wstratum[cond] = dat_proc$Wstratum[cond] + 100
# Bserostatus==1 & RAPDIAG==1: 201-270
cond = dat_proc$Bserostatus==1 & dat_proc$RAPDIAG=="POSITIVE"
dat_proc$tps.stratum[cond] = dat_proc$tps.stratum[cond] + 200
dat_proc$Wstratum[cond] = dat_proc$Wstratum[cond] + 200

mytable(dat_proc$tps.stratum)
mytable(dat_proc$Wstratum)


with(subset(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Trialstage==1), table.prop(!is.na(Day43bindSpike), Bserostatus))
with(subset(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Trialstage==1), table.prop(!is.na(Day43bindSpike), Bserostatus))

with(subset(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Trialstage==2), table.prop(!is.na(Day43bindSpike), Bserostatus))
with(subset(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Trialstage==2), table.prop(!is.na(Day43bindSpike), Bserostatus))

with(subset(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Trialstage==1), table.prop(!is.na(Day43bindSpike), RAPDIAG))
with(subset(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Trialstage==1), table.prop(!is.na(Day43bindSpike), RAPDIAG))

with(subset(dat_proc, EventIndPrimaryD22==1 & Trt==0 & Trialstage==2), table.prop(!is.na(Day43bindSpike), RAPDIAG))
with(subset(dat_proc, EventIndPrimaryD22==1 & Trt==1 & Trialstage==2), table.prop(!is.na(Day43bindSpike), RAPDIAG))

table(dat_proc$tps.stratum) 
table(dat_proc$Wstratum)


# all strata for cases ends with 10
# there are 12 of them: Stage (2) x Trt (2) x neg/pos (3)
tmp=unique(dat_proc$Wstratum)
sort(tmp[tmp%%10==0])

}



###############################################################################
# 4. Define ph1, ph2, and weights
# Note that Wstratum may have NA if any variables to form strata has NA

for (tp in timepoints) {
  dat_proc[["ph1.D"%.%tp]] = with(dat_proc, 
                                  Perprotocol==1
                                  & get("EarlyinfectionD"%.%tp)==0
                                  & get("EventTimePrimaryD"%.%tp) >= 7)
}

# for (tp in timepoints) {
#   dat_proc[["ph1.D"%.%tp%.%".st2"]] = with(dat_proc, 
#                                   Perprotocol==1
#                                   & get("EarlyinfectionD"%.%tp)==0
#                                   & get("EventTimePrimaryD"%.%tp) >= 7
#                                   & Trialstage==2)
# }


{# bAb, any nAb
  
# baseline 
dat_proc$baseline.bAb = with(dat_proc, 
                             !is.na(BbindSpike) | 
                               !is.na(BbindSpike_beta) | 
                               !is.na(BbindSpike_alpha) | 
                               !is.na(BbindSpike_gamma) | 
                               !is.na(BbindSpike_delta1) | 
                               !is.na(BbindSpike_delta2) | 
                               !is.na(BbindSpike_delta3) | 
                               !is.na(BbindSpike_omicron)) 
# D43 
dat_proc$D43.bAb = with(dat_proc, 
                        !is.na(Day43bindSpike) | 
                          !is.na(Day43bindSpike_beta) | 
                          !is.na(Day43bindSpike_alpha) | 
                          !is.na(Day43bindSpike_gamma) | 
                          !is.na(Day43bindSpike_delta1) | 
                          !is.na(Day43bindSpike_delta2) | 
                          !is.na(Day43bindSpike_delta3) | 
                          !is.na(Day43bindSpike_omicron)) 
# D22 
dat_proc$D22.bAb = with(dat_proc, 
                        !is.na(Day22bindSpike) | 
                          !is.na(Day22bindSpike_beta) | 
                          !is.na(Day22bindSpike_alpha) | 
                          !is.na(Day22bindSpike_gamma) | 
                          !is.na(Day22bindSpike_delta1) | 
                          !is.na(Day22bindSpike_delta2) | 
                          !is.na(Day22bindSpike_delta3) | 
                          !is.na(Day22bindSpike_omicron)) 

dat_proc[["TwophasesampIndbAb"]] = dat_proc$baseline.bAb & dat_proc$D22.bAb & dat_proc$D43.bAb
}



{# nAb, any variant nAb. don't include ancestral b/c there is a batch 0 for ancestral which we don't want to include

# baseline 
dat_proc$baseline.nAb = with(dat_proc, 
                               !is.na(Bpseudoneutid50_B.1.351) | 
                               !is.na(Bpseudoneutid50_BA.1) | 
                               !is.na(Bpseudoneutid50_BA.2) | 
                               !is.na(Bpseudoneutid50_BA.4.5)) 
# D43 
dat_proc$D43.nAb = with(dat_proc, 
                          !is.na(Day43pseudoneutid50_B.1.351) | 
                          !is.na(Day43pseudoneutid50_BA.1) | 
                          !is.na(Day43pseudoneutid50_BA.2) | 
                          !is.na(Day43pseudoneutid50_BA.4.5)) 
# D22 
dat_proc$D22.nAb = with(dat_proc, 
                          !is.na(Day22pseudoneutid50_B.1.351) | 
                          !is.na(Day22pseudoneutid50_BA.1) | 
                          !is.na(Day22pseudoneutid50_BA.2) | 
                          !is.na(Day22pseudoneutid50_BA.4.5)) 

# require availability of data at all 3 time points
dat_proc[["TwophasesampIndnAb"]] = dat_proc$baseline.nAb & dat_proc$D43.nAb & dat_proc$D22.nAb
}

# force TwophasesampIndbAb to be the intersection of TwophasesampIndbAb and TwophasesampIndnAb because very few ptids with bAb don't have nAb
# this way for multivariate analysis, we just need to use TwophasesampIndbAb and there is no need to define another set of variables
dat_proc[["TwophasesampIndbAb"]] = dat_proc[["TwophasesampIndbAb"]] & dat_proc[["TwophasesampIndnAb"]]



{
  
  # generate SubcohortInd for bAb and nAb separately
  dat_proc$SubcohortIndbAb=0
  dat_proc$SubcohortIndnAb=0
  
  # all non-cases with markers are included
  dat_proc$SubcohortIndbAb[dat_proc$EventIndFirstInfectionD1==0 & dat_proc$TwophasesampIndbAb==1]=1
  dat_proc$SubcohortIndnAb[dat_proc$EventIndFirstInfectionD1==0 & dat_proc$TwophasesampIndnAb==1]=1
  
  # pick controls using bAb and use for both SubcohortIndbAb and SubcohortIndnAb
  tab=mytable(dat_proc$TwophasesampIndbAb, dat_proc$tps.stratum, dat_proc$EventIndFirstInfectionD1)[,,1]
  px=tab[2,]/(tab[1,]+tab[2,])
  px
  
  tab=mytable(dat_proc$TwophasesampIndbAb, dat_proc$tps.stratum, dat_proc$EventIndFirstInfectionD1)[,,2]
  tab
  n.cases.to.sample=round((tab[1,]+tab[2,])*px)
  
  # sample cases to include stratum by stratum based on bAb
  ii = sort(unique(dat_proc$tps.stratum))
  ii
  picks=c()
  for (i in ii) {
    picks=c(picks, sample(subset(dat_proc, TwophasesampIndbAb==1 & tps.stratum==i & EventIndFirstInfectionD1==1, Ptid, drop=T))[1:n.cases.to.sample[i%.%""]])
  }
  
  dat_proc$SubcohortIndbAb[dat_proc$Ptid %in% picks]=1
  dat_proc$SubcohortIndnAb[dat_proc$Ptid %in% picks]=1
  
  
}

# a helper function
get.strata.merge.to = function(strata.to.merge) {
  # check if sorted
  if(!all(strata.to.merge==sort(strata.to.merge))) stop("strata.to.merge needs to be sorted")
  
  # first try merging senior and non-senior strata, which differ by 4
  # luckily, this is all we have to do, i.e. stop is never triggered
  strata.merge.to=sapply (strata.to.merge, function(x) {
    if (x %% 10 <= 4) {
      if ((x+4) %in% strata.to.merge) {
        stop("both present: ", x, " and ", x+4, "\n")
      } else {
        x+4
      }
      
    } else {
      if ((x-4) %in% strata.to.merge) {
        stop("both present: ", x, " and ", x-4, "\n")
      } else {
        x-4
      }
    }
  })
  
  strata.merge.to
}


#### Define ph2 and wt variables 


# use a single set of strata for bAb and nAb, or correlates and immuno, for sensitivity (stage 2) and main
if (TRUE) {

  # use bAb instead of nAb because there are less bAb samples
  Ab="bAb"
  tp=43 # use D43 for this. D22 will also likely be fine
  dat_proc[["ph2.D"%.%tp%.%"."%.%Ab]] = dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampInd"%.%Ab]]
  wts_table <- with(dat_proc[dat_proc[["ph1.D"%.%tp]]==1, ], table(Wstratum, get("ph2.D"%.%tp%.%"."%.%Ab)))
  strata.to.merge.1 = sort(as.integer(rownames(wts_table[wts_table[,2]==0, ,drop=F])))
  print(strata.to.merge.1)
  
  # sensitivity study in stage 2
  tp=43
  dat_proc[["ph2.D"%.%tp%.%".st2.nAb.sen"]] = dat_proc[["ph1.D"%.%tp]] & dat_proc$Trialstage==2 & dat_proc[["TwophasesampIndnAb"]] & dat_proc$nAbBatch==2
  wts_table <- with(dat_proc[dat_proc[["ph1.D"%.%tp]]==1 & dat_proc$Trialstage==2, ], 
                    table(Wstratum, get("ph2.D"%.%tp%.%".st2.nAb.sen")))
  strata.to.merge.2 = sort(as.integer(rownames(wts_table[wts_table[,2]==0, ,drop=F])))
  print(strata.to.merge.2)
    
  strata.to.merge = sort(unique(c(strata.to.merge.1, strata.to.merge.2)))
  strata.merge.to = get.strata.merge.to (strata.to.merge)
  print(strata.merge.to)
  
  # merge Wstratum and tps.stratum. Note that there are no case strata to merge
  for (i in 1:length(strata.to.merge)) {
    dat_proc$Wstratum[dat_proc$Wstratum==strata.to.merge[i]] = strata.merge.to[i]
    dat_proc$tps.stratum[dat_proc$tps.stratum==strata.to.merge[i]] = strata.merge.to[i]
  }
  
  
  # compute weights
  
  # cor weights
  for (k in 1:2) { #1: bAb; 2: nAb
    if (k==1) Ab="bAb"; if (k==2) Ab="nAb"
    # create weights
    for (tp in timepoints) {
      dat_proc[["ph2.D"%.%tp%.%"."%.%Ab]] = dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampInd"%.%Ab]]
      dat_proc = add.wt(dat_proc, 
                        ph1="ph1.D"%.%tp, 
                        ph2="ph2.D"%.%tp%.%"."%.%Ab, 
                        Wstratum="Wstratum", 
                        wt="wt.D"%.%tp%.%"."%.%Ab, verbose=F) 
    }
  }  
  
  # stage 2 sensitivity weights for nAb
  for (tp in timepoints) {
    dat_proc[["ph1.D"%.%tp%.%".st2"]]         = dat_proc[["ph1.D"%.%tp]] & dat_proc$Trialstage==2
    dat_proc[["ph2.D"%.%tp%.%".st2.nAb.sen"]] = dat_proc[["ph1.D"%.%tp]] & dat_proc$Trialstage==2 & dat_proc[["TwophasesampIndnAb"]] & dat_proc$nAbBatch==2
    dat_proc = add.wt(dat_proc, 
                      ph1="ph1.D"%.%tp%.%".st2", 
                      ph2="ph2.D"%.%tp%.%".st2.nAb.sen", 
                      Wstratum="Wstratum", 
                      wt="wt.D"%.%tp%.%".st2.nAb.sen", verbose=F) 
  }
  
  # immuno
  dat_proc[["ph1.immuno"]] = with(dat_proc, Perprotocol==1 & EarlyinfectionD43==0)
  
  dat_proc[["ph2.immuno.nAb"]] = dat_proc$ph1.D43 & dat_proc$SubcohortIndnAb 
  dat_proc = add.wt(dat_proc, 
                    ph1="ph1.immuno", 
                    ph2="ph2.immuno.nAb", 
                    Wstratum="tps.stratum", 
                    wt="wt.immuno.nAb", verbose=F) 
  
  dat_proc[["ph2.immuno.bAb"]] = dat_proc$ph1.D43 & dat_proc$SubcohortIndbAb
  dat_proc = add.wt(dat_proc, 
                    ph1="ph1.immuno", 
                    ph2="ph2.immuno.bAb", 
                    Wstratum="tps.stratum", 
                    wt="wt.immuno.bAb", verbose=F) 
  

    
} else {
    
    
  dat_proc$Wstratum.bAb = dat_proc$Wstratum
  dat_proc$Wstratum.nAb = dat_proc$Wstratum
  
  for (k in 1:2) { #1: bAb; 2: nAb
    
    if (k==1) Ab="bAb"; if (k==2) Ab="nAb"
    
    # need to merge Wstratum? 
    tp=43 # use D43 for this. D22 will also likely be fine
    dat_proc[["ph2.D"%.%tp%.%"."%.%Ab]] = dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampInd"%.%Ab]]
    wts_table <- with(dat_proc[dat_proc[["ph1.D"%.%tp]]==1, ], table(Wstratum, get("ph2.D"%.%tp%.%"."%.%Ab)))
    strata.to.merge = sort(as.integer(rownames(wts_table[wts_table[,2]==0, ,drop=F])))
    print(strata.to.merge)
    strata.merge.to = get.strata.merge.to (strata.to.merge)
    print(strata.merge.to)
  
    # perform merging
    for (i in 1:length(strata.to.merge)) {
      if (k==1) {
        dat_proc$Wstratum.bAb[dat_proc$Wstratum.bAb==strata.to.merge[i]] = strata.merge.to[i]
      } else {
        dat_proc$Wstratum.nAb[dat_proc$Wstratum.nAb==strata.to.merge[i]] = strata.merge.to[i]
      }
    }
    
    # create weights
    for (tp in timepoints) {
      dat_proc[["ph2.D"%.%tp%.%"."%.%Ab]] = dat_proc[["ph1.D"%.%tp]] & dat_proc[["TwophasesampInd"%.%Ab]]
      dat_proc = add.wt(dat_proc, 
                        ph1="ph1.D"%.%tp, 
                        ph2="ph2.D"%.%tp%.%"."%.%Ab, 
                        Wstratum="Wstratum."%.%Ab, 
                        wt="wt.D"%.%tp%.%"."%.%Ab, verbose=F) 
    }
  }
  
  
  {#### (2) Define weights for sensitivity study in stage 2 non-naive using batch 2 nAb data only 
    
  # need to merge Wstratum? 
  tp=43
  dat_proc[["ph2.D"%.%tp%.%".st2.nAb.sen"]] = dat_proc[["ph1.D"%.%tp]] & dat_proc$Trialstage==2 & dat_proc[["TwophasesampIndnAb"]] & dat_proc$nAbBatch==2
  wts_table <- with(dat_proc[dat_proc[["ph1.D"%.%tp]]==1 & dat_proc$Trialstage==2, ], 
                    table(Wstratum, get("ph2.D"%.%tp%.%".st2.nAb.sen")))
  strata.to.merge = sort(as.integer(rownames(wts_table[wts_table[,2]==0, ,drop=F])))
  print(strata.to.merge)
  strata.merge.to = get.strata.merge.to (strata.to.merge)
  print(strata.merge.to)
  
  # yes, perform merging
  dat_proc$Wstratum.st2.nAb.sen = dat_proc$Wstratum
  for (i in 1:length(strata.to.merge)) {
    dat_proc$Wstratum.st2.nAb.sen[dat_proc$Wstratum.st2.nAb.sen==strata.to.merge[i]] = strata.merge.to[i]
  }
  
  # create weights
  for (tp in timepoints) {
    dat_proc[["ph1.D"%.%tp%.%".st2"]]         = dat_proc[["ph1.D"%.%tp]] & dat_proc$Trialstage==2
    dat_proc[["ph2.D"%.%tp%.%".st2.nAb.sen"]] = dat_proc[["ph1.D"%.%tp]] & dat_proc$Trialstage==2 & dat_proc[["TwophasesampIndnAb"]] & dat_proc$nAbBatch==2
    dat_proc = add.wt(dat_proc, 
                      ph1="ph1.D"%.%tp%.%".st2", 
                      ph2="ph2.D"%.%tp%.%".st2.nAb.sen", 
                      Wstratum="Wstratum.st2.nAb.sen", 
                      wt="wt.D"%.%tp%.%".st2.nAb.sen", verbose=F) 
  }
  
  }
  
  # remove Wstratum because it is not used
  dat_proc$Wstratum = NULL
  
  
  
  {
  #### Define weight computation variables for immuno analysis
    
  # create weights for bAb and nAb separately
  tp=43 # use D43 for immuno by our convention
  
  # need to merge stratum? 
  # use SubcohortIndbAb, but same results when using SubcohortIndnAb
  wts_table <- with(dat_proc[dat_proc$ph1.D43==1, ], table(tps.stratum, SubcohortIndbAb))
  strata.to.merge = sort(as.integer(rownames(wts_table[wts_table[,2]==0, ,drop=F])))
  print(strata.to.merge)
  
  # yes, perform merging
  strata.merge.to = get.strata.merge.to (strata.to.merge)
  dat_proc$tps.stratum.immuno = dat_proc$tps.stratum
  for (i in 1:length(strata.to.merge)) {
    dat_proc$tps.stratum.immuno[dat_proc$tps.stratum.immuno==strata.to.merge[i]] = strata.merge.to[i]
  }
  
  
  dat_proc[["ph1.immuno"]] = with(dat_proc, Perprotocol==1 & get("EarlyinfectionD"%.%tp)==0)
  
  
  dat_proc[["ph2.immuno.nAb"]] = dat_proc$ph1.D43 & dat_proc$SubcohortIndnAb 
  dat_proc = add.wt(dat_proc, 
                    ph1="ph1.immuno", 
                    ph2="ph2.immuno.nAb", 
                    Wstratum="tps.stratum.immuno", 
                    wt="wt.immuno.nAb", verbose=F) 
  
  dat_proc[["ph2.immuno.bAb"]] = dat_proc$ph1.D43 & dat_proc$SubcohortIndbAb
  dat_proc = add.wt(dat_proc, 
                    ph1="ph1.immuno", 
                    ph2="ph2.immuno.bAb", 
                    Wstratum="tps.stratum.immuno", 
                    wt="wt.immuno.bAb", verbose=F) 
  
  }
  
  
}






###############################################################################
# 5. impute missing biomarkers in ph2 (assay imputation)
#     separately within stage 1 and 2, vaccine and placebo, baseline pos and neg

n.imp <- 1

for (step in 1:2) { 
  # step 1: bAb, stage 1 and 2
  # step 2: nAb, stage 1 and 2

  for (tp in c("B","Day22","Day43")) {
    # impute each time point separately
  
    # step=1; tp="Day22"
    
    if (step==1) {
      # bAb
      if (tp=="B")     kp = dat_proc$baseline.bAb
      if (tp=="Day22") kp = dat_proc$D22.bAb
      if (tp=="Day43") kp = dat_proc$D43.bAb
      imp.markers=tp %.% bAb
      
    } else if (step==2) {
      # nAb
      if (tp=="B")     kp = dat_proc$baseline.nAb 
      if (tp=="Day22") kp = dat_proc$D22.nAb 
      if (tp=="Day43") kp = dat_proc$D43.nAb 
      imp.markers=tp %.% nAb
      
    } else stop ("wrong step")

    dat.tmp.impute <- dat_proc[kp,]
    
    for (trt in unique(dat.tmp.impute$Trt)) {
      for (sero in unique(dat.tmp.impute$Bserostatus)) {  
        # note that has to use dat.tmp.impute$Trialstage
        for (stage in unique(dat.tmp.impute$Trialstage)) {
          # trt=1; sero=0; stage=1
          
          imp <- dat.tmp.impute %>% dplyr::filter(Trt == trt & Bserostatus==sero & Trialstage==stage) %>% select(all_of(imp.markers))         
          if(any(is.na(imp))) {
            # if there is no variability, fill in NA with constant values
            for (a in names(imp)) {
              if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
            }            
            # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
            imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
            dat.tmp.impute[dat.tmp.impute$Trt == trt & dat.tmp.impute$Bserostatus == sero & dat.tmp.impute$Trialstage==stage, 
                           imp.markers] <- mice::complete(imp, action = 1)
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
    dat_proc[kp==1, imp.markers] <-
      dat.tmp.impute[imp.markers][match(dat_proc[kp==1, "Ptid"], dat.tmp.impute$Ptid), ]
    
    assertthat::assert_that(
      all(complete.cases(dat_proc[kp == 1, imp.markers])),
      msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
    )
    
  }    
}

  
  
###############################################################################
# 6. transformation of the markers


# # batch normalization/batch correction
# 
# nassays=c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1", "pseudoneutid50_BA.2", "pseudoneutid50_BA.4.5")
# 
# # save copies of markers that need to be corrected: nAb
# for (a in nassays) {
# for (t in c("Day43", "Day22", "B")) {
#   dat_proc[, t%.%a%.%"_save"] = dat_proc[, t%.%a] 
# }
# }    
#     
# chat = cbind(st1=c(0.1, 0.2, 0.25, 0.45, 0.25), st2=c(0.2, 0.35, 0.35, 0.3, 0.25))
# rownames(chat)=nassays
#   
# for (st in 1:2) {
# for (a in nassays) {
# for (t in c("Day43", "Day22", "B")) {
#   # shift by \hat{c}
#   kp = dat_proc$Trialstage==st & dat_proc$Bserostatus==1 & dat_proc$nAbBatch==1 & !is.na(dat_proc$nAbBatch)
#   dat_proc[kp, t%.%a] = dat_proc[kp, t%.%a] - chat[a, st]
#   # censoring by lod
#   dat_proc[kp, t%.%a] = ifelse (dat_proc[kp, t%.%a] < log10(llods[a]),
#                                 log10(llods[a]/2),
#                                 dat_proc[kp, t%.%a])
# }
# }
# }



###############################################################################
# 7. add mdw scores

# use D43, stage 2, nnaive, vaccine to derive weights for use in all scenarios

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
  





###############################################################################
# 8. add delta for dat_proc
# assuming data has been censored at the lower limit. Thus, there is no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta

assays1=assays
# keep mdw, because same weights are used for different time points, thus score of delta is same as delta of score

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



###############################################################################
# 9. add discrete/trichotomized markers

for (trt in c(0,1)) {
  for (sero in c(0,1)) {  
    # put 43 in front such that that way Bmarkercat will be defined for more ptids
    for (tp in c("43","22")) {
      
      # bAb + mdw
      for (stage in c(1,2)) {
        myprint(trt, sero, tp, stage)
        dat_proc$tmp = with(dat_proc, Trialstage==stage & Trt==trt & Bserostatus==sero & get("ph2.D"%.%tp%.%".bAb")) 
        dat_proc = add.trichotomized.markers (dat_proc, 
                                              c("Day"%.%tp%.%bAb.1, "Delta"%.%tp%.%"overB"%.%bAb.1, "B"%.%bAb.1), 
                                              ph2.col.name="tmp", 
                                              wt.col.name="wt.D"%.%tp%.%".bAb", verbose=T)
        dat_proc$tmp = NULL
      }
      
      # nAb + mdw
      for (stage in c(1,2)) {
        dat_proc$tmp = with(dat_proc, Trialstage==stage & Trt==trt & Bserostatus==sero & get("ph2.D"%.%tp%.%".nAb")) 
        dat_proc = add.trichotomized.markers (dat_proc, 
                                              c("Day"%.%tp%.%nAb.1, "Delta"%.%tp%.%"overB"%.%nAb.1, "B"%.%nAb.1), 
                                              ph2.col.name="tmp", 
                                              wt.col.name="wt.D"%.%tp%.%".nAb", verbose=T)
        dat_proc$tmp = NULL
      }
      
    }
  }
}




###############################################################################
# 10. impute covariates if necessary
# do this last so as not to change earlier values
###############################################################################

n.imp <- 1
dat.tmp.impute <- dat_proc

imp.markers=c("FOI", "risk_score", "HighRiskInd", "Sex", "Age", "BMI")

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
###############################################################################





###############################################################################
# digest check
###############################################################################

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         vat08_combined = "0d9bb40b71e5ab1a8496342c2d68b862", 
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))



print("run time: "%.%format(Sys.time()-begin, digits=1))
