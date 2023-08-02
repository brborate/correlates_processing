renv::activate(here::here())

Sys.setenv(TRIAL = "moderna_boost")
TRIAL=Sys.getenv("TRIAL")

library(here)
library(dplyr)
library(kyotil)
library(mice)

config <- config::get(config = TRIAL)


###############################################################################
# combine stage1 analysis-ready dataset and stage 2 mapped dataset

# read stage 2 mapped data with risk score made from make riskscore_analysis
load('riskscore_baseline/output/moderna_boost/inputFile_with_riskscore.rda')
dat_stage2_mapped = inputFile_with_riskscore
dat_stage2_mapped$naive = 1-dat_stage2_mapped$nnaive
# write to a csv file for Dean
mywrite.csv(dat_stage2_mapped, file=paste0("/trials/covpn/p3001/analysis/mapping_immune_correlates/Part_C_Unblinded_Phase_Data/adata/COVID_Moderna_stage2_", format(Sys.Date(), "%Y%m%d"), "_withRiskScores"))
print(paste0("write /trials/covpn/p3001/analysis/mapping_immune_correlates/Part_C_Unblinded_Phase_Data/adata/COVID_Moderna_stage2_", format(Sys.Date(), "%Y%m%d"), "_withRiskScores"))

# read stage1 analysis ready data
dat_stage1_adata = read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_with_riskscore.csv")
# remove two risk-related columns that are from stage 1 since we are only need one risk-related column from stage 2
dat_stage1_adata= subset(dat_stage1_adata, select=-c(Riskscorecohortflag,standardized_risk_score))

# dat_stage2_mapped has about 14K rows while dat_stage1_adata has about 29K rows
# there are 15 ptids in stage 2 mapped data that are not in stage 1 mapped data
# and there are an additional 104 ptids in stage 2 mapped data that are not in stage 1 adata
# the exact reasons are impossible to track down
ptids2mapped=dat_stage2_mapped$Ptid
ptids1adata=dat_stage1_adata$Ptid
ptids.2minus1 =    sort(ptids2mapped[!ptids2mapped %in% ptids1adata])
summary(ptids.2minus1)


# read stage1 mapped data
dat_stage1_mapped=read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdata_correlates_originaldatafromModerna_v1.0_Oct28_2021.csv")
ptids1mapped=dat_stage1_mapped$Ptid
ptids.2minus1raw = sort(ptids2mapped[!ptids2mapped %in% ptids1mapped])
ptids.stage1.rawonly=setdiff(ptids.2minus1, ptids.2minus1raw)
summary(ptids.2minus1raw)
summary(ptids.stage1.rawonly)


# merge stage2 mapped data and stage1 adata to create stage 2 adata
# keep all rows
# for the duplicated columns, take from stage2 mapped data 
dat_stage2 = merge(dat_stage2_mapped, dat_stage1_adata, by="Ptid", all=T, suffixes=c("",".y"))
dat_stage2 = dat_stage2[,!endsWith(names(dat_stage2),".y")]

# taking from stage 2 for the duplicate columns has the consequence that the subjects in stage 1 but not stage 2 would have NA for the duplicated columns
# so we need to redefine some variables, e.g. 
dat_stage2$Senior = ifelse(dat_stage2$Age>=65, 1, 0)

# demo.stratum
tmp = with(dat_stage2, ifelse (URMforsubcohortsampling==1, ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3)), 3+ifelse(Senior, 1, ifelse(HighRiskInd == 1, 2, 3))))
# demo.stratum.new has no NAs
with(dat_stage2, table(tmp, demo.stratum, useNA = "ifany"))
dat_stage2$demo.stratum = tmp

# there are 253 sampled ptids for marker assays
sum(subset(dat_stage2, T, Stage2SamplingInd), na.rm=T)



###############################################################################
# define ph1.BD29

# remove cases that are not EventIndPrimaryOmicronBD29 but are EventIndOmicronBD29
#   to focus on primary cases and not include the cases based on the more relaxed criterion
# if we want to study EventIndOmicronBD29, then comment this line out, uncomment the next line, and make a new dataset
# missingness: there are 7 ptids that are 0 in EventIndPrimaryOmicronBD29 but NA in EventIndOmicronBD29
# these 7 are removed in this step, if they are not, they will also be removed by EventTimeOmicronBD29 >= 7 since they have EventTimeOmicronBD29=0
dat_stage2$ph1.BD29 = with(dat_stage2, !(EventIndOmicronBD29==1 & EventIndPrimaryOmicronBD29==0))
# dat_stage2$ph1.BD29 = T

with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29, useNA="ifany"))
with(subset(dat_stage2, ph1.BD29 & naive==1), table(Trt, EventIndOmicronBD29))
with(dat_stage2, table(is.na(EventIndOmicronBD29), is.na(EventTimeOmicronBD29)))
sum(subset(dat_stage2, ph1.BD29, Stage2SamplingInd), na.rm=T)


# remove ptids that is NA in the three bucket variables: Trt, naive and time period
# Wstratum is made up of the demo variables, CalendarBD1Interval, naive, trt
# controls with missing Wstratum won't be sampled, hence not part of ph1
# cases having missing Wstratum due to  may still be sampled, hence part of ph1
#     if cases have missing demo variables, we want to impute them so that we can assign weights, otherwise it gets too complicated
#     imputation is done for controls without missing demo and all cases. controls are included to improve imputation performance

# # there are two subjects, "US3252458" "US3702017", who have NA as URMforsubcohortsampling in stage 1 data, but not in stage 2 data, from Moderna
# # the two ptids are white with unknown hispanic ethnicity. since by convention these are assigned MinorityInd 0, we assign their URMforsubcohortsampling 0 as well
# dat_stage1_adata[dat_stage1_adata$Ptid %in% c("US3252458", "US3702017"),c("URMforsubcohortsampling","demo.stratum")]
# dat_stage2[dat_stage2$Ptid %in% c("US3252458", "US3702017"),c("URMforsubcohortsampling","demo.stratum")]

dat_stage2$ph1.BD29 = with(dat_stage2,
                             !is.na(naive) & !is.na(Trt) &
                             !is.na(Perprotocol) & !is.na(BDPerprotocol) &
                             !is.na(NumberdaysBD1toBD29) &
                             !is.na(CalendarBD1Interval) &
                             !is.na(EventTimeOmicronBD29)
                            )

# lose 1 ph2 due to na in EventTimeOmicronBD29
# EventTimeOmicronBD29 missingness and NumberdaysBD1toBD29 missingness are concordant 
#   and they are the same participants who missed the BD29 visit
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29, useNA="ifany"))
with(subset(dat_stage2, ph1.BD29 & naive==1), table(Trt, EventIndOmicronBD29))
sum(subset(dat_stage2, ph1.BD29, Stage2SamplingInd), na.rm=T)


# filter by BDPerprotocolIncludeSeroPos, which is based on BDPerprotocol but keep sero+ at BD baseline
# lose 2 cases and 10 controls from this step in the nnaive population
# Perprotocol is part of BDPerprotocol def
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & BDPerprotocolIncludeSeroPos)

with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))
sum(subset(dat_stage2, ph1.BD29, Stage2SamplingInd), na.rm=T)


# require not censored and no evidence of infection from BD1 to BD7
# no loss of ph1 samples from this in the nnaive population
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & EventTimeOmicronBD29 >= 7)

with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))
sum(subset(dat_stage2, ph1.BD29, Stage2SamplingInd), na.rm=T)


# interval bt BD1 and BD29 has to be [19,45] days. 45 is chosen because 45 and 49 lead to the same number of samples
# lose 0 cases and 6 controls from this step in the nnaive population
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & NumberdaysBD1toBD29 >= 19 & NumberdaysBD1toBD29 <= 45)

with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))
sum(subset(dat_stage2, ph1.BD29, Stage2SamplingInd), na.rm=T)


# filter out HIV positive ptids. this includes a highly influential ptid US3632155
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & !HIVinfection)

with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))
sum(subset(dat_stage2, ph1.BD29, Stage2SamplingInd), na.rm=T)
nrow(subset(dat_stage2, ph1.BD29))


# controls should not be NA in the demo vars for stratification, it does not matter for cases since we will impute
# lose none
demo.var=c("HighRiskInd", "URMforsubcohortsampling", "Senior")
dat_stage2$ph1.BD29 = dat_stage2$ph1.BD29 & (complete.cases(dat_stage2[demo.var]) | dat_stage2$EventIndOmicronBD29==1)

with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))
sum(subset(dat_stage2, ph1.BD29, Stage2SamplingInd), na.rm=T)
nrow(subset(dat_stage2, ph1.BD29))
  
  # impute demo variables if there are missingness
  imp.markers = demo.var # imp.markers may be a superset of demo.var to improve imputation performance
  dat.tmp.impute <- subset(dat_stage2, ph1.BD29) # only seek to impute within ph1 samples
  imp <- dat.tmp.impute %>%  select(all_of(imp.markers))         
  stopifnot(0==sum(is.na(imp)))
  
  # # the code for imputation if there is missingness
  # if(any(is.na(imp))) {
  #   n.imp <- 1
  #   
  #   # imputation. diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  #   imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  #   dat.tmp.impute[, imp.markers] <- mice::complete(imp, action = 1)
  # 
  #   # missing covariates imputed properly?
  #   assertthat::assert_that(
  #     all(complete.cases(dat.tmp.impute[, imp.markers])),
  #     msg = "missing covariates imputed properly?"
  #   )    
  #   
  #   # merged imputed values in dat.tmp.impute with dat_stage2
  #   dat_stage2[dat_stage2$ph1.BD29, imp.markers] <- dat.tmp.impute[imp.markers][match(dat_stage2[dat_stage2$ph1.BD29, "Ptid"], dat.tmp.impute$Ptid), ]
  #   
  #   # imputed values of missing covariates merged properly for all individuals in ph1?
  #   assertthat::assert_that(
  #     all(complete.cases(dat_stage2[dat_stage2$ph1.BD29, imp.markers])),
  #     msg = "imputed values of missing covariates merged properly for all individuals?"
  #   )
  # }


# don't allow NA in ph1 riskscore, but there may be NA in the whole dataset, probably due to missing data covariates
stopifnot(all(!is.na(dat_stage2[dat_stage2$ph1.BD29,"risk_score"])))

## a transient solution to missing risk score
# if (any(is.na(dat_stage2[dat_stage2$ph1.BD29,"risk_score"]))) {
#   print("impute risk score")
#   dat_stage2[is.na(dat_stage2$risk_score), "risk_score"]=mean(dat_stage2$risk_score, na.rm=T)
#   dat_stage2[is.na(dat_stage2$standardized_risk_score), "standardized_risk_score"]=mean(dat_stage2$standardized_risk_score, na.rm=T)
# }


# need to impute regression covariates?
if(any(is.na(subset(dat_stage2, ph1.BD29, select=c(MinorityInd, HighRiskInd, Senior))))) {
  stop("need to immpute missing regression covariates")
}

# by construction, there should be no NA. Check to make sure
stopifnot(!any(is.na(dat_stage2$ph1.BD29)))



###############################################################################
#### define ph2

# This step can be done before Wstratum is defined because Wstratum missingness does not affect ph1 anymore

# marker data missing pattern
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD1bindSpike), is.na(BD29bindSpike)))
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD1bindSpike_BA.1), is.na(BD29bindSpike_BA.1)))
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD1bindRBD), is.na(BD29bindRBD)))
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD1pseudoneutid50), is.na(BD29pseudoneutid50)))
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD29bindSpike), is.na(BD29pseudoneutid50)))
# 
# with(subset(dat_stage2, ph1.BD29 & xor(is.na(BD1pseudoneutid50), is.na(BD29pseudoneutid50))), cbind(naive, Trt))
# 
# mypairs(subset(dat_stage2, select=c(BD1bindSpike, BD29bindSpike, BD1bindSpike_BA.1, BD29bindSpike_BA.1, 
#                                     BD1bindRBD, BD29bindRBD, BD1pseudoneutid50, BD29pseudoneutid50)))

must_have_assays <- c("pseudoneutid50")

# require ph2 to have either BD1 or BD29 ID50
dat_stage2$ph2.BD29 = with(dat_stage2, ph1.BD29 & (!is.na(BD1pseudoneutid50) | !is.na(BD29pseudoneutid50)) )

with(subset(dat_stage2, ph2.BD29 & naive==0), table(Trt, EventIndOmicronBD29))
with(subset(dat_stage2, ph2.BD29 & naive==1), table(Trt, EventIndOmicronBD29))
nrow(subset(dat_stage2, ph1.BD29 & Stage2SamplingInd))
nrow(subset(dat_stage2, ph2.BD29))


# DD1 is available in a different subset of people from BD29
dat_stage2$ph2.DD1 = dat_stage2$ph1.BD29 & complete.cases(dat_stage2[,"DD1"%.%must_have_assays])      



###############################################################################
#### define Wstratum and compute weights
###############################################################################

# tps.stratum, which does not depend on EventInd, is no longer needed

n.demo = max(dat_stage2$demo.stratum,na.rm=T) 

# sampling_bucket: 0-31
dat_stage2$sampling_bucket = with(dat_stage2, 
                                  strtoi(paste0(
                                    Trt, 
                                    naive,
                                    EventIndOmicronBD29,
                                    dec_to_binary(CalendarBD1Interval-1, 2)
                                  ), base = 2)
)

dat_stage2$Wstratum = with(dat_stage2, demo.stratum + sampling_bucket * n.demo)

# there should not be na in Wstraum
stopifnot (0 == sum(with(dat_stage2, ph1.BD29 & is.na(Wstratum))))


# collapse Wstratum if there are empty ph2 strata

# sampling_bucket_formergingstrata is used in the second collapsing step
dat_stage2$sampling_bucket_formergingstrata = with(dat_stage2, 
                                                   strtoi(paste0(
                                                     Trt, 
                                                     naive,
                                                     EventIndOmicronBD29
                                                   ), base = 2))

dat.ph1.tmp=subset(dat_stage2, ph1.BD29, select=c(Ptid, Trt, naive, sampling_bucket, ph2.BD29, Wstratum, CalendarBD1Interval, sampling_bucket_formergingstrata))
dat.ph1.tmp$ph2 = dat.ph1.tmp$ph2.BD29

with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))
with(subset(dat_stage2, ph1.BD29 & naive==1), table(Trt, EventIndOmicronBD29))

# the algorithm is implemented in kyotil::cove.boost.collapse.strata because it is also needed by reporting3
dat.ph1.tmp2 = cove.boost.collapse.strata (dat.ph1.tmp, n.demo)

# replace dat_stage2 Wstratum with dat.ph1.tmp2$Wstratum
dat_stage2[dat_stage2$ph1.BD29, "Wstratum"] <- 
  dat.ph1.tmp2$Wstratum[match(dat_stage2[dat_stage2$ph1.BD29, "Ptid"], dat.ph1.tmp2$Ptid)]

# sanity checks
# 1. there should be no overlap in vacc and plac
tab=with(subset(dat_stage2, ph1.BD29), table(Wstratum, Trt))
stopifnot(! any(tab[,1]>0 & tab[,2]>0) )
# 2. there should be no overlap in naive and nnaive
tab=with(subset(dat_stage2, ph1.BD29), table(Wstratum, naive))
stopifnot(! any(tab[,1]>0 & tab[,2]>0) )


# there are too few ph1 ptids in the nnaive and original placebo population
# to stabilize weights, collapse over all calendarbd1interval and demo.stratum
# i.e. merge all cases into one Wstratum and all controls into one Wstratum
select=with(dat_stage2, ph1.BD29 & naive==0 & Trt==0 & EventIndOmicronBD29==0)
dat_stage2[select, "Wstratum"] = min(dat_stage2[select, "Wstratum"])
select=with(dat_stage2, ph1.BD29 & naive==0 & Trt==0 & EventIndOmicronBD29==1)
dat_stage2[select, "Wstratum"] = min(dat_stage2[select, "Wstratum"])

with(dat_stage2, table(Wstratum[EventIndPrimaryOmicronBD29==1 & naive==0 & ph1.BD29]))
with(dat_stage2, table(Wstratum[EventIndPrimaryOmicronBD29==1 & naive==0 & ph2.BD29]))
with(dat_stage2, table(Wstratum[EventIndPrimaryOmicronBD29==0 & naive==0 & ph1.BD29]))
with(dat_stage2, table(Wstratum[EventIndPrimaryOmicronBD29==0 & naive==0 & ph2.BD29]))
with(subset(dat_stage2, EventIndPrimaryOmicronBD29==0 & naive==0 & ph1.BD29), table(CalendarBD1Interval, Wstratum, Trt))
with(subset(dat_stage2, EventIndPrimaryOmicronBD29==0 & naive==1 & ph1.BD29), table(CalendarBD1Interval, Wstratum, Trt))


# compute inverse probability sampling weights wt.BD29
wts_table <- with(subset(dat_stage2, ph1.BD29), table(Wstratum, ph2.BD29))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_stage2$wt.BD29 = ifelse(dat_stage2$ph1.BD29, wts_norm[dat_stage2$Wstratum %.% ""], NA)
assertthat::assert_that(
  all(!is.na(subset(dat_stage2, ph1.BD29 & !is.na(Wstratum))[["wt.BD29"]])),
  msg = "missing wt.BD29")


# compute inverse probability sampling weights wt.DD1
wts_table <- with(subset(dat_stage2, ph1.BD29 & EventIndOmicronBD29), table(Wstratum, ph2.DD1))
if (any(wts_table[,2]==0)) {
  # there are empty ph2 cells when computing wt.DD1, need to compute a second collapsed Wstratum
  dat.ph1.tmp=subset(dat_stage2, ph1.BD29 & EventIndOmicronBD29, 
    select=c(Ptid, Trt, naive, sampling_bucket, ph2.DD1, Wstratum, CalendarBD1Interval, sampling_bucket_formergingstrata))
  dat.ph1.tmp$ph2 = dat.ph1.tmp$ph2.DD1
  
  # collapse Wstratum 
  dat.ph1.tmp2 = cove.boost.collapse.strata (dat.ph1.tmp, n.demo)
  # with(dat.ph1.tmp2, table(Wstratum, ph2.DD1))
  
  # replace dat_stage2 Wstratum with dat.ph1.tmp2$Wstratum
  dat_stage2$WstratumDD1 = NA
  dat_stage2[dat_stage2$ph1.BD29 & dat_stage2$EventIndOmicronBD29, "WstratumDD1"] <- 
    dat.ph1.tmp2$Wstratum[match(dat_stage2[dat_stage2$ph1.BD29 & dat_stage2$EventIndOmicronBD29, "Ptid"], dat.ph1.tmp2$Ptid)]
  
} else {
  dat_stage2$WstratumDD1 = dat_stage2$Wstratum
  
}
wts_table <- with(subset(dat_stage2, ph1.BD29 & EventIndOmicronBD29), table(WstratumDD1, ph2.DD1))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_stage2$wt.DD1 = ifelse(dat_stage2$ph1.BD29 & dat_stage2$EventIndOmicronBD29, wts_norm[dat_stage2$WstratumDD1 %.% ""], NA)
assertthat::assert_that(
  all(!is.na(subset(dat_stage2, ph1.BD29 & EventIndOmicronBD29 & !is.na(WstratumDD1))[["wt.DD1"]])),
  msg = "missing wt.DD1")



###############################################################################
#### impute assay values 
###############################################################################

n.imp <- 1
dat.tmp.impute <- subset(dat_stage2, ph2.BD29)

assay_metadata = read.csv(config$assay_metadata)

imp.markers=c(outer(c("BD1", "BD29"), assay_metadata$assay, "%.%"))

# impute for naive and nnaive separately
for (naive.status in 0:1) {    
  #summary(subset(dat.tmp.impute, Trt == 1 & Bserostatus==0)[imp.markers])      
  imp <- dat.tmp.impute %>% dplyr::filter(naive == naive.status) %>% select(all_of(imp.markers))         
  if(any(is.na(imp))) {
    # if there is no variability, fill in NA with constant values
    for (a in names(imp)) {
      if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
    }            
    # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
    imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
    dat.tmp.impute[dat.tmp.impute$naive == naive.status , imp.markers] <- mice::complete(imp, action = 1)
  }                
}

# missing markers imputed properly?
assertthat::assert_that(
  all(complete.cases(dat.tmp.impute[, imp.markers])),
  msg = "missing markers imputed properly?"
)    

# populate dat_stage2 imp.markers with the imputed values
dat_stage2[dat_stage2$ph2.BD29==1, imp.markers] <- 
  dat.tmp.impute[imp.markers][match(dat_stage2[dat_stage2$ph2.BD29==1, "Ptid"], dat.tmp.impute$Ptid), ]

# imputed values of missing markers merged properly for all individuals in the two phase sample?
assertthat::assert_that(
  all(complete.cases(dat_stage2[dat_stage2$ph2.BD29 == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)


###############################################################################
# using the same conversion factor for stage 1 analyses to make variables comparable to stage 1 analyses
# units will still be labeled as AU/ml
###############################################################################

# stage 1 conversion factors
# convf=c(bindSpike=0.0090, bindRBD=0.0272, bindN=0.0024, pseudoneutid50=0.242, pseudoneutid80=1.502)   

# ID50 in the mapped data have not been converted

for (a in assay_metadata$assay) {
  for (t in c("BD1","BD29","DD1") ) {
    convf = ifelse(
      # ID50 assays:
      startsWith(a,"pseudoneutid50"), 0.242/1.04, 
      # other assays:
      1) 
    dat_stage2[[t %.% a]] <- dat_stage2[[t %.% a]] + log10(convf)
  }
}


###############################################################################
# define delta for dat_stage2
###############################################################################

# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta

# assay limits in metadata has already been converted

tmp=list()
for (a in assay_metadata$assay) {
  uloq=assay_metadata$uloq[assay_metadata$assay==a]
  for (t in c("BD1", "BD29") ) {
    tmp[[t %.% a]] <- ifelse(dat_stage2[[t %.% a]] > log10(uloq), log10(uloq), dat_stage2[[t %.% a]])
  }
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

for (tp in 29) {
  dat_stage2["DeltaBD"%.%tp%.%"overBD1" %.% assay_metadata$assay] <- tmp["BD"%.%tp %.% assay_metadata$assay] - tmp["BD1" %.% assay_metadata$assay]
}   



###############################################################################
# digest check
###############################################################################

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
  tmp = switch(attr(config, "config"),
               moderna_boost = "bb2f964025ec347ec3f7970682542484",
               NA)    
  if (!is.na(tmp)) assertthat::assert_that(digest(dat_stage2[order(names(dat_stage2))])==tmp, msg = "failed make_dat_stage2 digest check. new digest "%.%digest(dat_stage2[order(names(dat_stage2))]))    
}


# save
data_name = paste0(attr(config, "config"), "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")
write.csv(dat_stage2, file = here("data_clean", data_name), row.names=F)
