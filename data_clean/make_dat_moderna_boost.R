#why stage 2 have ptids than stage 1

Sys.setenv(TRIAL = "moderna_boost")
TRIAL=Sys.getenv("TRIAL")

library(here)
renv::activate(here::here())

library(dplyr)
library(kyotil)
library(mice)

config <- config::get(config = TRIAL)


###############################################################################
# combine stage1 analysis-ready dataset and stage 2 mapped dataset

# read stage 2 mapped data 
dat_stage2_mapped = read.csv(config$mapped_data)
if (colnames(dat_stage2_mapped)[1]=="Subjectid")  colnames(dat_stage2_mapped)[1] <- "Ptid" else stop("the first column is unexpectedly not Subjectid")
dat_stage2_mapped$naive = 1-dat_stage2_mapped$nnaive

# read stage1 analysis ready dataset 
dat_stage1_adata = read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_with_riskscore.csv")
# use new risk score, which 99.5% correlated with the old one and is derived for all ptids, including baseline pos
dat_risk_score = read.csv("/trials/covpn/p3001/analysis/correlates/Part_C_Unblinded_Phase_Data/adata/inputFile_with_riskscore.csv")
dat_stage1_adata$risk_score              = dat_risk_score$risk_score             [match(dat_stage1_adata$Ptid, dat_risk_score$Ptid)]
dat_stage1_adata$Riskscorecohortflag     = dat_risk_score$Riskscorecohortflag    [match(dat_stage1_adata$Ptid, dat_risk_score$Ptid)]
dat_stage1_adata$standardized_risk_score = dat_risk_score$standardized_risk_score[match(dat_stage1_adata$Ptid, dat_risk_score$Ptid)]


#
ptids2=dat_stage2_mapped$Ptid
ptids1=dat_stage1_adata$Ptid
dat_stage1_mapped=read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdata_correlates_originaldatafromModerna_v1.0_Oct28_2021.csv")
ptids1.raw=dat_stage1_mapped$Ptid

ptids.2minus1 = sort(ptids2[!ptids2 %in% ptids1])
ptids.2minus1raw = sort(ptids2[!ptids2 %in% ptids1.raw])
ptids.stage1.rawonly=setdiff(ptids.2minus1, ptids.2minus1raw)

subset(dat_stage2_mapped, Ptid %in% ptids.2minus1, c(Ptid, Bserostatus, EventTimePrimaryD1, EventTimePrimaryD29, EventTimePrimaryD57))

subset(dat_stage1_mapped, Ptid %in% ptids.stage1.rawonly, c(Ptid, Bserostatus, EventTimePrimaryD1, EventTimePrimaryD29, EventTimePrimaryD57))

subset(dat_stage1_mapped, Ptid=="US3002290", Perprotocol)
subset(dat_stage2_mapped, Ptid=="US3002290", Perprotocol)

subset(dat_stage1_adata, Ptid=="US3002290")
subset(dat_stage1_mapped, Ptid=="US3002290")




# merge stage2 mapped data and stage1 analysis-ready data, keep all rows
# dat_stage2_mapped has about 14K rows while dat_stage1_adata has about 29K rows. \
# dat_stage2 keeps all rows from dat_stage1_adata
dat_stage2 = merge(dat_stage2_mapped, dat_stage1_adata, by="Ptid", all=T, suffixes=c("",".y"))
#setdiff(names(dat_stage2_mapped), names(dat_stage1_adata))
# for the duplicated columns, take from stage2 mapped data because stage 1 analysis ready data misses a few ptids, e.g. US3632155, for reasons that are impossible to track down
# this has the consequence that the subjects in stage 1 but not stage 2 would have NA for the duplicated columns
dat_stage2 = dat_stage2[,!endsWith(names(dat_stage2),".y")]
# subset(dat_stage2_mapped, Ptid=="US3632155")
# subset(dat_stage1_adata, Ptid=="US3632155")
# redefine Senior because some ptids are missing in stage 1 and Senior was not defined in stage 2 mappe data
dat_stage2$Senior = ifelse(dat_stage2$Age>=65, 1, 0)



###############################################################################
# define ph1.BD29

dat_stage2$ph1.BD29=T
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29, useNA="ifany"))
with(subset(dat_stage2, ph1.BD29 & naive==1), table(Trt, EventIndOmicronBD29))
# 
with(dat_stage2, table(is.na(EventIndOmicronBD29), is.na(EventTimeOmicronBD29)))

# not NA in the three bucket variables: Trt, naive and time period
# EventTimeOmicronBD29 missingness and NumberdaysBD1toBD29 missingness are concordant and they are participants who missed the BD29 visit.
dat_stage2$ph1.BD29 = with(dat_stage2,
                             !is.na(naive) & !is.na(Trt) &
                             !is.na(Perprotocol) & !is.na(BDPerprotocol) &
                             !is.na(NumberdaysBD1toBD29) &
                             !is.na(CalendarBD1Interval) &
                             !is.na(EventTimeOmicronBD29)
                            )
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29, useNA="ifany"))
with(subset(dat_stage2, ph1.BD29 & naive==1), table(Trt, EventIndOmicronBD29))

# EventIndPrimaryOmicronBD29 usea a more stringent case def. Our main focus is EventIndOmicronBD29
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndPrimaryOmicronBD29))
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))

# Perprotocol
# no loss of ph1 samples from this in the nnaive population
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & Perprotocol)
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))

# BDPerprotocolIncludeSeroPos, which is based on BDPerprotocol but keep sero+ at BD baseline
# lose 2 cases and 10 controls from this step in the nnaive population
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & BDPerprotocolIncludeSeroPos)
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))

# not censored and no evidence of infection from BD1 to BD7
# no loss of ph1 samples from this in the nnaive population
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & EventTimeOmicronBD29 >= 7)
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))

# interval bt BD1 and BD29 has to be [19,49] days
# lose 0 cases and 6 controls from this step in the nnaive population
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & NumberdaysBD1toBD29 >= 19 & NumberdaysBD1toBD29 <= 49)
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))

# controls should not be NA in the demo vars for stratification, it does not matter for cases since we will impute
# lose 0 cases and 1 control from this step in the nnaive population
demo.var=c("HighRiskInd", "URMforsubcohortsampling", "Senior")
dat_stage2$ph1.BD29 = dat_stage2$ph1.BD29 & (complete.cases(dat_stage2[demo.var]) | dat_stage2$EventIndOmicronBD29==1)
with(subset(dat_stage2, ph1.BD29 & naive==0), table(Trt, EventIndOmicronBD29))


# Wstratum is made up of the demo variables, CalendarBD1Interval, naive, trt
# controls with missing Wstratum won't be sampled, hence not part of ph1
# cases having missing Wstratum due to  may still be sampled, hence part of ph1
#     if cases have missing demo variables, we want to impute them so that we can assign weights, otherwise it gets too complicated
#     imputation is done for controls without missing demo and all cases. controls are included to improve imputation performance

# there are two subjects with missing URMforsubcohortsampling
# subset(dat_stage2, ph1.BD29 & is.na(URMforsubcohortsampling))
# with(subset(dat_stage2, ph1.BD29 & !EventIndPrimaryOmicronBD29), table(Trt, naive))
# with(subset(dat_stage2, ph1.BD29 & EventIndPrimaryOmicronBD29), table(Trt, naive))
# with(subset(dat_stage2, ph1.BD29 & is.na(URMforsubcohortsampling) & EventIndPrimaryOmicronBD29), table(Trt, naive))
# with(subset(dat_stage2, ph1.BD29 & EventIndPrimaryOmicronBD29 & !is.na(BD29bindSpike)), table(Trt, naive))
# with(subset(dat_stage2, ph1.BD29 & is.na(URMforsubcohortsampling) & EventIndPrimaryOmicronBD29), cbind(Trt, naive, MinorityInd))

# there are two subjects, "US3252458" "US3702017", who have NA as URMforsubcohortsampling in stage 1 data, but not in stage 2 data, from Moderna
# the two ptids are white with unknown hispanic ethnicity. since by convention these are assigned MinorityInd 0, we assign their URMforsubcohortsampling 0 as well
dat_stage1_adata[dat_stage1_adata$Ptid %in% c("US3252458", "US3702017"),c("URMforsubcohortsampling","demo.stratum")]
dat_stage2[dat_stage2$Ptid %in% c("US3252458", "US3702017"),c("URMforsubcohortsampling","demo.stratum")]

# since we use stage 2 data for duplicate columns, we just need to fix the demo.stratum in the merged dataset
dat_stage2[dat_stage2$Ptid %in% c("US3252458", "US3702017"),"demo.stratum"]=3+ifelse(tmp$Senior, 1, ifelse(tmp$HighRiskInd == 1, 2, 3))

# impute demo variables if there are missingness
imp.markers = demo.var # imp.markers may be a superset of demo.var to improve imputation performance
dat.tmp.impute <- subset(dat_stage2, ph1.BD29) # only seek to impute within ph1 samples
imp <- dat.tmp.impute %>%  select(all_of(imp.markers))         

# after assigning URMforsubcohortsampling, there should be no missingness
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

#stopifnot(all(!is.na(dat_stage2[,"risk_score"]))) # there are some NA in risk score in the whole dataset, probably due to missing data covariates
stopifnot(all(!is.na(dat_stage2[dat_stage2$ph1.BD29,"risk_score"])))

# need to impute regression covariates?
if(any(is.na(subset(dat_stage2, ph1.BD29, select=c(MinorityInd, HighRiskInd, risk_score))))) {
  stop("need to immpute missing regression covariates")
}

# by construction, there should be no NA. Check to make sure
stopifnot(!any(is.na(dat_stage2$ph1.BD29)))



###############################################################################
#### define ph2
# This can be done before Wstratum is defined because Wstratum missingness does not affect ph1 anymore

must_have_assays <- c("pseudoneutid50")

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

dat_stage2$ph2.BD29 = with(dat_stage2, ph1.BD29 & (!is.na(BD1pseudoneutid50) | !is.na(BD29pseudoneutid50)) )
dat_stage2$ph2.BD29 = with(dat_stage2, ph1.BD29 & (!is.na(BD29pseudoneutid50) & !is.na(BD29bindSpike)) )
with(subset(dat_stage2, ph2.BD29 & naive==0), table(Trt, EventIndOmicronBD29))
with(subset(dat_stage2, ph2.BD29 & naive==1), table(Trt, EventIndOmicronBD29))
nrow(subset(dat_stage2, ph2.BD29))


# DD1 may be available for a different subset of people than BD29
dat_stage2$ph2.DD1 = dat_stage2$ph1.BD29 & complete.cases(dat_stage2[,"DD1"%.%must_have_assays])      



###############################################################################
#### define Wstratum, which is used to compute inverse sampling prob weights
###############################################################################

# demo.stratum is defined in the same way as in the original COVE correlates, 
#     but if there are NAs in demo.stratum in cases, we want to impute them
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


# there should not be na in demo.
# stopifnot (0 == sum(with(dat_stage2, ph1.BD29 & is.na(Wstratum))))
dat_stage2$ph1.BD29 = dat_stage2$ph1.BD29 & !is.na(dat_stage2$Wstratum)


###############################################################################
#### collapse sparse Wstratum and compute weights

# the code for collapsing strata is wrapped in the function kyotil::cove.boost.collapse.strata 
#   because it is needed in bootstrapping code in reporting3 repo

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


# adjust Wstratum
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


# now compute inverse probability sampling weights wt.BD29
wts_table <- with(subset(dat_stage2, ph1.BD29), table(Wstratum, ph2.BD29))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_stage2$wt.BD29 = ifelse(dat_stage2$ph1.BD29, wts_norm[dat_stage2$Wstratum %.% ""], NA)
assertthat::assert_that(
  all(!is.na(subset(dat_stage2, ph1.BD29 & !is.na(Wstratum))[["wt.BD29"]])),
  msg = "missing wt.BD29")


# now compute inverse probability sampling weights wt.DD1
# we assume that there will be no empty cells in the following table. If there are, we may need to re-collapse strata.
wts_table <- with(subset(dat_stage2, ph1.BD29 & EventIndOmicronBD29), table(Wstratum, ph2.DD1))
if (any(wts_table[,2]==0)) {
  # there are empty ph2 cells when computing wt.DD1, need to compute a second collapsed Wstratum
  dat.ph1.tmp=subset(dat_stage2, ph1.BD29 & EventIndOmicronBD29, 
    select=c(Ptid, Trt, naive, sampling_bucket, ph2.DD1, Wstratum, CalendarBD1Interval, sampling_bucket_formergingstrata))
  dat.ph1.tmp$ph2 = dat.ph1.tmp$ph2.DD1
  
  # adjust Wstratum
  dat.ph1.tmp2 = cove.boost.collapse.strata (dat.ph1.tmp, n.demo)
  # with(dat.ph1.tmp2, table(Wstratum, ph2.DD1))
  
  # replace dat_stage2 Wstratum with dat.ph1.tmp2$Wstratum
  dat_stage2$WstratumDD1 = NA
  dat_stage2[dat_stage2$ph1.BD29 & dat_stage2$EventIndOmicronBD29, "WstratumDD1"] <- 
    dat.ph1.tmp2$Wstratum[match(dat_stage2[dat_stage2$ph1.BD29 & dat_stage2$EventIndOmicronBD29, "Ptid"], dat.ph1.tmp2$Ptid)]
  
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

# this depends on missingness pattern
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD1bindSpike), is.na(BD29bindSpike)))
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD1bindSpike_BA.1), is.na(BD29bindSpike_BA.1)))
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD1bindRBD), is.na(BD29bindRBD)))
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD1pseudoneutid50), is.na(BD29pseudoneutid50)))
# with(subset(dat_stage2, ph1.BD29), table(is.na(BD29bindSpike), is.na(BD29pseudoneutid50)))
# 
# mypairs(subset(dat_stage2, select=c(BD1bindSpike, BD29bindSpike, BD1bindSpike_BA.1, BD29bindSpike_BA.1, 
#                                     BD1bindRBD, BD29bindRBD, BD1pseudoneutid50, BD29pseudoneutid50)))


n.imp <- 1
dat.tmp.impute <- subset(dat_stage2, ph2.BD29)

assay_metadata = read.csv(config$assay_metadata)

imp.markers=c(outer(c("BD1", "BD29"), assay_metadata$assay, "%.%"))

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

# populate dat_proc imp.markers with the imputed values
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
               moderna_boost = "82a44fade6369d21d0b13b6a93d1cb9c",
               NA)    
  if (!is.na(tmp)) assertthat::assert_that(digest(dat_stage2[order(names(dat_stage2))])==tmp, msg = "failed make_dat_stage2 digest check. new digest "%.%digest(dat_stage2[order(names(dat_stage2))]))    
}




# save
data_name = paste0(attr(config, "config"), "_data_processed.csv")
write.csv(dat_stage2, file = here("data_clean", data_name), row.names=F)
