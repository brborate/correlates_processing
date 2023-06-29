Sys.setenv(TRIAL = "moderna_boost")
TRIAL=Sys.getenv("TRIAL")

library(here)
renv::activate(here::here())

library(dplyr)
library(kyotil)
library(mice)

###############################################################################
#### combine stage1 analysis-ready dataset and stage 2 mapped dataset
###############################################################################


config <- config::get(config = TRIAL)
# this line makes config elements variables in the global scope, it makes coding easier but kind of dangerous
#for(opt in names(config)) eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))

# read stage 2 mapped data 
dat_raw = read.csv(config$mapped_data)
if (colnames(dat_raw)[1]=="Subjectid")  colnames(dat_raw)[1] <- "Ptid" else stop("the first column is unexpectedly not Subjectid")

# read stage1 analysis ready dataset 
dat_stage1 = read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_with_riskscore.csv")

# use new risk score, which 99.5% correlated with the old one and is derived for all ptids, including baseline pos
dat_risk_score = read.csv("/trials/covpn/p3001/analysis/correlates/Part_C_Unblinded_Phase_Data/adata/inputFile_with_riskscore.csv")
dat_stage1$risk_score              = dat_risk_score$risk_score             [match(dat_stage1$Ptid, dat_risk_score$Ptid)]
dat_stage1$Riskscorecohortflag     = dat_risk_score$Riskscorecohortflag    [match(dat_stage1$Ptid, dat_risk_score$Ptid)]
dat_stage1$standardized_risk_score = dat_risk_score$standardized_risk_score[match(dat_stage1$Ptid, dat_risk_score$Ptid)]


#setdiff(names(dat_raw), names(dat_stage1))

# merge two files
# dat_raw has about 14K rows while dat_stage1 has about 29K rows
# dat_stage2 keeps all rows from dat_stage1
dat_stage2 = merge(dat_stage1, dat_raw, by="Ptid", all=T, suffixes=c("",".y"))

# remove columns ending in .y
dat_stage2 = dat_stage2[,!endsWith(names(dat_stage2),".y")]

# hack, need to be replaced with true definition
dat_stage2$naive = 1-dat_stage2$nnaive



###############################################################################
# define ph1 

# not NA in the three bucket variables: Trt, naive and time period
dat_stage2$ph1.BD29 = with(dat_stage2, !is.na(naive) & !is.na(CalendarBD1Interval) & !is.na(Trt) & !is.na(Perprotocol) & !is.na(EventTimeOmicronBD29))

# Perprotocol and BDPerprotocol are actually all 1 already in the mapped data stage 2 dataset
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & Perprotocol)

# not censored and no evidence of infection from BD1 to BD7, implemented with EventTimeOmicronBD29
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & EventTimeOmicronBD29 >= 7)


# interval bt BD1 and BD29 has to be [19,49] days
dat_stage2$ph1.BD29 = with(dat_stage2, ph1.BD29 & NumberdaysBD1toBD29 >= 19 & NumberdaysBD1toBD29 <= 49)


# controls should not be NA in the demo vars for stratification, it does not matter for cases since we will impute
demo.var=c("HighRiskInd", "URMforsubcohortsampling", "Senior")
dat_stage2$ph1.BD29 = dat_stage2$ph1.BD29 & (complete.cases(dat_stage2[demo.var]) | dat_stage2$EventIndOmicronBD29==1)

#stopifnot(all(!is.na(dat_stage2[,"risk_score"]))) # there are some NA in risk score in the whole dataset, probably due to missing data covariates
stopifnot(all(!is.na(dat_stage2[dat_stage2$ph1.BD29,"risk_score"])))

# with(subset(dat_stage2, is.na(risk_score)), table(Trt, Bserostatus))
# with(subset(dat_stage2), table(is.na(risk_score), Bserostatus, Trt))


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

# these two subjects ware white with unknown hispanic ethnicity
# since by convention these are assigned MinorityInd 0, we assign their URMforsubcohortsampling 0 as well
dat_stage2[with(dat_stage2, ph1.BD29 & is.na(URMforsubcohortsampling) & EventIndPrimaryOmicronBD29), "URMforsubcohortsampling"]=0


# impute demo variables if there are missingness
imp.markers = demo.var # imp.markers may be a superset of demo.var to improve imputation performance
dat.tmp.impute <- subset(dat_stage2, ph1.BD29) # only seek to impute within ph1 samples
imp <- dat.tmp.impute %>%  select(all_of(imp.markers))         

# after assigning URMforsubcohortsampling, there should be no missingness
stopifnot(0==sum(is.na(imp)))

# the code for imputation if there is missingness
if(any(is.na(imp))) {
  n.imp <- 1
  
  # imputation. diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  dat.tmp.impute[, imp.markers] <- mice::complete(imp, action = 1)

  # missing covariates imputed properly?
  assertthat::assert_that(
    all(complete.cases(dat.tmp.impute[, imp.markers])),
    msg = "missing covariates imputed properly?"
  )    
  
  # merged imputed values in dat.tmp.impute with dat_stage2
  dat_stage2[dat_stage2$ph1.BD29, imp.markers] <- dat.tmp.impute[imp.markers][match(dat_stage2[dat_stage2$ph1.BD29, "Ptid"], dat.tmp.impute$Ptid), ]
  
  # imputed values of missing covariates merged properly for all individuals in ph1?
  assertthat::assert_that(
    all(complete.cases(dat_stage2[dat_stage2$ph1.BD29, imp.markers])),
    msg = "imputed values of missing covariates merged properly for all individuals?"
  )
}

# by construction, there should be no NA. Check to make sure
stopifnot(!any(is.na(dat_stage2$ph1.BD29)))



###############################################################################
# impute regression covariates

# no missing data
if(any(is.na(subset(dat_stage2, ph1.BD29, select=c(MinorityInd, HighRiskInd, risk_score))))) {
  stop("need to immpute missing regression covariates")
}



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


# this was to amend ph1 definition, but it is not necessary anymore
# a check
# stopifnot (0 == sum(with(dat_stage2, ph1.BD29 & is.na(Wstratum))))
# It removes the controls with missing Wstratum since all cases with missing Wstratum are imputed
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
