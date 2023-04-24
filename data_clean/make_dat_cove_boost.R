library(here)
renv::activate(here::here())

library(dplyr)
library(kyotil)
library(mice)

###############################################################################
#### combine stage1 analysis-ready dataset and stage 2 mapped dataset
###############################################################################

config <- config::get(config = Sys.getenv("TRIAL"))
# this line makes config elements variables in the global scope, it makes coding easier but kind of dangerous
#for(opt in names(config)) eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))

# read stage1 analysis ready dataset 
dat_stage1 = read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_with_riskscore.csv")

# read stage 2 mapped data 
dat_raw = read.csv(config$mapped_data)
if (colnames(dat_raw)[1]=="Subjectid")  colnames(dat_raw)[1] <- "Ptid" else stop("the first column is unexpectedly not Subjectid")

#setdiff(names(dat_raw), names(dat_stage1))

# merge two files
dat_stage2 = merge(dat_stage1, dat_raw, by="Ptid", all=T, suffixes=c("",".y"))

# remove columns ending in .y
dat_stage2 = dat_stage2[,!endsWith(names(dat_stage2),".y")]

# hack, need to be replaced with true definition
dat_stage2$naive = 1-dat_stage2$Bserostatus



###############################################################################
# define ph1 

# not NA in the three bucket variables: Trt, naive and time period
dat_stage2$ph1.BD29 = !is.na(dat_stage2$naive) & !is.na(dat_stage2$CalendarBD1Interval) & !is.na(dat_stage2$Trt)

# not censored and no evidence of infection from BD1 to BD7
# hack alert: is it safe to use EventTimeOmicronBD29 for this?
dat_stage2$ph1.BD29 = dat_stage2$ph1.BD29 & !is.na(dat_stage2$EventTimeOmicronBD29) & dat_stage2$EventTimeOmicronBD29 >= 7

# controls should not be NA in the demo vars for stratification, it does not matter for cases since we will impute
demo.var=c("HighRiskInd", "URMforsubcohortsampling", "Senior")
dat_stage2$ph1.BD29 = dat_stage2$ph1.BD29 & (complete.cases(dat_stage2[demo.var]) | dat_stage2$EventIndOmicronBD29==1)

# Wstratum is made up of the demo variables, CalendarBD1Interval, naive, trt
# controls with missing Wstratum won't be sampled, hence not part of ph1
# cases having missing Wstratum due to  may still be sampled, hence part of ph1
#     if cases have missing demo variables, we want to impute them so that we can assign weights, otherwise it gets too complicated
#     imputation is done for controls without missing demo and all cases. controls are included to improve imputation performance

imp.markers = demo.var # imp.markers may be a superset of demo.var to improve imputation performance
dat.tmp.impute <- subset(dat_stage2, ph1.BD29) # only seek to impute within ph1 samples
imp <- dat.tmp.impute %>%  select(all_of(imp.markers))         

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
#### define ph2
# This can bde done before Wstratum is defined because Wstratum missingness does not affect ph1 anymore

must_have_assays <- c("bindSpike", "bindRBD")

dat_stage2$ph2.BD29 = dat_stage2$ph1.BD29 & complete.cases(dat_stage2[,c("BD1"%.%must_have_assays, "BD29"%.%must_have_assays)])      



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

# this was to amend ph1 definition. 
# It removes the controls with missing Wstratum since all cases with missing Wstratum are imputed
# But it is not necessary anymore
table(dat_stage2$ph1.BD29, !is.na(dat_stage2$Wstratum))
#dat_stage2$ph1.BD29 = dat_stage2$ph1.BD29 & !is.na(dat_stage2$Wstratum)



###############################################################################
#### collapse sparse Wstratum and compute weights

# the code for collapsing strata is wrapped in the function kyotil::cove.boost.collapse.strata 
#   because it is needed in bootstrapping code in reporting3 repo

dat_stage2$sampling_bucket_formergingstrata = with(dat_stage2, 
                                                strtoi(paste0(
                                                  Trt, 
                                                  naive,
                                                  EventIndOmicronBD29
                                                ), base = 2)
)

dat.ph1.tmp=subset(dat_stage2, ph1.BD29, 
          select=c(Ptid, sampling_bucket, ph2.BD29, Wstratum, CalendarBD1Interval, sampling_bucket_formergingstrata))
dat.ph1.tmp$ph2 = dat.ph1.tmp$ph2.BD29

dat.ph1.tmp = cove.boost.collapse.strata (dat.ph1.tmp, n.demo)

# replace dat_stage2 Wstratum with dat.ph1.tmp$Wstratum
dat_stage2[dat_stage2$ph1.BD29, "Wstratum"] <- 
    dat.ph1.tmp$Wstratum[match(dat_stage2[dat_stage2$ph1.BD29, "Ptid"], dat.ph1.tmp$Ptid)]

# sanity check, there should be no overlap in the two columns
tab=with(subset(dat_stage2, ph1.BD29), table(Wstratum, Trt))
stopifnot(! any(tab[,1]>0 & tab[,2]>0) )

# sanity check, there should be no overlap in the two columns
tab=with(subset(dat_stage2, ph1.BD29), table(Wstratum, naive))
stopifnot(! any(tab[,1]>0 & tab[,2]>0) )


# now compute inverse probability sampling weights
tp=29
tmp = with(dat_stage2, ph1.BD29)
wts_table <- with(dat_stage2[tmp,], table(Wstratum, get("ph2.BD"%.%tp)))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_stage2[["wt.BD"%.%tp]] = ifelse(dat_stage2$ph1.BD29, wts_norm[dat_stage2$Wstratum %.% ""], NA)

assertthat::assert_that(
  all(!is.na(subset(dat_stage2, tmp & !is.na(Wstratum))[["wt.BD"%.%tp]])),
  msg = "missing wt.BD for D analyses ph1 subjects")




###############################################################################
#### impute assay values
###############################################################################

# this depends on missingness pattern



###############################################################################
# define delta for dat_stage2
###############################################################################

# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta

assay_metadata = read.csv(config$assay_metadata)

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
#               moderna_boost = "34e297fd1a736f9320573ff1d2944904",
               NA)    
  if (!is.na(tmp)) assertthat::assert_that(digest(dat_stage2[order(names(dat_stage2))])==tmp, msg = "failed make_dat_stage2 digest check. new digest "%.%digest(dat_stage2[order(names(dat_stage2))]))    
}




# save
data_name = paste0(attr(config, "config"), "_data_processed.csv")
write.csv(dat_stage2, file = here("data_clean", data_name), row.names=F)
