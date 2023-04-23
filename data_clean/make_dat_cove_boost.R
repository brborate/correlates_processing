library(here)
renv::activate(here::here())

library(dplyr)
library(kyotil)
library(mice)

###############################################################################
#### bring stage1 analysis-ready dataset and stage 2 mapped dataset together
###############################################################################

config <- config::get(config = Sys.getenv("TRIAL"))
# this line makes config elements variables in the global scope, it makes coding easier but kind of dangerous
#for(opt in names(config)) eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))

# read stage1 analysis ready dataset 
dat_stage1 = read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_with_riskscore.csv")

# read stage 2 mapped data 
dat_raw = read.csv(config$mapped_data)
if (colnames(dat_raw)[1]=="Subjectid")  colnames(dat_raw)[1] <- "Ptid" else stop("the first column is unexpectedly not Subjectid")

setdiff(names(dat_raw), names(dat_stage1))

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



###############################################################################
#### define sampling stratum variables
###############################################################################

# demo.stratum is defined in the same way as in the original COVE correlates, 
#     but if there are NAs in demo.stratum in cases, we want to impute them
# tps.stratum has to be updated

# take d digits
dec_to_bin <- function(x,d) sapply(x, function(xx) ifelse(is.na(xx), NA, paste(rev(as.integer(intToBits(xx))[1:d]), collapse="")))

# sampling_bucket: 0-31
dat_stage2$sampling_bucket = with(dat_stage2, 
                                  strtoi(paste0(
                                    Trt, 
                                    naive,
                                    EventIndOmicronBD29,
                                    dec_to_bin(CalendarBD1Interval-1, 2)
                                  ), base = 2)
)

# used for defining tps stratum. case status is not used in this variable
dat_stage2$sampling_bucket_fortpsstratum = with(dat_stage2, 
                                                strtoi(paste0(
                                                  Trt, 
                                                  naive,
                                                  EventIndOmicronBD29,
                                                  dec_to_bin(CalendarBD1Interval-1,2) # two digits
                                                ), base = 2)
)

# used for step 2 of merging sampling strata. calendar period is not used 
dat_stage2$sampling_bucket_forstratamerging = with(dat_stage2, 
                                                strtoi(paste0(
                                                  Trt, 
                                                  naive,
                                                  EventIndOmicronBD29
                                                ), base = 2)
)

n.demo = max(dat_stage2$demo.stratum,na.rm=T) 

# Wstratum is used to compute sampling weights. 
# in case cohort, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+16. 
# Differs from tps stratum in that case is a separate stratum within each of the 16 groups 
#   defined by Trt, Bserostatus and CalendarBD1Interval
# in case control, there are as many case strata as control strata
dat_stage2$Wstratum = with(dat_stage2, demo.stratum + sampling_bucket * n.demo)

# tps stratum, used in tps regression
# [origina, crossover] x [Naive, non-naive] x [CalendarBD1Interval]
# 1 ~ 16*max(demo.stratum), 
dat_stage2$tps.stratum = with(dat_stage2, demo.stratum + sampling_bucket_fortpsstratum * n.demo)


# with(dat_stage2, table(tps.stratum))
# with(dat_stage2, table(Wstratum))
# with(dat_stage2, table(Wstratum, EventIndOmicronBD29, Trt))
# with(dat_stage2, table(sampling_bucket, Trt))
# with(dat_stage2, table(Wstratum, sampling_bucket))


# amend ph1 definition now that we have defined Wstratum
# this only removes the controls with missing Wstratum since all cases with missing Wstratum are imputed
dat_stage2$ph1.BD29 = dat_stage2$ph1.BD29 & !is.na(dat_stage2$Wstratum)

# by construction, there should be no NA. Check to make sure
stopifnot(!any(is.na(dat_stage2$ph1.BD29)))



###############################################################################
#### ph2

must_have_assays <- c("bindSpike", "bindRBD")

dat_stage2$ph2.BD29 = dat_stage2$ph1.BD29 & complete.cases(dat_stage2[,c("BD1"%.%must_have_assays, "BD29"%.%must_have_assays)])      

# collapse strata to deal with sparsity

# first, do it across demo strata within each sampling bucket 
sampling_buckets=0:max(dat_stage2$sampling_bucket, na.rm=T)
for (i in sampling_buckets) {
  
  select = dat_stage2$ph1.BD29 & dat_stage2$sampling_bucket==i 
  
  # make sure there are such ph1 samples
  if (sum(select, na.rm=T)>0) {
  
    # make sure there is at least 1 ph2 sample
    if (sum(dat_stage2[select,"ph2.BD29"], na.rm=T)>=1) {
      tab = with(dat_stage2[select,], table(Wstratum, ph2.BD29))
      # merging
      if (any(tab[,2]==0)) {
        dat_stage2[select,"Wstratum"] = min(dat_stage2[select,"Wstratum"], na.rm=T)
      }
    } else {
      # if there are no ph2 samples, will collapse in the next level
    }
  } # okay if there are no samples in the bucket
}


# second, do it across the 4 calendar periods
# merge a period with the next period if not the last, merge with the last period with the previous if needed

sampling_buckets.2=0:max(dat_stage2$sampling_bucket_forstratamerging, na.rm=T)
# need this in the merging process as it needs to be updated
dat_stage2$tmpCalendarBD1Interval=dat_stage2$CalendarBD1Interval

for (i in sampling_buckets.2) {
  select = dat_stage2$ph1.BD29 & dat_stage2$sampling_bucket_forstratamerging==i 
  
  # make sure there are such ph1 samples
  if (sum(select, na.rm=T)>0) {
    
    # make sure there is at least 1 ph2 sample
    if (sum(dat_stage2[select,"ph2.BD29"], na.rm=T)>=1) {
      tab = with(dat_stage2[select,], table(tmpCalendarBD1Interval, ph2.BD29))
      # merging
      (tab)
      jj = which(tab[,2]==0)
      for (j in jj) {
        select2 = select & dat_stage2$tmpCalendarBD1Interval==j
        # if the last is empty, set all to the first
        dat_stage2[select2,"Wstratum"] = dat_stage2[select2,"Wstratum"] + n.demo * ifelse(j<4, 1, -3)
        # also need to update sampling bucket so that step 3 can work
        dat_stage2[select2,"sampling_bucket"] = dat_stage2[select2,"sampling_bucket"] + ifelse(j<4, 1, -3)
        # need to update this so that select2 will pick these up
        dat_stage2$tmpCalendarBD1Interval[select2] = dat_stage2$tmpCalendarBD1Interval[select2] + 1
      }
    } else {
      stop("something wrong - no ph2 samples in sampling_buckets.2 "%.%i)
    }
  } else {
    # no ph1 samples in this sampling_buckets.2, probably case
  }
}

# third, one more time do it across demo strata within each sampling bucket 
#       because it is possible that collapsing buckets introduced empty demo strata
sampling_buckets=0:max(dat_stage2$sampling_bucket, na.rm=T)
for (i in sampling_buckets) {
  
  select = dat_stage2$ph1.BD29 & dat_stage2$sampling_bucket==i 
  
  # make sure there are such ph1 samples
  if (sum(select, na.rm=T)>0) {
    
    # make sure there is at least 1 ph2 sample
    if (sum(dat_stage2[select,"ph2.BD29"], na.rm=T)>=1) {
      tab = with(dat_stage2[select,], table(Wstratum, ph2.BD29))
      (tab)
      # merging
      if (any(tab[,2]==0)) {
        dat_stage2[select,"Wstratum"] = min(dat_stage2[select,"Wstratum"], na.rm=T)
      }
    } else {
      # if there are no ph2 samples, will collapse in the next level
    }
  } # okay if there are no samples in the bucket
}


# make sure no empty cells after all this
tab=with(subset(dat_stage2, ph1.BD29), table(Wstratum, ph2.BD29))
print(tab)
stopifnot(all(tab[,2]>0))

# sanity check, there should be no overlap in the two columns
tab=with(dat_stage2, table(Wstratum, Trt))
print(tab)
stopifnot(! any(tab[,1]>0 & tab[,2]>0) )

# sanity check, there should be no overlap in the two columns
tab=with(dat_stage2, table(Wstratum, naive))
print(tab)
stopifnot(! any(tab[,1]>0 & tab[,2]>0) )




###############################################################################
#### inverse probability sampling weights

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
