library(here)
renv::activate(here::here())


###############################################################################
#### bring stage1 analysis-ready dataset and stage 2 mapped dataset together
###############################################################################

config <- config::get(config = Sys.getenv("TRIAL"))
for(opt in names(config)){
  eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
}

data_name = paste0(attr(config, "config"), "_data_processed_with_riskscore.csv")

# read stage1 analysis ready dataset 
dat_stage1 = read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_with_riskscore.csv")

# read stage 2 mapped data 
dat_raw = read.csv(mapped_data)
if (colnames(dat_raw)[1]=="Subjectid")  colnames(dat_raw)[1] <- "Ptid" else stop("the first column is unexpectedly not Subjectid")

setdiff(names(dat_raw), names(dat_stage1))

# merge two files
dat_stage2 = merge(dat_stage1, dat_raw, by="Ptid", all=T, suffixes=c("",".y"))

# remove columns ending in .y
dat_stage2 = dat_stage2[,!endsWith(names(dat_stage2),".y")]


###############################################################################
#### define sampling stratum variables
###############################################################################

# demo.stratum has not changed from stage 1, but tps.stratum has to be updated

# sampling_bucket: 0-15
dat_stage2$sampling_bucket = with(dat_stage2, 
    strtoi(paste0(
      strtoi(paste0( 
        Trt, 
        Bserostatus # Bserostatus needs to be replaced by naive/non-naive in the real dataset
      ), base = 2), 
      CalendarBD1Interval-1
    ), base=4)
)

# tps stratum, used in tps regression
# [origina, crossover] x [Naive, non-naive] x [CalendarBD1Interval]
# 1 ~ 16*max(demo.stratum), 
dat_stage2 <- dat_stage2 %>% 
  mutate(
    tps.stratum = demo.stratum + sampling_bucket * max(demo.stratum,na.rm=T) 
  )

# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+16. 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the 16 groups defined by Trt, Bserostatus and CalendarBD1Interval
# A case will have a Wstratum even if its tps.stratum is NA

max.tps=max(dat_stage2$tps.stratum,na.rm=T)
dat_stage2$Wstratum = dat_stage2$tps.stratum

dat_stage2[!is.na(dat_stage2$EventIndOmicronBD29) & dat_stage2$EventIndOmicronBD29==1, "Wstratum"] = max.tps + 
  dat_stage2[!is.na(dat_stage2$EventIndOmicronBD29) & dat_stage2$EventIndOmicronBD29==1, "sampling_bucket"] + 1

with(dat_stage2, table(tps.stratum))
with(dat_stage2, table(Wstratum))
with(dat_stage2, table(Wstratum, EventIndOmicronBD29, Trt))

with(dat_stage2, table(sampling_bucket, Trt))

with(dat_stage2, table(Wstratum, sampling_bucket))



###############################################################################
#### inverse probability sampling weights
###############################################################################

#### Note that Wstratum may have NA if any variables to form strata has NA

# define must_have_assays for ph2 definition
must_have_assays <- c("bindSpike", "bindRBD")

# define ph1 as: (naive or non-naive) & (not censored and no evidence of infection from BD1 to BD7)
dat_stage2$ph1.BD29 = !is.na(dat_stage2$Bserostatus) & dat_stage2$EventTimePrimaryBD29 >= 7
# expect no NAs
stopifnot(!any(is.na(dat_stage2$ph1.BD29)))

# TwophasesampInd: be in ph1 &  have the necessary markers

# require baseline
dat_stage2[["TwophasesampIndBD29"]] = dat_stage2$ph1 & 
  complete.cases(dat_stage2[,c("BD1"%.%must_have_assays, "BD29"%.%must_have_assays)])      

dat_stage2[["ph2.BD"%.%tp]]= ifelse(dat_stage2$ph1, dat_stage2[["ph1.BD"%.%tp]] & dat_stage2[["TwophasesampIndBD"%.%tp]], NA)

# weights 
tp=29
tmp = with(dat_stage2, ph1)
wts_table <- with(dat_stage2[tmp,], table(Wstratum, get("TwophasesampIndBD"%.%tp)))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_stage2[["wt.BD"%.%tp]] = ifelse(dat_stage2$ph1, wts_norm[dat_stage2$Wstratum %.% ""], NA)

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

tmp=list()
for (a in assays.includeN) {
  for (t in c("B", paste0("Day", config$timepoints)) ) {
    tmp[[t %.% a]] <- ifelse(dat_stage2[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat_stage2[[t %.% a]])
  }
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

for (tp in rev(timepoints)) {
  dat_stage2["Delta"%.%tp%.%"overB" %.% assays.includeN] <- tmp["Day"%.%tp %.% assays.includeN] - tmp["B" %.% assays.includeN]
}   
if(two_marker_timepoints) {
  dat_stage2["Delta"%.%timepoints[2]%.%"over"%.%timepoints[1] %.% assays.includeN] <- tmp["Day"%.% timepoints[2]%.% assays.includeN] - tmp["Day"%.%timepoints[1] %.% assays.includeN]
}




###############################################################################
# digest check
###############################################################################

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
  tmp = switch(attr(config, "config"),
               moderna_boost = "34e297fd1a736f9320573ff1d2944904",
               NA)    
  if (!is.na(tmp)) assertthat::assert_that(digest(dat_proc[order(names(dat_proc))])==tmp, msg = "failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))]))    
}




# save
write.csv(dat_stage2, file = here("data_clean", data_name), row.names=F)
