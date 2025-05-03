Sys.setenv(TRIAL = "nextgen_mock")
source(here::here("_common.R"))
# no need to run renv::activate(here::here()) b/c .Rprofile exists

{
library(tidyverse)
library(glue)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)
library(here)
library(mdw)

begin=Sys.time()
}


########################################################################################################
# 1. read mapped data with risk score added

{
  # read mapped data
  dat_raw = read.csv(mapped_data)
  dat_proc = preprocess(dat_raw, study_name)   
  names(dat_proc)[[1]]="Ptid"

  assay_metadata=read.csv(config$assay_metadata)
  assays=assay_metadata$assay
  
  dat_proc$Track = as.factor(dat_proc$Track)
  
}


{
# # Youyi's quick fix
# kp = which(dat_proc$ph2.ps.immuno==0 & dat_proc$ph2.immuno==1 & dat_proc$ph2.D31_7==1 & dat_proc$COVIDIndD31_7toM12==1 & dat_proc$Track!="A")#[sample.int(268)[1:240]]
# dat_proc[kp,"ph2.immuno"]=0
# # dat_proc[kp,"ph2.D31_7"]=0
# 
# kp = which(dat_proc$ph2.ps.immuno==0 & dat_proc$ph2.ps.D31_7==1 & dat_proc$COVIDIndD31_7toM12==1)
# dat_proc[kp,"ph2.ps.immuno"]=1
# dat_proc[kp,"ph2.ps.D31_7"]=1
# 
# kp = which(dat_proc$ph2.ps.immuno==1 & dat_proc$ph2.immuno==0) # all cases
# dat_proc[kp, "ph2.ps.immuno"] = 0
# dat_proc[kp, "ph2.ps.D31_7"] = 0


mytable(dat_proc$Track, dat_proc$ph1.ps.D31_7)
mytable(dat_proc$Track, dat_proc$ph1.D31_7)
table.prop(dat_proc$ph2.immuno, dat_proc$Track)
table.prop(dat_proc$ph2.ps.immuno, dat_proc$Track)
mytable(dat_proc$ph2.immuno, dat_proc$ph2.ps.immuno)
mytable(dat_proc$ph2.immuno, dat_proc$ph2.D31_7, dat_proc$COVIDIndD31_7toM12)
mytable(dat_proc$ph2.ps.immuno, dat_proc$ph2.ps.D31_7, dat_proc$COVIDIndD31_7toM12)
mytable(dat_proc$ph2.ps.trackA, dat_proc$ph2.ps.immuno)

mytable(dat_proc$COVIDIndD31_7toM12, dat_proc$ph1.ps.D31_7)
}

###############################################################################
# 3. stratum variables
{
dat_proc$Bstratum = 1 # there are no demographics stratum for subcohort sampling

# demographics stratum are not used in making the mock data, but good to define it
dat_proc$demo.stratum = dat_proc$Senior + 2
dat_proc$demo.stratum[dat_proc$Senior==0 & dat_proc$HighRiskInd==1] = 1

# need one set for tcell markers and one set for antibody markers because the former cannot depend on Track

dat_proc <- dat_proc %>% mutate(tps.stratum = 3*(as.numeric(Track)-1)+demo.stratum)
dat_proc$tps.stratum[dat_proc$Trt==1] = dat_proc$tps.stratum[dat_proc$Trt==1] + 50
mytable(dat_proc$Trt, dat_proc$tps.stratum)

table.prop(dat_proc$ph2.immuno, dat_proc$demo.stratum)
table.prop(dat_proc$ph2.ps.immuno, dat_proc$demo.stratum)

dat_proc$Wstratum = dat_proc$tps.stratum
dat_proc$Wstratum[with(dat_proc, COVIDIndD31_7toM12==1 & Trt==0)]=99
dat_proc$Wstratum[with(dat_proc, COVIDIndD31_7toM12==1 & Trt==1)]=100 


# 
dat_proc <- dat_proc %>% mutate(tps.stratum.tcell = demo.stratum)
dat_proc$tps.stratum.tcell[dat_proc$Trt==1] = dat_proc$tps.stratum.tcell[dat_proc$Trt==1] + 50
mytable(dat_proc$Trt, dat_proc$tps.stratum.tcell)

dat_proc$Wstratum.tcell = dat_proc$tps.stratum.tcell
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD31_7toM12==1 & Trt==0)]=99
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD31_7toM12==1 & Trt==1)]=100 

mytable(dat_proc$ph1.ps.D31_7, dat_proc$Wstratum)
mytable(dat_proc$ph1.ps.D31_7, dat_proc$Wstratum.tcell)
mytable(dat_proc$ph1.ps.D31_7, dat_proc$Track)

mytable(dat_proc$Wstratum.tcell, dat_proc$Wstratum)
}


################################################################################
# 4. Define ph1, ph2, and weights
# Note that Wstratum may have NA if any variables to form strata has NA
{

tp="31_7"
dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp,    ph2="ph2.D"%.%tp,    Wstratum="Wstratum",       wt="wt.D"%.%tp, verbose=T) 
dat_proc = add.wt(dat_proc, ph1="ph1.ps.D"%.%tp, ph2="ph2.ps.D"%.%tp, Wstratum="Wstratum",       wt="wt.ps.D"%.%tp, verbose=T)
dat_proc = add.wt(dat_proc, ph1="ph1.ps.D"%.%tp, ph2="ph2.ps.D"%.%tp, Wstratum="Wstratum.tcell", wt="wt.ps.D"%.%tp%.%".tcell", verbose=T) 

}


###############################################################################
# 5. impute missing biomarkers in ph2 (assay imputation)
#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, each time point, but not Delta

{
tcellvv=assays[35:38]
abmarkers=assays[1:34]

n.imp=1
tp=31

imp.markers=c(outer(c("B", "Day"%.%tp), abmarkers, "%.%"))
dat.tmp.impute <- subset(dat_proc, get("ph2.D31_7") == 1)
imp <- dat.tmp.impute %>% select(all_of(imp.markers))
any(is.na(imp))

imp.markers=c(outer(c("B", "Day"%.%tp), tcellvv, "%.%"))
dat.tmp.impute <- subset(dat_proc, get("ph2.ps.D31_7") == 1)
imp <- dat.tmp.impute %>% select(all_of(imp.markers))
any(is.na(imp))


}




###############################################################################
# 6. transformation of the markers
# create S-stimulated markers by summing up S1 and S2 on the anti log scale


###############################################################################
# 7. add mdw scores for nAb
{
}


###############################################################################
# 8. add fold change markers
# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta
{
assays1=assays # here, mdw scores delta are computed as weighted average of delta, not as delta of mdw

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

}


###############################################################################
# 9. add discrete/trichotomized markers
{
dat_proc$tmp = with(dat_proc, ph2.D31_7)
assays1 = abmarkers
all.markers = c("B"%.%assays1, "Day31"%.%assays1, "Delta31overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D31_7", verbose=F)
dat_proc$tmp = NULL

dat_proc$tmp = with(dat_proc, ph2.ps.D31_7)
assays1 = tcellvv
all.markers = c("B"%.%assays1, "Day31"%.%assays1, "Delta31overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.ps.D31_7.tcell", verbose=F)
dat_proc$tmp = NULL
}

###############################################################################
# 10. impute covariates if necessary



###############################################################################
# special handling 




###############################################################################
# digest check

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         nextgen_mock = "a6331ba037f48bebf626e4c961689f37",
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, 
      msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))



print("run time: "%.%format(Sys.time()-begin, digits=1))
