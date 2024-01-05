Sys.setenv(TRIAL = "id27hpv")

# no need to load _common.R

library(tidyverse)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)
library(here)
library(kyotil)
library(mdw)

begin=Sys.time()

TRIAL=Sys.getenv("TRIAL")

config <- config::get(config = TRIAL)
for(opt in names(config))eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))

# disabling lower level parallelization in favor of higher level of parallelization
# set parallelization in openBLAS and openMP
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1L)
stopifnot(blas_get_num_procs() == 1L)
omp_set_num_threads(1L)

verbose=Sys.getenv("VERBOSE")=="1"

assay_metadata = read.csv(paste0(dirname(attr(config,"file")),"/",config$assay_metadata))
assays=assay_metadata$assay

# created named lists for assay metadata for easier access, e.g. assay_labels_short["bindSpike"]
assay_labels=assay_metadata$assay_label; names(assay_labels)=assays
assay_labels_short=assay_metadata$assay_label_short; names(assay_labels_short)=assays
uloqs=assay_metadata$uloq; names(uloqs)=assays
lloqs=assay_metadata$lloq; names(lloqs)=assays
llods=assay_metadata$lod; names(llods)=assays

markers='M18'%.%assays


# read mapped data 
dat_proc=read.csv(mapped_data)

# 26 cases
table(dat_proc$EventIndPrimaryAnyHPV)

# 4 cases have negative SusceptibilityTimeM18
tmp=subset(dat_proc, EventIndPrimaryAnyHPV==1, c(EligibilityorinitialsamplingTimeM18, SusceptibilityTimeM18  )); 
# tmp[order(tmp[[3]]),]


###############################################################################
# define strata variables
###############################################################################


dat_proc$Senior = as.integer(dat_proc$Age>=15)
  
Bstratum.labels <- c(
  "Age > 14",
  "Age <= 14"
)

# Bstratum: randomization strata
dat_proc$Bstratum = with(dat_proc, ifelse(Senior, 2, 1))
    
names(Bstratum.labels) <- Bstratum.labels

dat_proc$demo.stratum = dat_proc$Bstratum

# tps stratum, used in tps regression
dat_proc <- dat_proc %>%
  mutate(
    tps.stratum = demo.stratum + (Trt-1) * max(demo.stratum,na.rm=T)
  )

with(dat_proc, table(tps.stratum, Trt, useNA='ifany'))


# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+4. 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# A case will have a Wstratum even if its tps.stratum is NA
# The case is defined using EventIndPrimaryD29

max.tps=max(dat_proc$tps.stratum,na.rm=T)
max.tps=max(dat_proc$tps.stratum,na.rm=T)
dat_proc$Wstratum = dat_proc$tps.stratum

dat_proc$Wstratum[with(dat_proc, EventIndPrimaryAnyHPV==1)]=max.tps+1

with(dat_proc, table(tps.stratum, EventIndPrimaryAnyHPV, useNA='ifany'))




###############################################################################
# observation-level weights for bAb
###############################################################################

# markers are all or none, no imputation is needed, pick an arbitrary marker as must have
must_have_assays <- c("bind_HPV6")
dat_proc[["TwophasesampIndM18"]] = complete.cases(dat_proc[,c("M18"%.%must_have_assays)])      
        

# weights 
# ph1
kp = with(dat_proc, Perprotocol == 1 & EligibilityorinitialsamplingTimeM18>0)
wts_table <- with(dat_proc[kp,], table(Wstratum, get("TwophasesampIndM18")))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc[["wt.M18"]] <- wts_norm[dat_proc$Wstratum %.% ""]
# the step above assigns weights for some subjects outside ph1. the next step makes them NA
dat_proc[["wt.M18"]] = ifelse(kp, dat_proc[["wt.M18"]], NA) 
dat_proc[["ph1.M18"]]=!is.na(dat_proc[["wt.M18"]])
dat_proc[["ph2.M18"]]=dat_proc[["ph1.M18"]] & dat_proc[["TwophasesampIndM18"]]

assertthat::assert_that(
    all(!is.na(subset(dat_proc, kp & !is.na(Wstratum))[["wt.M18"]])),
    msg = "missing wt.D for D analyses ph1 subjects")


# alternatively, use SusceptibilityTimeM18>0 to define ph1
# ph1
kp = with(dat_proc, Perprotocol == 1 & SusceptibilityTimeM18>0)
wts_table <- with(dat_proc[kp,], table(Wstratum, get("TwophasesampIndM18")))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc[["wt.M18.sus"]] <- wts_norm[dat_proc$Wstratum %.% ""]
dat_proc[["wt.M18.sus"]] = ifelse(kp, dat_proc[["wt.M18.sus"]], NA) 
dat_proc[["ph1.M18.sus"]]=!is.na(dat_proc[["wt.M18.sus"]])
dat_proc[["ph2.M18.sus"]]=dat_proc[["ph1.M18.sus"]] & dat_proc[["TwophasesampIndM18"]]

assertthat::assert_that(
  all(!is.na(subset(dat_proc, kp & !is.na(Wstratum))[["wt.M18.sus"]])),
  msg = "missing wt.D for D analyses ph1 subjects")




###############################################################################
# observation-level weights for nAb
###############################################################################

# markers are all or none, no imputation is needed, pick an arbitrary marker as must have
must_have_assays <- c("pseudoneutid50_HPV6")
dat_proc[["TwophasesampIndM18nAb"]] = complete.cases(dat_proc[,c("M18"%.%must_have_assays)])      


# weights 
# ph1
kp = with(dat_proc, Perprotocol == 1 & EligibilityorinitialsamplingTimeM18>0)
wts_table <- with(dat_proc[kp,], table(Wstratum, get("TwophasesampIndM18nAb")))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc[["wt.M18.nAb"]] <- wts_norm[dat_proc$Wstratum %.% ""]
# the step above assigns weights for some subjects outside ph1. the next step makes them NA
dat_proc[["wt.M18.nAb"]] = ifelse(kp, dat_proc[["wt.M18.nAb"]], NA) 
stopifnot(all(dat_proc[["ph1.M18"]]==!is.na(dat_proc[["wt.M18.nAb"]])))
dat_proc[["ph2.M18.nAb"]]=dat_proc[["ph1.M18"]] & dat_proc[["TwophasesampIndM18nAb"]]

assertthat::assert_that(
  all(!is.na(subset(dat_proc, kp & !is.na(Wstratum))[["wt.M18.nAb"]])),
  msg = "missing wt.D for D analyses ph1 subjects")


# alternatively, use SusceptibilityTimeM18>0 to define ph1
# ph1
kp = with(dat_proc, Perprotocol == 1 & SusceptibilityTimeM18>0)
wts_table <- with(dat_proc[kp,], table(Wstratum, get("TwophasesampIndM18nAb")))
wts_norm <- rowSums(wts_table) / wts_table[, 2]
dat_proc[["wt.M18.sus.nAb"]] <- wts_norm[dat_proc$Wstratum %.% ""]
dat_proc[["wt.M18.sus.nAb"]] = ifelse(kp, dat_proc[["wt.M18.sus.nAb"]], NA) 
stopifnot(all(dat_proc[["ph1.M18.sus"]]==!is.na(dat_proc[["wt.M18.sus.nAb"]])))
dat_proc[["ph2.M18.sus.nAb"]]=dat_proc[["ph1.M18.sus"]] & dat_proc[["TwophasesampIndM18nAb"]]

assertthat::assert_that(
  all(!is.na(subset(dat_proc, kp & !is.na(Wstratum))[["wt.M18.sus.nAb"]])),
  msg = "missing wt.D for D analyses ph1 subjects")



###############################################################################
# define mdw scores

bAb = assays[1:5]
t = 'M18'
mdw.weights=tryCatch({
  tree.weight(cor(dat_proc[, t%.%bAb], use='complete.obs'))
}, error = function(err) {
  print(err$message)
  rep(1/length(bAb), length(bAb))
})
print(mdw.weights)  
dat_proc[, t%.%'bind_mdw'] = scale(dat_proc[, t%.%bAb]) %*% mdw.weights
write.csv(mdw.weights, file = here("data_clean", "csv", TRIAL%.%"_mdw_weights_bAb.csv"))


nAb = assays[11:15]
t = 'M18'
mdw.weights=tryCatch({
  tree.weight(cor(dat_proc[, t%.%nAb], use='complete.obs'))
}, error = function(err) {
  print(err$message)
  rep(1/length(nAb), length(nAb))
})
print(mdw.weights)  
dat_proc[, t%.%'pseudoneutid50_mdw'] = scale(dat_proc[, t%.%nAb]) %*% mdw.weights
write.csv(mdw.weights, file = here("data_clean", "csv", TRIAL%.%"_mdw_weights_nAb.csv"))


###############################################################################
# impute covariates if necessary
# not needed b/c:
# age is present for everybody
# EligibilityorinitialsamplingTimeM18 is not imputable




###############################################################################
# digest check
###############################################################################

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         id27hpv = "332781220e74ea3a013dc44e7b85d91e", 
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, msg = "Warning: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%'  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")
write_csv(dat_proc, file = here("data_clean", "csv", data_name))

print("run time: "%.%format(Sys.time()-begin, digits=1))
