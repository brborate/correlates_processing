library(kyotil)
library(dplyr)
library(mice); set.seed(2022)

# this is the one used to make data P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.csv for the lvmn manuscript
#dat  =read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.0_Oct28_2021.csv") # analysis ready dataset containing bAb and PsV markers
# this is the regular data flow
dat  =read.csv("/home/yfong/correlates_processing/data_clean/moderna_real_data_processed_with_riskscore.csv") # analysis ready dataset containing bAb and PsV markers
dat.1=read.csv("/trials/covpn/p3001/download_data/Moderna COVE mRNA 1273P301_immune/P3001ModernaCOVEimmunemarkerdata_correlates_LVMN_2022-01-13.csv") # LVMN marker data

conv.f.lvmn50=0.276

for (t in c("B","Day29","Day57")) {
    x=dat.1[[t%.%"liveneutmn50"]]
    #x=ifelse (x=="<LOD", log10(82.1/2*conv.f.lvmn50), log10(as.numeric(x)*conv.f.lvmn50))
    # the following two lines do the same thing as the previous line but will not generate warnings about NA
    x[x=="<LOD"] =  82.1/2
    x=log10(as.numeric(x)*conv.f.lvmn50)
    dat[[t%.%"liveneutmn50"]]=x[match(dat$Ptid, dat.1$Ptid)]
}

with(dat[dat$Trt==1 & dat$ph2.D57,], table(!is.na(Bliveneutmn50), Wstratum))
#       Wstratum
#         13  14  15  16  17  18  27
#  FALSE   0   0   0   0   0   0   4
#  TRUE  171 143 188 168 146 189  32
with(dat[dat$Trt==1 & dat$ph2.D57,], table(!is.na(Day29liveneutmn50), Wstratum))
#       Wstratum
#         13  14  15  16  17  18  27
#  FALSE   4   5   9   2   6  11   5
#  TRUE  167 138 179 166 140 178  31
with(dat[dat$Trt==1 & dat$ph2.D57,], table(!is.na(Day57liveneutmn50), Wstratum))
#       Wstratum
#         13  14  15  16  17  18  27
#  FALSE   1   0   0   1   1   0   5
#  TRUE  170 143 188 167 145 189  31
  
with(dat[dat$Trt==1 & dat$ph2.D57,], table(!is.na(Bliveneutmn50), Wstratum, SubcohortInd))
with(dat[dat$Trt==1 & dat$ph2.D57,], table(!is.na(Day29liveneutmn50), Wstratum, SubcohortInd))
with(dat[dat$Trt==1 & dat$ph2.D57,], table(!is.na(Day57liveneutmn50), Wstratum, SubcohortInd))

# 5 cases missing lv mn50
subset(dat[dat$Trt==1 & dat$ph2.D57,], (is.na(Bliveneutmn50) | is.na(Day29liveneutmn50) | is.na(Day57liveneutmn50)) & Wstratum==27, select=c(Ptid, Bliveneutmn50, Day29liveneutmn50, Day57liveneutmn50))

# some also miss PsV ID50
dat.0=read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdata_correlates_originaldatafromModerna_v1.0_Oct28_2021.csv")
subset(dat.0, Ptid %in% c("US3002132", "US3042072", "US3382235", "US3482227", "US3732399", "US3742357"), select=c(Ptid, Bpseudoneutid50, Day29pseudoneutid50, Day57pseudoneutid50))



###############################################################################
# impute missing lvmt biomarkers in ph2
#     impute vaccine and placebo (baseline neg only, no baseline pos)
#     use bAb and PsV
#     use baseline, D29 and D57, but not Delta
###############################################################################

dat_proc=dat
n.imp <- 1

assays=c("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80", "liveneutmn50")
imp.markers=c(outer(c("B", "Day29", "Day57"), assays, "%.%"))
dat.tmp.impute <- subset(dat_proc, TwophasesampIndD57 == 1)
summary(subset(dat.tmp.impute, Trt == 1 & Bserostatus==0)[imp.markers])

for (trt in unique(dat_proc$Trt)) {      
  imp <- dat.tmp.impute %>%
    dplyr::filter(Trt == trt) %>%
    select(all_of(imp.markers)) 
    
  # deal with constant variables
  for (a in names(imp)) {
    if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
  }
    
  # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)
    
  dat.tmp.impute[dat.tmp.impute$Trt == trt, imp.markers] <- mice::complete(imp, action = 1)    
}

# missing markers imputed properly?
assertthat::assert_that(
    all(complete.cases(dat.tmp.impute[, imp.markers])),
    msg = "missing markers imputed properly?"
)    

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc$TwophasesampIndD57==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc$TwophasesampIndD57==1, "Ptid"], dat.tmp.impute$Ptid), ]

# imputed values of missing markers merged properly for all individuals in the two phase sample?
assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc$TwophasesampIndD57 == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)


###############################################################################
# impute again for TwophasesampIndD29
#     use baseline and D29

n.imp <- 1

assays=c("bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80", "liveneutmn50")
imp.markers=c(outer(c("B", "Day29"), assays, "%.%"))
dat.tmp.impute <- subset(dat_proc, TwophasesampIndD29 == 1)
summary(subset(dat.tmp.impute, Trt == 1)[imp.markers])

for (trt in unique(dat_proc$Trt)) {    
  imp <- dat.tmp.impute %>%
    dplyr::filter(Trt == trt) %>%
    select(all_of(imp.markers)) 
  
  # deal with constant variables  
  for (a in names(imp)) {
    if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
  }
  
  # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)
    
  dat.tmp.impute[dat.tmp.impute$Trt == trt, imp.markers] <- mice::complete(imp, action = 1)    
}

# missing markers imputed properly?
assertthat::assert_that(
    all(complete.cases(dat.tmp.impute[, imp.markers])),
    msg = "missing markers imputed properly for day 29?"
) 

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc$TwophasesampIndD29==1, imp.markers] <- dat.tmp.impute[imp.markers][match(dat_proc[dat_proc$TwophasesampIndD29==1, "Ptid"], dat.tmp.impute$Ptid), ]    

# imputed values of missing markers merged properly for all individuals in the two phase sample?
assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc$TwophasesampIndD29 == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)


###############################################################################
# define delta for dat_proc
###############################################################################

a="liveneutmn50"
lloq=159.79*0.276

# lloq censoring before taking delta
tmp=list()
for (t in c("B", "Day29", "Day57") ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(lloq), log10(lloq / 2), dat_proc[[t %.% a]])
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

dat_proc["Delta57overB"  %.% a] <- tmp["Day57" %.% a] - tmp["B"     %.% a]
dat_proc["Delta29overB"  %.% a] <- tmp["Day29" %.% a] - tmp["B"     %.% a]
dat_proc["Delta57over29" %.% a] <- tmp["Day57" %.% a] - tmp["Day29" %.% a]

 



###############################################################################
# save

#readr::write_csv(dat_proc, file = "/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.csv")
readr::write_csv(dat_proc, file = "/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_with_riskscore.csv")

corplot(dat_proc$Day57liveneutmn50[dat_proc$Trt==1], dat_proc$Day57pseudoneutid50[dat_proc$Trt==1], xlab="LV MN50", ylab="PsV ID50", add.diag=F)

summary(dat_proc$Day57liveneutmn50[dat_proc$Trt==1 & dat_proc$ph2.D57==1])
mean(dat_proc$Day57liveneutmn50[dat_proc$Trt==1 & dat_proc$ph2.D57==1]>1.614, na.rm=T)
summary(dat_proc$Day57pseudoneutid50[dat_proc$Trt==1 & dat_proc$ph2.D57==1])
mean(dat_proc$Day57pseudoneutid50[dat_proc$Trt==1 & dat_proc$ph2.D57==1]>0.084, na.rm=T)
