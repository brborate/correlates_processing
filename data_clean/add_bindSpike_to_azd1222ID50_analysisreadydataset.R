dat.bAb =read.csv("/trials/covpn/p3002/analysis/correlates/Part_A_Blinded_Phase_Data/adata/azd1222_bAb_data_processed_with_riskscore.csv") 

for (t in c("B","Day29","Day57")) {
    x=dat.bAb[[t%.%"bindSpike"]]
    dat_proc[[t%.%"bindSpike"]]=x[match(dat_proc$Ptid, dat.bAb$Ptid)]
}


with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D57,], table(!is.na(BbindSpike), Wstratum))
with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D57,], table(!is.na(Day29bindSpike), Wstratum))
with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D57,], table(!is.na(Day57bindSpike), Wstratum))
  
with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D57,], table(!is.na(BbindSpike), Wstratum, SubcohortInd))
with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D57,], table(!is.na(Day29bindSpike), Wstratum, SubcohortInd))
with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D57,], table(!is.na(Day57bindSpike), Wstratum, SubcohortInd))

# 0 cases missing bindSpike
subset(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D57,], (is.na(BbindSpike) | is.na(Day29bindSpike) | is.na(Day57bindSpike)) & Wstratum==27, select=c(Ptid, BbindSpike, Day29bindSpike, Day57bindSpike))




###############################################################################
# impute missing bindSpike biomarkers in ph2
#     impute vaccine and placebo, baseline neg and pos
#     use bAb and PsV
#     use baseline, D29 and D57, but not Delta
###############################################################################

set.seed(2022)

n.imp <- 1

assays=c("bindSpike", "pseudoneutid50")
imp.markers=c(outer(c("B", "Day29", "Day57"), assays, "%.%"))
dat.tmp.impute <- subset(dat_proc, TwophasesampIndD57 == 1)
summary(subset(dat.tmp.impute, Trt == 1 & Bserostatus==0)[imp.markers])

for (trt in unique(dat_proc$Trt)) {      
for (sero in unique(dat_proc$Bserostatus)) {    
  imp <- dat.tmp.impute %>%
    dplyr::filter(Trt == trt & Bserostatus==sero) %>%
    select(all_of(imp.markers)) 
    
  # deal with constant variables
  for (a in names(imp)) {
    if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
  }
    
  # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)
    
  dat.tmp.impute[dat.tmp.impute$Trt == trt & dat.tmp.impute$Bserostatus == sero, imp.markers] <- mice::complete(imp, action = 1)    
}
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

assays=c("bindSpike", "pseudoneutid50")
imp.markers=c(outer(c("B", "Day29"), assays, "%.%"))
dat.tmp.impute <- subset(dat_proc, TwophasesampIndD29 == 1)
summary(subset(dat.tmp.impute, Trt == 1)[imp.markers])

for (trt in unique(dat_proc$Trt)) {    
for (sero in unique(dat_proc$Bserostatus)) {    
  imp <- dat.tmp.impute %>%
    dplyr::filter(Trt == trt & Bserostatus==sero) %>%
    select(all_of(imp.markers)) 
  
  # deal with constant variables  
  for (a in names(imp)) {
    if (all(imp[[a]]==min(imp[[a]], na.rm=TRUE), na.rm=TRUE)) imp[[a]]=min(imp[[a]], na.rm=TRUE)
  }
  
  # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)
    
  dat.tmp.impute[dat.tmp.impute$Trt == trt & dat.tmp.impute$Bserostatus == sero, imp.markers] <- mice::complete(imp, action = 1)    
}
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

a="bindSpike"
lloq=62.8*0.0090 # 0.5652


# lloq censoring before taking delta
tmp=list()
for (t in c("B", "Day29", "Day57") ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(lloq), log10(lloq / 2), dat_proc[[t %.% a]])
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

dat_proc["Delta57overB"  %.% a] <- tmp["Day57" %.% a] - tmp["B"     %.% a]
dat_proc["Delta29overB"  %.% a] <- tmp["Day29" %.% a] - tmp["B"     %.% a]
dat_proc["Delta57over29" %.% a] <- tmp["Day57" %.% a] - tmp["Day29" %.% a]

 
#corplot(dat_proc$Day57bindSpike[dat_proc$Trt==1], dat_proc$Day57pseudoneutid50[dat_proc$Trt==1], xlab="LV MN50", ylab="PsV ID50", add.diag=F)
#summary(dat_proc$Day57bindSpike[dat_proc$Trt==1 & dat_proc$ph2.D57==1])
#mean(dat_proc$Day57bindSpike[dat_proc$Trt==1 & dat_proc$ph2.D57==1]>1.614, na.rm=T)
#summary(dat_proc$Day57pseudoneutid50[dat_proc$Trt==1 & dat_proc$ph2.D57==1])
#mean(dat_proc$Day57pseudoneutid50[dat_proc$Trt==1 & dat_proc$ph2.D57==1]>0.084, na.rm=T)
