# read LVMN marker data
dat.1=read.csv("/trials/covpn/p3004/analysis/mapping_immune_correlates/adata/COVID_Novavax_RBD_20221005.csv") 

for (t in c("B","Day35")) {
    x=dat.1[[t%.%"bindRBD"]]
    dat_proc[[t%.%"bindRBD"]]=x[match(dat_proc$Ptid, dat.1$Subjectid)]
}

with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D35,], table(!is.na(BbindRBD), Wstratum))

with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D35,], table(!is.na(Day35bindRBD), Wstratum))


with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D35,], table(!is.na(BbindRBD), Wstratum, SubcohortInd))

with(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D35,], table(!is.na(Day35bindRBD), Wstratum, SubcohortInd))

# 0 cases missing RBD
subset(dat_proc[dat_proc$Trt==1 & dat_proc$ph2.D35,], (is.na(BbindRBD) | is.na(Day35bindRBD) ) & Wstratum==43, select=c(Ptid, BbindRBD, Day35bindRBD))




###############################################################################
# impute missing RBD biomarkers in ph2
#     impute vaccine and placebo (we expect data to have baseline neg only or baseline pos only)
#     use bindSpike and PsV
#     use baseline and D35, but not Delta
###############################################################################

set.seed(2022)

n.imp <- 1

assays=c("bindSpike", "bindRBD", "pseudoneutid50")
imp.markers=c(outer(c("B", "Day35"), assays, "%.%"))
dat.tmp.impute <- subset(dat_proc, TwophasesampIndD35 == 1)
summary(subset(dat.tmp.impute, Trt == 1 & Bserostatus==0)[imp.markers])

for (trt in unique(dat_proc$Trt)) {      
for (sero in unique(dat_proc$Bserostatus)) {    
  imp <- dat.tmp.impute %>% dplyr::filter(Trt == trt & Bserostatus==sero) %>% select(all_of(imp.markers)) 
    
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
dat_proc[dat_proc$TwophasesampIndD35==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc$TwophasesampIndD35==1, "Ptid"], dat.tmp.impute$Ptid), ]

# imputed values of missing markers merged properly for all individuals in the two phase sample?
assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc$TwophasesampIndD35 == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)



###############################################################################
# define delta for dat_proc
###############################################################################

a="bindRBD"
lloq=1126.7*0.0272

# lloq censoring before taking delta
tmp=list()
for (t in c("B", "Day35") ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] < log10(lloq), log10(lloq / 2), dat_proc[[t %.% a]])
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame

dat_proc["Delta35overB"  %.% a] <- tmp["Day35" %.% a] - tmp["B"     %.% a]

 
