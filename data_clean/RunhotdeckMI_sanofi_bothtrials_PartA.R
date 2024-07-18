# Make hotdeck data sets for Sanofi trials (both trials) final part A.
#
# T:/covpn/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/code/RunhotdeckMI_sanofi_bothtrials_PartA.R
# Peter Gilbert
# October 24, 2023
# Updated Feb 23, 2024 to match more the finalized stratified marker sampling plan.
#
# Sanofi Final Part A data sets: 
# /trials/covpn/p3005/analysis/mapping_immune_correlates/adata/COVID_Sanofi_stage1&2_20231024.csv
# T:/covpn/p3005/analysis/mapping_immune_correlates/combined/adata/COVID_Sanofi_stage1&2_20231013.csv
# Youyi: this should be the data file with all of the genotype marks in it and none of the hotdeck variables in it

# Add 10 hotdeck variable columns 
# seq1.variant.hotdeck1, ...., seq1.variant.hotdeck10, etc.


begin=Sys.time()
print(begin)

################################################################################

#renv::activate(here::here()) # manages the R packages, optional

# get mapped_data through config. 
source(here::here("_common.R"))
# Alternatively, e.g.
# mapped_data = 'T:/covpn/p3005/analysis/mapping_immune_correlates/combined/adata/COVID_Sanofi_stage1&2_20231024.csv'
# devtools::install_github("CoVPN/copcor")
# library(copcor) # installed from github CoVPN/copcor, needed for hotdeckMI


outputfile_name = sub(".csv", "_hotdeck.csv", mapped_data)

# quit if the output file already exists
if (file.exists(outputfile_name)) quit()

# For Sanofi the viral load variable is not used in neighborhoods for hotdeck multiple imputation

################################################################################
# read mapped data, adds sieve data

dat_mapped=read.csv(mapped_data)
# dat_mapped$EventSubvariantOmicron is the Omicron subvariant.  If this variable is !NA, it means Omicron.
# Not relevant though, as Omicron subvariant is not used for imputations; 
# we use seq1.variant that takes possible values NA, "Delta", or "Omicron".

################################################################################
# Run hotdeck imputation

# From Sanofi correlates SAP:

#Hotdeck multiple imputation uses nearest 
#neighborhoods $S_{ki}$ of participants with SARS-CoV-2 sequence observed ($\delta_{ki}=1$) based on z-scores of the following variables: 
#calendar time of COVID-19 failure time and matched on 
#coursened geographic region defined as
#(Honduras, not Honduras for the Stage 1 trial; India, Mexico, Other/Else for the Stage 2 trial).
# So the 5 region variables are (TrialStage==1 & Country=="Honduras", TrialStage==1 & Country!="Honduras", 
# TrialStage==2 & Country=="India", TrialStage==2 & Country=="Mexico", TrialStage==2 & (Country!="India" & Country!="Mexico"))

# Update Feb 2024 to harmonize with the updated correlates SAP:
# Definition of countries for stratified sampling:
#For the Stage 1 trial, the country category was defined by the four levels
#(United States, Japan) (Honduras, Colombia) (Kenya, Ghana) (India, Nepal).
#For the Stage 2 trial, 
#the country category was defined by the three levels
#(Colombia, Mexico) (Kenya, Uganda, Ghana), (India, Nepal), i.e., by continent Latin America, Africa, SE Asia. 

# For hotdeck, doesn't make sense to combine U.S. and Japan, so 
# For Stage 1 use United States, Japan, (Honduras, Colombia) (Kenya, Ghana) (India, Nepal).
# For Stage 2, use (Colombia, Mexico) (Kenya, Uganda, Ghana), (India, Nepal), i.e., by continent Latin America, Africa, SE Asia. 



# Define the event indicator based on the COVID-19 primary endpoint regardless of lineage information
EventInd <- ifelse(dat_mapped[,"EventIndKnownLineageOmicronD1"]==1 | 
dat_mapped[,"EventIndMissingLineageD1"]==1 | 
dat_mapped[,"EventIndKnownLineageNonOmicronD1"]==1 ,1,0)

# this kp is used at the end as well, so don't redefine it
kp <- dat_mapped[,"Perprotocol"]==1 & !is.na(EventInd)  
dat <- dat_mapped[kp,]

X <- rep(1,nrow(dat))
Delta <- EventInd 

# Z1discrete variables (matches required): Coursened country as noted above, no missing values allowed
#"Colombia" = 1
#"Ghana" = 2
#"Honduras" = 3
#"India" = 4
#"Japan" = 5
#"Kenya" = 6
#"Nepal" = 7
#"United States" = 8
#"Mexico" = 9
#"Uganda" = 10
#"Ukraine" = 11

#Country <- dat[,"Country"]
#TrialStage <- dat[,"Trialstage"]
#courseregion <- rep(NA,nrow(dat))
#courseregion[TrialStage==1 & Country==3] <- 0
#courseregion[TrialStage==1 & Country!=3] <- 1
#courseregion[TrialStage==2 & Country==4] <- 2
#courseregion[TrialStage==2 & Country==9] <- 3
#courseregion[TrialStage==2 & (Country!=4 & Country!=9)] <- 4


# Update:
# For Stage 1 use United States, Japan, (Honduras, Colombia) (Kenya, Ghana) (India, Nepal).
# For Stage 2, use (Colombia, Mexico) (Kenya, Uganda, Ghana), (India, Nepal), i.e., by continent Latin America, Africa, SE Asia. 

Country <- dat[,"Country"]
TrialStage <- dat[,"Trialstage"]
courseregion <- rep(NA,nrow(dat))
courseregion[TrialStage==1 & Country==8] <- 0
courseregion[TrialStage==1 & Country==5] <- 1
courseregion[TrialStage==1 & (Country==1 | Country==3)] <- 2
courseregion[TrialStage==1 & (Country==6 | Country==2)] <- 3
courseregion[TrialStage==1 & (Country==4 | Country==7)] <- 4

courseregion[TrialStage==2 & (Country==1 | Country==9)] <- 5
courseregion[TrialStage==2 & (Country==6 | Country==2 | Country==10)] <- 6
courseregion[TrialStage==2 & (Country==4 | Country==7)] <- 7

Z1discrete <- as.matrix(cbind(dat[,"Trt"],courseregion))
# Turn off the effect of Z1scalar on the nearest neighbors, as Avscalar is much more relevant and should carry the weight
Z1scalar <- matrix(rep(1,nrow(dat)),ncol=1)
#Z2 <- dat[,"Day29pseudoneutid50"]
Z2 <- rep(NA,nrow(dat))  # make all values of Z2 NA so that the immune marker does not affect the neighborhoods for hotdeck imputations
                         # that is the hotdeck ignores the immune marker
epsilonz <- ifelse(!is.na(Z2),1,0)

V <- dat[,"seq1.variant"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"seq1.variant"])] <- 0
epsilonv[!is.na(dat[,"seq1.variant"]) & dat[,"seq1.variant"]=="MissingLineage"] <- 0
Avdiscrete <- matrix(rep(NA,length(epsilonv)),ncol=1)
# Note "EventTimeFirstInfectionDate" is the date of the COVID-19 primary endpoint ignoring lineage information;
# the primary endpoint of the study
Numberdaysfromfirstpersonenrolleduntilprimaryendpoint <- as.Date(dat[,"EventTimeFirstInfectionDate"]) - as.Date(dat[,"FirstEnrollmentDate"])
Avscalar <- matrix(Numberdaysfromfirstpersonenrolleduntilprimaryendpoint,ncol=1)  
epsilona <- Delta
epsilona[is.na(Numberdaysfromfirstpersonenrolleduntilprimaryendpoint)] <- 0
M <- 10
L <- 5

ansvariant <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)


# Appends new columns to dat_mapped to make a new dataset

# use the same kp defined at the beginning
newdat <- cbind(dat_mapped,matrix(rep(NA,nrow(dat_mapped)*M),ncol=M))
j <- 0
# add new columns 
for (i in 1:nrow(newdat)) {
  if (kp[i]) {
    j <- j+1
    newdat[i,(ncol(dat_mapped)+1):(ncol(dat_mapped)+M)] <- ansvariant[[1]][j,] }}

newcolnames <-  c("seq1.variant.hotdeck1","seq1.variant.hotdeck2","seq1.variant.hotdeck3","seq1.variant.hotdeck4", 
                  "seq1.variant.hotdeck5","seq1.variant.hotdeck6","seq1.variant.hotdeck7","seq1.variant.hotdeck8",
                  "seq1.variant.hotdeck9","seq1.variant.hotdeck10")

colnames(newdat) <- c(colnames(dat_mapped),newcolnames)
write.csv(newdat, file=outputfile_name, row.names = F)

print("run time: "%.%format(Sys.time()-begin, digits=1))

# Note: For data analysis with a given hotdeck iteration, say 1 (seq1.variant.hotdeck1), 
# the failure time variables (EventIndOmicronD22hotdeck1,EventTimeOmicronD22hotdeck1=min(T,C)) 
# are computed as follows (for D22 marker correlates analysis):
# EventIndOmicronD22hotdeck1 <- ifelse(!is.na(seq1.variant.hotdeck1) & seq1.variant.hotdeck1=="Omicron"
# & !is.na(EventIndPrimaryD22),1,0)
# EventTimeOmicronD22hotdeck1 <- ifelse(!is.na(seq1.variant.hotdeck1) & seq1.variant.hotdeck1=="Omicron"
# & !is.na(EventIndPrimaryD22),min(EventTimeKnownLineageOmicronD22,EventTimeMissingLineageD22),
# max(EventTimeKnownLineageNonOmicronD22,EventTimeMissingLineageD22))

# Similarly for D43 marker correlates analyses (replace D22 with D43)



# 
# #######################################################################################
# # Checks on the data set 
# # Note that in make_dat_proc.R, we define EventIndOmicronD##hotdeck## and EventTimeOmicronD##hotdeck## variables
# # EventIndOmicronD22hotdeck1 and EventIndOmicronD22hotdeck10 defined below match the two variables defined in make_dat_proc.R perfectly
# 
# 
# newdat$EventIndPrimaryD1  = newdat$EventIndFirstInfectionD1
# newdat$EventIndPrimaryD22 = newdat$EventIndFirstInfectionD22
# newdat$EventIndPrimaryD43 = newdat$EventIndFirstInfectionD43
# 
# EventIndOmicronD22hotdeck1 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck1"]) 
#                                      & newdat[,"seq1.variant.hotdeck1"]=="Omicron"
#                                      & !is.na(newdat[,"EventIndPrimaryD22"]),1,0)
# 
# EventIndOmicronD22hotdeck10 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck10"]) 
#                                       & newdat[,"seq1.variant.hotdeck10"]=="Omicron"
#                                       & !is.na(newdat[,"EventIndPrimaryD22"]),1,0)
# 
# # Need to also define for hotdeck 2, ..., 9
# 
# table(newdat[,"seq1.variant.hotdeck1"],EventIndOmicronD22hotdeck1,useNA="ifany")
# table(newdat[,"seq1.variant.hotdeck10"],EventIndOmicronD22hotdeck10,useNA="ifany")
# 
# # Checks against original failure indicators:
# table(newdat[,"seq1.variant.hotdeck1"],newdat[,"EventIndKnownLineageOmicronD1"],useNA="ifany")
# table(newdat[,"seq1.variant.hotdeck10"],newdat[,"EventIndKnownLineageOmicronD1"],useNA="ifany")
# table(newdat[,"seq1.variant.hotdeck1"],newdat[,"EventIndPrimaryD1"],useNA="ifany")  #All cases
# table(newdat[,"seq1.variant.hotdeck10"],newdat[,"EventIndPrimaryD1"],useNA="ifany") #All cases
# 
# table(newdat[,"EventIndKnownLineageOmicronD1"],newdat[,"EventIndPrimaryD1"],useNA="ifany")
# table(newdat[,"EventIndKnownLineageNonOmicronD1"],newdat[,"EventIndPrimaryD1"],useNA="ifany")
# table(newdat[,"EventIndMissingLineageD1"],newdat[,"EventIndPrimaryD1"],useNA="ifany")
# # Adds up to 703. checks out.
# 
# # From hotdeck1, compare calendar time of case between Omicron and non-Omicron
# 
# kp1 <- EventIndOmicronD22hotdeck1==1  # Omicron
# kp2 <- EventIndOmicronD22hotdeck1==0 & newdat[,"EventIndPrimaryD22"]==1  # Non-Omicron
# kp <- kp1 | kp2
# grp <- ifelse(kp1[kp],1,0)
# 
# boxplot(as.integer(Numberdaysfromfirstpersonenrolleduntilprimaryendpoint)[kp] ~ grp,
#         col = 'steelblue',
#         main = 'Hotdeck imputation 1',
#         xlab = 'Non-Omicron vs. Omicron',
#         ylab = 'No. days to COVID-19'
# )
# 
# # Repeat for hotdeck10
# kp1 <- EventIndOmicronD22hotdeck10==1  # Omicron
# kp2 <- EventIndOmicronD22hotdeck10==0 & newdat[,"EventIndPrimaryD22"]==1  # Non-Omicron
# kp <- kp1 | kp2
# grp <- ifelse(kp1[kp],1,0)
# 
# boxplot(as.integer(Numberdaysfromfirstpersonenrolleduntilprimaryendpoint)[kp] ~ grp,
#         col = 'steelblue',
#         main = 'Hotdeck imputation 10',
#         xlab = 'Non-Omicron vs. Omicron',
#         ylab = 'No. days to COVID-19'
# )
# 
# 
# # Checks for D43 correlates (repeat the above)
# EventIndOmicronD43hotdeck1 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck1"]) 
#                                      & newdat[,"seq1.variant.hotdeck1"]=="Omicron"
#                                      & !is.na(newdat[,"EventIndPrimaryD43"]),1,0)
# 
# EventIndOmicronD43hotdeck10 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck10"]) 
#                                       & newdat[,"seq1.variant.hotdeck10"]=="Omicron"
#                                       & !is.na(newdat[,"EventIndPrimaryD43"]),1,0)
# 
# # Need to also define for hotdeck 2, ..., 9
# 
# table(newdat[,"seq1.variant.hotdeck1"],EventIndOmicronD43hotdeck1,useNA="ifany")
# table(newdat[,"seq1.variant.hotdeck10"],EventIndOmicronD43hotdeck10,useNA="ifany")
# 
# # Checks against original failure indicators:
# table(newdat[,"seq1.variant.hotdeck1"],newdat[,"EventIndKnownLineageOmicronD1"],useNA="ifany")
# table(newdat[,"seq1.variant.hotdeck10"],newdat[,"EventIndKnownLineageOmicronD1"],useNA="ifany")
# table(newdat[,"seq1.variant.hotdeck1"],newdat[,"EventIndPrimaryD1"],useNA="ifany")  #All cases
# table(newdat[,"seq1.variant.hotdeck10"],newdat[,"EventIndPrimaryD1"],useNA="ifany") #All cases
# 
# table(newdat[,"EventIndKnownLineageOmicronD1"],newdat[,"EventIndPrimaryD1"],useNA="ifany")
# table(newdat[,"EventIndKnownLineageNonOmicronD1"],newdat[,"EventIndPrimaryD1"],useNA="ifany")
# table(newdat[,"EventIndMissingLineageD1"],newdat[,"EventIndPrimaryD1"],useNA="ifany")
# 
# # From hotdeck1, compare calendar time of case between Omicron and non-Omicron
# 
# kp1 <- EventIndOmicronD43hotdeck1==1  # Omicron
# kp2 <- EventIndOmicronD43hotdeck1==0 & newdat[,"EventIndPrimaryD43"]==1  # Non-Omicron
# kp <- kp1 | kp2
# grp <- ifelse(kp1[kp],1,0)
# 
# boxplot(as.integer(Numberdaysfromfirstpersonenrolleduntilprimaryendpoint)[kp] ~ grp,
#         col = 'steelblue',
#         main = 'Hotdeck imputation 1',
#         xlab = 'Non-Omicron vs. Omicron',
#         ylab = 'No. days to COVID-19'
# )
# 
# # Repeat for hotdeck10
# kp1 <- EventIndOmicronD43hotdeck10==1  # Omicron
# kp2 <- EventIndOmicronD43hotdeck10==0 & newdat[,"EventIndPrimaryD43"]==1  # Non-Omicron
# kp <- kp1 | kp2
# grp <- ifelse(kp1[kp],1,0)
# 
# boxplot(as.integer(Numberdaysfromfirstpersonenrolleduntilprimaryendpoint)[kp] ~ grp,
#         col = 'steelblue',
#         main = 'Hotdeck imputation 10',
#         xlab = 'Non-Omicron vs. Omicron',
#         ylab = 'No. days to COVID-19'
# )
# 
# # Everything looks good, a few more checks.
# 
# table(dat_mapped[,"EventIndPrimaryD43"],dat_mapped[,"EventIndPrimaryD22"],useNA="ifany")
# table(dat_mapped[,"EventIndPrimaryD43"],useNA="ifany")
# table(dat_mapped[,"EventIndPrimaryD22"],useNA="ifany")
# 
# table(newdat[,"EventIndPrimaryD43"],newdat[,"EventIndPrimaryD22"],useNA="ifany")
# table(newdat[,"EventIndKnownLineageOmicronD43"],newdat[,"EventIndKnownLineageOmicronD22"],useNA="ifany")
# table(newdat[,"EventIndKnownLineageNonOmicronD43"],newdat[,"EventIndKnownLineageNonOmicronD22"],useNA="ifany")
# table(newdat[,"EventIndMissingLineageD43"],newdat[,"EventIndMissingLineageD22"],useNA="ifany")
# 
# EventInd1 <- ifelse(dat_mapped[,"EventIndKnownLineageOmicronD1"]==1 | 
#                       dat_mapped[,"EventIndMissingLineageD1"]==1 | 
#                       dat_mapped[,"EventIndKnownLineageNonOmicronD1"]==1 ,1,0)
# 
# EventInd22 <- ifelse(dat_mapped[,"EventIndKnownLineageOmicronD22"]==1 | 
#                        dat_mapped[,"EventIndMissingLineageD22"]==1 | 
#                        dat_mapped[,"EventIndKnownLineageNonOmicronD22"]==1 ,1,0)
# 
# EventInd43 <- ifelse(dat_mapped[,"EventIndKnownLineageOmicronD43"]==1 | 
#                        dat_mapped[,"EventIndMissingLineageD43"]==1 | 
#                        dat_mapped[,"EventIndKnownLineageNonOmicronD43"]==1 ,1,0)
# 
# table(EventInd1,EventInd22,EventInd43,useNA="ifany")
# 
# table(EventInd1,newdat[,"EventIndPrimaryD1"],useNA="ifany")
# table(EventInd22,newdat[,"EventIndPrimaryD22"],useNA="ifany")
# table(EventInd43,newdat[,"EventIndPrimaryD43"],useNA="ifany")



# #################################################################################
# # Move the rest of the material to data processing outside of this hotdeck data 
# # set construction file
# 
# EventIndOmicronD22hotdeck1 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck1"]) 
#                                      & newdat[,"seq1.variant.hotdeck1"]=="Omicron"
#                                      & !is.na(newdat[,"EventIndPrimaryD22"]),1,0)
# 
# EventTimeOmicronD22hotdeck1 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck1"]) & newdat[,"seq1.variant.hotdeck1"]=="Omicron"
#                                       & !is.na(newdat[,"EventIndPrimaryD22"]),min(newdat[,"EventTimeKnownLineageOmicronD22"],
#                                                                                   newdat[,"EventTimeMissingLineageD22"]),
#                                       max(newdat[,"EventTimeKnownLineageNonOmicronD22"],newdat[,"EventTimeMissingLineageD22"]))
# 
# # For defining EventTimeOmicronD22hotdeck1, consider special cases with 2 COVID-19 endpoints
# # 1. For the 7 participants with first event a known-lineage Omicron case, 
# #    this first event is counted as an endpoint and the second event is ignored.
# #    The above code works b/c the first event is known Omicron
# # 2. For the 2 participants with first event a missing-lineage case in 2021 and 
# #    second event known-lineage Omicron, the first event is ignored and the second 
# #    event is counted as an endpoint.
# #    The above code fails b/c the missing lineage event time would be used.
# 
# kp2 <- (PENDING, this flag is TRUE if one of the 2 participants and FALSE otherwise)
# EventTimeOmicronD22hotdeck1[kp2] <- newdat[,"EventTimeKnownLineageOmicronD22"]
# 
# # 3. For the 9 participants with first event a missing-lineage case in 2022 and second 
# #    event known-lineage Omicron, the first event is counted as an endpoint 
# #    (with hotdeck imputation employed as is generally done) and the second event is ignored.
# #    The above code works b/c EventTimeMissingLineageD22 corresponds to the first event.
# 
# 
# EventIndOmicronD22hotdeck10 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck10"]) 
#                                       & newdat[,"seq1.variant.hotdeck10"]=="Omicron"
#                                       & !is.na(newdat[,"EventIndPrimaryD22"]),1,0)
# 
# EventTimeOmicronD22hotdeck10 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck10"]) & newdat[,"seq1.variant.hotdeck10"]=="Omicron"
#                                        & !is.na(newdat[,"EventIndPrimaryD22"]),min(newdat[,"EventTimeKnownLineageOmicronD22"],
#                                                                                    newdat[,"EventTimeMissingLineageD22"]),
#                                        max(newdat[,"EventTimeKnownLineageNonOmicronD22"],newdat[,"EventTimeMissingLineageD22"]))
# 
# # Checks for D43 correlates (repeat the above)
# EventIndOmicronD43hotdeck1 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck1"]) 
#                                      & newdat[,"seq1.variant.hotdeck1"]=="Omicron"
#                                      & !is.na(newdat[,"EventIndPrimaryD43"]),1,0)
# 
# EventTimeOmicronD43hotdeck1 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck1"]) & newdat[,"seq1.variant.hotdeck1"]=="Omicron"
#                                       & !is.na(newdat[,"EventIndPrimaryD43"]),min(newdat[,"EventTimeKnownLineageOmicronD43"],
#                                                                                   newdat[,"EventTimeMissingLineageD43"]),
#                                       max(newdat[,"EventTimeKnownLineageNonOmicronD43"],newdat[,"EventTimeMissingLineageD43"]))
# 
# EventIndOmicronD43hotdeck10 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck10"]) 
#                                       & newdat[,"seq1.variant.hotdeck10"]=="Omicron"
#                                       & !is.na(newdat[,"EventIndPrimaryD43"]),1,0)
# 
# EventTimeOmicronD43hotdeck10 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck10"]) & newdat[,"seq1.variant.hotdeck10"]=="Omicron"
#                                        & !is.na(newdat[,"EventIndPrimaryD43"]),min(newdat[,"EventTimeKnownLineageOmicronD43"],
#                                                                                    newdat[,"EventTimeMissingLineageD43"]),
#                                        max(newdat[,"EventTimeKnownLineageNonOmicronD43"],newdat[,"EventTimeMissingLineageD43"]))
# 
# # For completeness, compute the survival variables for the D1 markers
# EventIndOmicronD1hotdeck1 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck1"]) 
#                                     & newdat[,"seq1.variant.hotdeck1"]=="Omicron"
#                                     & !is.na(newdat[,"EventIndPrimaryD1"]),1,0)
# 
# EventTimeOmicronD1hotdeck1 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck1"]) & newdat[,"seq1.variant.hotdeck1"]=="Omicron"
#                                      & !is.na(newdat[,"EventIndPrimaryD1"]),min(newdat[,"EventTimeKnownLineageOmicronD1"],
#                                                                                 newdat[,"EventTimeMissingLineageD1"]),
#                                      max(newdat[,"EventTimeKnownLineageNonOmicronD1"],newdat[,"EventTimeMissingLineageD1"]))
# 
# EventIndOmicronD1hotdeck10 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck10"]) 
#                                      & newdat[,"seq1.variant.hotdeck10"]=="Omicron"
#                                      & !is.na(newdat[,"EventIndPrimaryD1"]),1,0)
# 
# EventTimeOmicronD1hotdeck10 <- ifelse(!is.na(newdat[,"seq1.variant.hotdeck10"]) & newdat[,"seq1.variant.hotdeck10"]=="Omicron"
#                                       & !is.na(newdat[,"EventIndPrimaryD1"]),min(newdat[,"EventTimeKnownLineageOmicronD1"],
#                                                                                  newdat[,"EventTimeMissingLineageD1"]),
#                                       max(newdat[,"EventTimeKnownLineageNonOmicronD1"],newdat[,"EventTimeMissingLineageD1"]))
