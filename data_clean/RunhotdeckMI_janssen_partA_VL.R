# Make hotdeck data sets for Janssen ENSEMBLE final part A.
# The Day 29 pseudoneutralization ID50 titer marker is used.  

# ENSEMBLE Final Part A data set: 
# T:\covpn\p3003\analysis\correlates\Part_A_Blinded_Phase_Data\adata\janssen_la_partA_data_processed_with_riskscore.csv
# Youyi: this should be the data file with all of the genotype marks in it and none of the hotdeck variables in it

# For LA, with a total of 14 marks it adds 10*14 = 140 columns:
# seq1.spike.weighted.hamming.hotdeck1, ..., hotdeck10   seq1.s1.weighted.hamming.hotdeck1, ..., hotdeck10
# seq1.rbd.weighted.hamming.hotdeck1, ..., hotdeck10   seq1.ntd.weighted.hamming.hotdeck1, ..., hotdeck10
# seq1.variant.hotdeck1, ...., seq1.variant.hotdeck10, etc.
#
# Variable names of the 14 marks (13 distances and one the variant/lineage)

# Variable names of the 13 distances that are included for la:
#seq1.spike.weighted.hamming
#seq1.s1.weighted.hamming
#seq1.rbd.weighted.hamming
#seq1.ntd.weighted.hamming
#seq1.variant
# Distances with FWER p < 0.05 in the Magaret et al. ENSEMBLE sieve manuscript
#dms.seq1.RBD_antibody_escape_score
#dms.seq1.RBD_antibody_escape_score_cluster2
#dms.seq1.RBD_antibody_escape_score_cluster6
#dms.seq1.RBD_antibody_escape_score_cluster7
#dms.seq1.RBD_antibody_escape_score_cluster8
#pdb.seq1.mhrp.ab.dist.RBD4
#pdb.seq1.mhrp.ab.dist.RBD7
#pdb.seq1.mhrp.ab.dist.RBD8
#pdb.seq1.mhrp.ab.dist.NTD13

# source('T:/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/code/HotdeckMultipleImputationRfunction_v5.R')

# dat <- read.csv("T:/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/janssen_pooled_partA_data_processed_with_riskscore.csv")
# Youyi: this should be the data file with all of the genotype marks in it and none of the hotdeck variables in it, although the code
# still works with hotdeck variables in it

dat = dat_proc

kp <- dat[,"Bserostatus"]==0 & dat[,"Perprotocol"]==1 & dat[,"ph1.D29"]
dat <- dat[kp,]
X <- rep(1,nrow(dat))
Delta <- ifelse(dat[,"EventIndPrimaryIncludeNotMolecConfirmedD29"]==1 & !is.na(dat[,"seq1.log10vl"]),1,0)

# Z1discrete variables (matches required): Coursened region (4 levels LA Columbia, LA not Columbia, NA, SA)
# USA =0, ARG  = 1, BRA = 2, CHL = 3,COL = 4, MEX = 5,PER = 6, ZAF = 7, no missing values allowed

region <- dat[,"Country"]
courseregion <- rep(NA,nrow(dat))
courseregion[region==0] <- 0
courseregion[region==1 | region==2 | region==3 | region==5 | region==6] <- 1
courseregion[region==4] <- 2
courseregion[region==7] <- 3
Z1discrete <- as.matrix(cbind(dat[,"Trt"],courseregion))
# Turn off the effect of Z1scalar on the nearest neighbors, as Avscalar is much more relevant and should carry the weight
Z1scalar <- matrix(rep(1,nrow(dat)),ncol=1)
Z2 <- dat[,"Day29pseudoneutid50"]
epsilonz <- ifelse(dat[,"TwophasesampIndD29"],1,0)
epsilonz[dat[,"Trt"]==0] <- 0
Z2[epsilonz==0] <- NA

V <- dat[,"seq1.spike.weighted.hamming"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"seq1.spike.weighted.hamming"])] <- 0
Avdiscrete <- matrix(rep(NA,length(epsilonv)),ncol=1)
Avscalar <- as.matrix(dat[,c("Numberdaysfromfirstpersonenrolleduntilprimaryendpoint","seq1.log10vl")])
epsilona <- Delta
epsilona[is.na(dat[,"Numberdaysfromfirstpersonenrolleduntilprimaryendpoint"])] <- 0
epsilona[is.na(dat[,"seq1.log10vl"])] <- 0
M <- 10
L <- 5

anspooledspikedists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

# Repeat for S1 distance
V <- dat[,"seq1.s1.weighted.hamming"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"seq1.s1.weighted.hamming"])] <- 0
anspooledS1dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

# Repeat for RBD distance
V <- dat[,"seq1.rbd.weighted.hamming"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"seq1.rbd.weighted.hamming"])] <- 0
anspooledRBDdists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

# Repeat for NTD distance
V <- dat[,"seq1.ntd.weighted.hamming"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"seq1.ntd.weighted.hamming"])] <- 0
anspooledNTDdists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

# Repeat for variant
V <- dat[,"seq1.variant"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"seq1.variant"])] <- 0
anspooledvariant <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

# Repeat for the new epitope escape distances:

V <- dat[,"dms.seq1.RBD_antibody_escape_score"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"dms.seq1.RBD_antibody_escape_score"])] <- 0
anspooledDMSdists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

V <- dat[,"dms.seq1.RBD_antibody_escape_score_cluster2"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"dms.seq1.RBD_antibody_escape_score_cluster2"])] <- 0
anspooledDMS2dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

V <- dat[,"dms.seq1.RBD_antibody_escape_score_cluster6"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"dms.seq1.RBD_antibody_escape_score_cluster6"])] <- 0
anspooledDMS6dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

V <- dat[,"dms.seq1.RBD_antibody_escape_score_cluster7"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"dms.seq1.RBD_antibody_escape_score_cluster7"])] <- 0
anspooledDMS7dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

V <- dat[,"dms.seq1.RBD_antibody_escape_score_cluster8"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"dms.seq1.RBD_antibody_escape_score_cluster8"])] <- 0
anspooledDMS8dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

V <- dat[,"pdb.seq1.mhrp.ab.dist.RBD4"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"pdb.seq1.mhrp.ab.dist.RBD4"])] <- 0
anspooledPDB4dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

V <- dat[,"pdb.seq1.mhrp.ab.dist.RBD7"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"pdb.seq1.mhrp.ab.dist.RBD7"])] <- 0
anspooledPDB7dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

V <- dat[,"pdb.seq1.mhrp.ab.dist.RBD8"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"pdb.seq1.mhrp.ab.dist.RBD8"])] <- 0
anspooledPDB8dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

V <- dat[,"pdb.seq1.mhrp.ab.dist.NTD13"]
epsilonv <- rep(1,length(V))
epsilonv[is.na(dat[,"pdb.seq1.mhrp.ab.dist.NTD13"])] <- 0
anspooledPDB13dists <- hotdeckMI(X,Delta,Z1discrete,Z1scalar,Z2,epsilonz,V,epsilonv,Avdiscrete,Avscalar,epsilona,M,L)

#####
# Write out the data set:

dat <- dat_proc
# Youyi, this is the line where the columns with hotdeck variables need to be deleted.
kp <- dat[,"Bserostatus"]==0 & dat[,"Perprotocol"]==1 & dat[,"ph1.D29"]
newdat <- cbind(dat,matrix(rep(NA,nrow(dat)*140),ncol=140))
j <- 0
for (i in 1:nrow(newdat)) {
if (kp[i]) {
j <- j+1
newdat[i,(ncol(dat)+1):(ncol(dat)+140)] <- c(anspooledspikedists[[1]][j,],anspooledS1dists[[1]][j,],anspooledRBDdists[[1]][j,],
anspooledNTDdists[[1]][j,],anspooledvariant[[1]][j,],anspooledDMSdists[[1]][j,],anspooledDMS2dists[[1]][j,],anspooledDMS6dists[[1]][j,],
anspooledDMS7dists[[1]][j,],anspooledDMS8dists[[1]][j,],anspooledPDB4dists[[1]][j,],anspooledPDB7dists[[1]][j,],anspooledPDB8dists[[1]][j,],
anspooledPDB13dists[[1]][j,]) }}


newcolnames <- c("seq1.spike.weighted.hamming.hotdeck1","seq1.spike.weighted.hamming.hotdeck2","seq1.spike.weighted.hamming.hotdeck3",
                 "seq1.spike.weighted.hamming.hotdeck4","seq1.spike.weighted.hamming.hotdeck5","seq1.spike.weighted.hamming.hotdeck6",
                 "seq1.spike.weighted.hamming.hotdeck7","seq1.spike.weighted.hamming.hotdeck8","seq1.spike.weighted.hamming.hotdeck9",
                 "seq1.spike.weighted.hamming.hotdeck10",
                 "seq1.s1.weighted.hamming.hotdeck1","seq1.s1.weighted.hamming.hotdeck2","seq1.s1.weighted.hamming.hotdeck3",
                 "seq1.s1.weighted.hamming.hotdeck4","seq1.s1.weighted.hamming.hotdeck5","seq1.s1.weighted.hamming.hotdeck6",
                 "seq1.s1.weighted.hamming.hotdeck7","seq1.s1.weighted.hamming.hotdeck8","seq1.s1.weighted.hamming.hotdeck9",
                 "seq1.s1.weighted.hamming.hotdeck10",
                 "seq1.rbd.weighted.hamming.hotdeck1","seq1.rbd.weighted.hamming.hotdeck2","seq1.rbd.weighted.hamming.hotdeck3",
                 "seq1.rbd.weighted.hamming.hotdeck4","seq1.rbd.weighted.hamming.hotdeck5","seq1.rbd.weighted.hamming.hotdeck6",
                 "seq1.rbd.weighted.hamming.hotdeck7","seq1.rbd.weighted.hamming.hotdeck8","seq1.rbd.weighted.hamming.hotdeck9",
                 "seq1.rbd.weighted.hamming.hotdeck10",
                 "seq1.ntd.weighted.hamming.hotdeck1","seq1.ntd.weighted.hamming.hotdeck2","seq1.ntd.weighted.hamming.hotdeck3",
                 "seq1.ntd.weighted.hamming.hotdeck4","seq1.ntd.weighted.hamming.hotdeck5","seq1.ntd.weighted.hamming.hotdeck6",
                 "seq1.ntd.weighted.hamming.hotdeck7","seq1.ntd.weighted.hamming.hotdeck8","seq1.ntd.weighted.hamming.hotdeck9",
                 "seq1.ntd.weighted.hamming.hotdeck10",
                 "seq1.variant.hotdeck1","seq1.variant.hotdeck2","seq1.variant.hotdeck3","seq1.variant.hotdeck4", 
                 "seq1.variant.hotdeck5","seq1.variant.hotdeck6","seq1.variant.hotdeck7","seq1.variant.hotdeck8",
                 "seq1.variant.hotdeck9","seq1.variant.hotdeck10",
                 "dms.seq1.RBD_antibody_escape_score.hotdeck1","dms.seq1.RBD_antibody_escape_score.hotdeck2",
                 "dms.seq1.RBD_antibody_escape_score.hotdeck3","dms.seq1.RBD_antibody_escape_score.hotdeck4",
                 "dms.seq1.RBD_antibody_escape_score.hotdeck5","dms.seq1.RBD_antibody_escape_score.hotdeck6",
                 "dms.seq1.RBD_antibody_escape_score.hotdeck7","dms.seq1.RBD_antibody_escape_score.hotdeck8",
                 "dms.seq1.RBD_antibody_escape_score.hotdeck9","dms.seq1.RBD_antibody_escape_score.hotdeck10",
                 "dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck1","dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck2",
                 "dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck3","dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck4",
                 "dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck5","dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck6",
                 "dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck7","dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck8",
                 "dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck9","dms.seq1.RBD_antibody_escape_score_cluster2.hotdeck10",
                 "dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck1","dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck2",
                 "dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck3","dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck4",
                 "dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck5","dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck6",
                 "dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck7","dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck8",
                 "dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck9","dms.seq1.RBD_antibody_escape_score_cluster6.hotdeck10",
                 "dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck1","dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck2",
                 "dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck3","dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck4",
                 "dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck5","dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck6",
                 "dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck7","dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck8",
                 "dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck9","dms.seq1.RBD_antibody_escape_score_cluster7.hotdeck10",
                 "dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck1","dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck2",
                 "dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck3","dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck4",
                 "dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck5","dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck6",
                 "dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck7","dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck8",
                 "dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck9","dms.seq1.RBD_antibody_escape_score_cluster8.hotdeck10",
                 "pdb.seq1.mhrp.ab.dist.RBD4.hotdeck1","pdb.seq1.mhrp.ab.dist.RBD4.hotdeck2",
                 "pdb.seq1.mhrp.ab.dist.RBD4.hotdeck3","pdb.seq1.mhrp.ab.dist.RBD4.hotdeck4",
                 "pdb.seq1.mhrp.ab.dist.RBD4.hotdeck5","pdb.seq1.mhrp.ab.dist.RBD4.hotdeck6",
                 "pdb.seq1.mhrp.ab.dist.RBD4.hotdeck7","pdb.seq1.mhrp.ab.dist.RBD4.hotdeck8",
                 "pdb.seq1.mhrp.ab.dist.RBD4.hotdeck9","pdb.seq1.mhrp.ab.dist.RBD4.hotdeck10",
                 "pdb.seq1.mhrp.ab.dist.RBD7.hotdeck1","pdb.seq1.mhrp.ab.dist.RBD7.hotdeck2",
                 "pdb.seq1.mhrp.ab.dist.RBD7.hotdeck3","pdb.seq1.mhrp.ab.dist.RBD7.hotdeck4",
                 "pdb.seq1.mhrp.ab.dist.RBD7.hotdeck5","pdb.seq1.mhrp.ab.dist.RBD7.hotdeck6",
                 "pdb.seq1.mhrp.ab.dist.RBD7.hotdeck7","pdb.seq1.mhrp.ab.dist.RBD7.hotdeck8",
                 "pdb.seq1.mhrp.ab.dist.RBD7.hotdeck9","pdb.seq1.mhrp.ab.dist.RBD7.hotdeck10",
                 "pdb.seq1.mhrp.ab.dist.RBD8.hotdeck1","pdb.seq1.mhrp.ab.dist.RBD8.hotdeck2",
                 "pdb.seq1.mhrp.ab.dist.RBD8.hotdeck3","pdb.seq1.mhrp.ab.dist.RBD8.hotdeck4",
                 "pdb.seq1.mhrp.ab.dist.RBD8.hotdeck5","pdb.seq1.mhrp.ab.dist.RBD8.hotdeck6",
                 "pdb.seq1.mhrp.ab.dist.RBD8.hotdeck7","pdb.seq1.mhrp.ab.dist.RBD8.hotdeck8",
                 "pdb.seq1.mhrp.ab.dist.RBD8.hotdeck9","pdb.seq1.mhrp.ab.dist.RBD8.hotdeck10",
                 "pdb.seq1.mhrp.ab.dist..NTD13.hotdeck1","pdb.seq1.mhrp.ab.dist..NTD13.hotdeck2",
                 "pdb.seq1.mhrp.ab.dist..NTD13.hotdeck3","pdb.seq1.mhrp.ab.dist..NTD13.hotdeck4",
                 "pdb.seq1.mhrp.ab.dist..NTD13.hotdeck5","pdb.seq1.mhrp.ab.dist..NTD13.hotdeck6",
                 "pdb.seq1.mhrp.ab.dist..NTD13.hotdeck7","pdb.seq1.mhrp.ab.dist..NTD13.hotdeck8",
                 "pdb.seq1.mhrp.ab.dist..NTD13.hotdeck9","pdb.seq1.mhrp.ab.dist..NTD13.hotdeck10")

colnames(newdat) <- c(colnames(dat),newcolnames)
# write.csv(newdat, file="T:/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/janssen_pooled_partA_seq1_variant_hamming_hotdeckv7.csv")

dat_proc = newdat

# # Data checks:
# 
# mt <- read.csv("T:/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/janssen_pooled_partA_seq1_variant_hamming_hotdeckv7.csv")
# 
# kp <- mt[,"EventIndPrimaryIncludeNotMolecConfirmedD29"]==1 & !is.na(mt[,"seq1.log10vl"]) & mt[,"Bserostatus"]==0 & mt[,"Perprotocol"]==1 & mt[,"ph1.D29"]
# # 1046
# table(mt[kp,"seq1.variant"])           # 890 with observed data
# # no ""s
# table(mt[kp,"seq1.variant.hotdeck1"])  # 1046 with observed data
# # a handful have "" imputed!
# sum(table(mt[kp,"seq1.variant.hotdeck10"]))
# # Same story
# subset(mt, Ptid=="VAC31518COV3001-3009010")
