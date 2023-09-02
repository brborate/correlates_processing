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


begin=Sys.time()
print(begin)


################################################################################

renv::activate(here::here()) # manages the R packages, optional

# get mapped_data through config. 
source(here::here("_common.R")) 
# Alternatively, 
# mapped_data = '/trials/covpn/p3003/analysis/mapping_immune_correlates/adata/COVID_ENSEMBLE_PartAComplete_variant_mapped_20230809.csv'

# quit if the output file already exists
outputfile_name = sub(".csv", "_hotdeck.csv", mapped_data)
if (file.exists(outputfile_name)) quit()


################################################################################
# read mapped data, adds sieve data

dat_mapped=read.csv(mapped_data)

# add Spike physics-chemical weighted Hamming distance pertaining to the sequence that was obtained from the first chronological sample
dat_tmp = read.csv("/trials/covpn/p3003/analysis/post_covid/sieve/Part_A_Blinded_Phase_Data/adata/omnibus/cpn3003_sieve_cases_firstseq_v10a.csv")

dat_mapped$seq1.log10vl = dat_tmp$seq1.log10vl[match(dat_mapped$Subjectid, dat_tmp$USUBJID)]
dat_mapped$seq1.variant = dat_tmp$seq1.who.label[match(dat_mapped$Subjectid, dat_tmp$USUBJID)]

new.names = c("seq1.spike.weighted.hamming", "seq1.s1.weighted.hamming", "seq1.rbd.weighted.hamming", "seq1.ntd.weighted.hamming",
              "dms.seq1.RBD_antibody_escape_score", 'dms.seq1.RBD_antibody_escape_score_cluster2', 'dms.seq1.RBD_antibody_escape_score_cluster6', 'dms.seq1.RBD_antibody_escape_score_cluster7', 'dms.seq1.RBD_antibody_escape_score_cluster8', 
              'pdb.seq1.mhrp.ab.dist.RBD4', 'pdb.seq1.mhrp.ab.dist.RBD7', 'pdb.seq1.mhrp.ab.dist.RBD8', 
              'pdb.seq1.mhrp.ab.dist.NTD13')
old.names = c('seq1.hdist.zspace.spike', 'seq1.hdist.zspace.s1', 'seq1.hdist.zspace.rbd', 'seq1.hdist.zspace.ntd', 
              'seq1.RBD_antibody_escape_score', 'seq1.RBD_antibody_escape_score_cluster2', 'seq1.RBD_antibody_escape_score_cluster6', 'seq1.RBD_antibody_escape_score_cluster7', 'seq1.RBD_antibody_escape_score_cluster8',
              'seq1.mhrp.ab.dist.RBD4', 'seq1.mhrp.ab.dist.RBD7', 'seq1.mhrp.ab.dist.RBD8',
              'seq1.mhrp.ab.dist.NTD13' )
for (i in 1:length(new.names)) {
  dat_mapped[[new.names[i]]] = dat_tmp[[old.names[i]]][match(dat_mapped$Subjectid, dat_tmp$USUBJID)]
}



################################################################################
# run hotdeck imputation

library(copcor) # installed from github CoVPN/copcor, needed for hotdeckMI

# this kp is used at the end as well, so don't redefine it
kp <- dat_mapped[,"Bserostatus"]==0 & dat_mapped[,"Perprotocol"]==1 & !is.na(dat_mapped[,"EventIndPrimaryIncludeNotMolecConfirmedD29"])
dat <- dat_mapped[kp,]

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
epsilonz <- ifelse(!is.na(Z2),1,0) 
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



################################################################################
# Appends new columns to dat_mapped to make a new dataset

# use the same kp defined at the beginning
newdat <- cbind(dat_mapped,matrix(rep(NA,nrow(dat_mapped)*140),ncol=140))
j <- 0
# add new columns 
for (i in 1:nrow(newdat)) {
  if (kp[i]) {
    j <- j+1
    newdat[i,(ncol(dat_mapped)+1):(ncol(dat_mapped)+140)] <- c(anspooledspikedists[[1]][j,],anspooledS1dists[[1]][j,],anspooledRBDdists[[1]][j,],
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
                 "pdb.seq1.mhrp.ab.dist.NTD13.hotdeck1","pdb.seq1.mhrp.ab.dist.NTD13.hotdeck2",
                 "pdb.seq1.mhrp.ab.dist.NTD13.hotdeck3","pdb.seq1.mhrp.ab.dist.NTD13.hotdeck4",
                 "pdb.seq1.mhrp.ab.dist.NTD13.hotdeck5","pdb.seq1.mhrp.ab.dist.NTD13.hotdeck6",
                 "pdb.seq1.mhrp.ab.dist.NTD13.hotdeck7","pdb.seq1.mhrp.ab.dist.NTD13.hotdeck8",
                 "pdb.seq1.mhrp.ab.dist.NTD13.hotdeck9","pdb.seq1.mhrp.ab.dist.NTD13.hotdeck10")

colnames(newdat) <- c(colnames(dat_mapped),newcolnames)
write.csv(newdat, file=outputfile_name, row.names = F)

print("run time: "%.%format(Sys.time()-begin, digits=1))


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
# subset(mt, Subjectid=="VAC31518COV3001-3009010")

