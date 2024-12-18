dat <- read.csv("T:/covpn/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/adata/vat08_combined_data_processed_20240723.csv")
# t-tests for whether antibody markers differ by region

#######################################################################################
# Checking hotdeck imputations
#
# violin boxplot of times from first person enrolled until COVID-19 endpoint for
# Stage 1 D43 marker primary endpoints, stratified by Omicron vs. not Omicron

# Restrict to cases with a COVID-19 endpoint
keep <- dat$Trialstage==1 & dat$Bserostatus==1 & dat$ph2.D43.nAb==TRUE & !is.na(dat$seq1.variant.hotdeck10) & 
dat$EventTimeFirstInfectionD43 >=7 & dat$EventTimeFirstInfectionD22 <= 180 &
dat$EventIndFirstInfectionD43==1

dat <- dat[keep,]

# Install and load the vioplot package
install.packages("vioplot")
library(vioplot)

#calendartimetoCOVID <- dat$CalendarDateEnrollment + dat$EventTimeFirstInfectionD1
# matches the way it is done in the hotdeck code, except is one day short, so use the hotdeck code:
calendartimetoCOVID <- as.Date(dat[,"EventTimeFirstInfectionDate"]) - as.Date(dat[,"FirstEnrollmentDate"])
calendartimetoCOVID <- as.numeric(calendartimetoCOVID)


# Definition of having an event in the Runhotdeck code:
# EventInd <- ifelse(dat_mapped[,"EventIndKnownLineageOmicronD1"]==1 | 
# dat_mapped[,"EventIndMissingLineageD1"]==1 | 
# dat_mapped[,"EventIndKnownLineageNonOmicronD1"]==1 ,1,0)


group11 <- calendartimetoCOVID[dat$seq1.variant.hotdeck1=="Omicron"]
group12 <- calendartimetoCOVID[dat$seq1.variant.hotdeck1!="Omicron"]
group21 <- calendartimetoCOVID[dat$seq1.variant.hotdeck2=="Omicron"]
group22 <- calendartimetoCOVID[dat$seq1.variant.hotdeck2!="Omicron"]
group31 <- calendartimetoCOVID[dat$seq1.variant.hotdeck3=="Omicron"]
group32 <- calendartimetoCOVID[dat$seq1.variant.hotdeck3!="Omicron"]
group41 <- calendartimetoCOVID[dat$seq1.variant.hotdeck4=="Omicron"]
group42 <- calendartimetoCOVID[dat$seq1.variant.hotdeck4!="Omicron"]
group51 <- calendartimetoCOVID[dat$seq1.variant.hotdeck5=="Omicron"]
group52 <- calendartimetoCOVID[dat$seq1.variant.hotdeck5!="Omicron"]
group61 <- calendartimetoCOVID[dat$seq1.variant.hotdeck6=="Omicron"]
group62 <- calendartimetoCOVID[dat$seq1.variant.hotdeck6!="Omicron"]
group71 <- calendartimetoCOVID[dat$seq1.variant.hotdeck7=="Omicron"]
group72 <- calendartimetoCOVID[dat$seq1.variant.hotdeck7!="Omicron"]
group81 <- calendartimetoCOVID[dat$seq1.variant.hotdeck8=="Omicron"]
group82 <- calendartimetoCOVID[dat$seq1.variant.hotdeck8!="Omicron"]
group91 <- calendartimetoCOVID[dat$seq1.variant.hotdeck9=="Omicron"]
group92 <- calendartimetoCOVID[dat$seq1.variant.hotdeck9!="Omicron"]
group101 <- calendartimetoCOVID[dat$seq1.variant.hotdeck10=="Omicron"]
group102 <- calendartimetoCOVID[dat$seq1.variant.hotdeck10!="Omicron"]


min(c(group12,group22,group32,group42,group52,group62,group72,group82,group92,group102))
max(c(group12,group22,group32,group42,group52,group62,group72,group82,group92,group102))
# min = 147 days = October 20, 2021, max = 217 days = December 29, 2021


# Range of Not Omicron/Delta imputations: October 20 2021 to December 29, 2021

# Latest date of an actual Delta COVID-19: 
# Calculate the latest actual observed Delta COVID-19 onset:
mm <- max(calendartimetoCOVID[dat$seq1.variant!="MissingLineage" & dat$seq1.variant!="Omicron"])
# 217 days = December 29, 2021    # Fixed
dat$seq1.variant[calendartimetoCOVID==mm]
#[1] "Delta" "Delta" "Delta"

min(c(group11,group21,group31,group41,group51,group61,group71,group81,group91,group101))
max(c(group11,group21,group31,group41,group51,group61,group71,group81,group91,group101))
#min = 147 days = October 20, 2021, max = 272 days = February 22, 2022
# 2nd smallest is 191 days: November 14, 2021

# Range of Omicron imputations: October 20, 2021 to February 22, 2022

# Number D43 marker endpoints = 116


yy <- calendartimetoCOVID
xx1 <- ifelse(dat$seq1.variant.hotdeck1=="Omicron",1,0)
xx2 <- ifelse(dat$seq1.variant.hotdeck2=="Omicron",1,0)
xx3 <- ifelse(dat$seq1.variant.hotdeck3=="Omicron",1,0)
xx4 <- ifelse(dat$seq1.variant.hotdeck4=="Omicron",1,0)
xx5 <- ifelse(dat$seq1.variant.hotdeck5=="Omicron",1,0)
xx6 <- ifelse(dat$seq1.variant.hotdeck6=="Omicron",1,0)
xx7 <- ifelse(dat$seq1.variant.hotdeck7=="Omicron",1,0)
xx8 <- ifelse(dat$seq1.variant.hotdeck8=="Omicron",1,0)
xx9 <- ifelse(dat$seq1.variant.hotdeck9=="Omicron",1,0)
xx10 <- ifelse(dat$seq1.variant.hotdeck10=="Omicron",1,0)

table(xx1)
table(xx2)
table(xx3)
table(xx4)
table(xx5)
table(xx6)
table(xx7)
table(xx8)
table(xx9)
table(xx10)

# hotdeck1 is the modal result: 5 Delta, 111 Omicron
 sort(calendartimetoCOVID[xx1==1])
 sort(calendartimetoCOVID[xx1==1])
  [1] 191 197 205 211 212 214 218 218 219 219 221 221 222 222 222 222 222 222 223 223 223
 [22] 226 227 227 227 227 227 227 228 229 230 230 230 230 230 231 231 231 232 232 232 233
 [43] 233 233 233 233 235 235 235 235 236 236 236 236 236 236 236 237 238 239 239 239 239
 [64] 239 239 239 239 239 239 240 240 240 240 240 240 240 240 240 240 241 241 241 241 241
 [85] 241 242 242 242 243 244 244 245 245 245 245 246 247 247 247 247 248 249 249 250 250
[106] 259 261 266 271 272


sort(calendartimetoCOVID[xx7==1])
  [1] 191 197 205 207 211 212 214 218 218 219 219 221 221 222 222 222 222 222 222 223 223
 [22] 223 226 227 227 227 227 227 227 228 229 230 230 230 230 230 231 231 231 232 232 232
 [43] 233 233 233 233 233 235 235 235 235 236 236 236 236 236 236 236 237 238 239 239 239
 [64] 239 239 239 239 239 239 239 240 240 240 240 240 240 240 240 240 240 241 241 241 241
 [85] 241 241 242 242 242 243 244 244 245 245 245 245 246 247 247 247 247 248 249 249 250
[106] 250 259 261 266 271 272
# Resolve by moving to hotdeck 7, also with 5 Delta and 111 Omicron

1,2,4,5,6,7,10 candidate hotdecks because earliest Omicron event is >= 191 days
3,6,10 are 5 Delta, 111 Omicron

# So choose hotdeck 6 or 10

# 191 days is December 3, 2021

# Numbers with known lineage
table(dat$seq1.variant)
Delta MissingLineage        Omicron 
             4             58             54 

# So the L=5-nearest neighbors are always defined from the
# 54 known Omicron and 4 known Delta and filling in 58 cases with missing values.
# Sampling with replacement

# What are the calendar times of the 4 observed Delta viruses?

calendartimetoCOVID[dat$seq1.variant=="Delta"]
[1] 216 217 217 217

# What are the calendar times of the 58 with missing lineages?

sort(calendartimetoCOVID[dat$seq1.variant=="MissingLineage"])
 [1] 147 191 197 207 214 218 219 222 222 223 227 227 227 227 227 227 228 230 231 231 231 232 232 232 233 233 233 235 235 235 236 236
[33] 236 238 239 239 239 239 239 239 239 239 240 240 240 240 241 241 241 242 243 247 247 248 249 261 266 272

# What are the calendar times of the 54 with known Omicron?
sort(calendartimetoCOVID[dat$seq1.variant=="Omicron"])

 [1] 205 211 212 218 219 221 221 222 222 222 222 223 223 226 229 230 230 230 230 233 233 235 236 236 236 236 237 239 239 240 240 240
[33] 240 240 240 241 241 241 242 242 244 244 245 245 245 245 246 247 247 249 250 250 259 271

# So, the 147, 191, 197, 207, 214 missing ones all have closer known Omicron cases then the known Deltas.
# The 218 missing one would likely select Omicron.  The 214 missing one could select Delta or Omicron.
# Really only the 214, 218, 219 ones could potentially select Omicron
sort(calendartimetoCOVID[xx1==0])
sort(calendartimetoCOVID[xx2==0])
sort(calendartimetoCOVID[xx3==0])
sort(calendartimetoCOVID[xx4==0])
sort(calendartimetoCOVID[xx5==0])
sort(calendartimetoCOVID[xx6==0])
sort(calendartimetoCOVID[xx7==0])
sort(calendartimetoCOVID[xx8==0])
sort(calendartimetoCOVID[xx9==0])
sort(calendartimetoCOVID[xx10==0])

sort(calendartimetoCOVID[xx1==1])
sort(calendartimetoCOVID[xx2==1])
sort(calendartimetoCOVID[xx3==1])
sort(calendartimetoCOVID[xx4==1])
sort(calendartimetoCOVID[xx5==1])
sort(calendartimetoCOVID[xx6==1])
sort(calendartimetoCOVID[xx7==1])
sort(calendartimetoCOVID[xx8==1])
sort(calendartimetoCOVID[xx9==1])
sort(calendartimetoCOVID[xx10==1])


# Note: For a given subject i that is a case with missing lineage, a lineage from a nearest neighbor is randomly sampled,
but the values of this subject i's other variables are retained.  That is, subject i does NOT have its variable values 
(like calendar time from first person enrolled until COVID-19) changed to be those values from the selected nearest neighbor.

# Average number omicron
mean(c(table(xx1)[2],table(xx2)[2],table(xx3)[2],table(xx4)[2],table(xx5)[2],
table(xx6)[2],table(xx7)[2],table(xx8)[2],table(xx9)[2],table(xx10)[2]))
110.4
# So average Not Omicron = 5.6
# Range number Omicron
range(c(table(xx1)[2],table(xx2)[2],table(xx3)[2],table(xx4)[2],table(xx5)[2],
table(xx6)[2],table(xx7)[2],table(xx8)[2],table(xx9)[2],table(xx10)[2]))
109 to 112

# Percentage agreement across all 45 pairs

vect <- c(mean(xx1==xx2),mean(xx1==xx3),mean(xx1==xx4),mean(xx1==xx5),mean(xx1==xx6),mean(xx1==xx7),
mean(xx1==xx8),mean(xx1==xx9),mean(xx1==xx10),mean(xx2==xx3),mean(xx2==xx4),mean(xx2==xx5),
mean(xx2==xx6),mean(xx2==xx7),mean(xx2==xx8),mean(xx2==xx9),mean(xx2==xx10),mean(xx3==xx4),
mean(xx3==xx5),mean(xx3==xx6),mean(xx3==xx7),mean(xx3==xx8),mean(xx3==xx9),mean(xx3==xx10),
mean(xx4==xx5),mean(xx4==xx6),mean(xx4==xx7),mean(xx4==xx8),mean(xx4==xx9),mean(xx4==xx10),
mean(xx5==xx6),mean(xx5==xx7),mean(xx5==xx8),mean(xx5==xx9),mean(xx5==xx10),mean(xx6==xx7),
mean(xx6==xx8),mean(xx6==xx9),mean(xx6==xx10),mean(xx7==xx8),mean(xx7==xx9),mean(xx7==xx10),
mean(xx8==xx9),mean(xx8==xx10),mean(xx9==xx10))

mean(vect)
0.9842912 concordance in lineage calls averaging across all pairs of the 10 vectors (each of length 116) of lineage imputations

# Create side-by-side violin plots
pdf("T:/covpn/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/reports/TablesFigures/Stage1HotdeckcalendartimesFINAL.pdf")
par(mfrow=c(3,4),cex.axis=0.9,cex.lab=0.9,las=2,oma=c(3,3,6,3),mai = c(0.5, 0.1, 0.5, 0.1))
vioplot(yy ~ xx1, names = c("HD1 Omicron", "HD1 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx1, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 1")

vioplot(yy ~ xx2, names = c("HD2 Omicron", "HD2 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx2, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 2")

vioplot(yy ~ xx3, names = c("HD3 Omicron", "HD3 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx3, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 3")

vioplot(yy ~ xx4, names = c("HD4 Omicron", "HD4 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx4, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 4")

vioplot(yy ~ xx5, names = c("HD5 Omicron", "HD5 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx5, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 5")

vioplot(yy ~ xx6, names = c("HD6 Omicron", "HD6 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx6, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 6")

vioplot(yy ~ xx7, names = c("HD7 Omicron", "HD7 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx7, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 7")

vioplot(yy ~ xx8, names = c("HD8 Omicron", "HD8 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx8, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 8")

vioplot(yy ~ xx9, names = c("HD9 Omicron", "HD9 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy ~ xx9, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 9")

vioplot(yy ~ xx10, names = c("HD10 Omicron", "HD10 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,310),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22"))
stripchart(yy~xx10, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 10")

plot(1,1,type='n',axes=FALSE,xlab="",ylab="")
text(x=1,y=1.1, "Lineage Imputed:",cex=0.72)
legend(x="bottomright",c("Omicron", "Not Omicron"),fill=TRUE,col=c(2,1),
border=c(2,1),cex=0.722)

vioplot(vect,col=4,ylab="",xlab="",axes=FALSE,horizontal=TRUE,ylim=c(0.96,1.03),cex=0.9)
axis(1,at=c(0.97,0.98,0.99,1.0),labels=c("97%","98%","99%","100%"))
stripchart(vect, method = "jitter", jitter=0.19, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("% Agreement")

mtext("Stage 1 Non-Naive: Number Days post FPI to COVID-19 by Imputed Omicron vs. Not Omicron",
outer=T,cex=0.95,line=2,side=3,las=1)

dev.off()

# First person enrolled: (05/26/2021) 
# In the caption include the number of cases: I think 116.  Put lots of important details in the caption.
# and report the average number non-Omicron ... etc.


############################################
# Repeat hotdeck checks for the Stage 2 trial

# Restrict to cases with a COVID-19 endpoint  
keep <- dat$Trialstage==2 & dat$Bserostatus==1 & dat$ph2.D43.nAb==TRUE & !is.na(dat$seq1.variant.hotdeck10) & 
dat$EventTimeFirstInfectionD43 >=7 & dat$EventTimeFirstInfectionD22 <= 180 &
dat$EventIndFirstInfectionD43==1

dat <- dat[keep,]
calendartimetoCOVID <- as.Date(dat[,"EventTimeFirstInfectionDate"]) - as.Date("2021-05-26") # First enrollment date for the Stage 1 trial
calendartimetoCOVID <- as.numeric(calendartimetoCOVID)

group11 <- calendartimetoCOVID[dat$seq1.variant.hotdeck1=="Omicron"]
group12 <- calendartimetoCOVID[dat$seq1.variant.hotdeck1!="Omicron"]
group21 <- calendartimetoCOVID[dat$seq1.variant.hotdeck2=="Omicron"]
group22 <- calendartimetoCOVID[dat$seq1.variant.hotdeck2!="Omicron"]
group31 <- calendartimetoCOVID[dat$seq1.variant.hotdeck3=="Omicron"]
group32 <- calendartimetoCOVID[dat$seq1.variant.hotdeck3!="Omicron"]
group41 <- calendartimetoCOVID[dat$seq1.variant.hotdeck4=="Omicron"]
group42 <- calendartimetoCOVID[dat$seq1.variant.hotdeck4!="Omicron"]
group51 <- calendartimetoCOVID[dat$seq1.variant.hotdeck5=="Omicron"]
group52 <- calendartimetoCOVID[dat$seq1.variant.hotdeck5!="Omicron"]
group61 <- calendartimetoCOVID[dat$seq1.variant.hotdeck6=="Omicron"]
group62 <- calendartimetoCOVID[dat$seq1.variant.hotdeck6!="Omicron"]
group71 <- calendartimetoCOVID[dat$seq1.variant.hotdeck7=="Omicron"]
group72 <- calendartimetoCOVID[dat$seq1.variant.hotdeck7!="Omicron"]
group81 <- calendartimetoCOVID[dat$seq1.variant.hotdeck8=="Omicron"]
group82 <- calendartimetoCOVID[dat$seq1.variant.hotdeck8!="Omicron"]
group91 <- calendartimetoCOVID[dat$seq1.variant.hotdeck9=="Omicron"]
group92 <- calendartimetoCOVID[dat$seq1.variant.hotdeck9!="Omicron"]
group101 <- calendartimetoCOVID[dat$seq1.variant.hotdeck10=="Omicron"]
group102 <- calendartimetoCOVID[dat$seq1.variant.hotdeck10!="Omicron"]

min(c(group12,group22,group32,group42,group52,group62,group72,group82,group92,group102))
max(c(group12,group22,group32,group42,group52,group62,group72,group82,group92,group102))
# min = 239 days = January 20, 2022, max = 274 days = February 24, 2022

# Range of Not Omicron/Delta imputations: January 20, 2022 to February 24, 2022

# Latest date of an actual Delta COVID-19: 
# Calculate the latest actual observed Delta COVID-19 onset:
mm <- max(calendartimetoCOVID[dat$seq1.variant!="MissingLineage" & dat$seq1.variant!="Omicron"])
# 274 days = February 24, 2022
dat$seq1.variant[calendartimetoCOVID==mm]
#[1] "Delta" "Delta"

min(c(group11,group21,group31,group41,group51,group61,group71,group81,group91,group101))
max(c(group11,group21,group31,group41,group51,group61,group71,group81,group91,group101))
# min = 234 days = January 15, 2022, max = 435 days = August 4, 2022

# Range of Not Omicron/Delta imputations: January 15, 2022 to August 4, 2022

# First enrollee Stage 2: 2021-10-19 

yy <- calendartimetoCOVID
xx1 <- ifelse(dat$seq1.variant.hotdeck1=="Omicron",1,0)
xx2 <- ifelse(dat$seq1.variant.hotdeck2=="Omicron",1,0)
xx3 <- ifelse(dat$seq1.variant.hotdeck3=="Omicron",1,0)
xx4 <- ifelse(dat$seq1.variant.hotdeck4=="Omicron",1,0)
xx5 <- ifelse(dat$seq1.variant.hotdeck5=="Omicron",1,0)
xx6 <- ifelse(dat$seq1.variant.hotdeck6=="Omicron",1,0)
xx7 <- ifelse(dat$seq1.variant.hotdeck7=="Omicron",1,0)
xx8 <- ifelse(dat$seq1.variant.hotdeck8=="Omicron",1,0)
xx9 <- ifelse(dat$seq1.variant.hotdeck9=="Omicron",1,0)
xx10 <- ifelse(dat$seq1.variant.hotdeck10=="Omicron",1,0)

table(xx1)
table(xx2)
table(xx3)
table(xx4)
table(xx5)
table(xx6)
table(xx7)
table(xx8)
table(xx9)
table(xx10)

# hotdeck1 is modal: 4 Delta, 73 Omicron
# hotdeck 2 has 2 Delta, 75 Omicron
# hotdeck 10 is 2 Delta, 75 Omicron
# pick hotdeck10

# Which hotdeck to pick?  Look at calendar dates
keep3 <- xx1==0
dat[keep3,"EventTimeFirstInfectionDate"]

keep3 <- xx2==0
dat[keep3,"EventTimeFirstInfectionDate"]

keep3 <- xx7=0
dat[keep3,"EventTimeFirstInfectionDate"]

keep3 <- xx1==1
sort(dat[keep3,"EventTimeFirstInfectionDate"])

# Observed Delta times
keep3 <- dat$seq1.variant=="Delta"
sort(dat[keep3,"EventTimeFirstInfectionDate"])
> sort(dat[keep3,"EventTimeFirstInfectionDate"])
[1] "2022-02-24" "2022-02-24"

# Numbers with known lineage
table(dat$seq1.variant)

         Delta MissingLineage        Omicron 
             2             30             45 

# So the L=5-nearest neighbors are always defined from the
# 45 known Omicron and 2 known Delta and filling in 30 cases with missing values.
# Sampling with replacement

# Average number omicron
mean(c(table(xx1)[2],table(xx2)[2],table(xx3)[2],table(xx4)[2],table(xx5)[2],
table(xx6)[2],table(xx7)[2],table(xx8)[2],table(xx9)[2],table(xx10)[2]))
74.5
# So average Not Omicron = 2.5
# Range number Omicron
range(c(table(xx1)[2],table(xx2)[2],table(xx3)[2],table(xx4)[2],table(xx5)[2],
table(xx6)[2],table(xx7)[2],table(xx8)[2],table(xx9)[2],table(xx10)[2]))
74 to 75

# Percentage agreement across all 45 pairs

vect <- c(mean(xx1==xx2),mean(xx1==xx3),mean(xx1==xx4),mean(xx1==xx5),mean(xx1==xx6),mean(xx1==xx7),
mean(xx1==xx8),mean(xx1==xx9),mean(xx1==xx10),mean(xx2==xx3),mean(xx2==xx4),mean(xx2==xx5),
mean(xx2==xx6),mean(xx2==xx7),mean(xx2==xx8),mean(xx2==xx9),mean(xx2==xx10),mean(xx3==xx4),
mean(xx3==xx5),mean(xx3==xx6),mean(xx3==xx7),mean(xx3==xx8),mean(xx3==xx9),mean(xx3==xx10),
mean(xx4==xx5),mean(xx4==xx6),mean(xx4==xx7),mean(xx4==xx8),mean(xx4==xx9),mean(xx4==xx10),
mean(xx5==xx6),mean(xx5==xx7),mean(xx5==xx8),mean(xx5==xx9),mean(xx5==xx10),mean(xx6==xx7),
mean(xx6==xx8),mean(xx6==xx9),mean(xx6==xx10),mean(xx7==xx8),mean(xx7==xx9),mean(xx7==xx10),
mean(xx8==xx9),mean(xx8==xx10),mean(xx9==xx10))

mean(vect)
0.992785 concordance in lineage calls averaging across all pairs of the 10 vectors (each of length 116) of lineage imputations

# Create side-by-side violin plots
pdf("T:/covpn/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/reports/TablesFigures/Stage2HotdeckcalendartimesFINAL.pdf")
par(mfrow=c(3,4),cex.axis=0.9,cex.lab=0.9,las=2,oma=c(3,3,6,3),mai = c(0.5, 0.1, 0.5, 0.1))
vioplot(yy ~ xx1, names = c("HD1 Omicron", "HD1 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx1, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
text(unique(group72),0.7,"*",col="black")
text(unique(group72),1.3,"*",col="black")
title("Hotdeck 1")

vioplot(yy ~ xx2, names = c("HD2 Omicron", "HD2 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx2, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
text(unique(group22),0.7,"*",col="black")
text(unique(group22),1.3,"*",col="black")
title("Hotdeck 2")

vioplot(yy ~ xx3, names = c("HD3 Omicron", "HD3 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx3, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 3")

vioplot(yy ~ xx4, names = c("HD4 Omicron", "HD4 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx4, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
text(unique(group72),0.7,"*",col="black")
text(unique(group72),1.3,"*",col="black")
title("Hotdeck 4")

vioplot(yy ~ xx5, names = c("HD5 Omicron", "HD5 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx5, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 5")

vioplot(yy ~ xx6, names = c("HD6 Omicron", "HD6 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx6, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 6")

vioplot(yy ~ xx7, names = c("HD7 Omicron", "HD7 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx7, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
text(unique(group72),0.7,"*",col="black")
text(unique(group72),1.3,"*",col="black")
title("Hotdeck 7")

vioplot(yy ~ xx8, names = c("HD8 Omicron", "HD8 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx8, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
text(unique(group72),0.7,"*",col="black")
text(unique(group72),1.3,"*",col="black")
title("Hotdeck 8")

vioplot(yy ~ xx9, names = c("HD9 Omicron", "HD9 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy ~ xx9, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("Hotdeck 9")

vioplot(yy ~ xx10, names = c("HD10 Omicron", "HD10 Not-Omicron"), 
col=c(1,2),axes=FALSE,ylim=c(128,435),horizontal=TRUE,xlab="",ylab="")
axis(1,at=c(128,159,189,220, 251, 279, 310, 340, 371, 401, 432),
labels=c("Oct 1 21","Nov 1 21","Dec 1 21", "Jan 1 22", 
"Feb 1 22", "March 1 22", "April 1 22", "May 1 22","June 1 22","July 1 22","Aug 1 22"))
stripchart(yy~xx10, method = "jitter", jitter=0.155, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
text(unique(group72),0.7,"*",col="black")
text(unique(group72),1.3,"*",col="black")
title("Hotdeck 10")

plot(1,1,type='n',axes=FALSE,xlab="",ylab="")
text(x=1,y=1.1, "Lineage Imputed:",cex=0.72)
legend(x="bottomright",c("Omicron", "Not Omicron"),fill=TRUE,col=c(2,1),
border=c(2,1),cex=0.722)

vioplot(vect,col=4,ylab="",xlab="",axes=FALSE,horizontal=TRUE,ylim=c(0.96,1.03),cex=0.9)
axis(1,at=c(0.97,0.98,0.99,1.0),labels=c("97%","98%","99%","100%"))
stripchart(vect, method = "jitter", jitter=0.19, col = "yellow",
           vertical = FALSE, pch = 19, add = TRUE,cex=0.10)
title("% Agreement")

mtext("Stage 2 Non-Naive: Number Days post FPI to COVID-19 by Imputed Omicron vs. Not Omicron",
outer=T,cex=0.95,line=2,side=3,las=1)

dev.off()

