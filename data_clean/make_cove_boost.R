library(here)
renv::activate(here::here())


config <- config::get(config = Sys.getenv("TRIAL"))
for(opt in names(config)){
  eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
}

data_name = paste0(attr(config, "config"), "_data_processed_with_riskscore.csv")


# read stage 1analysis ready dataset 
dat_stage1 = read.csv("/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_with_riskscore.csv")

# read stage 2 mapped data 
dat_raw = read.csv(mapped_data)
if (colnames(dat_raw)[1]=="Subjectid")  colnames(dat_raw)[1] <- "Ptid" else stop("the first column is unexpectedly not Subjectid")

setdiff(names(dat_raw), names(dat_stage1))

# merge two files
dat_stage2 = merge(dat_stage1, dat_raw, by="Ptid", all=T, suffixes=c("",".y"))

# remove columns ending in .y
dat_stage2 = dat_stage2[,!endsWith(names(dat_stage2),".y")]


#### define sampling stratum variables

# demo.stratum from stage 1 study can be reused, but tps.stratum needs to be redefined 

# 0-15
dat_stage2$sampling_bucket = with(dat_stage2, 
    strtoi(paste0(
      strtoi(paste0( 
        Trt, 
        Bserostatus 
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


#### impute assay values



#### inverse probability sampling weights



# save
write.csv(dat_stage2, file = here("data_clean", data_name), row.names=F)
