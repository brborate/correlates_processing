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

# save
write.csv(dat_stage2, file = here("data_clean", data_name), row.names=F)
