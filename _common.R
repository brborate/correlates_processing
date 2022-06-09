library(kyotil)
library(methods)
library(dplyr)
library(digest)
set.seed(98109)
 

config <- config::get(config = Sys.getenv("TRIAL"))
for(opt in names(config)){
  eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
}
 
data_name = paste0(attr(config, "config"), "_data_processed_with_riskscore.csv")


# disabling lower level parallelization in favor of higher level of parallelization

# set parallelization in openBLAS and openMP
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1L)
stopifnot(blas_get_num_procs() == 1L)
omp_set_num_threads(1L)


verbose=Sys.getenv("VERBOSE")=="1"

names(assays)=assays # add names so that lapply results will have names

# if this flag is true, then the N IgG binding antibody is reported 
# in the immuno report (but is not analyzed in the cor or cop reports).
include_bindN <- !study_name %in% c("PREVENT19","AZD1222")


# For bAb, IU and BAU are the same thing
# all values on BAU or IU
# LOQ can not be NA, it is needed for computing delta
if(TRUE) {
    tmp=list(
        bindSpike=c(
            pos.cutoff=10.8424,
            LLOD = 0.3076,
            ULOD = 172226.2,
            LLOQ = 1.7968,
            ULOQ = 10155.95)
        ,
        bindRBD=c(
            pos.cutoff=14.0858,
            LLOD = 1.593648,
            ULOD = 223074,
            LLOQ = 3.4263,
            ULOQ = 16269.23)
        ,
        bindN=c( 
            pos.cutoff=23.4711,
            LLOD = 0.093744,
            ULOD = 52488,
            LLOQ = 4.4897,
            ULOQ = 574.6783)
        ,
        pseudoneutid50=c( 
            pos.cutoff=2.42,# as same lod
            LLOD = 2.42,
            ULOD = NA,
            LLOQ = 4.477,
            ULOQ = 10919)
        ,
        pseudoneutid80=c( 
            pos.cutoff=15.02,# as same lod
            LLOD = 15.02,
            ULOD = NA,
            LLOQ = 21.4786,
            ULOQ = 15368)
        ,
        liveneutmn50=c( 
            pos.cutoff=82.1*0.276,# as same lod
            LLOD = 82.11*0.276,
            ULOD = NA,
            LLOQ =  159.79*0.276,
            ULOQ = 11173.21*0.276)
    )
    
    pos.cutoffs=sapply(tmp, function(x) unname(x["pos.cutoff"]))
    llods=sapply(tmp, function(x) unname(x["LLOD"]))
    lloqs=sapply(tmp, function(x) unname(x["LLOQ"]))
    uloqs=sapply(tmp, function(x) unname(x["ULOQ"]))        
    
    if(study_name %in% c("COVE", "MockCOVE", "MockENSEMBLE")) {
        
        # nothing to do
        
    } else if(study_name=="ENSEMBLE") {
        
        if (contain(attr(config, "config"), "real")) {
        # EUA data
            
            # data less than pos cutoff is set to pos.cutoff/2
            llods["bindSpike"]=NA 
            lloqs["bindSpike"]=1.7968 
            uloqs["bindSpike"]=238.1165 
        
            # data less than pos cutoff is set to pos.cutoff/2
            llods["bindRBD"]=NA                 
            lloqs["bindRBD"]=3.4263                 
            uloqs["bindRBD"]=172.5755    
                    
            # data less than lloq is set to lloq/2
            llods["pseudoneutid50"]=NA  
            lloqs["pseudoneutid50"]=42*0.0653  #2.7426
            pos.cutoffs["pseudoneutid50"]=lloqs["pseudoneutid50"]
            uloqs["pseudoneutid50"]=9484*0.0653 # 619.3052
            
                # repeat for two synthetic markers that are adapted to SA and LA
                llods["pseudoneutid50sa"]=NA  
                lloqs["pseudoneutid50sa"]=42*0.0653  #2.7426
                pos.cutoffs["pseudoneutid50sa"]=lloqs["pseudoneutid50sa"]
                uloqs["pseudoneutid50sa"]=9484*0.0653 # 619.3052
        
                llods["pseudoneutid50la"]=NA  
                lloqs["pseudoneutid50la"]=42*0.0653  #2.7426
                pos.cutoffs["pseudoneutid50la"]=lloqs["pseudoneutid50la"]
                uloqs["pseudoneutid50la"]=9484*0.0653 # 619.3052
    
            # data less than lod is set to lod/2
            llods["ADCP"]=11.57
            lloqs["ADCP"]=8.87
            pos.cutoffs["ADCP"]=11.57# as same lod
            uloqs["ADCP"]=211.56
            
        } else if (contain(attr(config, "config"), "partA")) {
        # complete part A data
            
            # data less than pos cutoff is set to pos.cutoff/2
            llods["bindSpike"]=NA 
            lloqs["bindSpike"]=1.7968 
            uloqs["bindSpike"]=238.1165 
        
            # data less than pos cutoff is set to pos.cutoff/2
            llods["bindRBD"]=NA                 
            lloqs["bindRBD"]=3.4263                 
            uloqs["bindRBD"]=172.5755    
                    
            # data less than lloq is set to lloq/2
            llods["pseudoneutid50"]=NA  
            lloqs["pseudoneutid50"]=75*0.0653  #4.8975
            pos.cutoffs["pseudoneutid50"]=lloqs["pseudoneutid50"]
            uloqs["pseudoneutid50"]=12936*0.0653 # 844.7208
    
            # data less than lod is set to lod/2
            llods["ADCP"]=11.57
            lloqs["ADCP"]=8.87
            pos.cutoffs["ADCP"]=11.57# as same lod
            uloqs["ADCP"]=211.56
        }
        
    } else if(study_name=="PREVENT19") {
        # Novavax
        
        # data less than lloq is set to lloq/2 in the raw data
        llods["bindSpike"]=NA 
        lloqs["bindSpike"]=150.4*0.0090 # 1.3536
        uloqs["bindSpike"]=770464.6*0.0090 # 6934.181
        pos.cutoffs["bindSpike"]=10.8424 # use same as COVE
    
        # data less than lod is set to lod/2 in the raw data
        llods["pseudoneutid50"]=2.612 # 40 * 0.0653
        lloqs["pseudoneutid50"]=51*0.0653 # 3.3303
        uloqs["pseudoneutid50"]=127411*0.0653 # 8319.938
        pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
        
    } else if(study_name=="AZD1222") {
           
        # data less than lloq is set to lloq/2 in the raw data, Nexelis
        llods["bindSpike"]=NA 
        lloqs["bindSpike"]=62.8*0.0090 # 0.5652
        uloqs["bindSpike"]=238528.4*0.0090 # 2146.756
        pos.cutoffs["bindSpike"]=10.8424 # use same as COVE
    
        # data less than lod is set to lod/2
        llods["pseudoneutid50"]=2.612  
        lloqs["pseudoneutid50"]=56*0.0653 # 3.6568
        uloqs["pseudoneutid50"]=47806*0.0653 # 3121.732
        pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
        
    } else if(study_name=="VAT08m") {
        # Sanofi
           
        # data less than lod is set to lod/2
        llods["pseudoneutid50"]=2.612  
        lloqs["pseudoneutid50"]=95*0.0653 # 3.6568
        uloqs["pseudoneutid50"]=191429*0.0653 # 3121.732
        pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
        
    } else stop("unknown study_name 1")
    
    # llox is for plotting and can be either llod or lloq depending on trials
    lloxs=ifelse(config$llox_label=="LOD", llods[names(config$llox_label)], lloqs[names(config$llox_label)])
    
}


#assays_to_be_censored_at_uloq_cor <- c(
#  "bindSpike", "bindRBD", "pseudoneutid50", "pseudoneutid80"
#  # NOTE: the live neutralization marker will eventually be available
#  #"liveneutmn50"
#)

###############################################################################
# figure labels and titles for markers
###############################################################################

markers <- c(outer(times[which(times %in% c("B", paste0("Day", config$timepoints)))], 
                   assays, "%.%"))

# race labeling
labels.race <- c(
  "White", 
  "Black or African American",
  "Asian", 
  if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & startsWith(attr(config, "config"),"janssen_la")) "Indigenous South American" else "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander", 
  "Multiracial",
  if ((study_name=="COVE" | study_name=="MockCOVE")) "Other", 
  "Not reported and unknown"
)

# ethnicity labeling
labels.ethnicity <- c(
  "Hispanic or Latino", "Not Hispanic or Latino",
  "Not reported and unknown"
)

labels.assays.short <- c("Anti N IgG (BAU/ml)", 
                         "Anti Spike IgG (BAU/ml)", 
                         "Anti RBD IgG (BAU/ml)", 
                         "Pseudovirus-nAb cID50", 
                         "Pseudovirus-nAb cID80", 
                         "Live virus-nAb cMN50")
names(labels.assays.short) <- c("bindN",
  "bindSpike",
  "bindRBD",
  "pseudoneutid50",
  "pseudoneutid80",
  "liveneutmn50")

# hacky fix for tabular, since unclear who else is using
# the truncated labels.assays.short later
labels.assays.short.tabular <- labels.assays.short

labels.time <- c("Day 1", paste0("Day ", config$timepoints), paste0("D", config$timepoints, " fold-rise over D1"), "D57 fold-rise over D29")

names(labels.time) <- c("B", paste0("Day", config$timepoints), paste0("Delta", config$timepoints, "overB"), "Delta57over29")

# axis labeling
labels.axis <- outer(
  rep("", length(times)),
  labels.assays.short[assays],
  "%.%"
)
labels.axis <- as.data.frame(labels.axis)
rownames(labels.axis) <- times

labels.assays <- c("Binding Antibody to Spike", 
                   "Binding Antibody to RBD",
                   "PsV Neutralization 50% Titer",
                   "PsV Neutralization 80% Titer",
                   "WT LV Neutralization 50% Titer")

names(labels.assays) <- c("bindSpike", 
                          "bindRBD", 
                          "pseudoneutid50",
                          "pseudoneutid80",
                          "liveneutmn50")

# title labeling
labels.title <- outer(
  labels.assays[assays],
  ": " %.% c("Day 1", paste0("Day ", config$timepoints), paste0("D", config$timepoints, " fold-rise over D1"), "D57 fold-rise over D29"),
  paste0
)
labels.title <- as.data.frame(labels.title)
colnames(labels.title) <- times
# NOTE: hacky solution to deal with changes in the number of markers
rownames(labels.title)[seq_along(assays)] <- assays
labels.title <- as.data.frame(t(labels.title))

# creating short and long labels
labels.assays.short <- labels.axis[1, ]
labels.assays.long <- labels.title

# baseline stratum labeling
if (study_name=="COVE" | study_name=="MockCOVE") {
    Bstratum.labels <- c(
      "Age >= 65",
      "Age < 65, At risk",
      "Age < 65, Not at risk"
    )
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    Bstratum.labels <- c(
      "Age < 60, Not at risk",
      "Age < 60, At risk",
      "Age >= 60, Not at risk",
      "Age >= 60, At risk"
    )
    
} else if (study_name %in% c("PREVENT19","AZD1222")) {
    Bstratum.labels <- c(
      "Age >= 65",
      "Age < 65"
    )

} else if (study_name %in% c("VAT08m","VAT08m")) {
    Bstratum.labels <- c(
      "Age >= 60",
      "Age < 60"
    )

} else if (study_name=="HVTN705") {
    # do nothing

} else stop("unknown study_name 2")



# baseline stratum labeling
if (study_name=="COVE" | study_name=="MockCOVE") {
    demo.stratum.labels <- c(
      "Age >= 65, URM",
      "Age < 65, At risk, URM",
      "Age < 65, Not at risk, URM",
      "Age >= 65, White non-Hisp",
      "Age < 65, At risk, White non-Hisp",
      "Age < 65, Not at risk, White non-Hisp"
    )
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    demo.stratum.labels <- c(
      "US URM, Age 18-59, Not at risk",
      "US URM, Age 18-59, At risk",
      "US URM, Age >= 60, Not at risk",
      "US URM, Age >= 60, At risk",
      "US White non-Hisp, Age 18-59, Not at risk",
      "US White non-Hisp, Age 18-59, At risk",
      "US White non-Hisp, Age >= 60, Not at risk",
      "US White non-Hisp, Age >= 60, At risk",
      "Latin America, Age 18-59, Not at risk",
      "Latin America, Age 18-59, At risk",
      "Latin America, Age >= 60, Not at risk",
      "Latin America, Age >= 60, At risk",
      "South Africa, Age 18-59, Not at risk",
      "South Africa, Age 18-59, At risk",
      "South Africa, Age >= 60, Not at risk",
      "South Africa, Age >= 60, At risk"
    )
    
} else if (study_name=="PREVENT19") {
    demo.stratum.labels <- c(
      "US White non-Hisp, Age 18-64, Not at risk",
      "US White non-Hisp, Age 18-64, At risk",
      "US White non-Hisp, Age >= 65, Not at risk",
      "US White non-Hisp, Age >= 65, At risk",
      "US URM, Age 18-64, Not at risk",
      "US URM, Age 18-64, At risk",
      "US URM, Age >= 65, Not at risk",
      "US URM, Age >= 65, At risk",
      "Mexico, Age 18-64",
      "Mexico, Age >= 65"
    )

} else if (study_name=="AZD1222") {
    demo.stratum.labels <- c(
      "US White non-Hisp, Age 18-64",
      "US White non-Hisp, Age >= 65",
      "US URM, Age 18-64",
      "US URM, Age >= 65",
      "Non-US, Age 18-64",
      "Non-US, Age >= 65"
    )

} else if (study_name=="VAT08m") {
    demo.stratum.labels <- c(
      "US URM, Age 18-59, Not at risk",
      "US URM, Age 18-59, At risk",
      "US URM, Age >= 60, Not at risk",
      "US URM, Age >= 60, At risk",
      "US White non-Hisp, Age 18-59, Not at risk",
      "US White non-Hisp, Age 18-59, At risk",
      "US White non-Hisp, Age >= 60, Not at risk",
      "US White non-Hisp, Age >= 60, At risk",
      "Japan, Age 18-59, Not at risk",
      "Japan, Age 18-59, At risk",
      "Japan, Age >= 60, Not at risk",
      "Japan, Age >= 60, At risk"
    )

} else if (study_name=="HVTN705") {
    # do nothing

} else stop("unknown study_name 3")


labels.regions.ENSEMBLE =c("0"="Northern America", "1"="Latin America", "2"="Southern Africa")
regions.ENSEMBLE=0:2
names(regions.ENSEMBLE)=labels.regions.ENSEMBLE

labels.countries.ENSEMBLE=c("0"="United States", "1"="Argentina", "2"="Brazil", "3"="Chile", "4"="Columbia", "5"="Mexico", "6"="Peru", "7"="South Africa")
countries.ENSEMBLE=0:7
names(countries.ENSEMBLE)=labels.countries.ENSEMBLE


###############################################################################
# reproduciblity options
###############################################################################

# NOTE: used in appendix.Rmd to store digests of input raw/processed data files
# hash algorithm picked based on https://csrc.nist.gov/projects/hash-functions
hash_algorithm <- "sha256"


###############################################################################
# theme options
###############################################################################

# fixed knitr chunk options
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  out.width = "80%",
  out.extra = "",
  fig.pos = "H",
  fig.show = "hold",
  fig.align = "center",
  fig.width = 6,
  fig.asp = 0.618,
  fig.retina = 0.8,
  dpi = 600,
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)

# global options
options(
  digits = 6,
  #scipen = 999,
  dplyr.print_min = 6,
  dplyr.print_max = 6,
  crayon.enabled = FALSE,
  bookdown.clean_book = TRUE,
  knitr.kable.NA = "NA",
  repos = structure(c(CRAN = "https://cran.rstudio.com/"))
)

# no complaints from installation warnings
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# overwrite options by output type
if (knitr:::is_html_output()) {
  #options(width = 80)

  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}
if (knitr:::is_latex_output()) {
  #knitr::opts_chunk$set(width = 67)
  #options(width = 67)
  options(cli.unicode = TRUE)

  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}

# create and set global ggplot theme
# borrowed from https://github.com/tidymodels/TMwR/blob/master/_common.R
theme_transparent <- function(...) {
  # use black-white theme as base
  ret <- ggplot2::theme_bw(...)

  # modify with transparencies
  trans_rect <- ggplot2::element_rect(fill = "transparent", colour = NA)
  ret$panel.background  <- trans_rect
  ret$plot.background   <- trans_rect
  ret$legend.background <- trans_rect
  ret$legend.key        <- trans_rect

  # always have legend below
  ret$legend.position <- "bottom"
  return(ret)
}

library(ggplot2)
theme_set(theme_transparent())
theme_update(
  text = element_text(size = 25),
  axis.text.x = element_text(colour = "black", size = 30),
  axis.text.y = element_text(colour = "black", size = 30)
)

# custom ggsave function with updated defaults
ggsave_custom <- function(filename = default_name(plot),
                          height= 15, width = 21, ...) {
  ggsave(filename = filename, height = height, width = width, ...)
}




# extract assay from marker name such as Daytp1pseudoneutid80, Bpseudoneutid80
marker.name.to.assay=function(marker.name) {
    if(endsWith(marker.name, "bindSpike")) {
        "bindSpike"
    } else if(endsWith(marker.name, "bindRBD")) {
        "bindRBD"
    } else if(endsWith(marker.name, "bindN")) {
        "bindN"
    } else if(endsWith(marker.name, "pseudoneutid50")) {
        "pseudoneutid50"
    } else if(endsWith(marker.name, "pseudoneutid80")) {
        "pseudoneutid80"
    } else if(endsWith(marker.name, "liveneutmn50")) {
        "liveneutmn50"
    } else stop("marker.name.to.assay: wrong marker.name")
}


# x is the marker values
# assay is one of assays, e.g. pseudoneutid80
report.assay.values=function(x, assay){
    lars.quantiles=seq(0,1,length.out=30) [round(seq.int(1, 30, length.out = 10))]
    sens.quantiles=c(0.15, 0.85)
    # cannot have different lengths for different assays, otherwise downstream code may break
    fixed.values = log10(c("500"=500, "1000"=1000))#, "llod/2"=unname(llods[assay]/2))) # llod/2 may not be in the observed values
    out=sort(c(quantile(x, c(lars.quantiles,sens.quantiles), na.rm=TRUE), fixed.values))    
    out
    #out[!duplicated(out)] # unique strips away the names. But don't take out duplicates because 15% may be needed and because we may want the same number of values for each assay
}
#report.assay.values (dat.vac.seroneg[["Daytp1pseudoneutid80"]], "pseudoneutid80")


preprocess.for.risk.score=function(dat_raw, study_name) {
    dat_proc=dat_raw
    
    dat_proc=subset(dat_proc, !is.na(Bserostatus))
    
    # EventTimePrimaryIncludeNotMolecConfirmedD29 are the endpoint of interest and should be used to compute weights
    if(study_name=="ENSEMBLE") {
        dat_proc$EventTimePrimaryD29=dat_proc$EventTimePrimaryIncludeNotMolecConfirmedD29
        dat_proc$EventIndPrimaryD29 =dat_proc$EventIndPrimaryIncludeNotMolecConfirmedD29
        dat_proc$EventTimePrimaryD1 =dat_proc$EventTimePrimaryIncludeNotMolecConfirmedD1
        dat_proc$EventIndPrimaryD1  =dat_proc$EventIndPrimaryIncludeNotMolecConfirmedD1
    }
        
    for(tp in timepoints) dat_proc=dat_proc[!is.na(dat_proc[["EventTimePrimaryD"%.%tp]]), ]
    
    
    for(tp in timepoints) {
        dat_proc[["EarlyendpointD"%.%tp]] <- with(dat_proc, ifelse(get("EarlyinfectionD"%.%tp)==1 | (EventIndPrimaryD1==1 & EventTimePrimaryD1 < get("NumberdaysD1toD"%.%tp) + 7),1,0))
    }

    # ENSEMBLE only, since we are not using this variable to define Riskscorecohortflag and we are not doing D29start1 analyses for other trials
    if (study_name %in% c("MockENSEMBLE", "ENSEMBLE")) {
        dat_proc[["EarlyendpointD29start1"]]<- with(dat_proc, ifelse(get("EarlyinfectionD29start1")==1| (EventIndPrimaryD1==1 & EventTimePrimaryD1 < get("NumberdaysD1toD29") + 1),1,0))
    } else dat_proc$EarlyinfectionD29start1=dat_proc$EarlyinfectionD29 # this is not necessary, but it is kept here to make the hash checks for mock datasets happy
    
    
    # Indicator of membership in the cohort included in the analysis that defines the risk score in the placebo arm. It requires:
    # 1. baseline SARS-CoV-2 negative, 
    # 2. per-protocol, 
    # 3. no evidence of SARS-CoV-2 infection or right-censoring up to time point tinterm (2 dose) or tpeak (1 dose)
    # 4. lack of missing data on a certain set of baseline input variables (not enfored here because the developer of this script need not have knowledge of risk score requirements)
    # no NAs allowed. 
    if (study_name %in% c("MockCOVE", "COVE")) {
        # special case, redefined for backward compatibility
        dat_proc$Riskscorecohortflag <- with(dat_proc, ifelse(Bserostatus==0 & Perprotocol==1, 1, 0))

    } else if (study_name %in% c("ENSEMBLE", "MockENSEMBLE")){
        dat_proc$Riskscorecohortflag <-
          with(dat_proc, ifelse(Bserostatus==0 & Perprotocol==1 & get("EarlyendpointD"%.%timepoints[1]%.%"start1")==0 & get("EventTimePrimaryD"%.%timepoints[1])>=1, 1, 0))

    } else if (study_name == "PREVENT19") {
      dat_proc <- dat_proc %>%
        mutate(Riskscorecohortflag = ifelse(Bserostatus==0 & Perprotocol==1, 1, 0),
               RiskscoreAUCflag = case_when(Trt==0 & Bserostatus==0 & Perprotocol==1 ~ 1,
                                            Trt==1 & Bserostatus==0 & Perprotocol==1 & EarlyendpointD35==0 & EventTimePrimaryD35>=7 ~ 1,
                                            TRUE ~ 0))
    } else if (study_name == "AZD1222") {
        dat_proc <- dat_proc %>%
          mutate(Riskscorecohortflag = ifelse(Bserostatus==0 & Perprotocol==1, 1, 0),
                 RiskscoreAUCflag = case_when(Trt==0 & Bserostatus==0 & Perprotocol==1 ~ 1,
                                              Trt==1 & Bserostatus==0 & Perprotocol==1 & EarlyendpointD57==0 & EventTimePrimaryD57>=7 ~ 1,
                                              TRUE ~ 0))
    } else if (study_name == "VAT08m") {
        dat_proc <- dat_proc %>%
          mutate(Riskscorecohortflag = ifelse(Bserostatus==0 & Perprotocol==1, 1, 0),
                 RiskscoreAUCflag = case_when(Trt==0 & Bserostatus==0 & Perprotocol==1 ~ 1,
                                              Trt==1 & Bserostatus==0 & Perprotocol==1 & EarlyendpointD43==0 & EventTimePrimaryD43>=7 ~ 1,
                                              TRUE ~ 0))
    } else stop("unknown study_name 4")

    assertthat::assert_that(
        all(!is.na(dat_proc$Riskscorecohortflag)),
        msg = "missing Riskscorecohortflag")
    
    dat_proc    
}
