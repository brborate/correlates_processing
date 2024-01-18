# Add input variable definition column and comments column to dataset
# @param dataframe containing Variable Name of input variable used in risk score analysis
# @return dataframe with two new columns: Definition and Comments
get_defs_comments_riskVars <- function(data){
  if(study_name == "COVE" | study_name == "MockCOVE"){
    data <- data %>%
      mutate(Definition = case_when(
        `Variable Name` == "Age" ~ "Age at enrollment in years, between 18 and 85",
        `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male)",
        `Variable Name` == "BMI" ~ "BMI at enrollment (kg/m^2)",
        `Variable Name` == "MinorityInd" ~ "Baseline covariate underrepresented minority status (1=minority, 0=non-minority)",
        `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity = Hispanic (0 = Non-Hispanic)",
        `Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (0 = Non-Hispanic)",
        `Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (0 = Non-Hispanic)",
        `Variable Name` == "Black" ~ "Indicator race = Black (0 = White)",
        `Variable Name` == "Asian" ~ "Indicator race = Asian (0 = White)",
        `Variable Name` == "NatAmer" ~ "Indicator race = American Indian or Alaska Native (0 = White)",
        `Variable Name` == "PacIsl" ~ "Indicator race = Native Hawaiian or Other Pacific Islander (0 = White)",
        `Variable Name` == "WhiteNonHispanic" ~ "Indicator race = White or Caucasian (1 = White)",
        `Variable Name` == "Multiracial" ~ "Indicator race = Multiracial (0 = White)",
        `Variable Name` == "Other" ~ "Indicator race = Other (0 = White)",
        `Variable Name` == "Notreported" ~ "Indicator race = Not reported (0 = White)",
        `Variable Name` == "Unknown" ~ "Indicator race = unknown (0 = White)",
        `Variable Name` == "HighRiskInd" ~ "Baseline covariate high risk pre-existing condition (1=yes, 0=no)"
      ),
      Comments = "")
  }
  
  if(study_name == "ENSEMBLE" | study_name == "MockENSEMBLE"){
    data <- data %>%
      mutate(Definition = case_when(
        `Variable Name` == "Age" ~ "Age at enrollment in years (integer >= 18, NA=missing). Note that the randomization strata included Age 18-59 vs. Age >= 60.",
        `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male/undifferentiated/unknown",
        `Variable Name` == "BMI" ~ "BMI at enrollment (Ordered categorical 1,2, 3, 4, NA=missing); 1 = Underweight BMI < 18.5; 2 = Normal BMI 18.5 to < 25; 3 = Overweight BMI 25 to < 30; 4 = Obese BMI >= 30",
        `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity = Hispanic (1 = Hispanic, 0 = complement)",
        `Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (1 = Not reported, 0 = complement)",
        `Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (1 = Unknown, 0 = complement)",
        `Variable Name` == "Black" ~ "Indicator race = Black (1=Black, 0=complement)",
        `Variable Name` == "Asian" ~ "Indicator race = Asian (1=Asian, 0=complement)",
        `Variable Name` == "NatAmer" ~ "Indicator race = American Indian or Alaska Native (1=NatAmer, 0=complement)",
        `Variable Name` == "PacIsl" ~ "Indicator race = Native Hawaiian or Other Pacific Islander (1=PacIsl, 0=complement)",
        `Variable Name` == "Multiracial" ~ "Indicator race = Multiracial (1=Multiracial, 0=complement)",
        `Variable Name` == "Notreported" ~ "Indicator race = Not reported (1=Notreported, 0=complement)",
        `Variable Name` == "Unknown" ~ "Indicator race = unknown (1=Unknown, 0=complement)",
        `Variable Name` == "URMforsubcohortsampling" ~ "Indicator of under-represented minority (1=Yes, 0=No)",
        `Variable Name` == "HighRiskInd" ~ "Baseline covariate indicating >= 1 Co-existing conditions (1=yes, 0=no, NA=missing)",
        #`Variable Name` == "Country" ~ "Country of the study site of enrollment (0=United States, 1=Argentina,2=Brazil, 3=Chile,4=Columbia, 5=Mexico, 6=Peru, 7=South Africa)",
        #`Variable Name` == "Region" ~ "Major geographic region of the study site of enrollment (0=Northern America, 1=Latin America, 2=Southern Africa).",
        `Variable Name` == "HIVinfection" ~ "Indicator HIV infected at enrollment (1=infected, 0=not infected)",
        #`Variable Name` == "CalendarDateEnrollment" ~ "Date variable (used to control for calendar time trends in COVID incidence). Coded as number of days since first person enrolled until the ppt is enrolled.",
        `Variable Name` == "Country.X1" ~ "Indicator country = Argentina (1 = Argentina, 0 = complement)",
        `Variable Name` == "Country.X2" ~ "Indicator country = Brazil (1 = Brazil, 0 = complement)",
        `Variable Name` == "Country.X3" ~ "Indicator country = Chile (1 = Chile, 0 = complement)",
        `Variable Name` == "Country.X4" ~ "Indicator country = Columbia (1 = Columbia, 0 = complement)",
        `Variable Name` == "Country.X5" ~ "Indicator country = Mexico (1 = Mexico, 0 = complement)",
        `Variable Name` == "Country.X6" ~ "Indicator country = Peru (1 = Peru, 0 = complement)",
        `Variable Name` == "Country.X7" ~ "Indicator country = South Africa (1 = South Africa, 0 = complement)",
        `Variable Name` == "Region.X1" ~ "Indicator region = Latin America (1 = Latin America, 0 = complement)",
        `Variable Name` == "Region.X2" ~ "Indicator country = Southern Africa (1 = Southern Africa, 0 = complement)",
        `Variable Name` == "CalDtEnrollIND.X1" ~ "Indicator variable representing enrollment occurring between 4-8 weeks periods of first subject enrolled (1 = Enrollment between 4-8 weeks, 0 = complement).",
        `Variable Name` == "CalDtEnrollIND.X2" ~ "Indicator variable representing enrollment occurring between 8-12 weeks periods of first subject enrolled (1 = Enrollment between 8-12 weeks, 0 = complement).",
        `Variable Name` == "CalDtEnrollIND.X3" ~ "Indicator variable representing enrollment occurring between 12-16 weeks periods of first subject enrolled (1 = Enrollment between 12-16 weeks, 0 = complement).",
        `Variable Name` == "Region.X1.x.CalDtEnrollIND.X1" ~ "Interaction term between Region.X1 and CalDtEnrollIND.X1",
        `Variable Name` == "Region.X1.x.CalDtEnrollIND.X2" ~ "Interaction term between Region.X1 and CalDtEnrollIND.X2",
        `Variable Name` == "Region.X1.x.CalDtEnrollIND.X3" ~ "Interaction term between Region.X1 and CalDtEnrollIND.X3",
        `Variable Name` == "Region.X2.x.CalDtEnrollIND.X1" ~ "Interaction term between Region.X2 and CalDtEnrollIND.X1",
        `Variable Name` == "Region.X2.x.CalDtEnrollIND.X2" ~ "Interaction term between Region.X2 and CalDtEnrollIND.X2",
        `Variable Name` == "Region.X2.x.CalDtEnrollIND.X3" ~ "Interaction term between Region.X2 and CalDtEnrollIND.X3"
      ),
      Comments = "")
  }
  if(study_name == "PREVENT19"){
    data <- data %>%
      mutate(Definition = case_when(
        `Variable Name` == "Age" ~ "Age at enrollment in years",
        `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male)",
        `Variable Name` == "BMI" ~ "BMI at enrollment (kg/m^2)",
        #`Variable Name` == "MinorityInd" ~ "Baseline covariate underrepresented minority status (1=minority, 0=non-minority)",
        `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity = Hispanic or Latino (0 = Non-Hispanic/Non-Latino)",
        #`Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (0 = Non-Hispanic)",
        #`Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (0 = Non-Hispanic)",
        `Variable Name` == "Black" ~ "Indicator race = Black (0 = White)",
        `Variable Name` == "Asian" ~ "Indicator race = Asian (0 = White)",
        #`Variable Name` == "NatAmer" ~ "Indicator race = American Indian or Alaska Native (0 = White)",
        #`Variable Name` == "PacIsl" ~ "Indicator race = Native Hawaiian or Other Pacific Islander (0 = White)",
        #`Variable Name` == "WhiteNonHispanic" ~ "Indicator race = White or Caucasian (1 = White)",
        #`Variable Name` == "Multiracial" ~ "Indicator race = Multiracial (0 = White)",
        #`Variable Name` == "Other" ~ "Indicator race = Other (0 = White)",
        #`Variable Name` == "Notreported" ~ "Indicator race = Not reported (0 = White)",
        #`Variable Name` == "Unknown" ~ "Indicator race = unknown (0 = White)",
        `Variable Name` == "HighRiskInd" ~ "Baseline covariate high risk pre-existing condition (1=yes, 0=no)",
        `Variable Name` == "Height" ~ "Height at baseline (cm)",
        `Variable Name` == "Weight" ~ "Weight at baseline (kg)"
      ),
      Comments = "")
  }
  
  if(study_name == "AZD1222"){
    data <- data %>%
      mutate(Definition = case_when(
        `Variable Name` == "Age" ~ "Age at enrollment in years",
        `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male)",
        `Variable Name` == "BMI" ~ "BMI at enrollment (kg/m^2)",
        #`Variable Name` == "MinorityInd" ~ "Baseline covariate underrepresented minority status (1=minority, 0=non-minority)",
        `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity = Hispanic or Latino (0 = Non-Hispanic/Non-Latino)",
        #`Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (0 = Non-Hispanic)",
        #`Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (0 = Non-Hispanic)",
        `Variable Name` == "Black" ~ "Indicator race = Black (0 = White)",
        #`Variable Name` == "Asian" ~ "Indicator race = Asian (0 = White)",
        `Variable Name` == "NatAmer" ~ "Indicator race = American Indian or Alaska Native (0 = White)",
        #`Variable Name` == "PacIsl" ~ "Indicator race = Native Hawaiian or Other Pacific Islander (0 = White)",
        #`Variable Name` == "WhiteNonHispanic" ~ "Indicator race = White or Caucasian (1 = White)",
        `Variable Name` == "Multiracial" ~ "Indicator race = Multiracial (0 = White)",
        #`Variable Name` == "Other" ~ "Indicator race = Other (0 = White)",
        #`Variable Name` == "Notreported" ~ "Indicator race = Not reported (0 = White)",
        #`Variable Name` == "Unknown" ~ "Indicator race = unknown (0 = White)",
        #`Variable Name` == "HighRiskInd" ~ "Baseline covariate high risk pre-existing condition (1=yes, 0=no)"
        `Variable Name` == "Country.X1" ~ "Indicator country = Peru? ",
        `Variable Name` == "Country.X2" ~ "Indicator country = USA"
      ),
      Comments = "")
  }
  
  if(study_name %in% c("VAT08m", "VAT08")){
    data <- data %>%
      mutate(Definition = case_when(
        `Variable Name` == "Age" ~ "Age at enrollment in years",
        `Variable Name` == "pooled.age.grp" ~ "Pooled age group",
        `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male)",
        `Variable Name` == "BMI" ~ "BMI at enrollment (kg/m^2)",
        #`Variable Name` == "MinorityInd" ~ "Baseline covariate underrepresented minority status (1=minority, 0=non-minority)",
        `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity = Hispanic or Latino (0 = Non-Hispanic/Non-Latino)",
        #`Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (0 = Non-Hispanic)",
        #`Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (0 = Non-Hispanic)",
        `Variable Name` == "Black" ~ "Indicator race = Black (0 = White)",
        `Variable Name` == "Asian" ~ "Indicator race = Asian (0 = White)",
        `Variable Name` == "NatAmer" ~ "Indicator race = American Indian or Alaska Native (0 = White)",
        #`Variable Name` == "PacIsl" ~ "Indicator race = Native Hawaiian or Other Pacific Islander (0 = White)",
        #`Variable Name` == "WhiteNonHispanic" ~ "Indicator race = White or Caucasian (1 = White)",
        #`Variable Name` == "Multiracial" ~ "Indicator race = Multiracial (0 = White)",
        #`Variable Name` == "Other" ~ "Indicator race = Other (0 = White)",
        #`Variable Name` == "Notreported" ~ "Indicator race = Not reported (0 = White)",
        #`Variable Name` == "Unknown" ~ "Indicator race = unknown (0 = White)",
        `Variable Name` == "URMforsubcohortsampling" ~ "URMforsubcohortsampling = ",
        `Variable Name` == "HighRiskInd" ~ "Baseline covariate high risk pre-existing condition (1=yes, 0=no)",
        #`Variable Name` == "HIVinfection" ~ "HIV infection (1=yes, 0=no)"
        #`Variable Name` == "Country.X1" ~ "Indicator country =  ",
        `Variable Name` == "Country.X2" ~ "Indicator country = ",
        `Variable Name` == "Country.X3" ~ "Indicator country = ",
        `Variable Name` == "Country.X4" ~ "Indicator country = ",
        `Variable Name` == "Country.X5" ~ "Indicator country = ",
        #`Variable Name` == "Country.X6" ~ "Indicator country = ",
        `Variable Name` == "Country.X7" ~ "Indicator country = ",
        `Variable Name` == "CalDtEnrollIND.X1" ~ "Indicator Calendar Date Enrollment = ",
        `Variable Name` == "CalDtEnrollIND.X2" ~ "Indicator Calendar Date Enrollment = ",
        `Variable Name` == "CalDtEnrollIND.X3" ~ "Indicator Calendar Date Enrollment = ",
        `Variable Name` == "CalDtEnrollIND.X4" ~ "Indicator Calendar Date Enrollment = "
        #`Variable Name` == "CalDtEnrollIND.X5" ~ "Indicator Calendar Date Enrollment = "
        #`Variable Name` == "CalDtEnrollIND.X6" ~ "Indicator Calendar Date Enrollment = "
      ),
      Comments = "")
  }
  
  if(Sys.getenv("TRIAL") == "covail"){
    data <- data %>%
      mutate(Definition = case_when(
        `Variable Name` == "Age" ~ "Age at enrollment in years",
        `Variable Name` == "Age65C" ~ "Binary indicator for Age >= 65 (1 if Age >=65, 0 otherwise)",
        `Variable Name` == "Sex" ~ "Sex assigned at birth (1=female, 0=male)",
        #`Variable Name` == "BMI" ~ "BMI at enrollment (kg/m^2)",
        #`Variable Name` == "MinorityInd" ~ "Baseline covariate underrepresented minority status (1=minority, 0=non-minority)",
        `Variable Name` == "EthnicityHispanic" ~ "Indicator ethnicity is Hispanic or Latino (1 if Hispanic or Latino, 0 otherwise)",
        #`Variable Name` == "EthnicityNotreported" ~ "Indicator ethnicity = Not reported (0 = Non-Hispanic)",
        #`Variable Name` == "EthnicityUnknown" ~ "Indicator ethnicity = Unknown (0 = Non-Hispanic)",
        `Variable Name` == "Black" ~ "Indicator race is Black (1 if Black, 0 otherwise)",
        `Variable Name` == "Asian" ~ "Indicator race is Asian (1 if Asian, 0 otherwise)",
        `Variable Name` == "NatAmer" ~ "Indicator race is American Indian or Alaska Native (1 if NatAmer, 0 otherwise)",
        `Variable Name` == "PacIsl" ~ "Indicator race is Native Hawaiian or Other Pacific Islander (1 if PacIsl, 0 otherwise)",
        #`Variable Name` == "WhiteNonHispanic" ~ "Indicator race = White or Caucasian (1 = White)",
        `Variable Name` == "Multiracial" ~ "Indicator race is Multiracial (1 if Multiracial, 0 otherwise)",
        #`Variable Name` == "Other" ~ "Indicator race = Other (0 = White)",
        #`Variable Name` == "Notreported" ~ "Indicator race = Not reported (0 = White)",
        `Variable Name` == "Unknown" ~ "Indicator race is Unknown (1 if Unknown, 0 otherwise)",
        #`Variable Name` == "URMforsubcohortsampling" ~ "URMforsubcohortsampling = ",
        #`Variable Name` == "HighRiskInd" ~ "Baseline covariate high risk pre-existing condition (1=yes, 0=no)",
        #`Variable Name` == "HIVinfection" ~ "HIV infection (1=yes, 0=no)"
        #`Variable Name` == "Country.X1" ~ "Indicator country =  ",
        `Variable Name` == "pre.study.booster.until.studydose1.day" ~ "Number of days from last prior vaccination until enrollment",
        `Variable Name` == "pre.study.booster.until.studydose1.ind" ~ "Binary indicator that the number of days from last prior vaccination until enrollment is greater than the median value (1 if yes, 0 otherwise)",
        `Variable Name` == "primary.booster.type.J.J" ~ "Indicator that primary booster type is J, J (Reference for primary.booster.type)",
        `Variable Name` == "primary.booster.type.J.M" ~ "Indicator that primary booster type is J, M (1 if yes, 0 otherwise)",
        `Variable Name` == "primary.booster.type.J.P" ~ "Indicator that primary booster type is J, P (1 if yes, 0 otherwise)",
        `Variable Name` == "primary.booster.type.M.M" ~ "Indicator that primary booster type is M, M (1 if yes, 0 otherwise)",
        `Variable Name` == "primary.booster.type.M.P" ~ "Indicator that primary booster type is M, P (1 if yes, 0 otherwise)",
        `Variable Name` == "primary.booster.type.P.M" ~ "Indicator that primary booster type is P, M (1 if yes, 0 otherwise)",
        `Variable Name` == "primary.booster.type.P.P" ~ "Indicator that primary booster type is P, P (1 if yes, 0 otherwise)"
      ),
      Comments = "")
  }
  
  return(data)
}
