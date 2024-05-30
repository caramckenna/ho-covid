## ------------------------------------------------------------------------
## Name:         2.0_ho_covid_analysis_match_vax_status.R
## Author:       Cara McKenna
## Date Created: 03/06/2023
##
## Purpose:      Preform exact matching followed by propensity score matching.
##               
## Notes:        none
## ------------------------------------------------------------------------

library(tidyverse)
library(survival)
library(ipw)
library(optmatch)
library(lubridate)
library(MatchIt)
library(naniar)
library(lmtest)
library(RItools)
library(sandwich)
library(broom)
library(knitr)

filesfx <- "vax_status_2024-05-08"

loc <- "/data/tide/projects/ho-covid/annals_resubmission/add_addtl_vars_dr5_sdm/sensitivity_analyses/vax_status/"

dataLoc <- paste0(loc,"/data")

epCombF <-"/data/tide/projects/ho-covid/annals_resubmission/add_addtl_vars_dr/data/NEWS2_episode_combined_2024-03-04.Rds"

# need this to add startDTS (forgot to add to input dataset)
epComb <- readRDS(epCombF) %>% 
  select(patientID,episodeID,startDTS)

  covid_data_simulation <- readRDS("/data/tide/projects/ho-covid/annals_resubmission/add_addtl_vars_dr4_sdm/data/ho_covid_ar_av_dr4_sdm_2024-04-28.Rds") %>% 
  # limit to locations from original analysis; April 2023 was last month of universal covid testing
  filter(atrib_location %in% c("BWH","BWF","NSMC","MGH","NWH") & admitDT < "2023-05-01") %>% 
  left_join(.,epComb,by=c("patientID","episodeID")) %>% 
  mutate(hospID               = paste(patientID,episodeID,sep = "-"),
         # move these to ho_covid
         # hospLOS = as.numeric(difftime(hospitalDischargeDTS,hospitalAdmitDTS,units = "days")),
         # round up to nearist whole number (this will include day of discharge)
         LOS = (as.numeric(date(hospitalDischargeDTS) - date(startDTS))) +1,
         covVaxEver = ifelse(!is.na(lastVaxDT), 1, 0),
         #  date%m-%months(4) rolls back to the the first 'real date'
         #  subtracting 4 months from 2022-03-31 was producing NA
         covVaxLast4months = if_else(!is.na(lastVaxDT) & lastVaxDT >= date%m-%months(4), 1, 0),
         # hospital=word(mrnSource,1),
         dischargeLocation2 = ifelse(str_detect(dischargeLocation,"SLM"),"NSMC", word(dischargeLocation,1)),
         serviceGroup2 = case_when(
           toupper(serviceGroup) == "HOSPICE"    ~ "Oncology",
           toupper(serviceGroup) == "OTHER"      ~ "Medicine",
           TRUE                              ~ serviceGroup
         ),
         # recode missing as NA
         raceEthnicity=ifelse(raceEthnicity=="Missing",NA,raceEthnicity),
         gender = ifelse(gender=="Missing",NA,gender),
         serviceGroup = ifelse(serviceGroup=="",NA,serviceGroup),
         serviceGroup2 = ifelse(serviceGroup2=="",NA,serviceGroup2),
         AdmitDay = ceiling(as.numeric(hospitalAdmitDTS)/86400),
         encDeathDerived = ifelse(!is.na(deathDTS_derived) & deathDTS_derived >= hospitalAdmitDTS &
                                    deathDTS_derived <= hospitalDischargeDTS,1,0)) %>% 
  rename(highestO2device=maxO2device, covFlagMinDT=minCovFlagDT, covFlagMaxDT=maxCovFlagDT) %>% 
  mutate_at(vars(starts_with(c("max","min","med","elix"))), as.numeric) %>% 
  # combine some comorbidities & oxygen device groups
  mutate(elix_CANCER = if_else(elix_CANCER_LEUK==1 | elix_CANCER_LYMPH==1 | elix_CANCER_METS==1 |
                   elix_CANCER_NSITU==1 | elix_CANCER_SOLID==1, 1, 0),
         elix_LIVER  = if_else(elix_LIVER_MLD==1 | elix_LIVER_SEV==1, 1, 0),
         elix_NEURO  = if_else(elix_NEURO_MOVT==1 | elix_NEURO_OTH==1 | elix_NEURO_SEIZ==1 | elix_PARALYSIS==1, 1, 0),
         elix_DIAB   = if_else(elix_DIAB_CX==1 | elix_DIAB_UNCX==1, 1, 0),
         elix_RENLFL = if_else(elix_RENLFL_MOD==1 | elix_RENLFL_SEV==1, 1, 0),
         elix_DRUGALC = if_else(elix_ALCOHOL==1 | elix_DRUG_ABUSE==1, 1, 0),
         o2Device_hf_bipap_vent = if_else(o2Device_hf==1| o2Device_bipap==1 | o2Device_vent==1, 1, 0),
         
         maxO2device_none = if_else(maxO2deviceN==0 | maxO2deviceN==1, 1, 0),
         maxO2device_nc =    if_else(maxO2deviceN==2, 1, 0),
         maxO2device_am =    if_else(maxO2deviceN >2 & maxO2deviceN < 7, 1, 0),
         maxO2device_bvent =    if_else(maxO2deviceN >=7, 1, 0),
         
         maxO2deviceN_clp = case_when(
           maxO2deviceN==0 | maxO2deviceN==1  ~ 0,
           maxO2deviceN==2                    ~ 1,
           maxO2deviceN >2 & maxO2deviceN < 7 ~ 2,
           maxO2deviceN >=7                    ~3)
         
         ) %>% 
  replace_na(list(ICUalgo=0)) 


summary(covid_data_simulation)
covid_data_simulation %>% distinct(atrib_location)
covid_data_simulation %>% select(admitDT,dischargeDT) %>% summary()
covid_data_simulation %>% filter(serviceGroup==""|is.na(serviceGroup)) %>% distinct(patientServiceDSC)
nrow(covid_data_simulation)
covid_data_simulation %>% select(starts_with("maxO2device")) %>% summary()
covid_data_simulation %>% distinct(maxO2deviceN,maxO2device_none,maxO2device_nc,maxO2device_am,maxO2device_bvent) %>% 
  arrange(maxO2deviceN)
# get max date available in dataset
# will use later to filter dataset when looking at 30-day readmission
# max_admitDT <- max(covid_data_simulation$admitDT)

excl <- covid_data_simulation %>% filter(toupper(serviceGroup)=="EXCLUDE") %>% 
  distinct(hospID)
nrow(excl)

nrow(distinct(covid_data_simulation,hospID))

# remove hospitalizations with serviceGroup=="exclude"
covid_data_simulation <- anti_join(covid_data_simulation,excl,by="hospID")

nrow(distinct(covid_data_simulation,hospID))

#######################
#Generate fake data####
#######################
# set.seed(1729)
# covid_data_simulation = data.frame(hospID = rep(1:1000, sample(20:30, 1000, replace = T)))
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(hospitalDay = row_number())
# covid_data_simulation$covidFlag = rbinom(nrow(covid_data_simulation), 1, 6/10)
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(hospital = sample(c("A", "B", "C", "D", "E"), 1))
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(service = sample(c("alpha", "beta", "gamma", "delta", "epsilon"), 1))
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(sex = sample(c("Female", "Male"), 1))
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(race = sample(c("White", "Black", "Asian", "Native"), 1))
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(age = sample(c(20:60), 1))
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(vital = rnorm(n()))
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(hospitalAdmitDTS = as.POSIXct(1607405580 + round(rnorm(1, 0, 1e6)), origin = "1970-01-01", tz = "UTC"))
# covid_data_simulation$ICU = rbinom(nrow(covid_data_simulation), 1, 1/5)

#################################################################################
#Function to find position of last element in a length n sequence of running 1's
#If no such sequence exists, this function returns 0's for the entire sequence
#Therefore, in order to both satisfy (1) COVID-19 flag >= hospital day 5 and
#COVID-19 flag active >= 4 days, this function must return >=8
#################################################################################
# g <- function(signal, n = 4) {
#   r <- rle(signal)
#   end = cumsum(r$lengths)
#   start = c(1, lag(end)[-1] + 1)
#   idx = (end - start + 1 >= n)&(r$values == 1)
#   
#   start_idx = head(start[idx][end[idx] >= 8], 1)
#   end_idx = head(end[idx][end[idx] >= 8], 1)
#   
#   case = numeric(length(signal))
#   if (length(end_idx) > 0) {
#     case[max(start_idx, 5)] = 1
#   }
#   return(case)
# }

#Define the `case` that results from logic regarding covidFlag
# covid_data_simulation = covid_data_simulation %>% group_by(hospID) %>% mutate(case = g(covidFlag))

#####################################
# Remove Community Acquired Covid
#####################################
covidFlagd1tod4 <- covid_data_simulation %>%
  filter(hospitalDay<5) %>%
  group_by(hospID) %>%
  mutate(covidFlagSumd1tod4=sum(covidFlag)) %>%
  distinct(hospID,covidFlagSumd1tod4) %>%
  ungroup()

# create dataset of hospitalizations with covid flag active before day 5
comCovid <- covidFlagd1tod4 %>%
  filter(covidFlagSumd1tod4>0)


# remove hospitalizations with covid flag active before day 5
# from email Re: FW: HO-Covid analysis on 9/23/2022
covid_data_simulation <- anti_join(covid_data_simulation,comCovid,by="hospID")

# define cases 
# Activation of COVID-19 flag on hospital day 5 or greater
# COVID-19 flag remains active for ≥4 days
covid_data_simulation <- covid_data_simulation %>% 
  group_by(hospID) %>% 
  arrange(date) %>% 
  mutate(covidFlagCumSum=cumsum(covidFlag),
         # COVID-19 flag remains active for ≥4 days includes post-discharge days
         case=if_else(covidFlagCumSum==1 & hospitalDay>=5 & tot_covflg_days>=4,1,0)) %>% 
  ungroup()



# moving LOS filter here so we can count total hospitalization before filtering
covid_data_simulation %>% distinct(patientID,episodeID) %>% nrow()

cdat <- covid_data_simulation %>% distinct(patientID,episodeID,case,omicron) %>% 
  group_by(patientID,episodeID) %>% 
  filter(case==max(case))

nrow(cdat)

cdat %>% group_by(omicron,case) %>% tally()

covid_data_simulation <- covid_data_simulation %>% filter(LOS >= 5)

# N cases
cases <- covid_data_simulation %>% filter(case==1)  %>% 
  distinct(hospID,hospitalDay,tot_covflg_days)


summary(cases)
nrow(cases)




#################################################################################
# Fill in missing labs & vitals 
#################################################################################
covid_data_simulation2 <- covid_data_simulation %>% 
  # indicators if values were imputed
  mutate(maxALT_imp=ifelse(is.na(maxALT),1,0),
         maxTBILI_imp=ifelse(is.na(maxTBILI),1,0),
         minALB_imp=ifelse(is.na(minALB),1,0),
         minSODIUM_imp=ifelse(is.na(minSODIUM),1,0),
         maxCREAT_imp=ifelse(is.na(maxCREAT),1,0),
         minPLT_imp=ifelse(is.na(minPLT),1,0),
         maxWBC_imp=ifelse(is.na(maxWBC),1,0),
         minHCT_imp=ifelse(is.na(minHCT),1,0),
         maxGLU_imp=ifelse(is.na(maxGLU),1,0),
         maxTemp_imp=ifelse(is.na(maxTemp),1,0),
         medRR_imp=ifelse(is.na(medRR),1,0),
         medSBP_imp=ifelse(is.na(medSBP),1,0),
         medDBP_imp=ifelse(is.na(medDBP),1,0)) %>% 
  group_by(hospID) %>% 
  arrange(hospitalDay) %>% 
  # fill in last known value from encounter
  fill(c("maxALT","maxTBILI","minALB","minSODIUM","maxCREAT","minPLT","maxWBC","minHCT","maxGLU","maxTemp","medRR","medSBP","medDBP")) %>% 
  # if still missing impute as normal
  mutate(maxALT_imp_nrm=ifelse(is.na(maxALT),1,0),
         maxTBILI_imp_nrm=ifelse(is.na(maxTBILI),1,0),
         minALB_imp_nrm=ifelse(is.na(minALB),1,0),
         minSODIUM_imp_nrm=ifelse(is.na(minSODIUM),1,0),
         maxCREAT_imp_nrm=ifelse(is.na(maxCREAT),1,0),
         minPLT_imp_nrm=ifelse(is.na(minPLT),1,0),
         maxWBC_imp_nrm=ifelse(is.na(maxWBC),1,0),
         minHCT_imp_nrm=ifelse(is.na(minHCT),1,0),
         maxGLU_imp_nrm=ifelse(is.na(maxGLU),1,0),
         maxALT=ifelse(is.na(maxALT),30,maxALT),
         maxTBILI=ifelse(is.na(maxTBILI),1.0,maxTBILI),
         minALB=ifelse(is.na(minALB),4,minALB),
         minSODIUM=ifelse(is.na(minSODIUM),140,minSODIUM),
         maxCREAT=ifelse(is.na(maxCREAT),1.0,maxCREAT),
         minPLT=ifelse(is.na(minPLT),200,minPLT),
         maxWBC=ifelse(is.na(maxWBC),8,maxWBC),
         minHCT=ifelse(is.na(minHCT),36,minHCT),
         maxGLU=ifelse(is.na(maxGLU),100,maxGLU))
         
saveRDS(covid_data_simulation2,
        file.path(dataLoc,paste0("covid_data_simulation_",filesfx,".Rds")))

# For vitals only look back 1 day prior
# covid_data_simulation_first_case3 <- covid_data_simulation_first_case2 %>% 
#   group_by(hospID) %>% 
#   arrange(hospitalDay) %>% 
#   mutate(maxTemp = ifelse(is.na(maxTemp),lag(maxTemp),maxTemp),
#          medRR = ifelse(is.na(medRR),lag(medRR),medRR),
#          medSBP = ifelse(is.na(medSBP),lag(medSBP),medSBP),
#          medDBP = ifelse(is.na(medDBP),lag(medDBP),medDBP))


# create match day -2 variables
# 4/3/2024 these _m2 variables are not correct once the case day is reset below
# covid_data_simulation3 <- covid_data_simulation2 %>% 
#   group_by(hospID) %>% 
#   mutate_at(c("ICUalgo","serviceGroup2","maxTemp","medRR",
#               "medSBP","medDBP","maxALT","maxTBILI","minALB",
#               "minSODIUM","maxCREAT","minPLT","maxWBC",
#               "minHCT","maxGLU","maxO2deviceN"), list(m2 = ~lag(.,n=2))) %>%
#   mutate(
#     #Define a hospitalDay0 that's hospitalDay - 1, for purposes of defining 
#     # (start, stop) when modeling survival probabilities for propensity scores
#     hospitalDay0 = hospitalDay - 1) %>% 
#   ungroup()


# Redefine case_m2 as case (first day covidFlag=1) to case - 2 days
# Mike wants all matching based on first covid day minus 2
covid_data_simulation3 <- covid_data_simulation2 %>% 
  group_by(hospID) %>% 
  mutate(case_m2=lead(case, n = 2, default = 0),
         #Define a hospitalDay0 that's hospitalDay - 1, for purposes of defining 
         # (start, stop) when modeling survival probabilities for propensity scores
         hospitalDay0 = hospitalDay - 1)

# Subset data to the case_m2 indicator for each individual (once case always case)
# This removes all days AFTER case_m2=1
# So does not include actual first covid day
# All code after this, case means case - 2 days
covid_data_simulation_first_case <- covid_data_simulation3 %>% 
  group_by(hospID) %>%
  filter(cumall(lag(case_m2, n = 1, default = 0) != 1)) %>% 
  # for simplicity, renaming case_m2=case here so that remainder of code does not have to change
  mutate(case=case_m2)

## Study Population
# N Hospitalizations
nrow(distinct(covid_data_simulation,hospID))

# Min Admit Date
min(covid_data_simulation$admitDT)

# Max Admit Date
max(covid_data_simulation$admitDT)

# N Hospital Days
nrow(covid_data_simulation)

# N Cases
nrow(cases)


# check missing data
# create dataset with only vars used in propensity score
covid_data_simulation_first_caseC <- covid_data_simulation_first_case  %>% ungroup() %>% 
  select(encAge , raceEthnicity , ICUalgo , gender , serviceGroup2,
           atrib_location , maxO2deviceN ,  maxTemp , medRR , medSBP , medDBP ,
           maxWBC , minHCT , minPLT , maxCREAT , minALB , minSODIUM , maxGLU , maxALT , maxTBILI , 
           covVaxCumSum , covVaxEver, covVaxLast4months, elix_ANEMDEF , elix_AUTOIMMUNE ,
           elix_CANCER , elix_CBVD , elix_DEMENTIA , elix_DIAB , elix_HF , elix_LIVER ,
           elix_LUNG_CHRONIC , elix_NEURO , elix_OBESE , elix_RENLFL , elix_DRUGALC ,
           elix_WGHTLOSS  , elix_index_mortality,starts_with("ccsr_cat_prefix"),
         case)

miss_var_summary(covid_data_simulation_first_caseC)

## Check Missing Data
# The propensity model will drop any observations with missing data.

### Overall missing data ###
covid_data_simulation_first_caseC %>% miss_var_summary()

covid_data_simulation_first_caseC %>% filter(case==1) %>% miss_var_summary()


# Remove episodes with missing data
missDataEps <- covid_data_simulation_first_case %>% filter(is.na(raceEthnicity) | is.na(serviceGroup2) | is.na(medRR) | 
                                                              is.na(maxTemp) | is.na(medDBP) | is.na(medSBP) | 
                                                              is.na(gender)) %>% 
  distinct(hospID)
  

# covid_data_simulation_first_case4 <- covid_data_simulation_first_case3 %>% filter(!is.na(raceEthnicity) & !is.na(serviceGroup) & !is.na(medRR) & 
#              !is.na(maxTemp) & !is.na(medDBP) & !is.na(medSBP) & !is.na(gender) & !is.na(serviceGroup_m2))
# return all rows from x without a match in y
nrow(distinct(covid_data_simulation_first_case,hospID))
nrow(missDataEps)
covid_data_simulation_first_case2 <- anti_join(covid_data_simulation_first_case,missDataEps)
nrow(distinct(covid_data_simulation_first_case2 ,hospID))


## Survival Propensity Score
#Compute survival propensity score for case
#Include ALL covariates (both exact matching covariates (e.g. hospital, service) and other fuzzy-matching covariates (e.g. sex, race, age))
# 4/4/2024 replaced maxO2deviceN with individual variables
PS_mod = coxph(Surv(hospitalDay0, hospitalDay, case) ~ encAge + raceEthnicity + ICUalgo + gender +
                 serviceGroup2 + atrib_location + maxTemp + medRR + medSBP + medDBP +
                 maxWBC + minHCT + minPLT + maxCREAT + minALB + minSODIUM + maxGLU + maxALT + maxTBILI +
                 elix_ANEMDEF + elix_AUTOIMMUNE +
                 elix_CANCER + elix_CBVD + elix_DEMENTIA + elix_DIAB + elix_HF + elix_LIVER +
                 elix_LUNG_CHRONIC + elix_NEURO + elix_OBESE + elix_RENLFL + elix_DRUGALC +
                 elix_WGHTLOSS + elix_index_mortality +
                 ccsr_cat_prefix_DIG+ccsr_cat_prefix_MUS+ccsr_cat_prefix_MBD+ccsr_cat_prefix_CIR+
                 ccsr_cat_prefix_END+ccsr_cat_prefix_INF+ccsr_cat_prefix_GEN+ccsr_cat_prefix_NVS+
                 ccsr_cat_prefix_INJ+ccsr_cat_prefix_NEO+ccsr_cat_prefix_SYM+ccsr_cat_prefix_RSP+
                 ccsr_cat_prefix_BLD+ccsr_cat_prefix_DEN+ccsr_cat_prefix_SKN+ccsr_cat_prefix_PRG+
                 # 3/18/2024 removing ccsr_cat_prefix_EAR bc no cases =1 
                 ccsr_cat_prefix_EXT+ccsr_cat_prefix_EYE+ccsr_cat_prefix_FAC+ccsr_cat_prefix_MAL+
                 maxO2deviceN_clp + covVaxEver + covVaxLast4months + covVaxCumSum,
               data = covid_data_simulation_first_case2)


#Because of how coxph() handles time-varying covariates, we need to 
#perform some algebra magic to get a survival probability
cumhaz = stepfun(y = c(0, basehaz(PS_mod)$hazard), x =c(basehaz(PS_mod)$time), 0)

covid_data_simulation_first_case5 <- covid_data_simulation_first_case2
# covid_data_simulation_first_case5 <- 
#   covid_data_simulation_first_case5 %>% 
#   mutate(hospital=word(mrnSource,1))


covid_data_simulation_first_case5$basecumhaz = cumhaz(covid_data_simulation_first_case5$hospitalDay)
covid_data_simulation_first_case5 = covid_data_simulation_first_case5 %>% group_by(hospID) %>% mutate(basehaz = c(0, diff(basecumhaz)))
covid_data_simulation_first_case5$haz = covid_data_simulation_first_case5$basehaz*predict(PS_mod, newdata = covid_data_simulation_first_case5, type = "risk")
covid_data_simulation_first_case5 = covid_data_simulation_first_case5 %>% group_by(hospID) %>% mutate(cumhaz = cumsum(haz))
#covid_data_simulation$haz2 = predict(PS_mod, newdata = covid_data_simulation, type = "expected")
covid_data_simulation_first_case5$survival = exp(-covid_data_simulation_first_case5$cumhaz)
#covid_data_simulation$survival2 = predict(PS_mod, newdata = covid_data_simulation, type = "survival")



#Convert hospital days (which is relative to admission date) to an absolute hospital day (measured in terms of days since 1970-01-01, which is arbirary)
covid_data_simulation_first_case5$day = ceiling(as.numeric(covid_data_simulation_first_case5$hospitalAdmitDTS)/86400) + covid_data_simulation_first_case5$hospitalDay0
min_day = 1
max_day = 365
covid_data_simulation_risk_set = as.data.frame(covid_data_simulation_first_case5)
covid_data_simulation_matched = NULL

clustid=1

#Risk set matching via sequence matching
for (i in min_day:max_day) {
  #Current risk set strata
  current_strata = covid_data_simulation_risk_set[covid_data_simulation_risk_set$hospitalDay == i,]
  #if(sum(current_strata$case) !=0){ #If there is at least 1 case in this set, perform matching
  if((sum(current_strata$case) !=0)&(sum(current_strata$case==0) > 0)){ #If there is at least 1 case in this set, perform matching
    #Hope to match 1:3 cases to controls; if not feasible, pairmatch() throws an error, and we'll forgo exact matching
    # print(i)
    #Variable to put a huge penalty if match day exceeds control LOS
    LOS.dist = outer(current_strata$hospitalDay[current_strata$case == 1], current_strata$LOS[current_strata$case == 0], ">")
    rownames(LOS.dist) = rownames(current_strata[current_strata$case == 1,])
    colnames(LOS.dist) = rownames(current_strata[current_strata$case == 0,])
    
    
    #Define the difference in Admission dates; will be used to put penalty if this is greater than 14
    AdmitDay.dist = match_on(case ~ AdmitDay, data=current_strata, method = "euclidean")
    
    #Define difference in age; will be used to exclude if ages differ by more than 2 years
    Age.dist = match_on(case ~ encAge, data=current_strata, method = "euclidean")
    
    #PS distances
    PS.dist = match_on(case ~ survival, data=current_strata, method = "euclidean")
    
    ##############
    #UPDATED: Tom#
    #elix_index_mortality matching
    elix_index_mortality.dist = match_on(case ~ elix_index_mortality, data=current_strata, method = "euclidean")
    ##############
    
    #case ~ hospital + service can be expanded to include other things (ICU, etc)
    #Cara: expand the `exactMatch` covariate list as see fit
    # first find exact matches
    # then perform caliper matching on LOS, admission day, and age
    # last perform propensity matching after all previous criteria satisfied
    
    ##############
    #UPDATED: Tom#
    matched = match_on(PS.dist + caliper(PS.dist, 0.2*sd(PS.dist)) + caliper(LOS.dist, 0.9) + caliper(AdmitDay.dist, 90), 
                       data = current_strata, 
                       within = exactMatch(case ~ atrib_location + serviceGroup2 + ICUalgo, data = current_strata), method = "euclidean")
    ##############
    res = try(match <- pairmatch(matched, data = current_strata, controls = 2), silent = T)
    
    failed = "try-error" %in% class(res)
    
    #If first-pass matching fails, try second-pass matching with looser conditions
    if (failed){
      ##############
      #UPDATED: Tom#
      matched = match_on(PS.dist + caliper(PS.dist, 0.2*sd(PS.dist)) + caliper(LOS.dist, 0.9) + caliper(AdmitDay.dist, 90), 
                         data = current_strata, 
                         within = exactMatch(case ~ atrib_location + serviceGroup2 + ICUalgo, data = current_strata), method = "euclidean")
      ##############
      
      res = try(match <- pairmatch(matched, data = current_strata, controls = 2), silent = T)
      
      #If second-pass matching fails, try third-pass matching with looser conditions
      if ("try-error" %in% class(res)) {
        ##############
        #UPDATED: Tom#
        matched = match_on(PS.dist + caliper(PS.dist, 0.2*sd(PS.dist)) + caliper(LOS.dist, 0.9) + caliper(AdmitDay.dist, 90), 
                           data = current_strata, 
                           within = exactMatch(case ~ atrib_location + serviceGroup2 + ICUalgo, data = current_strata),
                           method = "euclidean")
        ##############
        
        res = try(match <- pairmatch(matched, data = current_strata, controls = 2), silent = T)
        
        #If all this fails, let's drop this case
        if ("try-error" %in% class(res)) {
          match = NA
        }
      }
    }
    
    #Concatenate current risk set strata with matches
    covid_data_temp <- cbind(current_strata, matches = match)
    
    #Remove controls which were not matched
    covid_data_temp2 = covid_data_temp[!is.na(covid_data_temp$matches),]
    
    #Make unique matchid
    if (nrow(covid_data_temp2) > 0) {
      covid_data_temp2$matchid = as.integer(covid_data_temp2$matches) + (clustid - 1)
      clustid = clustid + max(as.integer(covid_data_temp2$matches))  
    }
    
    
    
    #Add these exposed/not exposed matches to existing matched data
    covid_data_simulation_matched <- rbind(covid_data_simulation_matched, covid_data_temp2)
    
    
    #Remove these exposed from risk set for next iteration (unexposed can remain for further matches)
    covid_data_temp3 = covid_data_temp2[covid_data_temp2$case == 1,]
    covid_data_simulation_risk_set = covid_data_simulation_risk_set[!(covid_data_simulation_risk_set$patientID %in% covid_data_temp3$patientID),]

    
  }
}

# dim(covid_data_simulation_matched)
# sum(covid_data_simulation_matched$case==1)
# sum(covid_data_simulation_matched$case==0)
# 
# sum(covid_data_simulation_first_case5$case==1)

# keep matched case controls in the same grouping for tables & analysis
covid_data_simulation_matched2 <- covid_data_simulation_matched %>% 
  mutate(cov_grp=ifelse(case==1,omicron,NA)) %>% 
  group_by(matchid) %>% 
  arrange(desc(case)) %>% 
  fill(cov_grp) %>% 
  ungroup()


saveRDS(covid_data_simulation_matched2,file.path(dataLoc,paste0("covid_data_simulation_matched_",filesfx,".Rds")))



nrow(distinct(covid_data_simulation_matched2,matchid))
nrow(distinct(covid_data_simulation_matched2,matches))

