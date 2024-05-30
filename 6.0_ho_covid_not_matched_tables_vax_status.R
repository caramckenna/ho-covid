## ------------------------------------------------------------------------
## Name:         6.0_ho_covid_not_matched_tables_vax_status.R
## Author:       Cara McKenna
## Date Created: 03/06/2023
##
## Purpose:      Compare the matched HO-COVID cases vs the unmatched HO-COVID patients.
##               
## Notes:        none
## ------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(labelled)
library(tableone)
library(cobalt)
filesfx <- "vax_status_2024-05-29-1"

dataLoc <- "/data/tide/projects/ho-covid/annals_resubmission/vax_status/data"
outLoc <- "/data/tide/projects/ho-covid/annals_resubmission/vax_status/output"

covid_data_simulation <- readRDS(file.path(dataLoc,"covid_data_simulation_vax_status_2024-05-08.Rds"))

matchedDSN <- "covid_data_simulation_matched_vax_status_2024-05-08.Rds"
covid_data_simulation_matched <- readRDS(file.path(dataLoc,matchedDSN))


covidDatIDs <- covid_data_simulation %>% distinct(hospID)
matchedIDs <- covid_data_simulation_matched %>% distinct(hospID)
notmatched <- anti_join(covidDatIDs,matchedIDs)
nrow(notmatched)

# Redefine case (first day covidFlag=1) to case - 2 days to match the matched cases
covidDat2 <- inner_join(covid_data_simulation,notmatched) %>% 
  group_by(hospID) %>% 
  mutate(case_m2=lead(case, n = 2, default = 0),
         #Define a hospitalDay0 that's hospitalDay - 1, for purposes of defining 
         # (start, stop) when modeling survival probabilities for propensity scores
         hospitalDay0 = hospitalDay - 1)

# Subset data to the case_m2 indicator for each individual (once case always case)
# This removes all days AFTER case_m2=1
# So does not include actual first covid day
nm_case <- covidDat2 %>% 
  group_by(hospID) %>%
  filter(cumall(lag(case_m2, n = 1, default = 0) != 1)) %>% 
  # for simplicity, renaming case_m2=case here so that remainder of code does not have to change
  mutate(case=case_m2,
         matched=0) %>% 
  filter(case==1) 

nrow(nm_case)

mat_case <- covid_data_simulation_matched %>% 
  filter(case==1) %>% 
  mutate(matched=1)

nrow(mat_case)

cases <- bind_rows(nm_case,mat_case)


# Create new events after match -------------------------------------------

cases2 <- cases %>% 
  select(hospID,matchid,matchday=hospitalDay) %>% 
  left_join(.,covid_data_simulation,by=c("hospID"),relationship = "many-to-many") %>% 
  arrange(hospID,hospitalDay)

new_event_after_match <- function(event) {
  
  matched_new_event <- cases2 %>% 
    group_by(hospID,matchid) %>% 
    # cannot have event on match day
    filter(matchday==hospitalDay & {{event}}!=1) %>% 
    distinct(hospID,matchid,matchday) %>% 
    left_join(.,covid_data_simulation,by=c("hospID"),relationship = "many-to-many") %>%
    filter(hospitalDay > matchday & {{event}}==1) %>% 
    group_by(hospID,matchid) %>% 
    arrange(hospID,hospitalDay) %>% 
    filter(row_number()==1) %>% 
    distinct(hospID,matchid,hospitalDay,{{event}}) %>% 
    ungroup()
  
  return(matched_new_event)
  
}


### New admission to ICU After Match Day ###

newICU <- new_event_after_match(ICUalgo) %>% 
  select(hospID,matchid,after_match_new_icu=ICUalgo,after_match_new_icu_T=hospitalDay)

### New Need for High Flow O2 After Match Day ###

newHF <- new_event_after_match(o2Device_hf) %>% 
  select(hospID,matchid,after_match_new_hf=o2Device_hf,after_match_new_hf_T=hospitalDay)


### New Need for BIPAP After Match Day ###

newBIPAP <- new_event_after_match(o2Device_bipap) %>% 
  select(hospID,matchid,after_match_new_bipap=o2Device_bipap,after_match_new_bipap_T=hospitalDay)


### New Need for Mechanical Ventilation After Match Day ###

newVent  <- new_event_after_match(o2Device_vent) %>% 
  select(hospID,matchid,after_match_new_vent=o2Device_vent,after_match_new_vent_T=hospitalDay)


### New Need for High Flow O2, BIPAP, Or Mechanical Ventilation After Match Day ###

newO2 <- new_event_after_match(o2Device_hf_bipap_vent) %>% 
  select(hospID,matchid,after_match_new_hfbipapvent=o2Device_hf_bipap_vent,after_match_new_hfbipapvent_T=hospitalDay)


max_admitDT <- max(covid_data_simulation$admitDT)

### 30-day Readmission (counting from match day) ###
readmit30day <- cases %>% 
  # only look at hospitalization that could have data 30 days post discharge
  filter(dischargeDT<=max_admitDT-30) %>% 
  mutate(nextAdmitDT = date(nextAdmitDTS),
         match_to_next_admit_days = as.numeric((nextAdmitDT - date)+1),
         readmit30day = if_else(!is.na(match_to_next_admit_days) & match_to_next_admit_days <= 30, 1, 0))

readmit30day <- readmit30day %>% 
  ungroup() %>% 
  select(hospID,matchid,readmit30day) %>% 
  distinct()


### 30-day Mortality (counting from match day) ###
death30day <- cases %>% 
  # only look at hospitalization that could have data 30 days post discharge
  filter(dischargeDT<=max_admitDT-30) %>% 
  mutate(deathDT_derived = date(deathDTS_derived),
         match_to_death_days = as.numeric((deathDT_derived - date) +1),
         death30day = if_else(!is.na(match_to_death_days) & match_to_death_days <= 30, 1, 0))

death30day <- death30day %>% 
  ungroup() %>% 
  select(hospID,matchid,death30day) %>% 
  distinct()

### 30d hospital-free days (counting from match day) ###
hospitalfree30day <- cases %>%
  mutate(nextAdmitDT = date(nextAdmitDTS),
         match_to_next_admit_days = as.numeric((nextAdmitDT - date)+1),
         matchLOS = (LOS-hospitalDay)+1,
         hospitalfree30day_ = 30-matchLOS,
         hospitalfree30day = case_when(
           encDeathDerived==1 ~ 0,
           hospitalfree30day_ <0 ~ 0,
           match_to_next_admit_days <=30 ~ as.numeric((nextAdmitDT-dischargeDT)-1),
           TRUE ~ hospitalfree30day_
         ))

hospitalfree30day <- hospitalfree30day %>% 
  ungroup() %>% 
  select(hospID,matchid,hospitalfree30day) %>% 
  distinct()


new_events <- full_join(newICU,newHF,by=c("hospID","matchid")) %>% 
  full_join(.,newBIPAP,by=c("hospID","matchid")) %>% 
  full_join(.,newVent,by=c("hospID","matchid")) %>%
  full_join(.,newO2,by=c("hospID","matchid")) %>%
  full_join(.,readmit30day,by=c("hospID","matchid")) %>% 
  full_join(.,death30day,by=c("hospID","matchid")) %>%
  full_join(.,hospitalfree30day,by=c("hospID","matchid")) %>% 
  select(-ends_with("_T"))


# Dataset for table -------------------------------------------------------

# Add flag for imputed labs
# implab <- covid_data_simulation %>% select(hospID,date,ends_with("_imp_nrm"))
# nrow(cases)
# cases <- left_join(cases,implab,by=c("hospID","date"))
# nrow(cases)

cases_t <- left_join(cases,new_events,by=c("hospID","matchid")) %>% 
  mutate(ttmatch_days= as.numeric((date-admitDT)+1),
         matchtodischarge_days= as.numeric(dischargeDT-date),
         
         deathDT_derived = date(deathDTS_derived),
         match_to_death_days = as.numeric((deathDT_derived - date) +1),
         death30day = if_else(!is.na(match_to_death_days) & match_to_death_days <= 30, 1, 0),
         
         nextAdmitDT = date(nextAdmitDTS),
         match_to_next_admit_days = as.numeric((nextAdmitDT - date)+1),
         readmit30day = if_else(!is.na(match_to_next_admit_days) & match_to_next_admit_days <= 30, 1, 0),
         race_asian=if_else(!is.na(raceEthnicity) & raceEthnicity=="Asian",1,0),
         race_black=if_else(!is.na(raceEthnicity) & raceEthnicity=="Black",1,0),
         race_hispanic=if_else(!is.na(raceEthnicity) & raceEthnicity=="Hispanic",1,0),
         race_other=if_else(!is.na(raceEthnicity) & raceEthnicity=="Other",1,0),
         race_two=if_else(!is.na(raceEthnicity) & raceEthnicity=="Two or More",1,0),
         race_white=if_else(!is.na(raceEthnicity) & raceEthnicity=="White",1,0),
         race_miss=if_else(raceEthnicity==""|is.na(raceEthnicity),1,0),
         
         maxO2deviceN = factor(maxO2deviceN, levels = c(0,1,2,3,4,5,6,7,8), 
                               labels = c("None and SpO2 >= 95%","None and SpO2 < 95%","Nasal Cannula","Simple Mask",
                                          "Oxygen Conserving Device","Non-rebreather","High Flow","BIPAP","Ventilator")),
         
         ccsr_cat_prefix_OTHER = if_else(ccsr_cat_prefix_SYM==1 | ccsr_cat_prefix_DEN==1 | ccsr_cat_prefix_EXT==1 |
                                           ccsr_cat_prefix_EYE==1 | ccsr_cat_prefix_FAC==1 | ccsr_cat_prefix_EAR==1 |
                                           ccsr_cat_prefix_MAL==1, 1, 0),
         
         serviceGroup2_miss = if_else(serviceGroup2==""|is.na(serviceGroup2),1,0),
         
         matchtodischarge_days_copy = matchtodischarge_days,
         LOS_copy                    = LOS,
         hospitalfree30day_copy = hospitalfree30day,
         ttmatch_days_copy = ttmatch_days
  ) %>% ungroup() %>% 
  mutate_at(vars(starts_with("after_match_new"),c("readmit30day","death30day")), ~replace_na(., 0))

race_vars <- cases_t %>% select(starts_with("race_")) %>% 
  colnames()

allvars <- c("atrib_location","gender","encAge","raceEthnicity",race_vars,"LOS","LOS_copy","encDeathDerived",
             "dischargeDispositionGroup2","prevDischarge90Days","facility_admission",
             #Comormidities
             "elix_AIDS","elix_ALCOHOL","elix_ANEMDEF","elix_AUTOIMMUNE","elix_BLDLOSS","elix_CANCER_LEUK",
             "elix_CANCER_LYMPH","elix_CANCER_METS","elix_CANCER_NSITU","elix_CANCER_SOLID","elix_COAG","elix_DEMENTIA",
             "elix_DEPRESS","elix_DIAB_CX","elix_DIAB_UNCX","elix_DRUG_ABUSE","elix_HF","elix_HTN_CX","elix_HTN_UNCX",
             "elix_LIVER_MLD","elix_LIVER_SEV","elix_LUNG_CHRONIC","elix_NEURO_MOVT","elix_NEURO_OTH","elix_NEURO_SEIZ",
             "elix_OBESE","elix_PARALYSIS","elix_PERIVASC","elix_PSYCHOSES","elix_PULMCIRC","elix_RENLFL_MOD",
             "elix_RENLFL_SEV","elix_THYROID_HYPO","elix_THYROID_OTH","elix_ULCER_PEPTIC",
             "elix_VALVE","elix_WGHTLOSS","elix_CBVD",
             "elix_CANCER","elix_LIVER","elix_NEURO","elix_DIAB","elix_RENLFL","elix_DRUGALC",
             "elix_index_mortality",
             "ccsr_cat_prefix_DIG","ccsr_cat_prefix_MUS","ccsr_cat_prefix_MBD","ccsr_cat_prefix_CIR","ccsr_cat_prefix_END","ccsr_cat_prefix_INF","ccsr_cat_prefix_GEN","ccsr_cat_prefix_NVS",
             "ccsr_cat_prefix_INJ","ccsr_cat_prefix_NEO","ccsr_cat_prefix_SYM","ccsr_cat_prefix_RSP","ccsr_cat_prefix_BLD","ccsr_cat_prefix_DEN","ccsr_cat_prefix_SKN","ccsr_cat_prefix_PRG",
             "ccsr_cat_prefix_EXT","ccsr_cat_prefix_EYE","ccsr_cat_prefix_FAC","ccsr_cat_prefix_EAR","ccsr_cat_prefix_MAL","ccsr_cat_prefix_OTHER")

catvars <- c("atrib_location","gender","raceEthnicity",race_vars,"encDeathDerived",
             "dischargeDispositionGroup2","prevDischarge90Days","facility_admission",
             "elix_AIDS","elix_ALCOHOL","elix_ANEMDEF","elix_AUTOIMMUNE","elix_BLDLOSS","elix_CANCER_LEUK",
             "elix_CANCER_LYMPH","elix_CANCER_METS","elix_CANCER_NSITU","elix_CANCER_SOLID","elix_COAG","elix_DEMENTIA",
             "elix_DEPRESS","elix_DIAB_CX","elix_DIAB_UNCX","elix_DRUG_ABUSE","elix_HF","elix_HTN_CX","elix_HTN_UNCX",
             "elix_LIVER_MLD","elix_LIVER_SEV","elix_LUNG_CHRONIC","elix_NEURO_MOVT","elix_NEURO_OTH","elix_NEURO_SEIZ",
             "elix_OBESE","elix_PARALYSIS","elix_PERIVASC","elix_PSYCHOSES","elix_PULMCIRC","elix_RENLFL_MOD",
             "elix_RENLFL_SEV","elix_THYROID_HYPO","elix_THYROID_OTH","elix_ULCER_PEPTIC",
             "elix_VALVE","elix_WGHTLOSS","elix_CBVD","elix_CANCER","elix_LIVER","elix_NEURO","elix_DIAB","elix_RENLFL","elix_DRUGALC",
             "ccsr_cat_prefix_DIG","ccsr_cat_prefix_MUS","ccsr_cat_prefix_MBD","ccsr_cat_prefix_CIR","ccsr_cat_prefix_END","ccsr_cat_prefix_INF","ccsr_cat_prefix_GEN","ccsr_cat_prefix_NVS",
             "ccsr_cat_prefix_INJ","ccsr_cat_prefix_NEO","ccsr_cat_prefix_SYM","ccsr_cat_prefix_RSP","ccsr_cat_prefix_BLD","ccsr_cat_prefix_DEN","ccsr_cat_prefix_SKN","ccsr_cat_prefix_PRG",
             "ccsr_cat_prefix_EXT","ccsr_cat_prefix_EYE","ccsr_cat_prefix_FAC","ccsr_cat_prefix_EAR","ccsr_cat_prefix_MAL","ccsr_cat_prefix_OTHER")


levelvars <- c("atrib_location","raceEthnicity","dischargeDispositionGroup2")


var_label_list <- list(atrib_location = "Site",
                       gender = "Gender",
                       encAge = "Age at Admission",
                       raceEthnicity = "Race/Ethnicity", 
                       LOS = "Hospital Length of Stay in Days",
                       encDeathDerived = "Death During Hospitalization",
                       dischargeDispositionGroup2 = "Discharge Disposition",
                       prevDischarge90Days = "MGB hospitalization with discharge within 90 days prior to index admission date",
                       facility_admission = "Admission From Facility",
                       elix_AIDS="Elixhauser: AIDS",
                       elix_ALCOHOL="Elixhauser: Alcohol abuse",
                       elix_ANEMDEF="Elixhauser: Deficiency anemias",
                       elix_AUTOIMMUNE="Elixhauser: Autoimmune conditions",
                       elix_BLDLOSS="Elixhauser: Chronic blood loss anemia",
                       elix_CANCER_LEUK="Elixhauser: Leukemia",
                       elix_CANCER_LYMPH="Elixhauser: Lymphoma",
                       elix_CANCER_METS="Elixhauser: Metastatic cancer",
                       elix_CANCER_NSITU="Elixhauser: Solid tumor without metastasis, in situ",
                       elix_CANCER_SOLID="Elixhauser: Solid tumor without metastasis, malignant",
                       elix_COAG="Elixhauser: Coagulopathy",
                       elix_DEMENTIA="Elixhauser: Dementia",
                       elix_DEPRESS="Elixhauser: Depression",
                       elix_DIAB_CX="Elixhauser: Diabetes with chronic complications",
                       elix_DIAB_UNCX="Elixhauser: Diabetes without chronic complications",
                       elix_DRUG_ABUSE="Elixhauser: Drug abuse",
                       elix_HF="Elixhauser: Heart failure",
                       elix_HTN_CX="Elixhauser: Hypertension, complicated",
                       elix_HTN_UNCX="Elixhauser: Hypertension, uncomplicated",
                       elix_LIVER_MLD="Elixhauser: Liver disease, mild",
                       elix_LIVER_SEV="Elixhauser: Liver disease, moderate to severe",
                       elix_LUNG_CHRONIC="Elixhauser: Chronic pulmonary disease",
                       elix_NEURO_MOVT="Elixhauser: Neurological disorders affecting movement",
                       elix_NEURO_OTH="Elixhauser: Other neurological disorders",
                       elix_NEURO_SEIZ="Elixhauser: Seizures and epilepsy",
                       elix_OBESE="Elixhauser: Obesity",
                       elix_PARALYSIS="Elixhauser: Paralysis",
                       elix_PERIVASC="Elixhauser: Peripheral vascular disease",
                       elix_PSYCHOSES="Elixhauser: Psychoses",
                       elix_PULMCIRC="Elixhauser: Pulmonary circulation disease",
                       elix_RENLFL_MOD="Elixhauser: Renal failure, moderate",
                       elix_RENLFL_SEV="Elixhauser: Renal failure, severe",
                       elix_THYROID_HYPO="Elixhauser: Hypothyroidism",
                       elix_THYROID_OTH="Elixhauser: Other thyroid disorders",
                       elix_ULCER_PEPTIC="Elixhauser: Peptic ulcer with bleeding",
                       elix_VALVE="Elixhauser: Valvular disease",
                       elix_WGHTLOSS="Elixhauser: Weight loss",
                       elix_CBVD="Elixhauser: Cerebrovascular disease",
                       elix_CANCER="Elixhauser: Combined Cancer",
                       elix_LIVER="Elixhauser: Combined Liver",
                       elix_NEURO="Elixhauser: Combined Neurological & Paralysis",
                       elix_DIAB="Elixhauser: Combined Diabetes",
                       elix_RENLFL="Elixhauser: Combined Renal failure",
                       elix_DRUGALC="Elixhauser: Combined Drug Abuse & Alcohol",
                       elix_index_mortality="Elixhauser: Index for the risk of in-hospital mortality")

after_match_vars <- new_events %>% select(starts_with("after_match_")) %>% 
  colnames()

impnrmvars <- cases_t %>% select(ends_with("_imp_nrm")) %>% 
  colnames()

show_median <- cases_t %>% select(ends_with("_copy")) %>% 
  colnames()

allvars2 <- c(allvars,"maxO2deviceN","maxO2device_none","maxO2device_nc","maxO2device_am","maxO2device_bvent","maxTemp","medRR","medSBP","medDBP","maxALT","maxTBILI","minALB","minSODIUM","maxCREAT",
              "minPLT","maxWBC","minHCT","maxGLU","ICUalgo","serviceGroup2","ttmatch_days","ttmatch_days_copy","matchtodischarge_days","matchtodischarge_days_copy",
              "covVaxEver","covVaxLast4months","covVaxCumSum",after_match_vars,"death30day","readmit30day","hospitalfree30day","hospitalfree30day_copy",
              impnrmvars)

catvars2 <- c(catvars,"maxO2deviceN","maxO2device_none","maxO2device_nc","maxO2device_am","maxO2device_bvent","ICUalgo","serviceGroup2",
              "covVaxEver","covVaxLast4months",after_match_vars,"death30day","readmit30day",impnrmvars)

var_label_list2 <- list(
  maxTemp   = "Max Temperature On Match Day Minus 2",
  medRR     = "Median RR On Match Day Minus 2",
  medSBP    = "Median SBP On Match Day Minus 2",
  medDBP    = "Median DBP On Match Day Minus 2",
  maxALT    = "Max ALT On Match Day Minus 2",
  maxTBILI  = "Max Bilirubin On Match Day Minus 2",
  maxCREAT  = "Max Creatinine  On Match Day Minus 2",
  maxGLU    = "Max Glucose On Match Day Minus 2",
  maxWBC    = "Max WBC On Match Day Minus 2",
  minALB    = "Min Albumin On Match Day Minus 2",
  minSODIUM = "Min Sodium On Match Day Minus 2",
  minPLT    = "Min Platelets On Match Day Minus 2",
  minHCT    = "Min Hematocrit On Match Day Minus 2",
  ICUalgo  = "ICU Status on Match Day Minus 2",
  serviceGroup2 = "Service Group on Match Day Minus 2",
  ttmatch_days     = "Days from Hospital Admission to Match Day",
  matchtodischarge_days = "Days from Match Day to Hospital Discharge",
  covVaxEver = "Received SARS-CoV-2 Vaccines Ever Prior to Match Day",
  covVaxLast4months = paste0("Received SARS-CoV-2 Vaccine in Last Four Months (counting from match day)"),
  covVaxCumSum = "Total SARS-CoV-2 Vaccines Prior to Match Day")

var_label_list <- c(var_label_list, var_label_list2) 

var_label(cases_t) <- var_label_list


# Create Balance Table ----------------------------------------------------

# Mike wants SMDs for each value of categorical variables
# this produces the same numbers as codiing each component as a 0/1 variable
create_bal_table <- function(indat,varlist){
  
  #Subset covariates (continuous, binary, categorical) for which we want to compute SMDs for
  #Do not include case/control status here
  x  <-  indat %>% select(all_of(varlist))
  
  #Subset the case/control status here
  mat  <-  indat$matched
  
  #Not relevant for us, but bal.tab() also requires a weight argument. Make dummy weights of all 1
  weights  <-  rep(1, length(mat))
  
  #Now get SMDs. Compare the continuous ones reported here with those produced from tableone to make sure we have consistency
  y <- bal.tab(x,
               treat = mat,
               weights = weights, binary = "std", continuous = "std")
  
  # make dataframe for output
  # the Diff.Adj is equal (in absolute value) to the SMDs from tableone
  smd <- y$Balance %>%
    mutate(SMD = round(abs(Diff.Adj),digits = 3))
  
  return(smd)
}



# Make tables -------------------------------------------------------------

smd_vars_mat <- cases_t %>% select(atrib_location,raceEthnicity,dischargeDispositionGroup2,maxO2deviceN,serviceGroup2) %>% colnames()

# Pre-Omicron
cases_po <- cases_t %>% filter(omicron==0)

cases_tab_po <- CreateTableOne(vars=allvars2,data=cases_po,
                                 factorVars=catvars2,
                                 strata="matched")

cases_tab_po_print <- print(cases_tab_po,cramVars=levelvars,varLabels = TRUE,smd = TRUE, 
                              quote = FALSE, noSpaces = TRUE, printToggle = FALSE, nonnormal = show_median) 

## Save to a CSV file
write.csv(cases_tab_po_print, file = file.path(outLoc,
                                                 paste0("cases_po_",filesfx,".csv")))

# matched pre-omicron balance table
cases_baltab_po <- create_bal_table(cases_po,smd_vars_mat)
write.csv(cases_baltab_po, file = file.path(outLoc,paste0("cases_baltab_po_",filesfx,".csv")))


# Omicron
cases_o <- cases_t %>% filter(omicron==1)

cases_tab_o <- CreateTableOne(vars=allvars2,data=cases_o,
                                factorVars=catvars2,
                                strata="matched")

cases_tab_o_print <- print(cases_tab_o,cramVars=levelvars,varLabels = TRUE,smd = TRUE, 
                             quote = FALSE, noSpaces = TRUE, printToggle = FALSE, nonnormal = show_median) 

## Save to a CSV file
write.csv(cases_tab_o_print, file = file.path(outLoc,
                                                paste0("cases_o_",filesfx,".csv")))

# matched pre-omicron balance table
cases_baltab_o <- create_bal_table(cases_o,smd_vars_mat)
write.csv(cases_baltab_o, file = file.path(outLoc,paste0("cases_baltab_o_",filesfx,".csv")))

