## ------------------------------------------------------------------------
## Name:         4.0_ho_covid_analysis_tables_vax_status.R
## Author:       Cara McKenna
## Date Created: 03/06/2023
##
## Purpose:      Output results table in csv files.
##               
## Notes:        none
## ------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(labelled)
library(tableone)
library(cobalt)
# library(data.table)
# library(readxl)

outfx <- "vax_status_2024-05-29"

loc <- "/data/tide/projects/ho-covid/annals_resubmission/vax_status"

dataLoc <- paste0(loc,"/data")
outLoc <- paste0(loc,"/output")

covidDat <- readRDS(file.path(dataLoc,paste0("covid_data_simulation_vax_status_2024-05-08.Rds")))

day5Dat <- readRDS("/data/tide/projects/ho-covid/annals_resubmission/vax_status_clp0.1/data/ho_covid_day5_2024-05-23.Rds")

after_match_new_events <- readRDS(file.path(dataLoc,
                                            paste0("ho_covid_after_match_new_events_vax_status_2024-05-22.Rds"))) %>% 
  select(-cov_grp,-case)

matchedDSN <- paste0("covid_data_simulation_matched_vax_status_2024-05-08.Rds")

## Pre-Matching Hospitalization Level Summary
# create dataset for table - one row per hospID
covid_data_simulation_t <- covidDat %>% 
  group_by(hospID) %>% 
  mutate(case=max(case)) %>% 
  ungroup() %>% 
  select(hospID,case,omicron,encAge,gender,raceEthnicity,atrib_location,atrib_location,encDeathDerived,
         dischargeDispositionGroup2,prevDischarge90Days,LOS,facility_admission,
         starts_with("elix_"),starts_with("ccsr_")) %>% 
  distinct() %>% 
  mutate(race_asian=if_else(!is.na(raceEthnicity) & raceEthnicity=="Asian",1,0),
         race_black=if_else(!is.na(raceEthnicity) & raceEthnicity=="Black",1,0),
         race_hispanic=if_else(!is.na(raceEthnicity) & raceEthnicity=="Hispanic",1,0),
         race_other=if_else(!is.na(raceEthnicity) & raceEthnicity=="Other",1,0),
         race_two=if_else(!is.na(raceEthnicity) & raceEthnicity=="Two or More",1,0),
         race_white=if_else(!is.na(raceEthnicity) & raceEthnicity=="White",1,0),
         race_miss=if_else(raceEthnicity==""|is.na(raceEthnicity),1,0),
         
         ccsr_cat_prefix_OTHER = if_else(ccsr_cat_prefix_SYM==1 | ccsr_cat_prefix_DEN==1 | ccsr_cat_prefix_EXT==1 |
                                           ccsr_cat_prefix_EYE==1 | ccsr_cat_prefix_FAC==1 | ccsr_cat_prefix_EAR==1 |
                                           ccsr_cat_prefix_MAL==1, 1, 0)
  )

covid_data_simulation_t <- left_join(covid_data_simulation_t,day5Dat,by="hospID") %>% 
  mutate(
    maxO2device_none_day5 = if_else(maxO2deviceN_day5==0 | maxO2deviceN_day5==1, 1, 0),
    maxO2device_nc_day5 =    if_else(maxO2deviceN_day5==2, 1, 0),
    maxO2device_am_day5 =    if_else(maxO2deviceN_day5 >2 & maxO2deviceN_day5 < 7, 1, 0),
    maxO2device_bvent_day5 =    if_else(maxO2deviceN_day5 >=7, 1, 0),
    
    maxO2deviceN_day5 = factor(maxO2deviceN_day5, levels = c(0,1,2,3,4,5,6,7,8), 
                               labels = c("None and SpO2 >= 95%","None and SpO2 < 95%","Nasal Cannula","Simple Mask",
                                          "Oxygen Conserving Device","Non-rebreather","High Flow","BIPAP","Ventilator")),
    
    serviceGroup2_day5_miss = if_else(serviceGroup2_day5==""|is.na(serviceGroup2_day5),1,0),
    
    day5_to_discharge_days_copy = day5_to_discharge_days,
    LOS_copy                    = LOS,
    hospitalfree30day_day5_copy = hospitalfree30day_day5
  )

nrow(covid_data_simulation_t)
nrow(distinct(covidDat,hospID))

# get imputed labs
covid_data_simulation_day5 <- covidDat %>% 
  filter(hospitalDay==5) %>%
  select(hospID,ends_with("imp_nrm")) %>% 
  rename_with(~paste0(.,"_day5"), -c("hospID"))

covid_data_simulation_t <- left_join(covid_data_simulation_t,covid_data_simulation_day5,by="hospID")

nrow(covid_data_simulation_t)
nrow(distinct(covidDat,hospID))


race_vars <- covid_data_simulation_t %>% select(starts_with("race_")) %>% 
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

impnrmvars <- covid_data_simulation_t %>% select(ends_with("_imp_nrm_day5")) %>% 
  colnames()

allvars_unmat <-  c(allvars,"maxO2deviceN_day5","maxO2device_none_day5","maxO2device_nc_day5","maxO2device_am_day5","maxO2device_bvent_day5","maxTemp_day5","medRR_day5","medSBP_day5","medDBP_day5","maxALT_day5","maxTBILI_day5","minALB_day5","minSODIUM_day5",
                    "maxCREAT_day5","minPLT_day5","maxWBC_day5","minHCT_day5","maxGLU_day5","ICUalgo_day5","serviceGroup2_day5","serviceGroup2_day5_miss",
                    "day5_to_discharge_days","day5_to_discharge_days_copy","covVaxEver_day5","covVaxLast4months_day5","covVaxCumSum_day5",
                    "death30day_day5","readmit30day_day5","hospitalfree30day_day5","hospitalfree30day_day5_copy",impnrmvars,
                    "after_day5_new_icu","after_day5_new_hf","after_day5_new_bipap","after_day5_new_vent","after_day5_new_hfbipapvent")

catvars_unmat <- c(catvars,"maxO2deviceN_day5","maxO2device_none_day5","maxO2device_nc_day5","maxO2device_am_day5","maxO2device_bvent_day5","ICUalgo_day5","serviceGroup2_day5","serviceGroup2_day5_miss","covVaxEver_day5","covVaxLast4months_day5",
                   "death30day_day5","readmit30day_day5",impnrmvars,
                   "after_day5_new_icu","after_day5_new_hf","after_day5_new_bipap","after_day5_new_vent","after_day5_new_hfbipapvent")

show_median <- covid_data_simulation_t %>% select(ends_with("_copy")) %>% 
  colnames()

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

var_label(covid_data_simulation_t) <- var_label_list


# Create Balance Table ----------------------------------------------------

# Mike wants SMDs for each value of categorical variables
# this produces the same numbers as codiing each component as a 0/1 variable
create_bal_table <- function(indat,varlist){
  
  #Subset covariates (continuous, binary, categorical) for which we want to compute SMDs for
  #Do not include case/control status here
  x  <-  indat %>% select(all_of(varlist))
  
  #Subset the case/control status here
  case  <-  indat$case
  
  #Not relevant for us, but bal.tab() also requires a weight argument. Make dummy weights of all 1
  weights  <-  rep(1, length(case))
  
  #Now get SMDs. Compare the continuous ones reported here with those produced from tableone to make sure we have consistency
  y <- bal.tab(x,
               treat = case,
               weights = weights, binary = "std", continuous = "std")
  
  # make dataframe for output
  # the Diff.Adj is equal (in absolute value) to the SMDs from tableone
  smd <- y$Balance %>%
    mutate(SMD = round(abs(Diff.Adj),digits = 3))
  
  return(smd)
}

smd_vars_unmat <- covid_data_simulation_t %>% select(atrib_location,raceEthnicity,dischargeDispositionGroup2,maxO2deviceN_day5,serviceGroup2_day5) %>% colnames()

# Unmatched data  ---------------------------
# pre-omicron
covid_data_simulation_t_po <- covid_data_simulation_t %>% filter(omicron==0)
unmattab_po <- CreateTableOne(vars=allvars_unmat,data=covid_data_simulation_t_po,factorVars=catvars_unmat,
                              strata=c("case"))

unmat_tab_po <- print(unmattab_po,cramVars=levelvars,varLabels = TRUE,smd = TRUE, 
                      quote = FALSE, noSpaces = TRUE, printToggle = FALSE, nonnormal = show_median) 

## Save to a CSV file
write.csv(unmat_tab_po, file = file.path(outLoc,paste0("unmatched_po_",outfx,".csv")))

# unmatched pre-omicron balance table
unmat_baltab_po <- create_bal_table(covid_data_simulation_t_po,smd_vars_unmat)
write.csv(unmat_baltab_po, file = file.path(outLoc,paste0("unmatched_baltab_po_",outfx,".csv")))

# omicron
covid_data_simulation_t_o <- covid_data_simulation_t %>% filter(omicron==1)
unmattab_o <- CreateTableOne(vars=allvars_unmat,data=covid_data_simulation_t_o,factorVars=catvars_unmat,
                             strata=c("case"))


unmat_tab_o <- print(unmattab_o,cramVars=levelvars,varLabels = TRUE,smd = TRUE, 
                     quote = FALSE, noSpaces = TRUE, printToggle = FALSE, nonnormal = show_median) 

## Save to a CSV file
write.csv(unmat_tab_o, file = file.path(outLoc,paste0("unmatched_o_",outfx,".csv")))

# unmatched omicron balance table
unmat_baltab_o <- create_bal_table(covid_data_simulation_t_o,smd_vars_unmat)
write.csv(unmat_baltab_o, file = file.path(outLoc,paste0("unmatched_baltab_o_",outfx,".csv")))



############################################################################
# Matched
############################################################################
covid_data_simulation_matched <- readRDS(file.path(dataLoc,matchedDSN))

# fileDT <- sub("data/covid_data_simulation_matched_", "", matchedDSN)
# fileDT <- sub(".Rds", "", fileDT)  

after_match_vars <- after_match_new_events %>% select(starts_with("after_match_")) %>% 
  select(-ends_with("_T")) %>% 
  colnames()

# Add flag for imputed labs
# implab <- covidDat %>% select(hospID,date,ends_with("_imp_nrm"))
# nrow(covid_data_simulation_matched)
# 
# covid_data_simulation_matched <- left_join(covid_data_simulation_matched,implab,by=c("hospID","date"))
# nrow(covid_data_simulation_matched)

impnrmvars <- covid_data_simulation_matched %>% select(ends_with("_imp_nrm")) %>% 
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


covid_data_simulation_matched_t <- covid_data_simulation_matched %>% 
  left_join(.,after_match_new_events,by=c("hospID","matchid")) %>% 
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
         
         matchtodischarge_days_copy = matchtodischarge_days,
         LOS_copy                    = LOS,
         hospitalfree30day_copy = hospitalfree30day,
         ttmatch_days_copy = ttmatch_days
         
  )

show_median <- covid_data_simulation_matched_t %>% select(ends_with("_copy")) %>% 
  colnames()

var_label(covid_data_simulation_matched_t) <- var_label_list


smd_vars_mat <- covid_data_simulation_matched_t %>% select(atrib_location,raceEthnicity,dischargeDispositionGroup2,maxO2deviceN,serviceGroup2) %>% colnames()


covid_data_simulation_matched_po <- covid_data_simulation_matched_t %>% filter(cov_grp==0)

matched_tab_po <- CreateTableOne(vars=allvars2,data=covid_data_simulation_matched_po,
                                 factorVars=catvars2,
                                 strata="case")

matched_tab_po_print <- print(matched_tab_po,cramVars=levelvars,varLabels = TRUE,smd = TRUE, 
                              quote = FALSE, noSpaces = TRUE, printToggle = FALSE, nonnormal = show_median) 

## Save to a CSV file
write.csv(matched_tab_po_print, file = file.path(outLoc,
                                                 paste0("matched_po_",outfx,".csv")))

# matched pre-omicron balance table
mat_baltab_po <- create_bal_table(covid_data_simulation_matched_po,smd_vars_mat)
write.csv(mat_baltab_po, file = file.path(outLoc,paste0("matched_baltab_po_",outfx,".csv")))


covid_data_simulation_matched_o <- covid_data_simulation_matched_t %>% filter(cov_grp==1)

matched_tab_o <- CreateTableOne(vars=allvars2,data=covid_data_simulation_matched_o,
                                factorVars=catvars2,
                                strata="case")

matched_tab_o_print <- print(matched_tab_o,cramVars=levelvars,varLabels = TRUE,smd = TRUE, 
                             quote = FALSE, noSpaces = TRUE, printToggle = FALSE, nonnormal = show_median) 

## Save to a CSV file
write.csv(matched_tab_o_print, file = file.path(outLoc,
                                                paste0("matched_o_",outfx,".csv")))

# matched pre-omicron balance table
mat_baltab_o <- create_bal_table(covid_data_simulation_matched_o,smd_vars_mat)
write.csv(mat_baltab_o, file = file.path(outLoc,paste0("matched_baltab_o_",outfx,".csv")))
