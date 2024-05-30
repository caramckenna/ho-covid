## ------------------------------------------------------------------------
## Name:         3.0_ho_covid_day5_vax_status.R
## Author:       Cara McKenna
## Date Created: 11/17/2023
##
## Purpose:      Create dataset with hospital day 5 to use for controls in summary tables
##               
## Notes:        none
## ------------------------------------------------------------------------

library(tidyverse)

loc <- "/data/tide/projects/ho-covid/annals_resubmission/add_addtl_vars_dr5_sdm/sensitivity_analyses/vax_status"

dataLoc <- paste0(loc,"/data")

hocovid <- readRDS(file.path(dataLoc,paste0("covid_data_simulation_vax_status_2024-05-08.Rds")))

hocovid_day5_a <- hocovid %>% ungroup() %>% 
  filter(hospitalDay==5) %>%
  select(patientID,episodeID,date,dischargeDT,serviceGroup,maxTemp,medRR,medSBP,medDBP,maxWBC,minHCT,minPLT,maxCREAT,
         minSODIUM,maxGLU,maxALT,maxTBILI,minALB,ICUalgo,lastVaxDT,covVaxCumSum,deathDTS_derived,nextAdmitDTS,
         maxO2deviceN,LOS,encDeathDerived,o2Device_hf,o2Device_bipap,o2Device_vent,o2Device_hf_bipap_vent) %>% 
  rename_with(~paste0(.,"_day5"), -c("patientID", "episodeID", "dischargeDT", "deathDTS_derived"))


hocovid_day5_b <- hocovid_day5_a %>% 
  mutate(hospID = paste(patientID,episodeID,sep = "-"),
         serviceGroup2_day5 = case_when(
           toupper(serviceGroup_day5) == "HOSPICE"    ~ "Oncology",
           toupper(serviceGroup_day5) == "OTHER"      ~ "Medicine",
           TRUE                              ~ serviceGroup_day5
         ),
         covVaxEver_day5 = ifelse(!is.na(lastVaxDT_day5), 1, 0),
         #  date%m-%months(4) rolls back to the the first 'real date'
         #  subtracting 4 months from 2022-03-31 was producing NA
         covVaxLast4months_day5 = if_else(!is.na(lastVaxDT_day5) & lastVaxDT_day5 >= date_day5%m-%months(4), 1, 0),
         
         day5_to_discharge_days = as.numeric(dischargeDT-date_day5),
         deathDT_derived = date(deathDTS_derived),
         day5_to_death_days = as.numeric((deathDT_derived - date_day5) +1),
         death30day_day5 = if_else(!is.na(day5_to_death_days) & day5_to_death_days <= 30, 1, 0),
         
         nextAdmitDT_day5 = date(nextAdmitDTS_day5),
         day5_to_next_admit_days = as.numeric((nextAdmitDT_day5 - date_day5)+1),
         readmit30day_day5 = if_else(!is.na(day5_to_next_admit_days) & day5_to_next_admit_days <= 30, 1, 0),
         
         match_to_next_admit_days = as.numeric((nextAdmitDT_day5 - date_day5)+1),
         matchLOS = (LOS_day5-5)+1,
         hospitalfree30day_ = 30-matchLOS,
         hospitalfree30day_day5 = case_when(
           encDeathDerived_day5==1 ~ 0,
           hospitalfree30day_ <0 ~ 0,
           match_to_next_admit_days <=30 ~ as.numeric((nextAdmitDT_day5-dischargeDT)-1),
           TRUE ~ hospitalfree30day_
         )) %>% 
  mutate_at(vars(starts_with(c("min","max","med"))), ~ as.numeric(.))



saveRDS(hocovid_day5_b,file.path(dataLoc,paste0("ho_covid_day5_",today(),".Rds")))
