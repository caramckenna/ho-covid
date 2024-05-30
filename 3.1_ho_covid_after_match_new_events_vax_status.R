## ------------------------------------------------------------------------
## Name:         3.1_ho_covid_after_match_new_events_vax_status.R
## Author:       Cara McKenna
## Date Created: 11/30/2023
##
## Purpose:      Create dataset for new events that can be added to summary table.
##               
## Notes:        none
## ------------------------------------------------------------------------
library(tidyverse)

filesfx <- "vax_status_2024-05-08"
outfx <- "vax_status_2024-05-22"

loc <- "/data/tide/projects/ho-covid/annals_resubmission/add_addtl_vars_dr5_sdm/sensitivity_analyses/vax_status"

dataLoc <- paste0(loc,"/data")

covid_data_simulation <- readRDS(file.path(dataLoc,paste0("covid_data_simulation_",filesfx,".Rds")))

matchedF <- paste0("covid_data_simulation_matched_",filesfx,".Rds")

covid_data_simulation_matched <- readRDS(file.path(dataLoc,matchedF)) %>% 
  mutate(ttmatch_days= as.numeric((date-admitDT)+1),
         matchtodischarge_days= as.numeric(dischargeDT-date))

# N matches
nrow(distinct(covid_data_simulation_matched,matchid))

covid_data_simulation_matched2 <- covid_data_simulation_matched %>% 
  select(hospID,matchid,matchday=hospitalDay) %>% 
  left_join(.,covid_data_simulation,by=c("hospID"),relationship = "many-to-many") %>% 
  arrange(hospID,hospitalDay)


new_event_after_match <- function(event) {
  
  matched_new_event <- covid_data_simulation_matched2 %>% 
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
readmit30day <- covid_data_simulation_matched %>% 
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
death30day <- covid_data_simulation_matched %>% 
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
hospitalfree30day <- covid_data_simulation_matched %>%
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


new_events <- covid_data_simulation_matched %>% distinct(hospID,matchid,case,cov_grp,LOS) %>% 
  left_join(.,newICU,by=c("hospID","matchid")) %>% 
  left_join(.,newHF,by=c("hospID","matchid")) %>% 
  left_join(.,newBIPAP,by=c("hospID","matchid")) %>% 
  left_join(.,newVent,by=c("hospID","matchid")) %>%
  left_join(.,newO2,by=c("hospID","matchid")) %>%
  left_join(.,readmit30day,by=c("hospID","matchid")) %>% 
  left_join(.,death30day,by=c("hospID","matchid")) %>%
  left_join(.,hospitalfree30day,by=c("hospID","matchid")) %>% 
  # if patient does not have new event, set _T to LOS
  mutate_at(vars(ends_with("_T")), ~coalesce(., LOS)) %>% 
  mutate_at(vars(starts_with("after_match_new"),c("readmit30day","death30day")), ~replace_na(., 0)) %>% 
  select(-LOS)




saveRDS(new_events,paste0(file.path(dataLoc,
                                    paste0("ho_covid_after_match_new_events_",outfx,".Rds"))))
