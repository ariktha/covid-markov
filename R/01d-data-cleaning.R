library(here)
library(tidyverse)

rm(list = ls())

# Read data ---------------------------------------------------------------

load(here("data", "preprocess", "post-c-data.RData"))

# More data cleaning and processing ---------------------------------------

## Demographics -----------------------------------------------------------

# Use `demographics_and_event_table` to create a table of demographics
# Collapse levels for factor variables and set reference levels

dem <- demographics_and_event_table %>%
  dplyr::select(deid_enc_id, deid_pat_id, age, sex, race, ethnicity, language, insurance_type, smoking,
                BMI, end_in_death, LOS) %>%
  mutate(bmi_cat = case_when(BMI < 18.5 ~ "Underweight", BMI < 25 ~ "Healthy weight", BMI < 30 ~ "Overweight",
                             BMI > 30 ~ "Obese", TRUE ~ NA),
         language = case_when(language %in% c("English", "Spanish") ~ language, TRUE ~ "Other"),
         race = case_when(race %in% c("Unknown", "Declined", "Unknown/Declined", "NULL", "American Indian or Alaska Native", "Other") ~ "Others/Unknown",
                          race %in% c("Native Hawaiian", "Native Hawaiian or Other Pacific Islander", "Other Pacific Islander") ~ "Pacific Islanders",
                          TRUE ~ race),
         ethnicity = case_when(ethnicity %in% c("Unknown", "Declined", "Unknown/Declined", "NULL") ~ "Unknown",
                               TRUE ~ ethnicity),
         smoking = case_when(smoking %in% c("Current Every Day Smoker", "Current Some Day Smoker", "Heavy Tobacco Smoker", 
                                            "Light Tobacco Smoker", "Smoker, Current Status Unknown") ~ "Current Smoker",
                             smoking %in% c("Former Smoker") ~ "Former Smoker",
                             smoking %in% c("Never Smoker", "Passive Smoke Exposure - Never Smoker") ~ "Non Smoker",
                             TRUE ~ "Unknown"),
         insurance_type = case_when(insurance_type == "01 - Medicare" ~ "Medicare",
                                    insurance_type == "02 - Medi-Cal" ~ "Medi-Cal",
                                    insurance_type == "03 - Private Coverage" ~ "Private Coverage",
                                    TRUE ~ "Other")) %>%
  mutate(across(c(sex, bmi_cat, language, race, ethnicity, smoking, insurance_type, end_in_death), as.factor)) %>%
  mutate(across(c(age, LOS, BMI), as.numeric))

colSums(is.na(dem))

age_quintiles <- quantile(dem$age, probs = seq(0, 1, 0.2))

dem <- dem %>%
  mutate(
    age_cat = case_when(
      age <= age_quintiles[2] ~ "Q1",
      age > age_quintiles[2] & age <= age_quintiles[3] ~ "Q2",
      age > age_quintiles[3] & age <= age_quintiles[4] ~ "Q3",
      age > age_quintiles[4] & age <= age_quintiles[5] ~ "Q4",
      age > age_quintiles[5] ~ "Q5"
    ),
    age_cat = factor(age_cat, levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))
  )

## COVID treatment and vaccination ----------------------------------------

# Read list of meds that count as corticosteroids
# List is from Daniel over email on 18 Feb 2026

corticosteroids <- read_csv(here("data", "corticosteroids.csv"))

# Use `med_admin_table` to create a table of COVID treatment

tx_steroids <- med_admin_table %>% 
  dplyr::filter(name %in% corticosteroids$meds & mar_action %in% c("Given", "New Bag")) %>%
  mutate(date = date_fn(taken_time),
         corticosteroids = TRUE) %>%
  dplyr::select(deid_enc_id, date, corticosteroids) %>%
  distinct()

tx_remdesivir <- med_admin_table %>% 
  dplyr::filter(grepl("REMDESIVIR", name, fixed = TRUE) & mar_action %in% c("Given", "New Bag")) %>%
  mutate(date = date_fn(taken_time),
         remdesivir = TRUE) %>%
  dplyr::select(deid_enc_id, date, remdesivir) %>%
  distinct()

cov_tx <- full_join(tx_steroids, tx_remdesivir, by = c("deid_enc_id", "date")) %>%
  replace_na(list(corticosteroids = FALSE, remdesivir = FALSE)) %>%
  mutate(COVID_tx = case_when(corticosteroids & remdesivir ~ "Both",
                              corticosteroids & !remdesivir ~ "Corticosteroids Only",
                              !corticosteroids & remdesivir ~ "Remdesivir Only",
                              TRUE ~ "No Treatment")) %>%
  dplyr::select(deid_enc_id, date, COVID_tx)

table(cov_tx$COVID_tx)

saveRDS(cov_tx, here("data", "pt_covid_tx.rds"))

# Use `covid_vaccinations` to create a table of COVID vaccination status

cov_vax <- demographics_and_event_table %>% 
  dplyr::select(deid_pat_id, deid_enc_id) %>% 
  distinct() %>% 
  left_join(covid_vaccinations, by = c("deid_pat_id"), relationship = "many-to-many") %>% 
  drop_na() %>%
  group_by(deid_pat_id) %>% arrange(immune_date) %>%
  slice_head(n = 1) %>% ungroup() %>%
  dplyr::select(deid_enc_id, immune_date) %>%
  right_join(dm_covid_stg %>% dplyr::select(deid_enc_id, date) %>% distinct(), by = "deid_enc_id") %>%
  filter(immune_date <= date) %>%
  dplyr::select(deid_enc_id, date) %>%
  mutate(vax = 1)

## Comorbidities ----------------------------------------------------------

# Use `all_encounter_diagnoses` to create a table of comorbidities

comorb <- comorbidity::comorbidity(all_encounter_diagnoses, id = "deid_enc_id", code = "current_icd10_list", 
                                   map = "charlson_icd10_quan", assign0 = FALSE)
comorb$cci_score <- as.numeric(comorbidity::score(comorb, weights = NULL, assign0 = FALSE))
comorb$cci_fct <- factor(comorb$cci_score, ordered = T, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

comorb <- comorb %>% dplyr::select(deid_enc_id, cci_score, cci_fct, chf, cpd) %>%
  mutate(chf = if_else(chf == 1, TRUE, FALSE),
         cpd = if_else(cpd == 1, TRUE, FALSE)) %>%
  mutate(across(c(chf, cpd), as.factor)) %>%
  mutate(cci_cat = case_when(cci_score == 0 ~ "None",
                             cci_score %in% 1:2 ~ "Mild",
                             cci_score %in% 3:7 ~ "Moderate/Severe",
                             TRUE ~ as.character(NA))) %>%
  mutate(cci_cat = factor(cci_cat, levels = c("None", "Mild", "Moderate/Severe"), ordered = TRUE))

### COPD ------------------------------------------------------------

## ICD-10 (codes from J41 to J44) considered COPD 
### source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6805625/
### “In 34 of 38 studies … ICD-10 (codes from J41 to J44) coding were used as one part of the identification process”
#### J41: Simple and mucopurulent chronic bronchitis
#### J42: Unspecified chronic bronchitis
#### J43: Emphysema
#### J44: Other chronic obstructive pulmonary disease
#### -	J44.0: Chronic obstructive pulmonary disease with (acute) lower respiratory infection
#### -	J44.1: Chronic obstructive pulmonary disease with (acute) exacerbation
#### -	J44.9: Chronic obstructive pulmonary disease, unspecified

copd_codes <- c("J41", "J42", "J43", "J44")

# cpd in the comorb table is how the comorbidity package codes chronic pulmonary disease

all_copd <- all_encounter_diagnoses %>% 
  mutate(copd = ifelse(substring(current_icd10_list, 1, 3) %in% copd_codes, 1, 0)) %>%
  dplyr::select(deid_enc_id, copd) %>% distinct() %>% group_by(deid_enc_id) %>%
  summarise(copd_sum = sum(copd)) %>% mutate(copd = ifelse(copd_sum, 1, 0)) %>% dplyr::select(deid_enc_id, copd)

comorb <- left_join(comorb, all_copd, by = "deid_enc_id")

# COPD is a form of CPD, so presence of CPD doesn't necessarily indicate COPD; but absence of CPD indicates absence of COPD
# We are more interested in COPD, so we will use the presence of COPD as the primary indicator
table(comorb$copd, comorb$cpd)
comorb <- comorb %>% dplyr::select(-cpd)

### DNR status ------------------------------------------------------------

dnr <- all_encounter_diagnoses %>% dplyr::filter(dx_group == "DO NOT RESUSCITATE STATUS") %>% 
  dplyr::select(deid_enc_id, present_on_admit) %>% 
  mutate(dnr_on_admit = factor(as.character(present_on_admit)),
         DNR = TRUE) %>% dplyr::select(deid_enc_id, DNR, dnr_on_admit)

## COVID-19 diagnosis -----------------------------------------------------

cov_dx <- all_encounter_diagnoses %>% 
  left_join(diag_classification, join_by(DX_NAME == `Primary Diagnosis`)) %>%
  mutate(dx_covid = ifelse(DX_NAME %in% c("COVID-19", "Pneumonia due to coronavirus disease 2019"), 1, 0)) %>%
  mutate(cov_diag_strict = ifelse(dx_covid == 1 | cov_diag_strict == 1, 1, 0),
         cov_diag_medium = ifelse(dx_covid == 1 | cov_diag_medium == 1, 1, 0)) %>%
  replace_na(list(cov_diag_strict = 0, cov_diag_medium = 0)) %>%
  dplyr::select(c(deid_enc_id, dx_covid, cov_diag_strict, cov_diag_medium)) %>%
  group_by(deid_enc_id) %>%
  summarise(across(c(dx_covid, cov_diag_strict, cov_diag_medium), sum), .groups = "keep") %>%
  ungroup() %>% mutate(across(c(dx_covid, cov_diag_strict, cov_diag_medium), as.logical))

table(cov_dx$dx_covid)
table(cov_dx$dx_covid, cov_dx$cov_diag_strict)
table(cov_dx$dx_covid, cov_dx$cov_diag_medium)

## Summaries of pt-day stages data -----------------------------------------

dm_cov_summ <- dm_covid_stg %>% arrange(desc(DaysSinceEntry)) %>% 
  slice_head(n = 1) %>% ungroup() %>% 
  mutate(abs_state = ifelse(stage == 10, "Death", "Recovered"), LOS = DaysSinceEntry-1) %>% 
  dplyr::select(deid_enc_id, abs_state, LOS)

dm_cov_stg_summ <- dm_covid_stg %>%
  mutate(state = case_when(stage %in% c(4:6) ~ 1,
                           stage %in% c(7:9) ~ 2,
                           stage == 10 ~ 3)) %>% 
  group_by(deid_enc_id, state) %>% 
  summarise(DaysInState = n(), .groups = "keep") %>% ungroup() %>%
  dplyr::filter(state == 2) %>% rename(DaysInSevere = DaysInState) %>% dplyr::select(-state) %>%
  mutate(EverSevere = TRUE)

dm_cov_summ <- dm_cov_summ %>% left_join(dm_cov_stg_summ, by = "deid_enc_id") %>%
  replace_na(list(EverSevere = FALSE, DaysInSevere = 0))

# Merge and save tables ---------------------------------------------------

dem <- dem %>% 
  left_join(comorb, by = "deid_enc_id") %>% 
  left_join(dnr, by = "deid_enc_id") %>% 
  left_join(cov_dx, by = "deid_enc_id") %>%
  replace_na(list(copd = 0, DNR = FALSE, dnr_on_admit = "No",
                  dx_covid = FALSE, cov_diag_strict = FALSE, cov_diag_medium = FALSE))

colSums(is.na(dem))

library(mice)

mice_imputed <- mice(dem, formulas = list(BMI ~ sex + age + cci_score + race + ethnicity)) 
dem_imp <- complete(mice_imputed)

dm_cov_stg <- dm_covid_stg %>% 
  ungroup() %>%
  left_join(cov_vax, by = c("deid_enc_id", "date")) %>% 
  left_join(cov_tx, by = c("deid_enc_id", "date")) %>% 
  dplyr::select(-death_date) %>%
  replace_na(list(vax = 0, COVID_tx = "No Treatment")) %>%
  mutate(COVID_tx = factor(COVID_tx)) %>%
  mutate(COVID_tx = relevel(COVID_tx, ref = "No Treatment"))

saveRDS(dem_imp, here("data", "pt_demographics.rds"))
saveRDS(dm_cov_stg, here("data", "pt_stage_data.rds"))
saveRDS(dm_cov_summ, here("data", "pt_stage_summary.rds"))
