library(here)

## Setup

source(here("scripts", "00-config.R"))
source(here("scripts", "00-functions.R"))


# Read in data ------------------------------------------------------------

# Get a list of all .rds files in the folder
rds_files <- list.files(path = here("data"), pattern = "\\.rds$", full.names = TRUE)

# Loop through each .rds file and read it into a separate data table
for (file in rds_files) {
  # Create a variable name based on the file name (without extension)
  var_name <- tools::file_path_sans_ext(basename(file))
  # Read the .rds file into a data.table and assign it to a variable
  assign(paste0(var_name, "_raw"), readRDS(file), envir = .GlobalEnv)
}

# Basic data cleaning -----------------------------------------------------

## Assign appropriate column types

all_encounter_diagnoses <- all_encounter_diagnoses_raw %>%
  type_convert(na = c("", "NA", "NULL"), trim_ws = TRUE, guess_integer = TRUE) %>%
  mutate(across(c(present_on_admit, severity_of_dx, co_morbidity, cc_or_mcc), as.factor))

# If admit_time is missing and ed_dispo is "Admit", then using ed_disp_time as admit_time
# Otherwise using encounter_time as admit_time

demographics_and_event_table <- demographics_and_event_table_raw %>%
  mutate(across(ends_with("_time"), date_fn)) %>%
  mutate(admit_time = case_when(is.na(admit_time) & ED_dispo == "Admit" ~ ED_disp_time,
                                is.na(admit_time) ~ encounter_start_time,
                                TRUE ~ admit_time)) %>%
  rename(race = pat_race, LOS = hospital_LOS, covid_pos_date = covid_pos_time,
         admit_date = admit_time, discharge_date = discharge_time, death_date = death_time) %>%
  dplyr::filter(deid_enc_id %in% encounter_ids_raw) %>% 
  mutate(discharge_date = case_when(is.na(death_date) ~ discharge_date, TRUE ~ death_date),
         LOS = as.integer(LOS))

med_admin_table <- med_admin_table_raw %>% 
  type_convert(na = c("", "NA", "NULL"), trim_ws = TRUE, guess_integer = TRUE) %>%
  mutate(across(c(vasopressor_flag, COVID_tx_flag), as.factor))

flowsheet_data_table <- flowsheet_data_table_raw %>% 
  type_convert(na = c("", "NA", "NULL"), trim_ws = TRUE, guess_integer = TRUE)

covid_vaccinations <- covid_vaccinations_raw %>% 
  mutate(immune_date = date_fn(immune_date), dose_num = as.numeric(dose_num)) %>%
  dplyr::select(-manufacturer)

diag_classification <- primary_diagnoses_SB_raw %>% dplyr::select(-count) %>%
  rename(cov_diag_strict = Strict, cov_diag_medium = Medium)

## Check missingness

colSums(is.na(all_encounter_diagnoses))
colSums(is.na(demographics_and_event_table))
colSums(is.na(med_admin_table))
colSums(is.na(flowsheet_data_table))
colSums(is.na(covid_vaccinations))

# Assign COVID severity stage to each patient-day -------------------------

# Create `dm` table (reformat and select necessary variables from `demographics_and_event_table`)
# A dataset where each row represent a single patient-day
## Parse date from admit/discharge times and create new variable that spans from admission to discharge for each encounter

dm <- demographics_and_event_table %>%
  group_by(deid_enc_id) %>% 
  mutate(day_date = list(seq(from = min(admit_date), to = max(discharge_date), by = "day"))) %>% 
  unnest(day_date) %>% relocate(day_date, .after = 2) %>%  
  rename(date = day_date) %>% group_by(deid_enc_id) %>% 
  mutate(DaysSinceEntry = difftime(date, min(date), units = "days"),
         DaysSincePos = difftime(date, covid_pos_date, units = "days")) %>% ungroup() 

# Create `med` table (reformat and select necessary variables from `med_admin_table`)

med <- med_admin_table %>% 
  mutate(date = as.Date(taken_time, format = "%Y-%m-%d"),
         infusion_rate = as.numeric(infusion_rate)) %>% 
  dplyr::select(deid_enc_id, taken_time, date, vasopressor_flag, infusion_rate)

# Create `fs` table (reformat and select necessary variables from `flowsheet_data_table`)

fs <- flowsheet_data_table %>% 
  dplyr::filter(clinical_concept %nin% c("Respirations", "Urine_output", "Blood_pressure", "Pulse", "Temperature")) %>%
  mutate(date = as.Date(recorded_time, format = "%Y-%m-%d"), 
         meas_val_chr = meas_value,                 # create new var `meas_val_chr` to keep meas_value as character type
         meas_val_num = as.numeric(meas_value)) %>%   
  dplyr::select(deid_enc_id, flo_meas_name, meas_val_num, meas_val_chr, recorded_time, date, clinical_concept)

# Create table to indicate whether each patient received vasopressor each day
# -   Create new variable `vp_yn` to indicate whether patient received VP at timestamp
# -   If `vasopressor_flag`=Yes & `infusion_rate`\>0, then `vp_yn`=1
# -   Group by ID and date, then create variable `VP` to indicate whether pt received VP at any timestamp during day (take max value of `vp_yn`)
# -   Remove duplicate rows

vp.tab <- med %>% 
  mutate(vp_yn = ifelse((vasopressor_flag == "Yes" & infusion_rate > 0), 1, 0)) %>%
  group_by(deid_enc_id, date) %>%
  mutate(VP = max(vp_yn, na.rm = T)) %>%
  distinct(deid_enc_id, date, VP)

# Create table to indicate whether each patient received ECMO each day
# -   Create new variable `ecmo_yn` to indicate whether patient received ECMO at timestamp
# -   If `clinical_concept`='ECMO_settings', then `ecmo_yn`=1
# -   Group by ID and date, then create variable `ECMO` to indicate whether pt received ECMO at any timestamp during day (take max value of `ecmo_yn`)
# -   Remove duplicate rows

ecmo.tab <- fs %>% 
  mutate(ecmo_yn = ifelse(clinical_concept == "ECMO_settings", 1, 0)) %>%
  group_by(deid_enc_id, date) %>%
  mutate(ECMO = max(ecmo_yn, na.rm = T)) %>%
  distinct(deid_enc_id, date, ECMO)

# Create table to indicate whether patient received CRRT each day
# -   Create new variable `crrt_yn` to indicate whether patient received CRRT at timestamp
# -   If `clinical_concept`="HD_UF_CRRT_settings" & `flo_meas_name` %in% c("R UFR TRANSCRIBED CRRT IP_CD_UCSF", "R CRRT BLOOD FLOW RATE"), then `crrt_yn`=1
# -   Group by ID and date, then create variable `CRRT` to indicate whether pt received CRRT at any timestamp during day (take max value of `crrt_yn`)
# -   Remove duplicate rows

crrt.tab <- fs %>% 
  mutate(crrt_yn = ifelse((clinical_concept == "HD_UF_CRRT_settings" &
                             flo_meas_name %in% c("R UFR TRANSCRIBED CRRT IP_CD_UCSF", "R CRRT BLOOD FLOW RATE")), 1, 0)) %>%  
  group_by(deid_enc_id, date) %>% 
  mutate(CRRT = max(crrt_yn, na.rm = T)) %>%
  distinct(deid_enc_id, date, CRRT)

# Create table to indicate whether patient received support from NIV device each day
# -   Create new variable `niv_yn` to indicate whether patient received support from an NIV device at timestamp
# -   If `clinical_concept`="NIV_settings", then `niv_yn`=1
# -   Group by ID and date, then create variable `NIV_PER_DAY` to indicate number of times that an NIV device was recorded during day
# -   Remove duplicate rows

niv.tab <- fs %>% 
  mutate(niv_yn = ifelse(clinical_concept == "NIV_settings", 1, 0)) %>% 
  group_by(deid_enc_id, date) %>% 
  mutate(NIV = max(niv_yn),
         NIV_PER_DAY = sum(clinical_concept == "NIV_settings")) %>% 
  distinct(deid_enc_id, date, NIV, NIV_PER_DAY)

# Create table to indicate whether patient received HD each day
# -   Create new variable `hd_yn` to indicate whether patient received support from an HD device or HD settings at timestamp
# -   If `clinical_concept`="HD_UF_CRRT_settings" & `flo_meas_name` %in% c("R HD BLOOD FLOW RATE", "R HD ULTRAFILTRATION RATE"), then `hd_yn`=1
# -   Group by ID and date, then create variable `HD` to indicate whether pt received HD support at any timestamp during day (take max value of `hd_yn`)
# -   Remove duplicate rows

hd.tab <- fs %>% 
  mutate(hd_yn = ifelse((clinical_concept == "HD_UF_CRRT_settings" & 
                           flo_meas_name %in% c("R HD BLOOD FLOW RATE", "R HD ULTRAFILTRATION RATE")), 1, 0)) %>%
  group_by(deid_enc_id, date) %>% 
  mutate(HD = max(hd_yn, na.rm = T)) %>% 
  distinct(deid_enc_id, date, HD)

# Create table to indicate whether patient was intubated each day
# -   Create new variable `intub_yn` to indicate whether patient was intubated at timestamp
# -   If `clinical_concept`="Intubation_settings" , then `intub_yn`=1
# -   Group by ID and date, then create variable `INTUB` to indicate whether pt was intubated at any timestamp during day (take max value of `intub_yn`)
# -   Remove duplicate rows

intub.tab <- fs %>% mutate(intub_yn = ifelse(clinical_concept == "Intubation_settings", 1, 0)) %>% 
  group_by(deid_enc_id, date) %>% mutate(INTUB = max(intub_yn, na.rm = T)) %>% distinct(deid_enc_id, date, INTUB)

# Filter data for SpO2 and FiO2 and merge tables to calculate SF
# -   Filter data for SpO2 and FiO2
# -   Merge tables and calculate the pulse oximetric saturation SpO2/FiO2 (SF) ratio, a non-invasive surrogate measure for PF ratio
# -   `SF` = (SpO2/FiO2) \* 100

# Filter data for SpO2 values and assign to table `spo2.tab`
spo2.tab <- fs %>% filter(clinical_concept == "SpO2") %>% 
  dplyr::select(deid_enc_id, recorded_time, date, meas_val_num)

# Filter data for FiO2 values and assign to table `fio2.tab`
fio2.tab <- fs %>% filter(flo_meas_name == "R FIO2") %>% 
  dplyr::select(deid_enc_id, recorded_time, date, meas_val_num)

# Merge SpO2 and FiO2 tables and calculate SF 
SF.tab <- fio2.tab %>% left_join(spo2.tab, by = c("deid_enc_id", "recorded_time"), suffix = c(".fi", ".sp")) %>% 
  mutate(SF = (meas_val_num.sp / meas_val_num.fi * 100)) %>% dplyr::select(deid_enc_id, recorded_time, date.fi, SF) %>% 
  rename(date = date.fi)

# Create table to indicate whether patient experienced SpO2/FiO2\<200 at same timestamp, and assign to new variable `SF_LT_200`
# -   Using `SF.tab` created above, filter for patients that have a SF ratio of less than 200
# -   Create new variable `sf_lt_200_yn` to indicate whether patient had SF ratio \< 200 at timestamp
# -   Group by ID and date, then create variable `SF_LT_200` to indicate whether pt was intubated at any timestamp during day (take max value of `sf_lt_200_yn`)
# -   Remove duplicate rows

SF200.tab <- SF.tab %>% mutate(sf_lt_200_yn = case_when(SF < 200 ~ 1,
                                                        is.na(SF) ~ 0,
                                                        TRUE ~ 0)) %>% 
  group_by(deid_enc_id, date) %>% mutate(SF_LT_200 = max(sf_lt_200_yn, na.rm = T)) %>% 
  distinct(deid_enc_id, date, SF_LT_200)

# Assign all O2 devices to respiratory support categories

# Read in `o2_dev_names.csv` which lists all O2 devices and their respective O2 support category. Categories include:
# -   NONE: `dev_cat`="no_dev"
# -   SIMPLE: `dev_cat`="simple_dev"
# -   Non-Invasive Ventilation (NIV): `dev_cat`="NIV_dev"
# -   Invasive Ventilation (IV): `dev_cat`="iv_dev"
# -   CPAP: `dev_cat`="cpap_dev"
# -   NC (Nasal cannula and high-flow nasal cannula): `nc_dev`=1

o2_dev_names <- o2_dev_names_raw
no_dev_names <- filter(o2_dev_names, dev_cat == "no_dev")[[1]]
simple_dev_names <- filter(o2_dev_names, dev_cat == "simple_dev")[[1]]
niv_dev_names <- filter(o2_dev_names, dev_cat == "NIV_dev")[[1]]
iv_dev_names <- filter(o2_dev_names, dev_cat == "iv_dev")[[1]]
cpap_dev_names <- filter(o2_dev_names, dev_cat == "cpap_dev")[[1]]
nc_dev_names <- filter(o2_dev_names, nc_dev == 1)[[1]]

# Create table of O2 devices categorized by respiratory support type to indicate whether patient received respiratory support, and if yes, which type
# 
# Each column represents a respiratory support category (as coded in previous step): 
# `NONE`, `SIMPLE`, `NIV`, `IV`, `CPAP`, and `NC`. You can think of these as different types of O2 devices, varying in severity.

device.tab <- fs %>% filter(clinical_concept == "O2_device") %>%
  mutate(NOdev_yn = ifelse(meas_val_chr %in% no_dev_names, 1, 0),
         SIMPLEdev_yn = ifelse(meas_val_chr %in% simple_dev_names, 1, 0),
         NIVdev_yn = ifelse(meas_val_chr %in% niv_dev_names, 1, 0),
         IVdev_yn = ifelse(meas_val_chr %in% iv_dev_names, 1, 0),
         CPAPdev_yn = ifelse(meas_val_chr %in% cpap_dev_names, 1, 0),
         NCdev_yn = ifelse(meas_val_chr %in% nc_dev_names, 1, 0)) %>% 
  group_by(deid_enc_id, date) %>% 
  mutate(NODEV = max(NOdev_yn, na.rm = T),
         SIMPLEDEV = max(SIMPLEdev_yn, na.rm = T),
         SIMPLE_PER_DAY = sum(meas_val_chr %in% simple_dev_names),
         NIVDEV = max(NIVdev_yn, na.rm = T),
         IVDEV = max(IVdev_yn, na.rm = T),
         CPAPDEV = max(CPAPdev_yn, na.rm = T),
         NCDEV = max(NCdev_yn, na.rm = T),
         NC_PER_DAY = sum(meas_val_chr %in% nc_dev_names)) %>% 
  distinct(deid_enc_id, date, NODEV, SIMPLEDEV, SIMPLE_PER_DAY, NIVDEV, IVDEV, CPAPDEV, NCDEV, NC_PER_DAY)

# Create table to indicate whether patients received any O2, low O2, or high O2
# -   Categorize O2 usage into any, low, and high:
# -   **O2:** If `meas_val_num`\>0, then `O2`=1
# -   **LowO2:** If `meas_val_num`\>0 & `meas_val_num`\<=12, then `LowO2`=1
# -   **HighO2:** If `meas_val_num`\>12, then `HighO2`=1

# Create vector for O2 flow rate names
O2_names <- c("R OXYGEN FLOW RATE", "R OXYGEN FLOW RATE 2 IP_CD_UCSF")

# Create table with columns 'O2', 'LowO2', and 'HighO2'
O2.tab <- fs %>% 
  mutate(O2 = ifelse(flo_meas_name %in% O2_names & meas_val_num > 0, 1, 0),
         LowO2 = ifelse(flo_meas_name %in% O2_names & meas_val_num > 0 & meas_val_num <= 12, 1, 0),
         HighO2 = ifelse(flo_meas_name %in% O2_names & meas_val_num > 12, 1, 0)) %>% group_by(deid_enc_id, date) %>% 
  mutate(O2 = max(O2, na.rm = T),
         LowO2 = max(LowO2, na.rm = T), 
         HighO2 = max(HighO2, na.rm = T)) %>% distinct(deid_enc_id, date, O2, LowO2, HighO2)

# Merge `device.tab` and `O2.tab` to create table to indicate whether patient received support from NC and high O2
# -   Merge `device.tab` and `O2.tab` to create `O2_DEV.tab`
# -   Create variable `NC_GT_12` to indicate whether patient received support from nasal cannula device AND had high O2 flow rate
# -   If `NCDEV`=1 & `HighO2`=1, then `NC_GT_12`=1

# Left_join on 'ID' and 'date' columns
O2_DEV.tab <- O2.tab %>%
  left_join(device.tab, by = c("deid_enc_id", "date")) %>% 
  mutate(NC_GT_12 = ifelse(NCDEV == 1 & HighO2 == 1, 1, 0))

# Merge all tables to aggregate binary clinical scenarios and O2 device classifications

bin_clin_scen_df <- dm %>% 
  left_join(fs, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>% 
  left_join(vp.tab, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>%
  left_join(ecmo.tab, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>%
  left_join(crrt.tab, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>%
  left_join(niv.tab, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>%
  left_join(hd.tab, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>%
  left_join(intub.tab, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>%
  left_join(SF200.tab, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>%
  left_join(O2_DEV.tab, by = c("deid_enc_id", "date"), suffix = c("", ".dupcol")) %>%
  dplyr::select(-ends_with(".dupcol")) %>% distinct() %>%
  mutate(death = if_else(date == death_date, 1, 0, missing = NA))

# Get COVID stages based on patients' binary clinical scenarios and usage of O2 devices (and their respiratory support category)

bin_clin_scen_df <- bin_clin_scen_df %>% 
  mutate(stage = case_when(death==1 ~ 10,
                           ECMO==1  ~ 9,
                           SF_LT_200==1 & (INTUB==1 | IVDEV==1) & (VP==1 | CRRT==1) ~ 9, 
                           VP==1 | CRRT==1 ~ 8,
                           SF_LT_200==1 & (INTUB==1 | IVDEV==1) ~ 8,
                           SF_LT_200==0 & (INTUB==1 | IVDEV==1) ~ 7,
                           (NIVDEV==1 & CPAPDEV==0) & IVDEV==0 & 
                             ((NIV_PER_DAY > NC_PER_DAY) | (NC_PER_DAY>0 & NC_GT_12==1)) ~ 6,
                           NC_PER_DAY>0 & NC_GT_12==1 ~ 6,
                           SIMPLEDEV==1 & (LowO2==1 | HighO2==1 | SIMPLE_PER_DAY>1) ~ 5,
                           NC_PER_DAY>0 & NC_GT_12==0 ~ 5,
                           LowO2==1 & SIMPLEDEV==0 & NIVDEV==0 & IVDEV==0 & CPAPDEV==0 & NCDEV==0 & 
                             INTUB==0 ~ 5,
                           TRUE ~ 4))

# Clean table to just indicate the stage of COVID-19 for each patient-day

bin_clin_scen_df <- bin_clin_scen_df %>% 
  relocate(c("DaysSinceEntry", "DaysSincePos", "stage", "death"), .after = "end_in_death") %>%
  dplyr::select(deid_enc_id, date, DaysSinceEntry, DaysSincePos, stage, death_date) %>% distinct() %>%
  mutate(DaysSinceEntry = as.numeric(DaysSinceEntry), DaysSincePos = as.numeric(DaysSincePos))

# Add stage 11 to denote "recovered" for those who leave the model before death

dm_covid_stg <- bin_clin_scen_df %>% group_by(deid_enc_id) %>%
  do ({                        
    if (max(.$stage) == 10) {    
      dm_covid_stg <- .
    } else {                 
      new <- tibble(deid_enc_id = .$deid_enc_id[1],
                    date = max(.$date) + 1,
                    stage = 11, 
                    DaysSinceEntry = max(.$DaysSinceEntry) + 1,
                    DaysSincePos = max(.$DaysSincePos) + 1,
                    death_date = NA)
      dm_covid_stg <- rbind(., new)  
    }
    dm_covid_stg
  })

rm(list = ls(pattern = "\\.tab$"))
rm(list = ls(pattern = "\\_names$"))

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
  
## COVID treatment and vaccination ----------------------------------------

# Use `med_admin_table` to create a table of COVID treatment

cov_tx <- med_admin_table %>% dplyr::select(deid_enc_id, taken_time, name, COVID_tx_flag, mar_action) %>%
  dplyr::filter(COVID_tx_flag == "Yes" & grepl("Given", mar_action)) %>%
  mutate(date = date_fn(taken_time)) %>% 
  mutate(COVID_tx_name = case_when(grepl("PLACEBO", name, fixed = TRUE) ~ NA,
                                   grepl("REMDESIVIR", name, fixed = TRUE) ~ "Remdesivir",
                                   grepl("PREDNISONE", name, fixed = TRUE) ~ "Prednisone",
                                   grepl("METHYLPREDNISOLONE", name, fixed = TRUE) ~ "Methylprednisolone",
                                   grepl("DEXAMETHASONE", name, fixed = TRUE) ~ "Dexamethasone",
                                   grepl("BAMLANIVIMAB", name, fixed = TRUE) ~ "Bamlanvimab",
                                   TRUE ~ NA)) %>% dplyr::select(deid_enc_id, date, COVID_tx_name) %>% distinct() %>% 
  distinct() %>% group_by(deid_enc_id, date) %>% mutate(n = n()) %>% rowwise() %>%
  mutate(COVID_tx_name = ifelse(n>1, "Multiple", COVID_tx_name)) %>%
  ungroup() %>% distinct() %>% dplyr::select(deid_enc_id, date, COVID_tx_name)

# Use `covid_vaccinations` to create a table of COVID vaccination status

cov_vax <- demographics_and_event_table %>% dplyr::select(deid_pat_id, deid_enc_id, admit_date) %>% distinct() %>% 
  left_join(covid_vaccinations, by = c("deid_pat_id"), relationship = "many-to-many") %>% 
  dplyr::filter(immune_date < admit_date) %>% group_by(deid_enc_id) %>% arrange(desc(dose_num)) %>% 
  slice_head(n = 1) %>% ungroup() %>% dplyr::select(deid_enc_id, dose_num) %>%
  rename(COVID_vax_doses = dose_num) %>% mutate(COVID_vax = TRUE)


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
                             cci_score %in% 3:5 ~ "Moderate/Severe",
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

dem <- dem %>% left_join(cov_vax, by = "deid_enc_id") %>% left_join(comorb, by = "deid_enc_id") %>% 
  left_join(dnr, by = "deid_enc_id") %>% left_join(cov_dx, by = "deid_enc_id") %>%
  replace_na(list(COVID_vax_doses = 0, COVID_vax = FALSE, copd = 0, DNR = FALSE, dnr_on_admit = "No",
                  dx_covid = FALSE, cov_diag_strict = FALSE, cov_diag_medium = FALSE))

dm_covid_stg <- dm_covid_stg %>% left_join(cov_tx, by = c("deid_enc_id", "date")) %>% dplyr::select(-death_date)

saveRDS(dem, here("data", "pt_demographics.rds"))
saveRDS(dm_covid_stg, here("data", "pt_enc_staging.rds"))
saveRDS(dm_cov_summ, here("data", "pt_stage_summary.rds"))
