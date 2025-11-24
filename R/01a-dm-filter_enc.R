library(here)

## Setup

source(here("R", "00-config.R"))
source(here("R", "00a-fns-data_preprocess.R"))

## Read raw data files

dm_raw <- read_delim(here("data", "raw", "demographics_and_event_table_4-14-22.csv"), na = c("", "NA", "NULL"), quote = "", delim = "|")

# Check for missing data
colSums(is.na(dm_raw))

# Check rows where age is NA -- 2 rows and it's not real data; removing them
dm_raw[which(is.na(dm_raw$age)),]
dm_raw <- dm_raw %>% drop_na(age)

# 54404 encounters (each row corresponds to a unique encounter); 38433 patients
nrow(dm_raw)
length(unique(dm_raw$deid_enc_id))
length(unique(dm_raw$deid_pat_id))

# 81.2% of encounters (44175) are inpatient -- Only keeping these
# Remove encounters with LOS < 1 day
tbl_fn(dm_raw$enc_class)
dm_inp <- dm_raw %>% filter(enc_class == "Inpatient") %>% 
  mutate(across(ends_with("_time"), date_fn)) %>%
  mutate(admit_time = case_when(is.na(admit_time) & ED_dispo == "Admit" ~ ED_disp_time,
                                is.na(admit_time) ~ encounter_start_time,
                                TRUE ~ admit_time)) %>% 
  mutate(los = difftime(discharge_time, admit_time, units = "days")) %>%
  dplyr::filter(los > 1)

colSums(is.na(dm_inp))

# 52.05% of encounters have a COVID test; 4.24% of those tests are positive
# Somehow, 0.81% of encounters without a COVID test are positive
# Probably because the same COVID test has been assigned to all encounters of a patient
tbl_fn(dm_inp$covid_tested)
tbl_fn(dm_inp[which(dm_inp$covid_tested == "Yes"),]$covid_pos)
tbl_fn(dm_inp[which(dm_inp$covid_tested == "No"),]$covid_pos)

# Format columns needed for filtering
dm_inp <- dm_inp %>% 
  dplyr::select(starts_with("deid"), "admit_time") %>%
  mutate(across(ends_with("_time"), date_fn))

colSums(is.na(dm_inp))

# Only keeping encounters with a positive COVID test
covid_pos_dates <- dm_raw %>% mutate(covid_pos_time = date_fn(covid_pos_time)) %>%
  filter(!is.na(covid_pos_time)) %>% 
  dplyr::select(deid_pat_id, covid_pos_time) %>% distinct()

# Check admit_time -- 0 rows with NA; replacing with encounter_start_time
## In cases where both are available, for most cases the difference is 0 but the range is from -1 to 15 days
colSums(is.na(dm_inp))

dm_inp <- dm_inp %>%
  right_join(covid_pos_dates, by = "deid_pat_id", relationship = "many-to-many") %>%
  mutate(pos_adm = as.numeric(covid_pos_time - admit_time)) %>%
  drop_na(deid_enc_id)

# Look at the distribution of days from admission to positive COVID test
temp <- dm_inp %>% dplyr::filter(abs(pos_adm) <= 15)
hist(temp$pos_adm, breaks = 31, main = "Days from admission to positive COVID test", xlab = "Days")

# Filter based on the days from admission to positive COVID test
dm_inp <- dm_inp %>% dplyr::filter(pos_adm > -14 & pos_adm < 7) %>%
  group_by(deid_enc_id) %>% arrange(abs(pos_adm)) %>% slice_head(n = 1) %>% ungroup() %>%
  group_by(deid_pat_id, covid_pos_time) %>% arrange(abs(pos_adm)) %>% slice_head(n = 1) %>% ungroup()

length(unique(dm_inp$deid_enc_id))
length(unique(dm_inp$deid_pat_id))

# Since there are more unique encounters than patients, there are multiple encounters per patient
# Checking the time between encounters -- Keeping both encounters if the difference is more than 30 days
temp <- dm_inp %>% group_by(deid_pat_id) %>% mutate(n = n()) %>% filter(n > 1) %>% arrange(admit_time) %>%
  mutate(diff = difftime(admit_time, lag(admit_time), units = "days"))

print(temp$diff)

# Removing the second encounter if the difference is less than 30 days
dm_inp <- dm_inp %>% group_by(deid_pat_id) %>% mutate(diff_encs = as.numeric(difftime(admit_time, lag(admit_time), units = "days"))) %>% 
  mutate(diff_encs = if_else(is.na(diff_encs), 0, diff_encs)) %>%
  dplyr::filter(!(diff_encs > 0 & diff_encs < 30)) %>% ungroup()

# Distribution of days from admission to positive COVID test 
summary(dm_inp$pos_adm)
boxplot(dm_inp$pos_adm, ylab = "Days", main = "Days from admission to positive COVID test")

enc_ids <- dm_inp$deid_enc_id
saveRDS(enc_ids, here("data", "preprocess", "encounter_ids.rds"))

pat_ids <- unique(dm_inp$deid_pat_id)
saveRDS(pat_ids, here("data", "preprocess", "patient_ids.rds"))
