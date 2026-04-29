# Setup -----

library(here)

source(here("R", "00-config.R"))
source(here("R", "00a-fns-data_preprocess.R"))

# Load raw data ----

dm_raw <- read_delim(
  here("data", "raw", "demographics_and_event_table_4-14-22.csv"),
  na = c("", "NA", "NULL"),
  quote = "",
  delim = "|"
)

# Check for missing data
colSums(is.na(dm_raw))

# Check rows where age is NA -- 2 rows and it's not real data; removing them
dm_raw[which(is.na(dm_raw$age)), ]
dm_raw <- dm_raw %>% drop_na(age)

# 54404 encounters (each row corresponds to a unique encounter); 38433 patients
nrow(dm_raw)
length(unique(dm_raw$deid_enc_id))
length(unique(dm_raw$deid_pat_id))

# Filter to inpatient encounters with LOS > 1 day ----

# 81.2% of encounters (44175) are inpatient -- Only keeping these
tbl_fn(dm_raw$enc_class)
dm_inp <- dm_raw %>% filter(enc_class == "Inpatient")

# Estimating LOS from admit_time and discharge_time
dm_inp <- dm_inp %>%
  mutate(across(ends_with("_time"), date_fn)) %>%
  mutate(
    admit_time = case_when(
      is.na(admit_time) & ED_dispo == "Admit" ~ ED_disp_time,
      is.na(admit_time) ~ encounter_start_time,
      TRUE ~ admit_time
    )
  ) %>%
  mutate(los = as.numeric(difftime(discharge_time, admit_time, units = "days")))

summary(dm_inp$los)
boxplot(dm_inp$los, ylab = "Days", main = "Length of Stay (LOS)")

# Remove encounters with LOS < 1 day: 43673 encounters remaining
dm_inp <- dm_inp %>% filter(los >= 1)

# Filter to COVID-positive encounters ----

# 52.05% (22730) encounters have a COVID test; 4.24% (964) of those tests are positive
# Somehow, 0.81% (168) of encounters without a COVID test are positive
# Probably because the same COVID test has been assigned to all encounters of a patient
tbl_fn(dm_inp$covid_tested)
tbl_fn(dm_inp[which(dm_inp$covid_tested == "Yes"), ]$covid_pos)
tbl_fn(dm_inp[which(dm_inp$covid_tested == "No"), ]$covid_pos)

# Getting a list of unique patients with a positive COVID test and their test dates
covid_pos_dates <- dm_inp %>% filter(covid_pos == "Yes") %>%
  dplyr::select(deid_pat_id, covid_pos_time) %>% distinct()

# Check admit_time -- 0 rows with NA
sum(is.na(dm_inp$admit_time))

dm_inp <- dm_inp %>%
  dplyr::select(-covid_pos_time) %>%
  left_join(covid_pos_dates, by = "deid_pat_id", relationship = "many-to-many") %>%
  mutate(pos_adm = as.numeric(covid_pos_time - admit_time)) %>%
  drop_na(deid_enc_id)

# Look at the distribution of days from admission to positive COVID test
summary(dm_inp$pos_adm)
boxplot(dm_inp$pos_adm, ylab = "Days", main = "Days from admission to positive COVID test")

# Filter based on the days from admission to positive COVID test
# Keep encounters where the positive COVID test is within 14 days before admission to 7 days after admission
# Keep only one encounter per positive test a patient has --
#     If a patient has multiple encounters with the same positive test, keep the encounter closest to the positive test

dm_pos <- dm_inp %>% dplyr::filter(pos_adm >= -14 &
                                     pos_adm <= 2) %>%
  group_by(deid_enc_id) %>% arrange(abs(pos_adm)) %>% slice_head(n = 1) %>% ungroup() %>%
  group_by(deid_pat_id, covid_pos_time) %>% arrange(abs(pos_adm)) %>% slice_head(n = 1) %>% ungroup()

length(unique(dm_pos$deid_enc_id))
length(unique(dm_pos$deid_pat_id))

# Distribution of days from admission to positive COVID test
summary(dm_pos$pos_adm)
hist(dm_pos$pos_adm, ylab = "Days", main = "Days from admission to positive COVID test")

enc_ids <- dm_pos$deid_enc_id
saveRDS(enc_ids, here("data", "preprocess", "encounter_ids.rds"))

pat_ids <- unique(dm_pos$deid_pat_id)
saveRDS(pat_ids, here("data", "preprocess", "patient_ids.rds"))
