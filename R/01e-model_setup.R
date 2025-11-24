library(here)
library(tidyverse)
library(msm)

rm(list = ls())

source(here("R", "00a-fns-data_preprocess.R"))

# Read data ---------------------------------------------------------------

dem_raw <- readRDS(here("data", "pt_demographics.rds"))
stg_raw <- readRDS(here("data", "pt_stage_data.rds"))
model_specs_raw <- read_csv(here("data", "model-specs.csv"))

# Process model specifications --------------------------------------------
model_specs <- model_specs_raw %>%
  mutate(stage = strsplit(as.character(stage), ", ")) %>% 
  unnest(stage) %>%
  mutate(hx_sev_bin = strsplit(as.character(hx_sev_bin), ", ")) %>% 
  unnest(hx_sev_bin) %>%
  mutate(stage = as.integer(stage), hx_sev_bin = as.integer(hx_sev_bin))

# Create patient staging data with history of severity ------
pt_stg <- stg_raw %>% 
  group_by(deid_enc_id) %>%
  mutate(severe_bin = ifelse(stage %in% 7:9, 1, 0)) %>%
  mutate(hx_sev_time = slider::slide_index_sum(
    x = severe_bin,
    i = date,
    before = lubridate::days(365)
  )) %>%
  mutate(hx_sev_bin = ifelse(hx_sev_time > 0, 1, 0)) %>%
  select(-severe_bin) %>% 
  ungroup()

# Join with model specifications and demographics ------
pt_stg <- pt_stg %>% 
  cross_join(tibble(model = unique(model_specs$model))) %>%
  left_join(model_specs, by = c("stage", "hx_sev_bin", "model")) %>%
  add_calendar_time(date_column = "date") %>%
  left_join(dem_raw, by = "deid_enc_id") %>%
  arrange(deid_enc_id, date)

# Create time in current stage variable ------

pt_stg <- pt_stg %>%
  arrange(deid_enc_id, model, DaysSinceEntry) %>%
  group_by(deid_enc_id, model) %>%
  mutate(
    # Detect state changes
    state_changed = state_num != lag(state_num, default = first(state_num)),
    # Create state episode ID
    state_episode = cumsum(state_changed),
    # Calculate days in current state
    time_in_current_state = DaysSinceEntry - first(DaysSinceEntry[state_episode == state_episode])
  ) %>%
  ungroup()

# Calculate crude rates ------
source(here("R", "00b-qmat_setup.R"))
crude_rates <- calc_crude_init_rates(pt_stg, qmat_list)

# Save processed data ------
saveRDS(pt_stg, here("data", "pt_stg.rds"), compress = FALSE)
saveRDS(crude_rates, here("data", "crude_rates.rds"))
