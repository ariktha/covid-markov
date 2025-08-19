rm(list = ls())

library(here)
library(tidyverse)
library(msm)

source(here("scripts", "00-config.R"))
source(here("scripts", "00-functions.R"))

# Load data
dem_raw <- readRDS(here("data", "pt_demographics.rds"))
stg_raw <- readRDS(here("data", "pt_enc_staging.rds"))
model_specs_raw <- read_csv(here("data", "model-specs.csv"))

model_specs <- model_specs_raw %>%
  mutate(stage = strsplit(as.character(stage), ", ")) %>% 
  unnest(stage) %>%
  mutate(hx_sev_bin = strsplit(as.character(hx_sev_bin), ", ")) %>% 
  unnest(hx_sev_bin) %>%
  mutate(stage = as.integer(stage), hx_sev_bin = as.integer(hx_sev_bin))

pt_stg <- stg_raw %>% group_by(deid_enc_id) %>%
  mutate(severe_bin = ifelse(stage %in% 7:9, 1, 0)) %>%
  mutate(hx_sev_time = slider::slide_index_sum(
    x = severe_bin,
    i = date,
    before = lubridate::days(365)
  )) %>%
  mutate(hx_sev_bin = ifelse(hx_sev_time > 0, 1, 0)) %>%
  dplyr::select(-c(severe_bin)) %>% ungroup()

pt_stg <- pt_stg %>% 
  cross_join(tibble(model = unique(model_specs$model))) %>%
  left_join(model_specs, by = c("stage", "hx_sev_bin", "model"))

saveRDS(pt_stg, here("data", "temp", "pt_stg.rds"))

qmat_list <- list(
  base_model = matrix(c(1, 1, 1, 1,
                        1, 1, 1, 1,
                        0, 0, 0, 0,
                        0, 0, 0, 0), 
                      nrow = 4, ncol = 4,
                      byrow = TRUE,
                      dimnames = list(c("M", "S", "D", "R"),
                                      c("M", "S", "D", "R"))),
  
  mod_2 = matrix(c(1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0), 
                 nrow = 5, ncol = 5,
                 byrow = TRUE,
                 dimnames = list(c("M1", "M2", "S", "D", "R"),
                                 c("M1", "M2", "S", "D", "R"))),
  
  mod_3 = matrix(c(1, 1, 1, 1, 1, 1, 
                   1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 
                   1, 1, 1, 1, 1, 1, 
                   0, 0, 0, 0, 0, 0, 
                   0, 0, 0, 0, 0, 0), 
                 nrow = 6, ncol = 6,
                 byrow = TRUE,
                 dimnames = list(c("M1", "M2", "M3", "S", "D", "R"),
                                 c("M1", "M2", "M3", "S", "D", "R"))),
  sev_2 = matrix(c(1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0), 
                 nrow = 5, ncol = 5,
                 byrow = TRUE,
                 dimnames = list(c("M", "S1", "S2", "D", "R"),
                                 c("M", "S1", "S2", "D", "R"))),
  reduced_trans = matrix(c(1, 1, 0, 1,
                           1, 1, 1, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0), 
                         nrow = 4, ncol = 4,
                         byrow = TRUE,
                         dimnames = list(c("M", "S", "D", "R"),
                                         c("M", "S", "D", "R"))),
  hx_sev = matrix(c(1, 0, 1, 1, 1,
                    0, 1, 1, 1, 1,
                    0, 1, 1, 1, 1,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1), 
                  nrow = 5, ncol = 5,
                  byrow = TRUE,
                  dimnames = list(c("M", "MS", "S", "D", "R"),
                                  c("M", "MS", "S", "D", "R")))
)

crude_rate <- calc_crude_init_rates(pt_stg, qmat_list)

saveRDS(qmat_list, here("data", "temp", "mod_qmat_list.rds"))
saveRDS(crude_rate, here("data", "temp", "mod_crude_rates.rds"))

# Fit the models

fitted_base_models <- fit_msm_models(pt_stg, crude_rate)
tidied_base_results <- tidy_msm_models(fitted_base_models)
tidied_base_pmats <- tidy_msm_pmats(fitted_base_models)

saveRDS(fitted_base_models, here("data", "temp", "models.rds"))
saveRDS(tidied_base_results, here("data", "temp", "mod_rates_tidy.rds"))
saveRDS(tidied_base_pmats, here("data", "temp", "mod_pmats_tidy.rds"))

