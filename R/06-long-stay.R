rm(list = ls())

library(here)
library(tidyverse)
library(msm)

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))

source(here("R", "00-qmat_setup.R"))
pt_stg <- readRDS(here("data", "temp", "pt_stg.rds"))

# Identify long-stay patients (>30 days)
long_stay_patients <- pt_stg %>%
  dplyr::filter(model == "base_model") %>%
  group_by(deid_enc_id) %>%
  summarise(max_los = max(DaysSinceEntry), .groups = "drop") %>%
  filter(max_los > 30) %>%
  pull(deid_enc_id)

# Excluded analysis (remove long-stay patients)
excluded_data <- pt_stg %>%
  filter(!deid_enc_id %in% long_stay_patients)

excluded_crude_rates <- calc_crude_init_rates(
  excluded_data, 
  qmat_list
)

excluded_models <- fit_msm_models(
  patient_data = excluded_data,
  crude_rates = excluded_crude_rates,
  covariates = NULL,
  mc.cores = n.cores
)

excluded_models_univar <- univariate_progressive_timeout(
  patient_data = excluded_data,
  crude_rates = list(base_model = excluded_crude_rates[["base_model"]]),
  covariates = key_covariates,
  n.cores = n.cores,
  timeout_vector = c(1, 2, 5, 10), 
  save_prefix = "excl_univar_progressive"
)

excluded_multivar <- multivariate_selection(
  patient_data = excluded_data,
  crude_rates = list(base_model = excluded_crude_rates[["base_model"]]),
  candidate_covariates = key_covariates,
  method = "forward",
  alpha_enter = 0.05,
  max_variables = length(key_covariates),
  n_cores = n.cores,
  save_intermediate = TRUE
)

# Truncated analysis (truncate at 30 days)
truncated_data <- pt_stg %>%
  filter(DaysSinceEntry <= 30)

truncated_crude_rates <- calc_crude_init_rates(
  truncated_data,
  qmat_list
)

truncated_models <- fit_msm_models(
  patient_data = truncated_data,
  crude_rates = list(base_model = truncated_crude_rates[["base_model"]]),
  covariates = NULL,
  mc.cores = n.cores
)

truncated_models_univar <- univariate_progressive_timeout(
  patient_data = truncated_data,
  crude_rates = list(base_model = truncated_crude_rates[["base_model"]]),
  covariates = key_covariates,
  n.cores = n.cores,
  timeout_vector = c(1, 2, 5, 10), 
  save_prefix = "trunc_univar_progressive"
)

truncated_multivar <- multivariate_selection(
  patient_data = truncated_data,
  crude_rates = list(base_model = truncated_crude_rates[["base_model"]]),
  candidate_covariates = key_covariates,
  method = "forward",
  alpha_enter = 0.05,
  max_variables = length(key_covariates),
  n_cores = n.cores,
  save_intermediate = TRUE
)

long_stay_sensitivity <- list(
  long_stay_patients = long_stay_patients,
  n_long_stay = length(long_stay_patients),
  excluded = list(
    data = excluded_data,
    crude_rates = excluded_crude_rates,
    models = excluded_models,
    univariate_models = excluded_models_univar,
    multivar_selection = excluded_multivar
  ),
  truncated = list(
    data = truncated_data,
    crude_rates = truncated_crude_rates,
    models = truncated_models,
    univariate_models = truncated_models_univar,
    multivar_selection = truncated_multivar
  )
)

saveRDS(long_stay_sensitivity, here("data", "temp", "long_stay_sensitivity.rds"))

long_stay_models <- c(
  excluded_models = excluded_models,
  excluded_univariate_models = excluded_models_univar,
  excluded_multivar = list(base_model = excluded_multivar),
  truncated_models = truncated_models,
  truncated_univariate_models = truncated_models_univar,
  truncated_multivar = list(base_model = truncated_multivar)
)

saveRDS(long_stay_models, here("data", "temp", "long_stay_models.rds"))

# long_stay_model_comp <- run_comprehensive_msm_analysis(
#   long_stay_models, pt_stg, crude_rates, mc.cores = n.cores, 
#   analysis_config = config_core_cov
# )
# 
# saveRDS(long_stay_model_comp, here("data", "temp", "long_stay_model_comp.rds"))
