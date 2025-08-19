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


# Data prep ---------------------------------------------------------------

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

# Fit base models --------------------------------------------------------

fitted_base_models <- fit_msm_models(pt_stg, crude_rate)
tidied_base_results <- tidy_msm_models(fitted_base_models)
tidied_base_pmats <- tidy_msm_pmats(fitted_base_models)

saveRDS(fitted_base_models, here("data", "temp", "models.rds"))
saveRDS(tidied_base_results, here("data", "temp", "mod_rates_tidy.rds"))
saveRDS(tidied_base_pmats, here("data", "temp", "mod_pmats_tidy.rds"))


# Covariate models --------------------------------------------------------

time_invariant_covariates <- c("age", "sex", "race", "ethnicity", "language", 
                               "insurance_type", "smoking", "BMI", "bmi_cat", 
                               "COVID_vax", "cci_score", "chf", "cci_cat", 
                               "copd", "dnr_on_admit")

continuous_covariates <- c("age", "BMI", "cci_score")

# Add transition trend categories
trend_types <- add_transition_trends(pt_stg_raw)

pt_stg_enhanced <- pt_stg %>%
  left_join(dem_raw, by = "deid_enc_id") %>%
  arrange(deid_enc_id, date) %>%
  group_by(deid_enc_id) %>%
  mutate(
    from_state = state,
    to_state = lead(state)
  ) %>%
  ungroup() %>%
  left_join(trend_types, by = c("from_state" = "from", "to_state" = "to"))


## Univariate models ------------------------------------------------------

# Fit univariate models with spline testing for continuous variables
covariate_models <- fit_covariate_models(
  patient_data = pt_stg_enhanced,
  crude_rates = crude_rates,
  covariates = time_invariant_covariates,
  continuous_vars = continuous_covariates,
  spline_df = 3
)

# Extract model statistics
covariate_stats <- extract_covariate_stats(covariate_models, fitted_base_models)

# Save results
saveRDS(covariate_models, here("data", "temp", "covariate_models.rds"))
saveRDS(covariate_stats, here("data", "temp", "covariate_stats.rds"))

## Transition-specific covariate effects -----------------------------------

# Select key covariates for transition-specific analysis
key_covariates <- c("age", "cci_cat")

transition_specific_models <- fit_transition_specific_models(
  patient_data = pt_stg_enhanced,
  crude_rates = crude_rates,
  covariates = key_covariates
)

# Extract statistics for transition-specific models
transition_stats <- extract_covariate_stats(transition_specific_models, fitted_base_models)

# Save results
saveRDS(transition_specific_models, here("data", "temp", "transition_specific_models.rds"))
saveRDS(transition_stats, here("data", "temp", "transition_stats.rds"))

## Multivariable model with forward selection -----------------------------

model_data <- pt_stg %>% filter(model == "base_model")

multivariate_models <- multivariate_selection(
  patient_data = model_data,
  crude_rates = setNames(list(crude_rates[["base_model"]]), "base_model"),
  covariates = time_invariant_covariates,
  method = "forward",
  alpha_enter = 0.05
)

if (!is.null(multivariate_models[[best_base_model]])) {
  selected_vars <- multivariate_models[[best_base_model]]$selected_covariates
  cat("Selected variables:", paste(selected_vars, collapse = ", "), "\n")
  cat("Final AIC:", round(multivariate_models[[best_base_model]]$final_aic, 1), "\n\n")
}

# Save multivariate model results

saveRDS(multivariate_models, here("data", "temp", "multivariate_models.rds"))

# Time-varying models -----------------------------------------------------

# Fit models with time-dependent transition rates
time_varying_models <- fit_time_varying_models(
  patient_data = pt_stg_enhanced,
  crude_rates = crude_rates,
  time_covariates = c("linear", "spline", "piecewise")
)

# Compare time-varying models to base models
time_varying_stats <- extract_covariate_stats(time_varying_models, fitted_base_models)

# Save results
saveRDS(time_varying_models, here("data", "temp", "time_varying_models.rds"))
saveRDS(time_varying_stats, here("data", "temp", "time_varying_stats.rds"))

# Sensitivity analysis for long-stay patients ----------------------------

# Perform sensitivity analysis for patients with stays > 30 days
sensitivity_results <- sensitivity_analysis_outliers(
  patient_data = pt_stg_enhanced,
  crude_rates = crude_rates,
  cutoff_day = 30
)

# Calculate model comparison metrics for sensitivity analysis
sensitivity_comparison <- list(
  excluded_vs_original = compare_same_structure_models(
    fitted_models = c(fitted_base_models, sensitivity_results$excluded$models),
    model_names = names(fitted_base_models)
  ),
  truncated_vs_original = compare_same_structure_models(
    fitted_models = c(fitted_base_models, sensitivity_results$truncated$models),
    model_names = names(fitted_base_models)
  )
)

# Save results
saveRDS(sensitivity_results, here("data", "temp", "sensitivity_results.rds"))
saveRDS(sensitivity_comparison, here("data", "temp", "sensitivity_comparison.rds"))


# Save all ----------------------------------------------------------------

# Create a comprehensive summary
all_models <- list(
  base_models = fitted_base_models,
  covariate_models = covariate_models,
  time_varying_models = time_varying_models,
  transition_specific_models = transition_specific_models,
  sensitivity_results = sensitivity_results,
  model_statistics = list(
    covariate_stats = covariate_stats,
    time_varying_stats = time_varying_stats,
    transition_stats = transition_stats
  ),
  data_info = list(
    n_patients = length(unique(pt_stg_enhanced$deid_enc_id)),
    n_observations = nrow(pt_stg_enhanced),
    covariates_tested = time_invariant_covariates,
    continuous_covariates = continuous_covariates,
    long_stay_patients = length(sensitivity_results$long_stay_patients)
  )
)

# Save comprehensive summary
saveRDS(all_models, here("data", "temp", "all_models.rds"))

