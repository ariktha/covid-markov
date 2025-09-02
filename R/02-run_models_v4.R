rm(list = ls())

library(here)
library(tidyverse)
library(msm)

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))

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

pt_stg <- pt_stg %>% 
  cross_join(tibble(model = unique(model_specs$model))) %>%
  left_join(model_specs, by = c("stage", "hx_sev_bin", "model")) %>%
  add_calendar_time(date_column = "date") %>%
  left_join(dem_raw, by = "deid_enc_id") %>%
  arrange(deid_enc_id, date)

# Q-matrix specifications for different model structures
qmat_list <- list(
  base_model = matrix(c(1, 1, 1, 1,
                        1, 1, 1, 1,
                        0, 0, 0, 0,
                        0, 0, 0, 0), 
                      nrow = 4, ncol = 4, byrow = TRUE,
                      dimnames = list(c("M", "S", "D", "R"),
                                      c("M", "S", "D", "R"))),
  
  mod_2 = matrix(c(1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0), 
                 nrow = 5, ncol = 5, byrow = TRUE,
                 dimnames = list(c("M1", "M2", "S", "D", "R"),
                                 c("M1", "M2", "S", "D", "R"))),
  
  mod_3 = matrix(c(1, 1, 1, 1, 1, 1, 
                   1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 
                   1, 1, 1, 1, 1, 1, 
                   0, 0, 0, 0, 0, 0, 
                   0, 0, 0, 0, 0, 0), 
                 nrow = 6, ncol = 6, byrow = TRUE,
                 dimnames = list(c("M1", "M2", "M3", "S", "D", "R"),
                                 c("M1", "M2", "M3", "S", "D", "R"))),
  
  sev_2 = matrix(c(1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0), 
                 nrow = 5, ncol = 5, byrow = TRUE,
                 dimnames = list(c("M", "S1", "S2", "D", "R"),
                                 c("M", "S1", "S2", "D", "R"))),
  
  reduced_trans = matrix(c(1, 1, 0, 1,
                           1, 1, 1, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0), 
                         nrow = 4, ncol = 4, byrow = TRUE,
                         dimnames = list(c("M", "S", "D", "R"),
                                         c("M", "S", "D", "R"))),
  
  hx_sev = matrix(c(1, 0, 1, 1, 1,
                    0, 1, 1, 1, 1,
                    0, 1, 1, 1, 1,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1), 
                  nrow = 5, ncol = 5, byrow = TRUE,
                  dimnames = list(c("M", "MS", "S", "D", "R"),
                                  c("M", "MS", "S", "D", "R")))
)

crude_rates <- calc_crude_init_rates(pt_stg, qmat_list)

# Covariate specifications
time_invariant_covariates <- c("age", "sex", "race", "ethnicity", "language", 
                               "insurance_type", "smoking", "BMI", "bmi_cat", 
                               "COVID_vax", "cci_score", "chf", "cci_cat", 
                               "copd", "dnr_on_admit")

continuous_covariates <- c("age", "BMI", "cci_score")
key_covariates <- c("age", "sex", "race", "ethnicity", "language", "smoking", 
                    "BMI", "cci_cat", "COVID_vax", "dnr_on_admit", "chf", "copd")

# Save initial data and specifications
saveRDS(pt_stg, here("data", "temp", "pt_stg.rds"))
saveRDS(qmat_list, here("data", "temp", "qmat_list.rds"))
saveRDS(crude_rates, here("data", "temp", "crude_rates.rds"))

cat("Data preparation complete. N patients:", length(unique(pt_stg$deid_enc_id)), "\n")
cat("N observations:", nrow(pt_stg), "\n")

# Base models (no covariates) --------------------------------------------

# crude_rates_sub <- list(base_model = crude_rates[["base_model"]],
#                         hx_sev = crude_rates[["hx_sev"]])
# 
# test_mods <- fit_msm_models(
#   patient_data = pt_stg, 
#   crude_rates = crude_rates_sub,
#   covariates = list(
#     "no_covariates" = NULL,
#     "age_only" = "age",
#     "sex_only" = "sex"
#   ),
#   mc.cores = 6
# )

# Skip the computationally expensive components
config_core <- list(
  prevalence = list(skip = TRUE),
  cv = list(skip = TRUE), 
  residuals = list(skip = TRUE)
)

# test_comp_core <- run_comprehensive_msm_analysis(
#   test_mods, pt_stg, crude_rates_sub, mc.cores = 6, analysis_config = config_core)
# 
# test_time_mods_hosp <- fit_time_varying_models(
#   patient_data = pt_stg,
#   crude_rates = crude_rates_sub,
#   time_term = "days_since_entry",
#   time_covariates = c("linear", "piecewise")
# )
# 
# test_time_mods_cal <- fit_time_varying_models(
#   patient_data = pt_stg,
#   crude_rates = crude_rates_sub,
#   time_term = "calendar_time",
#   time_covariates = c("linear", "piecewise")
# )

cat("\n=== FITTING BASE MODELS ===\n")
fitted_base_models <- fit_msm_models(
  patient_data = pt_stg, 
  crude_rates = crude_rates,
  covariates = NULL,
  mc.cores = n.cores
)

# Tidy base model results
base_models_comp <- run_comprehensive_msm_analysis(
  fitted_base_models, pt_stg, crude_rates, mc.cores = n.cores, 
  analysis_config = config_core
)

saveRDS(fitted_base_models, here("data", "temp", "fitted_base_models.rds"))
saveRDS(base_models_comp, here("data", "temp", "base_models_comp.rds"))

gc()

# Markov assumption testing ----------------------------------------------

cat("\n=== TESTING MARKOV ASSUMPTIONS ===\n")

# 1. Time in severe as covariate (duration dependence)
markov_time_severe <- fit_msm_models(
  patient_data = pt_stg,
  crude_rates = crude_rates$base_model,
  covariates = list("time_in_severe" = "hx_sev_time"),
  mc.cores = n.cores
)

# 2. History of severe state model (already in base models as "hx_sev")
# Extract for comparison
hx_sev_model <- fitted_base_models[["hx_sev"]]

# Comparison: base vs time-in-severe covariate vs history state
# Create properly nested structure for Markov comparison
markov_comparison_models <- list(
  base_no_covariates = list(
    "~ 1" = fitted_base_models[["base_model"]][["~ 1"]]
  ),
  base_with_time_severe = list(
    "~ hx_sev_time" = markov_time_severe[["base_model"]][["~ hx_sev_time"]]
  ),
  history_severe_structure = list(
    "~ 1" = fitted_base_models[["hx_sev"]][["~ 1"]]
  )
)

markov_model_comp <- run_comprehensive_msm_analysis(
  markov_comparison_models, pt_stg, list(base_model = crude_rates[["base_model"]]), 
  mc.cores = n.cores, analysis_config = config_core
)

saveRDS(markov_time_severe, here("data", "temp", "markov_time_severe.rds"))
saveRDS(markov_comparison_models, here("data", "temp", "markov_comparison_models.rds"))
saveRDS(markov_model_comp, here("data", "temp", "markov_model_comp.rds"))

gc()

# Time-varying transition rates ------------------------------------------

cat("\n=== FITTING TIME-VARYING MODELS ===\n")

# Hospital time effects
hospital_time_models <- fit_time_varying_models(
  patient_data = pt_stg,
  crude_rates = crude_rates,
  time_term = "days_since_entry",
  time_covariates = c("linear", "piecewise"),
  spline_df = 3,
  spline_type = "ns"
)

# Calendar time effects  
calendar_time_models <- fit_time_varying_models(
  patient_data = pt_stg,
  crude_rates = crude_rates,
  time_term = "calendar_time",
  time_covariates = c("linear", "piecewise"),
  spline_df = 3,
  spline_type = "ns"
)

# # COVID wave periods
# covid_waves <- list(
#   "wave1" = c(0, 120),
#   "wave2" = c(200, 350), 
#   "wave3" = c(450, 600)
# )
# 
# calendar_wave_models <- fit_calendar_time_models(
#   patient_data = pt_stg,
#   crude_rates = crude_rates,
#   wave_periods = covid_waves
# )

# Tidy time-varying results
hospital_time_comp <- run_comprehensive_msm_analysis(
  hospital_time_models, pt_stg, crude_rates, mc.cores = n.cores, 
  analysis_config = config_core
)

calendar_time_comp <- run_comprehensive_msm_analysis(
  calendar_time_models, pt_stg, crude_rates, mc.cores = n.cores, 
  analysis_config = config_core
)

# calendar_wave_summary <- tidy_msm_models(calendar_wave_models)

# # Time model comparisons
# time_varying_comparison <- compare_time_varying_models(
#   c(hospital_time_models, calendar_time_models)
# )

saveRDS(hospital_time_models, here("data", "temp", "hospital_time_models.rds"))
saveRDS(calendar_time_models, here("data", "temp", "calendar_time_models.rds"))
# saveRDS(calendar_wave_models, here("data", "temp", "calendar_wave_models.rds"))
saveRDS(hospital_time_comp, here("data", "temp", "hospital_time_comp.rds"))
saveRDS(calendar_time_comp, here("data", "temp", "calendar_time_comp.rds"))

gc()

# Covariate effects ------------------------------------------------------

cat("\n=== FITTING COVARIATE MODELS ===\n")

# Univariate models
univariate_models <- fit_msm_models(
  patient_data = pt_stg,
  crude_rates = crude_rates,
  covariates = setNames(as.list(key_covariates), paste0("univar_", key_covariates)),
  mc.cores = n.cores
)

# # Spline models for continuous covariates
# spline_models <- fit_msm_models(
#   patient_data = pt_stg,
#   crude_rates = crude_rates_sub,
#   covariates = setNames(rep(list(NULL), length(continuous_covariates)), paste0("spline_", continuous_covariates)),
#   spline_vars = setNames(as.list(continuous_covariates), paste0("spline_", continuous_covariates)),
#   spline_df = list(default = 3),
#   spline_type = list(default = "ns"),
#   mc.cores = n.cores
# )
# 
# # Test nonlinear relationships
# nonlinear_tests <- map(continuous_covariates, function(cov) {
#   test_nonlinear_relationship(
#     patient_data = pt_stg,
#     crude_rates = crude_rates,
#     covariate = cov,
#     spline_df = 3,
#     spline_type = "ns"
#   )
# })
# names(nonlinear_tests) <- continuous_covariates

# Transition-specific effects
transition_specific_models <- fit_transition_specific_models(
  patient_data = pt_stg,
  crude_rates = crude_rates,
  covariates = list("age" = c("age")),
  constraint_configs = list(
    "constant" = "default",
    "by_transition" = NULL
  )
)

# # Interaction models
# interaction_models <- fit_interaction_models(
#   patient_data = pt_stg,
#   crude_rates = crude_rates,
#   base_covariates = key_covariates,
#   interaction_specs = "all_pairwise",
#   include_main_effects = TRUE,
#   mc.cores = n.cores
# )

# Multivariate selection (base model only for efficiency)
base_model_data <- pt_stg %>% filter(model == "base_model")
multivariate_models <- multivariate_selection(
  patient_data = base_model_data,
  crude_rates = list(base_model = crude_rates[["base_model"]]),
  candidate_covariates = key_covariates,
  method = "forward",
  alpha_enter = 0.05
)

# Tidy covariate results
univariate_comp <- run_comprehensive_msm_analysis(
  univariate_models, pt_stg, crude_rates, mc.cores = n.cores, 
  analysis_config = config_core
)

# spline_summary <- tidy_msm_models(spline_models)
transition_specific_comp <- run_comprehensive_msm_analysis(
  transition_specific_models, pt_stg, crude_rates, mc.cores = n.cores, 
  analysis_config = config_core
)

# interaction_summary <- tidy_msm_models(interaction_models)

covariate_comparison <- compare_covariate_effects(
  c(univariate_models),
  base_formula = "~ 1"
)

# interaction_comparison <- compare_interaction_models(interaction_models)


saveRDS(univariate_models, here("data", "temp", "univariate_models.rds"))
saveRDS(univariate_comp, here("data", "temp", "univariate_comp.rds"))
# saveRDS(spline_models, here("data", "temp", "spline_models.rds"))
saveRDS(transition_specific_models, here("data", "temp", "transition_specific_models.rds"))
saveRDS(transition_specific_comp, here("data", "temp", "transition_specific_comp.rds"))
# saveRDS(interaction_models, here("data", "temp", "interaction_models.rds"))
saveRDS(multivariate_models, here("data", "temp", "multivariate_models.rds"))
# saveRDS(nonlinear_tests, here("data", "temp", "nonlinear_tests.rds"))
# saveRDS(covariate_comparison, here("data", "temp", "covariate_comparison.rds"))
# saveRDS(interaction_comparison, here("data", "temp", "interaction_comparison.rds"))

gc()

# # Predictive performance evaluation --------------------------------------
# 
# cat("\n=== EVALUATING PREDICTIVE PERFORMANCE ===\n")
# 
# # Cross-validation for all models
# cv_all_models <- c(
#   fitted_base_models,
#   markov_time_severe,
#   hospital_time_models,
#   calendar_time_models,
#   univariate_models,
#   transition_specific_models
#   # spline_models,
#   # interaction_models
# )
# 
# key_models_for_cv <- list(
#   base_model = fitted_base_models[["base_model"]],
#   best_univariate = univariate_models[["base_model"]]
#   # best_interaction = interaction_models[["main_effects"]]
# )
# 
# cv_results <- calculate_predictive_performance(
#   patient_data = base_model_data,
#   fitted_models = key_models_for_cv,
#   crude_rates = list(base_model = crude_rates[["base_model"]]),
#   k_folds = 5,
#   prediction_times = c(14, 30),
#   parallel = TRUE,
#   n_cores = n.cores
# )
# 
# saveRDS(cv_results, here("data", "temp", "cv_results.rds"))
# cat("Cross-validation completed for", length(key_models_for_cv), "models\n")
# 
# # State prevalence evaluation --------------------------------------------
# 
# cat("\n=== CALCULATING STATE PREVALENCES ===\n")
# 
# # Calculate prevalences for key models
# prevalence_time_points <- seq(1, 30, by = 1)
# 
# base_prevalences <- tidy_msm_prevalences(
#   fitted_base_models,
#   time_points = prevalence_time_points,
#   mc.cores = n.cores,
#   ci = TRUE
# )
# 
# covariate_prevalences <- tidy_msm_prevalences(
#   list(base_univariate = univariate_models[["base_model"]]),
#   time_points = prevalence_time_points,
#   mc.cores = n.cores,
#   ci = TRUE
# )
# 
# saveRDS(base_prevalences, here("data", "temp", "base_prevalences.rds"))
# saveRDS(covariate_prevalences, here("data", "temp", "covariate_prevalences.rds"))
# 
# cat("Prevalences calculated for", length(prevalence_time_points), "time points\n")
# 
# # Residual diagnostics ---------------------------------------------------
# 
# cat("\n=== CALCULATING RESIDUAL DIAGNOSTICS ===\n")
# 
# # Transition residuals for key models
# transition_residuals <- calculate_transition_residuals(
#   fitted_msm_models = key_models_for_cv,
#   patient_data = base_model_data,
#   residual_type = "deviance",
#   debug = FALSE
# )
# 
# saveRDS(transition_residuals, here("data", "temp", "transition_residuals.rds"))
# cat("Transition residuals calculated for", nrow(transition_residuals), "transitions\n")
# 
# Long-stay patient sensitivity analysis --------------------------------

cat("\n=== LONG-STAY PATIENT SENSITIVITY ANALYSIS ===\n")

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

excluded_models_univar <- fit_msm_models(
  patient_data = excluded_data,
  crude_rates = excluded_crude_rates,
  covariates = setNames(as.list(key_covariates), paste0("univar_", key_covariates)),
  mc.cores = n.cores
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
  crude_rates = truncated_crude_rates,
  covariates = NULL,
  mc.cores = n.cores
)

truncated_models_univar <- fit_msm_models(
  patient_data = truncated_data,
  crude_rates = truncated_crude_rates,
  covariates = setNames(as.list(key_covariates), paste0("univar_", key_covariates)),
  mc.cores = n.cores
)

long_stay_sensitivity <- list(
  long_stay_patients = long_stay_patients,
  n_long_stay = length(long_stay_patients),
  excluded = list(
    data = excluded_data,
    crude_rates = excluded_crude_rates,
    models = excluded_models,
    univariate_models = excluded_models_univar
  ),
  truncated = list(
    data = truncated_data,
    crude_rates = truncated_crude_rates,
    models = truncated_models,
    univariate_models = truncated_models_univar
  )
)

saveRDS(long_stay_sensitivity, here("data", "temp", "long_stay_sensitivity.rds"))

long_stay_models <- c(
  base_models = fitted_base_models,
  univariate_models = univariate_models,
  excluded_models = excluded_models,
  excluded_univariate_models = excluded_models_univar,
  truncated_models = truncated_models,
  truncated_univariate_models = truncated_models_univar
)

long_stay_model_comp <- run_comprehensive_msm_analysis(
  long_stay_models, pt_stg, crude_rates, mc.cores = n.cores, 
  analysis_config = config_core
)

saveRDS(long_stay_models, here("data", "temp", "long_stay_models.rds"))
saveRDS(long_stay_model_comp, here("data", "temp", "long_stay_model_comp.rds"))