rm(list = ls())

library(here)
library(tidyverse)
library(msm)

source(here("R", "00-config.R"))
source(here("R", "00-functions_redo.R"))

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

# Save initial data and specifications
saveRDS(pt_stg, here("data", "temp", "pt_stg.rds"), compress = FALSE)
saveRDS(qmat_list, here("data", "temp", "qmat_list.rds"))
saveRDS(crude_rates, here("data", "temp", "crude_rates.rds"))

cat("Data preparation complete. N patients:", length(unique(pt_stg$deid_enc_id)), "\n")
cat("N observations:", nrow(pt_stg), "\n")

# Base models (no covariates) --------------------------------------------

if(fit_no_cov_models){
  cat("\n=== Fitting base models (No covariates) ===\n")
  start_time <- Sys.time()
  fitted_base_models <- fit_msm_models(
    patient_data = pt_stg, 
    crude_rates = crude_rates,
    covariates = NULL,
    mc.cores = n.cores
  )
  end_time <- Sys.time()
  cat("Base models fitting time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
  saveRDS(fitted_base_models, here("data", "temp", "fitted_base_models.rds"))
  
} else {
  cat("\n=== Loading pre-fitted base models (No covariates) ===\n")
  fitted_base_models <- readRDS(here("data", "temp", "fitted_base_models.rds"))
  
}

if(comp_no_cov_models){
  cat("\n=== Compiling base model results ===\n")
  start_time <- Sys.time()
  base_models_comp <- run_comprehensive_msm_analysis(
    fitted_base_models, pt_stg, crude_rates, mc.cores = n.cores, 
    analysis_config = config_core
  )
  end_time <- Sys.time()
  cat("Base models results compilation time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
  saveRDS(base_models_comp, here("data", "temp", "base_models_comp.rds"))
  
} else {
  cat("\n=== Loading pre-compiled base model results ===\n")
  base_models_comp <- readRDS(here("data", "temp", "base_models_comp.rds"))
  
}

gc()

pred <- calculate_predictive_performance(
  patient_data = pt_stg,
  fitted_models = fitted_base_models,
  crude_rates = crude_rates,
  k_folds = 5,
  prediction_times = seq(1, 15, by = 2),
  parallel = TRUE,
  n_cores = n.cores
)

saveRDS(pred, here("data", "temp", "base_model_predictive_performance.rds"))

# prev <- tidy_msm_prevalences(
#   fitted_base_models,
#   time_points = c(1, 10),
#   mc.cores = n.cores,
#   use_approximation = TRUE
# )

# Markov assumption testing ----------------------------------------------

cat("\n=== TESTING MARKOV ASSUMPTIONS ===\n")

# 1a. Time in severe as covariate: linear term
markov_time_severe <- fit_msm_models(
  patient_data = pt_stg,
  crude_rates =  list(base_model = crude_rates[["base_model"]]),
  covariates = list("time_in_severe" = "hx_sev_time"),
  mc.cores = n.cores
)

# 1b. Time in severe as covariate: spline term
markov_time_severe_spline <- fit_spline_msm_models(
  patient_data = pt_stg,
  crude_rates =  list(base_model = crude_rates[["base_model"]]),
  covariates = list("time_in_severe_spline" = NULL),
  spline_vars = list("time_in_severe_spline" = "hx_sev_time"),
  spline_df = list("time_in_severe_spline" = 3),
  mc.cores = n.cores
)

# 1c. Time in severe as covariate: transition-specific covariate
markov_time_severe_ts <- fit_transition_specific_models(
  patient_data = pt_stg,
  crude_rates = list(base_model = crude_rates[["base_model"]]),
  covariates = list("time_in_severe_ts" = "hx_sev_time")
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
  base_with_time_severe_spline = list(
    "~ hx_sev_time_ns1 + hx_sev_time_ns2 + hx_sev_time_ns3" = markov_time_severe_spline[["base_model"]][["~ hx_sev_time_ns1 + hx_sev_time_ns2 + hx_sev_time_ns3"]]
  ),
  base_with_time_severe_ts = list(
    "~ hx_sev_time (transition-specific)" = markov_time_severe_ts[["base_model"]][["~ hx_sev_time"]]
  ),
  history_severe_structure = list(
    "~ 1" = fitted_base_models[["hx_sev"]][["~ 1"]]
  )
)

markov_model_comp <- run_comprehensive_msm_analysis(
  markov_comparison_models, 
  pt_stg, 
  list(base_model = crude_rates[["base_model"]]), 
  mc.cores = n.cores, 
  analysis_config = config_core
)

# markov_prev <- tidy_msm_prevalences(
#   markov_comparison_models,
#   time_points = seq(1, 30, by = 5),
#   mc.cores = n.cores,
#   use_approximation = TRUE
# )

# markov_pred <- calculate_predictive_performance(
#   patient_data = pt_stg,
#   fitted_models = markov_comparison_models,
#   crude_rates = list(base_model = crude_rates[["base_model"]],
#                      hx_sev = crude_rates[["hx_sev"]]),
#   k_folds = 5,
#   prediction_times = seq(1, 15, by = 2),
#   parallel = TRUE,
#   n_cores = n.cores
# )

# saveRDS(markov_time_severe, here("data", "temp", "markov_time_severe.rds"))
saveRDS(markov_comparison_models, here("data", "temp", "markov_comparison_models.rds"))
saveRDS(markov_model_comp, here("data", "temp", "markov_model_comp.rds"))
# saveRDS(markov_prev, here("data", "temp", "markov_prev.rds"))
# saveRDS(markov_pred, here("data", "temp", "markov_pred.rds"))

gc()

# Covariate effects ------------------------------------------------------

cat("\n=== FITTING COVARIATE MODELS ===\n")

# Covariate specifications

continuous_covariates <- c("age", "cci_score", "BMI")

key_covariates <- c("age", "sex", "race", "ethnicity", "language", 
                    "insurance_type", "smoking", "BMI", "bmi_cat", 
                    "COVID_vax", "cci_score", "chf", "cci_cat", 
                    "copd", "dnr_on_admit")

## Linear univariate models ---------------------------------------

if(fit_cov_models){
  
  univar_results <- univariate_progressive_timeout(
    patient_data = pt_stg,
    crude_rates = crude_rates,
    covariates = key_covariates,
    n.cores = n.cores,
    timeout_vector = c(1, 2, 5, 10), 
    save_prefix = "univar_progressive"
  )
  
} else {
  univar_results <- readRDS(here("data", "temp", "univar_progressive_full_results.rds"))
}

univariate_models <- univar_results$models
failed_summary <- univar_results$failed_models

# Tidy covariate results
univariate_comp <- run_comprehensive_msm_analysis(
  univariate_models, 
  pt_stg, 
  crude_rates, 
  mc.cores = n.cores, 
  analysis_config = config_core_cov
)

## Spline univariate models ---------------------------------------

if(fit_spline_models){
  
  cat("\n=== FITTING SPLINE MODELS ===\n")
  
  individual_spline_models <- list()
  for (covar in continuous_covariates) {
    cat("Fitting spline model for:", covar, "\n")
    
    cov_name <- paste0("spline_", covar)
    model_result <- tryCatch({
      withTimeout({
        fit_spline_msm_models(
          patient_data = pt_stg,
          crude_rates = list(base_model = crude_rates[["base_model"]]),
          covariates = setNames(list(NULL), cov_name),           # spline_age = NULL
          spline_vars = setNames(list(covar), cov_name),         # spline_age = "age"
          spline_df = list(default = 3),
          mc.cores = n.cores
        )
      }, timeout = 20 * 60)  # 20 minutes in seconds
    }, TimeoutException = function(e) {
      cat("WARNING: Model for", covar, "timed out after 20 minutes - skipping\n")
      return(NULL)
    }, error = function(e) {
      cat("ERROR: Model for", covar, "failed with error:", e$message, "- skipping\n")
      return(NULL)
    })
    
    # Only save if model fitting succeeded
    if (!is.null(model_result)) {
      # Save the individual model
      individual_spline_models[[covar]] <- model_result
      
      # Optional: Save to disk immediately
      saveRDS(model_result, file = here("data", "temp", paste0("spline_model_", covar, ".rds")))
      
      cat("Completed model for:", covar, "\n")
    } else {
      cat("Skipped model for:", covar, "\n")
    }
  }
  
  spline_models <- do.call(c, individual_spline_models)
  saveRDS(spline_models, here("data", "temp", "spline_models.rds"))
  
} else {
  cat("\n=== LOADING PRE-FITTED SPLINE MODELS ===\n")
  
  spline_models <- readRDS(here("data", "temp", "spline_models.rds"))
}

if(comp_spline_models){
  cat("\n=== COMPILING SPLINE MODEL RESULTS ===\n")
  
  spline_comp <- run_comprehensive_msm_analysis(
    spline_models, 
    pt_stg, 
    crude_rates = list(base_model = crude_rates[["base_model"]]), 
    mc.cores = n.cores, 
    analysis_config = config_core_cov
  )
  saveRDS(spline_comp, here("data", "temp", "spline_comp.rds"))
  
} else {
  cat("\n=== LOADING PRE-COMPILED SPLINE MODEL RESULTS ===\n")
  
  spline_comp <- readRDS(here("data", "temp", "spline_comp.rds"))
}

if(hr_spline_models){

} else {


}


## Transition-specific effects ---------------------------------------

transition_specific_models <- fit_transition_specific_models(
  patient_data = pt_stg,
  crude_rates = list(base_model = crude_rates[["base_model"]]),
  covariates = setNames(as.list(key_covariates), paste0("univar_", key_covariates))
)

transition_specific_comp <- run_comprehensive_msm_analysis(
  transition_specific_models, 
  pt_stg, 
  crude_rates = list(base_model = crude_rates[["base_model"]]), 
  mc.cores = n.cores, 
  analysis_config = config_core_cov
)

## Interaction models ---------------------------------------

# interaction_models <- fit_interaction_models(
#   patient_data = pt_stg,
#   crude_rates = crude_rates,
#   base_covariates = key_covariates,
#   interaction_specs = "all_pairwise",
#   include_main_effects = TRUE,
#   mc.cores = n.cores
# )

# interaction_comparison <- compare_interaction_models(interaction_models)

## Multivariable models ---------------------------------------

# Multivariate selection (base model only for efficiency)
base_model_data <- pt_stg %>% filter(model == "base_model")
multivariate_models <- multivariate_selection(
  patient_data = base_model_data,
  crude_rates = list(base_model = crude_rates[["base_model"]]),
  candidate_covariates = key_covariates,
  required_covariates = c("age", "BMI", "cci_cat"),
  method = "forward",
  alpha_enter = 0.05,
  max_variables = 4,
  # max_variables = length(key_covariates),
  n_cores = n.cores,
  save_intermediate = TRUE
)

temp_multivariable <- fit_msm_models(
  patient_data = base_model_data,
  crude_rates = list(base_model = crude_rates[["base_model"]]),
  covariates = list("multivariable" = c("age", "BMI", "cci_cat")),
  mc.cores = n.cores
)

temp_multivariable_comp <- run_comprehensive_msm_analysis(
  temp_multivariable, 
  base_model_data, 
  crude_rates = list(base_model = crude_rates[["base_model"]]), 
  mc.cores = n.cores, 
  analysis_config = config_core_cov
)

## Save all covariate results ---------------------------------------

saveRDS(temp_multivariable, here("data", "temp", "temp_multivariable.rds"))
saveRDS(temp_multivariable_comp, here("data", "temp", "temp_multivariable_comp.rds"))

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

# Time-varying transition rates ------------------------------------------

cat("\n=== FITTING TIME-VARYING MODELS ===\n")

# Hospital time effects
hospital_time_models <- fit_time_varying_models(
  patient_data = pt_stg,
  crude_rates = list(base_model = crude_rates[["base_model"]]),
  time_term = "days_since_entry",
  time_covariates = c("linear", "piecewise"),
  spline_df = 3,
  spline_type = "ns"
)

# Calendar time effects  
calendar_time_models <- fit_time_varying_models(
  patient_data = pt_stg,
  crude_rates = list(base_model = crude_rates[["base_model"]]),
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
  hospital_time_models, 
  pt_stg, 
  crude_rates = list(base_model = crude_rates[["base_model"]]), 
  mc.cores = n.cores, 
  analysis_config = config_core
)

calendar_time_comp <- run_comprehensive_msm_analysis(
  calendar_time_models, 
  pt_stg, 
  crude_rates = list(base_model = crude_rates[["base_model"]]), 
  mc.cores = n.cores, 
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


# Proof of concept --------------------------------------------------------

poc_rates <- list(
  base_model = crude_rates[["base_model"]],
  hx_sev = crude_rates[["hx_sev"]]
)

poc_models <- list(
  base_model = fitted_base_models[["base_model"]],
  hx_sev = fitted_base_models[["hx_sev"]]
)

poc_data <- pt_stg %>% filter(model %in% c("base_model", "hx_sev"))

poc_prev <- tidy_msm_prevalences(
  poc_models,
  time_points = seq(1, 12, by = 3),
  mc.cores = 1,
  ci = TRUE
)

poc_pred <- calculate_predictive_performance(
  patient_data = poc_data,
  fitted_models = poc_models,
  crude_rates = poc_rates,
  k_folds = 3,
  prediction_times = c(7, 14),
  parallel = TRUE,
  n_cores = 6
)

poc_residuals <- calculate_transition_residuals(
  fitted_msm_models = poc_models,
  patient_data = poc_data,
  residual_type = "deviance",
  debug = FALSE
)

saveRDS(poc_prev, here("data", "temp", "poc_prevalences.rds"))
saveRDS(poc_pred, here("data", "temp", "poc_predictive_performance.rds"))
saveRDS(poc_residuals, here("data", "temp", "poc_transition_residuals.rds"))

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
  crude_rates = list(base_model = excluded_crude_rates[["base_model"]]),
  covariates = NULL,
  mc.cores = n.cores
)

excluded_models_univar <- fit_msm_models(
  patient_data = excluded_data,
  crude_rates = list(base_model = excluded_crude_rates[["base_model"]]),
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
  crude_rates = list(base_model = truncated_crude_rates[["base_model"]]),
  covariates = NULL,
  mc.cores = n.cores
)

truncated_models_univar <- fit_msm_models(
  patient_data = truncated_data,
  crude_rates = list(base_model = truncated_crude_rates[["base_model"]]),
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
  analysis_config = config_core_cov
)

saveRDS(long_stay_models, here("data", "temp", "long_stay_models.rds"))
saveRDS(long_stay_model_comp, here("data", "temp", "long_stay_model_comp.rds"))
gc()
