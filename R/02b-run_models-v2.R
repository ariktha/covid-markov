# Enhanced model running script
# Save as R/02b-run_enhanced_models.R

rm(list = ls())

library(here)
library(tidyverse)
library(msm)
library(splines)
library(furrr)
library(parallel)

source(here("scripts", "00-config.R"))
source(here("scripts", "00-functions.R"))  # Now includes enhanced functions

# ============================================================================
# LOAD DATA AND BASE MODELS
# ============================================================================

# Load previously fitted base models
fitted_base_models <- readRDS(here("data", "temp", "models.rds"))
pt_stg <- readRDS(here("data", "temp", "pt_stg.rds"))
qmat_list <- readRDS(here("data", "temp", "mod_qmat_list.rds"))
crude_rates <- readRDS(here("data", "temp", "mod_crude_rates.rds"))

# Define covariate lists
time_invariant_covariates <- c("age", "sex", "race", "ethnicity", "language", 
                               "insurance_type", "smoking", "BMI", "bmi_cat", 
                               "COVID_vax_doses", "COVID_vax", "cci_score", 
                               "cci_fct", "chf", "cci_cat", "copd", "DNR", 
                               "dnr_on_admit")

continuous_covariates <- c("age", "BMI", "COVID_vax_doses", "cci_score")

# Add transition trend categories
trend_types <- add_transition_trends(pt_stg)
pt_stg_enhanced <- pt_stg %>%
  arrange(deid_enc_id, date) %>%
  group_by(deid_enc_id) %>%
  mutate(
    from_state = state,
    to_state = lead(state)
  ) %>%
  ungroup() %>%
  left_join(trend_types, by = c("from_state" = "from", "to_state" = "to"))

# ============================================================================
# COVARIATE MODELS
# ============================================================================

cat("Fitting univariate covariate models...\n")

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

cat("Completed univariate covariate models\n")

# ============================================================================
# TRANSITION-SPECIFIC EFFECTS
# ============================================================================

cat("Fitting transition-specific effect models...\n")

# Select key covariates for transition-specific analysis
key_covariates <- c("age", "sex", "cci_score", "COVID_vax")

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

cat("Completed transition-specific effect models\n")

# ============================================================================
# TIME-VARYING RATE MODELS
# ============================================================================

cat("Fitting time-varying rate models...\n")

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

cat("Completed time-varying rate models\n")

# ============================================================================
# PREDICTIVE PERFORMANCE ANALYSIS
# ============================================================================

cat("Calculating predictive performance...\n")

# Combine all models for performance comparison
all_models <- c(
  fitted_base_models,
  unlist(covariate_models, recursive = FALSE),
  time_varying_models
)

# Calculate predictive performance using cross-validation
predictive_performance <- calculate_predictive_performance(
  fitted_models = all_models,
  patient_data = pt_stg_enhanced,
  k_folds = 5
)

# Save results
saveRDS(predictive_performance, here("data", "temp", "predictive_performance.rds"))

cat("Completed predictive performance analysis\n")

# ============================================================================
# SENSITIVITY ANALYSIS FOR OUTLIERS
# ============================================================================

cat("Performing sensitivity analysis for long-stay patients...\n")

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

cat("Completed sensitivity analysis\n")

# ============================================================================
# MODEL COMPARISONS
# ============================================================================

cat("Performing comprehensive model comparisons...\n")

# Compare models with same structure (different covariates)
same_structure_comparisons <- list(
  base_models = compare_same_structure_models(fitted_base_models),
  covariate_models = covariate_stats,
  time_varying_models = time_varying_stats,
  transition_specific_models = transition_stats
)

# Compare models with different structures
model_pairs <- expand_grid(
  model1 = c("base_model", "mod_2", "mod_3"),
  model2 = c("sev_2", "hx_sev", "reduced_trans")
) %>%
  filter(model1 != model2)

different_structure_comparisons <- compare_different_structure_models(
  fitted_models = fitted_base_models,
  model_pairs = model_pairs,
  cores = n.cores
)

# Save comparison results
saveRDS(same_structure_comparisons, here("data", "temp", "same_structure_comparisons.rds"))
saveRDS(different_structure_comparisons, here("data", "temp", "different_structure_comparisons.rds"))

cat("Completed model comparisons\n")

# ============================================================================
# DIAGNOSTIC ANALYSIS
# ============================================================================

cat("Calculating diagnostic metrics...\n")

# Calculate deviance residuals for key models
diagnostic_results <- map(names(fitted_base_models), function(model_name) {
  model_data <- pt_stg_enhanced %>% filter(model == model_name)
  fitted_model <- fitted_base_models[[model_name]]
  
  deviance_data <- calculate_transition_deviance(fitted_model, model_data)
  
  if (!is.null(deviance_data)) {
    deviance_data$model <- model_name
  }
  
  return(deviance_data)
}) %>%
  bind_rows()

# Save diagnostic results
saveRDS(diagnostic_results, here("data", "temp", "diagnostic_results.rds"))

cat("Completed diagnostic analysis\n")

# ============================================================================
# SAVE COMPREHENSIVE RESULTS SUMMARY
# ============================================================================

# Create a comprehensive summary
comprehensive_summary <- list(
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
  model_comparisons = list(
    same_structure = same_structure_comparisons,
    different_structure = different_structure_comparisons,
    sensitivity = sensitivity_comparison
  ),
  predictive_performance = predictive_performance,
  diagnostics = diagnostic_results,
  data_info = list(
    n_patients = length(unique(pt_stg_enhanced$deid_enc_id)),
    n_observations = nrow(pt_stg_enhanced),
    covariates_tested = time_invariant_covariates,
    continuous_covariates = continuous_covariates,
    long_stay_patients = length(sensitivity_results$long_stay_patients)
  )
)

# Save comprehensive summary
saveRDS(comprehensive_summary, here("data", "temp", "comprehensive_model_results.rds"))

cat("\n=== MODEL FITTING COMPLETED ===\n")
cat("Results saved to data/temp/ directory\n")
cat("Key files:\n")
cat("- comprehensive_model_results.rds: Complete analysis results\n")
cat("- covariate_models.rds: Univariate covariate models\n")
cat("- time_varying_models.rds: Time-dependent models\n")
cat("- sensitivity_results.rds: Outlier sensitivity analysis\n")
cat("- predictive_performance.rds: Cross-validation results\n")
cat("- diagnostic_results.rds: Model diagnostic metrics\n")

# Print summary statistics
cat("\nSUMMARY STATISTICS:\n")
cat("Number of patients:", comprehensive_summary$data_info$n_patients, "\n")
cat("Number of observations:", comprehensive_summary$data_info$n_observations, "\n")
cat("Number of long-stay patients (>30 days):", comprehensive_summary$data_info$long_stay_patients, "\n")
cat("Covariates tested:", length(comprehensive_summary$data_info$covariates_tested), "\n")

# Print convergence summary for covariate models
convergence_summary <- covariate_stats %>%
  group_by(covariate_type) %>%
  summarise(
    n_models = n(),
    n_converged = sum(converged, na.rm = TRUE),
    convergence_rate = round(n_converged / n_models * 100, 1)
  )

cat("\nCOVARIATE MODEL CONVERGENCE:\n")
print(convergence_summary)