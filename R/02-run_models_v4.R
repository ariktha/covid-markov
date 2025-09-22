rm(list = ls())

library(here)
library(tidyverse)
library(msm)

source(here("R", "00-config.R"))
source(here("R", "00-functions.R"))


# Data prep ---------------------------------------------------------------

if(run_data_setup){
  cat("\n=== Prepping data ===\n")
  
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
  
  source(here("R", "00-qmat_setup.R"))
  crude_rates <- calc_crude_init_rates(pt_stg, qmat_list)
  
  # Save initial data and specifications
  saveRDS(pt_stg, here("data", "temp", "pt_stg.rds"), compress = FALSE)
  saveRDS(crude_rates, here("data", "temp", "crude_rates.rds"))
  rm(dem_raw, stg_raw, model_specs_raw, model_specs, qmat_list)
  
} else {
  cat("\n=== Loading data ===\n")
  pt_stg <- readRDS(here("data", "temp", "pt_stg.rds"))
  # source(here("R", "00-qmat_setup.R"))
  crude_rates <- readRDS(here("data", "temp", "crude_rates.rds"))
}

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
  
  start_time <- Sys.time()
  pred <- calculate_predictive_performance(
    patient_data = pt_stg,
    fitted_models = fitted_base_models,
    crude_rates = crude_rates,
    k_folds = 5,
    prediction_times = seq(1, 15, by = 2),
    parallel = TRUE,
    n_cores = n.cores
  )
  end_time <- Sys.time()
  cat("Base models predictive performance calculation time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
  
  saveRDS(base_models_comp, here("data", "temp", "base_models_comp.rds"))
  saveRDS(pred, here("data", "temp", "base_model_predictive_performance.rds"))
  
} else {
  cat("\n=== Loading pre-compiled base model results ===\n")
  base_models_comp <- readRDS(here("data", "temp", "base_models_comp.rds"))
  base_model_pred <- readRDS(here("data", "temp", "base_model_predictive_performance.rds"))
}

gc()

# prev <- tidy_msm_prevalences(
#   fitted_base_models,
#   time_points = c(1, 10),
#   mc.cores = n.cores,
#   use_approximation = TRUE
# )

# Markov assumption testing ----------------------------------------------

cat("\n=== TESTING MARKOV ASSUMPTIONS ===\n")

if(fit_markov_models){
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
      "~ hx_sev_time" = markov_time_severe_ts[["base_model"]][["~ hx_sev_time"]]
    ),
    history_severe_structure = list(
      "~ 1" = fitted_base_models[["hx_sev"]][["~ 1"]]
    )
  )
  
  saveRDS(markov_comparison_models, here("data", "temp", "markov_comparison_models.rds"))
} else {
  markov_comparison_models <- readRDS(here("data", "temp", "markov_comparison_models.rds"))
}

if(comp_markov_models){
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
  
  # Create mapping with descriptive names
  markov_models_mapped <- list(
    "base_no_covariates" = markov_comparison_models[["base_no_covariates"]],
    "base_with_time_severe" = markov_comparison_models[["base_with_time_severe"]], 
    "base_with_time_severe_spline" = markov_comparison_models[["base_with_time_severe_spline"]],
    "base_with_time_severe_ts" = markov_comparison_models[["base_with_time_severe_ts"]],
    "history_severe_structure" = markov_comparison_models[["history_severe_structure"]]
  )
  
  # Map to data structure names but keep original names as nested level
  markov_for_cv <- list(
    "base_model" = list(
      "base_no_covariates" = markov_models_mapped[["base_no_covariates"]],
      "base_with_time_severe" = markov_models_mapped[["base_with_time_severe"]],
      "base_with_time_severe_spline" = markov_models_mapped[["base_with_time_severe_spline"]],
      "base_with_time_severe_ts" = markov_models_mapped[["base_with_time_severe_ts"]]
    ),
    "hx_sev" = list(
      "history_severe_structure" = markov_models_mapped[["history_severe_structure"]]
    )
  )
  
  markov_cv_results <- calculate_predictive_performance(
    patient_data = pt_stg,
    fitted_models = markov_for_cv,
    crude_rates = crude_rates,
    k_folds = 5,
    prediction_times = seq(1, 15, by = 2),
    parallel = TRUE,
    n_cores = n.cores
  )
  
  markov_pred <- calculate_predictive_performance(
    patient_data = pt_stg,
    fitted_models = markov_comparison_models,
    crude_rates = list(base_model = crude_rates[["base_model"]],
                       hx_sev = crude_rates[["hx_sev"]]),
    k_folds = 5,
    prediction_times = seq(1, 15, by = 2),
    parallel = TRUE,
    n_cores = 10
  )
  
  saveRDS(markov_model_comp, here("data", "temp", "markov_model_comp.rds"))
  # saveRDS(markov_prev, here("data", "temp", "markov_prev.rds"))
  saveRDS(markov_pred, here("data", "temp", "markov_pred.rds"))
} else {
  markov_model_comp <- readRDS(here("data", "temp", "markov_model_comp.rds"))
  # markov_prev <- readRDS(here("data", "temp", "markov_prev.rds"))
  markov_pred <- readRDS(here("data", "temp", "markov_pred.rds"))
}

gc()

# Covariate effects ------------------------------------------------------

cat("\n=== Fitting covariate models ===\n")

## Linear univariate models ---------------------------------------

if(fit_univar_models){
  
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

if(comp_univar_models){
  cat("\n=== Compile univariate model results ===\n")
  univariate_comp <- run_comprehensive_msm_analysis(
    univariate_models, 
    pt_stg, 
    crude_rates, 
    mc.cores = n.cores, 
    analysis_config = config_core_cov
  )
  saveRDS(univariate_comp, here("data", "temp", "univariate_comp.rds"))
} else {
  cat("\n=== Load univariate model results ===\n")
  univariate_comp <- readRDS(here("data", "temp", "univariate_comp.rds"))
}

# Tidy covariate results

## Spline univariate models ---------------------------------------

if(fit_spline_models){
  
  cat("\n=== Fitting Spline Models ===\n")
  cat("Constrain: effect on all transitions same")
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
          mc.cores = n.cores,
          constraint = "unconstrained"
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
  
  cat("Unconstrained")
  individual_spline_models_cons <- list()
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
      individual_spline_models_cons[[covar]] <- model_result
      
      # Optional: Save to disk immediately
      saveRDS(model_result, file = here("data", "temp", paste0("cons_spline_model_", covar, ".rds")))
      
      cat("Completed model for:", covar, "\n")
    } else {
      cat("Skipped model for:", covar, "\n")
    }
  }
  
  spline_models_cons <- do.call(c, individual_spline_models_cons)
  saveRDS(spline_models_cons, here("data", "temp", "spline_models_cons.rds"))
  
  
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
  
  spline_comp_cons <- run_comprehensive_msm_analysis(
    spline_models_cons, 
    pt_stg, 
    crude_rates = list(base_model = crude_rates[["base_model"]]), 
    mc.cores = n.cores, 
    analysis_config = config_core_cov
  )
  saveRDS(spline_comp_cons, here("data", "temp", "spline_comp_uncons.rds"))
  
} else {
  cat("\n=== LOADING PRE-COMPILED SPLINE MODEL RESULTS ===\n")
  
  spline_comp <- readRDS(here("data", "temp", "spline_comp.rds"))
}

if(hr_spline_models){
  
  spline_hrs <- calculate_spline_hazard_ratios(
    spline_models,
    n_breaks = 5
  )
  saveRDS(spline_hrs, here("data", "temp", "spline_hrs.rds"))
  
  spline_effects <- extract_all_spline_effects(
    models_list = spline_models,
    mc.cores = n.cores
  )
  saveRDS(spline_effects, here("data", "temp", "spline_effects.rds"))
  
  spline_hrs_cons <- calculate_spline_hazard_ratios(
    spline_models_cons,
    n_breaks = 5
  )
  saveRDS(spline_hrs_cons, here("data", "temp", "spline_hrs_cons.rds"))
  
  spline_effects_cons <- extract_all_spline_effects(
    models_list = spline_models_cons,
    mc.cores = n.cores
  )
  saveRDS(spline_effects_cons, here("data", "temp", "spline_effects_cons.rds"))
  
} else {
  spline_hrs <- readRDS(here("data", "temp", "spline_hrs.rds"))
  spline_effects <- readRDS(here("data", "temp", "spline_effects.rds"))
}


## Transition-specific effects ---------------------------------------

if(fit_trans_models){
  cat("\n=== Fitting transition-specific models ===\n")
  transition_specific_models <- fit_transition_specific_models(
    patient_data = pt_stg,
    crude_rates = list(base_model = crude_rates[["base_model"]]),
    covariates = setNames(as.list(key_covariates), paste0("univar_", key_covariates))
  )
  saveRDS(transition_specific_models, here("data", "temp", "transition_specific_models.rds"))
} else {
  cat("\n=== Loading pre-fitted transition-specific models ===\n")
  transition_specific_models <- readRDS(here("data", "temp", "transition_specific_models.rds"))
}

if(comp_trans_models){
  cat("\n=== Compiling transition-specific model results ===\n")
  transition_specific_comp <- run_comprehensive_msm_analysis(
    transition_specific_models, 
    pt_stg, 
    crude_rates = list(base_model = crude_rates[["base_model"]]), 
    mc.cores = n.cores, 
    analysis_config = config_core_cov
  )
  saveRDS(transition_specific_comp, here("data", "temp", "transition_specific_comp.rds"))
} else {
  cat("\n=== Loading pre-compiled transition-specific model results ===\n")
  transition_specific_comp <- readRDS(here("data", "temp", "transition_specific_comp.rds"))
}

## Multivariable models ---------------------------------------

# Multivariate selection (base model only for efficiency)

if(fit_multivar_models){
  cat("\n=== Fitting multivariable models ===\n")
  multivariate_models <- multivariate_selection(
    patient_data = pt_stg %>% filter(model == "base_model"),
    crude_rates = list(base_model = crude_rates[["base_model"]]),
    candidate_covariates = key_covariates,
    # required_covariates = c("age", "BMI", "cci_cat"),
    method = "forward",
    alpha_enter = 0.05,
    # max_variables = 4,
    max_variables = length(key_covariates),
    n_cores = n.cores,
    save_intermediate = TRUE
  )
  saveRDS(multivariate_models, here("data", "temp", "multivariate_models.rds"))
} else {
  cat("\n=== Loading pre-fitted multivariable models ===\n")
  multivariate_models <- readRDS(here("data", "temp", "multivariate_models.rds"))
} 

if(comp_multivar_models){
  cat("\n=== Compiling multivariable model results ===\n")
  multivariate_comp <- run_comprehensive_msm_analysis(
    multivariate_models, 
    pt_stg %>% filter(model == "base_model"), 
    crude_rates = list(base_model = crude_rates[["base_model"]]), 
    mc.cores = n.cores, 
    analysis_config = config_core_cov
  )
  saveRDS(multivariate_comp, here("data", "temp", "multivariate_comp.rds"))
} else {
  cat("\n=== Loading pre-compiled multivariable model results ===\n")
  multivariate_comp <- readRDS(here("data", "temp", "multivariate_comp.rds"))
}

gc()

# Time-varying transition rates ------------------------------------------

if(fit_time_vary_models){
  cat("\n=== Running time-varying model fitting ===\n")
  # Hospital time effects
  hospital_time_models <- fit_time_varying_models(
    patient_data = pt_stg,
    crude_rates = crude_rates,
    time_term = "days_since_entry",
    time_covariates = c("linear", "piecewise", "spline"),
    spline_df = 3,
    piecewise_breakpoints = seq(0, max(pt_stg$DaysSinceEntry), by = 10)
  )
  
  # Calendar time effects  
  calendar_time_models <- fit_time_varying_models(
    patient_data = pt_stg,
    crude_rates = crude_rates,
    time_term = "calendar_time",
    time_covariates = c("linear", "piecewise", "spline"),
    spline_df = 3,
    piecewise_breakpoints = seq(0, max(pt_stg$CalendarTime), by = 90)
  )
  
  saveRDS(hospital_time_models, here("data", "temp", "hospital_time_models.rds"))
  saveRDS(calendar_time_models, here("data", "temp", "calendar_time_models.rds"))
} else {
  cat("\n=== Loading pre-fitted time-varying models ===\n")
  hospital_time_models <- readRDS(here("data", "temp", "hospital_time_models.rds"))
  calendar_time_models <- readRDS(here("data", "temp", "calendar_time_models.rds"))
}

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
# saveRDS(calendar_wave_models, here("data", "temp", "calendar_wave_models.rds"))

if(comp_time_vary_models){
  cat("\n=== Compiling time-varying model results ===\n")
  
  hospital_time_comp <- run_comprehensive_msm_analysis(
    hospital_time_models, 
    pt_stg, 
    crude_rates = crude_rates, 
    mc.cores = n.cores, 
    analysis_config = config_core
  )
  
  calendar_time_comp <- run_comprehensive_msm_analysis(
    calendar_time_models, 
    pt_stg, 
    crude_rates = crude_rates, 
    mc.cores = n.cores, 
    analysis_config = config_core
  )
  saveRDS(hospital_time_comp, here("data", "temp", "hospital_time_comp.rds"))
  saveRDS(calendar_time_comp, here("data", "temp", "calendar_time_comp.rds"))
  
  hospital_plot_data <- extract_time_tis_and_hrs(
    hospital_time_models,
    pt_stg,
    time_interval = 14,
    mc.cores = n.cores
  )
  
  saveRDS(hospital_plot_data, here("data", "temp", "hospital_time_plot_data.rds"))
  
  calendar_plot_data <- extract_time_tis_and_hrs(
    calendar_time_models,
    pt_stg,
    time_interval = 14,
    mc.cores = n.cores
  )
  
  saveRDS(calendar_plot_data, here("data", "temp", "calendar_time_plot_data.rds"))
  
} else {
  cat("\n=== Loading pre-compiled time-varying model results ===\n")
  hospital_time_comp <- readRDS(here("data", "temp", "hospital_time_comp.rds"))
  calendar_time_comp <- readRDS(here("data", "temp", "calendar_time_comp.rds"))
}

gc()

# Long-stay patient sensitivity analysis --------------------------------

if(run_long_stay_analysis){
  cat("\n=== Running long-stay patient sensitivity analysis ===\n")
  source(here("R", "06-long-stay.R"))
}

cat("\n=== Loading long-stay patient sensitivity analysis ===\n")

long_stay_sensitivity <- readRDS(here("data", "temp", "long_stay_sensitivity.rds"))
long_stay_models <- readRDS(here("data", "temp", "long_stay_models.rds"))
# long_stay_model_comp <- readRDS(here("data", "temp", "long_stay_model_comp.rds"))


# Graveyard ---------------------------------------------------------------

## Proof of concept --------------------------------------------------------

# poc_rates <- list(
#   base_model = crude_rates[["base_model"]],
#   hx_sev = crude_rates[["hx_sev"]]
# )
# 
# poc_models <- list(
#   base_model = fitted_base_models[["base_model"]],
#   hx_sev = fitted_base_models[["hx_sev"]]
# )
# 
# poc_data <- pt_stg %>% filter(model %in% c("base_model", "hx_sev"))
# 
# poc_prev <- tidy_msm_prevalences(
#   poc_models,
#   time_points = seq(1, 12, by = 3),
#   mc.cores = 1,
#   ci = TRUE
# )
# 
# poc_pred <- calculate_predictive_performance(
#   patient_data = poc_data,
#   fitted_models = poc_models,
#   crude_rates = poc_rates,
#   k_folds = 3,
#   prediction_times = c(7, 14),
#   parallel = TRUE,
#   n_cores = 6
# )
# 
# poc_residuals <- calculate_transition_residuals(
#   fitted_msm_models = poc_models,
#   patient_data = poc_data,
#   residual_type = "deviance",
#   debug = FALSE
# )
# 
# saveRDS(poc_prev, here("data", "temp", "poc_prevalences.rds"))
# saveRDS(poc_pred, here("data", "temp", "poc_predictive_performance.rds"))
# saveRDS(poc_residuals, here("data", "temp", "poc_transition_residuals.rds"))

## Predictive performance evaluation --------------------------------------
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
## State prevalence evaluation --------------------------------------------
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
## Residual diagnostics ---------------------------------------------------
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

