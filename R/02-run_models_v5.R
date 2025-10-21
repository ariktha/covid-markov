rm(list = ls())

library(here)
library(tidyverse)
library(msm)

source(here("R", "00-config.R"))
# source(here("R", "00-functions.R"))
source(here("R", "00a-fns-model_fitting.R"))
source(here("R", "00b-fns-model_eval.R"))
source(here("R", "00c-fns-predictive_performance.R"))

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

if(do_no_cov){
  
  if(fit_no_cov_models){
    
    cat("\n=== Fitting base models (No covariates) ===\n")
    start_time <- Sys.time()
    
    base_model_specs <- tibble(
      model_name = names(crude_rates),
      model_structure = names(crude_rates)
    )
    
    base_models <- fit_msm_models(
      patient_data = pt_stg, 
      crude_rates = crude_rates,
      model_specs = base_model_specs,
      covariates = NULL,
      spline_vars = NULL,
      time_varying = NULL,
      constraint = "transition_specific",
      mc.cores = n.cores
    )
    
    end_time <- Sys.time()
    cat("Base models fitting time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
    saveRDS(base_models, here("data", "temp", "base_models.rds"))
    
  } else {
    cat("\n=== Loading pre-fitted base models (No covariates) ===\n")
    base_models <- readRDS(here("data", "temp", "base_models.rds"))
  }

  if(comp_no_cov_models){
    cat("\n=== Compiling base model results ===\n")
    start_time <- Sys.time()
    base_models_comp <- run_comprehensive_msm_analysis(
      base_models, pt_stg, crude_rates, mc.cores = n.cores, 
      analysis_config = config_core
    )
    end_time <- Sys.time()
    cat("Base models results compilation time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
    
    saveRDS(base_models_comp, here("data", "temp", "base_models_comp.rds"))
    
    cat("Model structure comparison:\n")
    start_time <- Sys.time()
    base_model_structure_comp <- compare_across_structures(
      base_models,
      methods = c("draic", "drlcv"),
      mc.cores = n.cores
    )
    end_time <- Sys.time()
    cat("Base models structure comparison time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
    
    saveRDS(base_model_structure_comp, here("data", "temp", "base_model_structure_comp.rds"))
  }

  gc()
}

base_pred <- cv_msm_models(fitted_models = base_models, 
                           patient_data = pt_stg, 
                           k_folds = 5,
                           parallel = TRUE,
                           n_cores = n.cores)

# Check Markov assumption --------------------------------------------------

if(do_markov){
  
  if(fit_markov_models){
    
    markov_model_specs <- tibble(
      model_name = c("base_time_severe", "base_time_severe_spline", 
                     "base_time_severe_ts"),
      model_structure = c("base_model", "base_model", "base_model"),
      covariates = list(c("hx_sev_time"), NULL, c("hx_sev_time")),
      spline_vars = list(NULL, c("hx_sev_time"), NULL),
      spline_df = c(NA, 3, NA),
      spline_type = c(NA, "ns", NA),
      constraint = c("transition_specific", "transition_specific", 
                     "constant_across_transitions")
    )
    
    markov_comparison_models <- fit_msm_models(
      patient_data = pt_stg,
      crude_rates = crude_rates,
      model_specs = markov_model_specs,
      mc.cores = n.cores
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
    
    saveRDS(markov_model_comp, here("data", "temp", "markov_model_comp.rds"))
  } else {
    markov_model_comp <- readRDS(here("data", "temp", "markov_model_comp.rds"))
  }
  
  gc()
}


# Covariate effects -------------------------------------------------------

## Univariate models: Linear and spline covariate effects -----------------

if(do_univar){
  if(fit_univar_models){
    
    univar_results <- univariate_progressive_timeout(
      patient_data = subset(pt_stg, model == "base_model"),
      crude_rates = list(base_model = crude_rates[["base_model"]]),
      # covariates = c("age", "sex"),
      # spline_vars = c("age"),
      covariates = key_covariates,
      spline_vars = continuous_covariates,
      n.cores = n.cores,
      timeout_vector = univar_timeout_times,
      # timeout_vector = c(1),
      save_prefix = "univar_progressive"
    )
    
  } else {
    univar_results <- readRDS(here("data", "temp", "univar_progressive_full_results.rds"))
  }
  
  univariate_models <- univar_results$models
  failed_summary <- univar_results$failed_models
  
  if(refit_failed_univar_models & nrow(failed_summary) > 0){
    cat("\n=== Refitting failed univariate models with extended timeout ===\n")
    
    retry_results <- retry_failed_models(
      failed_results = univar_results,
      patient_data = pt_stg,
      crude_rates = crude_rates,
      n.cores = 10,
      save_prefix = "univar_retry"
    )
    
    # Combine successful models from both runs
    all_univar <- combine_univariate_results(univar_results, retry_results)
    saveRDS(all_univar, here("data", "temp", "univar_and_retried_full_results.rds"))
  }
  
  if(comp_univar_models){
    cat("\n=== Compile univariate model results ===\n")
    
    # Use combined results if available, otherwise use progressive timeout results
    if (exists("combined_results")) {
      cat("Using combined results (progressive timeout + retry)\n")
      univariate_models <- combined_results$models
    } else if (exists("univar_results")) {
      cat("Using progressive timeout results only\n")
      univariate_models <- univar_results$models
    } else {
      stop("Neither combined_results nor univar_results found in environment")
    }
    
    univariate_comp <- run_comprehensive_msm_analysis(
      univariate_models, 
      subset(pt_stg, model == "base_model"),
      crude_rates = list(base_model = crude_rates[["base_model"]]), 
      mc.cores = n.cores, 
      analysis_config = config_core
    )
    saveRDS(univariate_comp, here("data", "temp", "univariate_comp.rds"))
  } else {
    cat("\n=== Load univariate model results ===\n")
    univariate_comp <- readRDS(here("data", "temp", "univariate_comp.rds"))
  }
  
  gc()
}

## Univariate models: Spline covariate effects ------------------------------

## Univariate models: Constrained transition effects ------------------------

if(do_constrained_trans){
  
  if(fit_trans_models){
    
    const_univar_results <- univariate_progressive_timeout(
      patient_data = subset(pt_stg, model == "base_model"),
      crude_rates = list(base_model = crude_rates[["base_model"]]),
      # covariates = c("age", "sex"),
      # spline_vars = c("age"),
      covariates = key_covariates,
      spline_vars = continuous_covariates,
      n.cores = n.cores,
      timeout_vector = univar_timeout_times,
      constraint = "constant_across_transitions",
      # timeout_vector = c(1),
      save_prefix = "const_univar_progressive"
    )
    
  } else{
    const_univar_results <- readRDS(here("data", "temp", "const_univar_progressive_full_results.rds"))
  }
  
  if(comp_trans_models){

    const_univariate_comp <- run_comprehensive_msm_analysis(
      const_univar_results$models, 
      subset(pt_stg, model == "base_model"), 
      list(base_model = crude_rates[["base_model"]]), 
      mc.cores = n.cores, 
      analysis_config = config_core
    )
    saveRDS(const_univariate_comp, here("data", "temp", "const_univariate_comp.rds"))
    
  } else {
    const_univariate_comp <- readRDS(here("data", "temp", "const_univariate_comp.rds"))
  }
  
  gc()
}

## Multivariable models -----------------------------------------------------

# Time-inhomogenous models --------------------------------------------------

if(do_time_inhomogeneity){  
  
  if(fit_time_vary_models){
    cat("\n=== Fitting time-inhomogeneous models ===\n")
    
    time_vary_models <- run_all_time_models(
      patient_data = subset(pt_stg, model == "base_model"),
      crude_rates = list(base_model = crude_rates[["base_model"]]),
      constraint = "transition_specific",
      mc.cores = n.cores
    )
    
    saveRDS(time_vary_models, here("data", "temp", "time_vary_models.rds"))
    
  } else {
    cat("\n=== Loading pre-fitted time-inhomogeneous models ===\n")
    time_vary_models <- readRDS(here("data", "temp", "time_vary_models.rds"))
  }
  
  if(comp_time_vary_models){
    cat("\n=== Compiling time-inhomogeneous model results ===\n")
    
    time_vary_comp <- run_comprehensive_msm_analysis(
      time_vary_models, 
      subset(pt_stg, model == "base_model"), 
      list(base_model = crude_rates[["base_model"]]), 
      mc.cores = n.cores, 
      analysis_config = config_core
    )
    
    saveRDS(time_vary_comp, here("data", "temp", "time_vary_comp.rds"))
    
  } else {
    cat("\n=== Loading compiled time-inhomogeneous model results ===\n")
    time_vary_comp <- readRDS(here("data", "temp", "time_vary_comp.rds"))
  }

  gc()
}

# Long-stay patient sensitivity analysis --------------------------------

if(run_long_stay_analysis){
  cat("\n=== Running long-stay patient sensitivity analysis ===\n")
  source(here("R", "06-long-stay.R"))
}

# cat("\n=== Loading long-stay patient sensitivity analysis ===\n")
# 
# long_stay_sensitivity <- readRDS(here("data", "temp", "long_stay_sensitivity.rds"))
# long_stay_models <- readRDS(here("data", "temp", "long_stay_models.rds"))
# long_stay_model_comp <- readRDS(here("data", "temp", "long_stay_model_comp.rds"))

gc()
