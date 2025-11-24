rm(list = ls())

library(here)
library(tidyverse)
library(msm)

source(here("R", "00-config.R"))
source(here("R", "00c-fns-model_fitting.R"))
source(here("R", "00d-fns-model_eval.R"))
source(here("R", "00e-fns-predictive_performance.R"))

# Workflow config --------------------------------------------------------------

workflow_config <- list(

  # Which analysis sections to run
  sections = c(
    # "base"           # Base models (no covariates)
    "markov"         # Markov assumption check
    # "univar",         # Univariate covariate effects
    # "univar_const",   # Univariate with constant transitions
    # "univar_deep_dive"
    # "multivar"     # Multivariable models (not implemented)
    # "time_inhomog"    # Time-inhomogeneous models
    # "long_stay"     # Long-stay sensitivity analysis
  ),
  
  # Control for each section: fit, compile, compare, cv
  base = list(
    fit = FALSE,
    compile = FALSE,
    compare = FALSE,
    cv = FALSE
  ),
  
  markov = list(
    fit = TRUE,
    compile = TRUE,
    compare = TRUE,
    cv = TRUE
  ),
  
  univar = list(
    fit = FALSE,
    compile = FALSE,
    retry_failed = FALSE,  # Unlimited-time retry for failed models
    compare = FALSE,
    cv = FALSE
  ),
  
  univar_const = list(
    fit = FALSE,
    compile = FALSE,
    retry_failed = FALSE,
    compare = FALSE,
    cv = FALSE
  ),
  
  univar_deep_dive = list(
    fit = FALSE,
    compile = FALSE,
    compare = FALSE,
    cv = FALSE
  ),
  
  time_inhomog = list(
    fit = FALSE,
    compile = FALSE,
    compare = FALSE,
    cv = FALSE
  ),
  
  long_stay = list(
    fit = FALSE,
    compile = FALSE,
    compare = FALSE,
    cv = FALSE
  ),
  
  # Model comparison configuration
  comparisons = list(
    # Compare base models across structures
    base_structures = list(
      run = FALSE,
      methods = c("draic", "drlcv")
    ),
    
    # Compare markov models within base structure
    markov_within = list(
      run = TRUE,
      methods = c("lrt")  # Likelihood ratio test for nested models
    ),
    
    # Compare univar models to base
    univar_to_base = list(
      run = FALSE,
      methods = c("draic")
    )
  )
)


## Analysis configurations for run_comprehensive_msm_analysis ------------------

analysis_configs <- list(
  
  base = list(
    tidy_models = list(),
    qmatrix = list(covariates_list = NULL, 
                   mc.cores = n.cores,
                   ci_method = "bootstrap"),
    pmats = list(t_values = time_vec, 
                 covariates_list = NULL, 
                 mc.cores = n.cores,
                 ci_method = "bootstrap"),
    sojourns = list(covariates_list = NULL,
                    ci_method = "bootstrap"),
    prevalence = list(
      time_points = time_vec,
      covariates_list = NULL,
      ci = TRUE,
      ci_method = "bootstrap",
      use_approximation = FALSE,
      mc.cores = n.cores
    ),
    hazards = list(skip = TRUE),  # No covariates
    cv = list(skip = TRUE),  # Handled separately
    residuals = list(residual_type = "deviance", debug = FALSE),
    auto_covariates = FALSE
  ),
  
  markov = list(
    tidy_models = list(),
    qmatrix = list(covariates_list = NULL, mc.cores = n.cores,
                   ci_method = "bootstrap"),
    pmats = list(t_values = time_vec, covariates_list = NULL, mc.cores = n.cores,
                 ci_method = "bootstrap"),
    sojourns = list(covariates_list = NULL,
                    ci_method = "bootstrap"),
    prevalence = list(
      time_points = time_vec,
      covariates_list = NULL,
      ci = TRUE,
      # ci_method = "normal",
      use_approximation = FALSE,
      mc.cores = n.cores,
      ci_method = "normal"
    ),
    hazards = list(hazard_scale = 1),
    cv = list(skip = TRUE),
    residuals = list(residual_type = "deviance", debug = FALSE)
  ),
  
  univar = list(
    tidy_models = list(),
    qmatrix = list(
      # covariates_list = NULL, mc.cores = n.cores, ci_method = "bootstrap",
      skip = TRUE
      ),
    pmats = list(
      # t_values = time_vec, covariates_list = NULL, mc.cores = n.cores, ci_method = "bootstrap",
      skip = TRUE
      ),
    sojourns = list(
      # covariates_list = NULL, ci_method = "bootstrap",
      skip = TRUE
      ),
    prevalence = list(
      skip = TRUE
      ),  # Too many models
    hazards = list(
      hazard_scale = 1
      ),
    cv = list(
      skip = TRUE
      ),
    residuals = list(
      skip = TRUE
      )
  ),
  
  univar_const = list(
    tidy_models = list(),
    qmatrix = list(
      # covariates_list = NULL, mc.cores = n.cores, ci_method = "bootstrap",
      skip = TRUE
    ),
    pmats = list(
      # t_values = time_vec, covariates_list = NULL, mc.cores = n.cores, ci_method = "bootstrap",
      skip = TRUE
    ),
    sojourns = list(
      # covariates_list = NULL, ci_method = "bootstrap",
      skip = TRUE
    ),
    prevalence = list(
      skip = TRUE
    ),  # Too many models
    hazards = list(
      hazard_scale = 1
    ),
    cv = list(
      skip = TRUE
    ),
    residuals = list(
      skip = TRUE
    )
  ),
  
  univar_deep_dive = list(
    tidy_models = list(),
    qmatrix = list(
      covariates_list = NULL, mc.cores = n.cores, ci_method = "bootstrap"
      # skip = TRUE
    ),
    pmats = list(
      t_values = time_vec, covariates_list = NULL, mc.cores = n.cores, ci_method = "bootstrap"
      # skip = TRUE
    ),
    sojourns = list(
      covariates_list = NULL, ci_method = "bootstrap"
      # skip = TRUE
    ),
    prevalence = list(
      skip = TRUE
    ),  # Too many models
    hazards = list(
      hazard_scale = 1
    ),
    cv = list(
      skip = TRUE
    ),
    residuals = list(
      # skip = TRUE
    )
  ),
  
  time_inhomog = list(
    tidy_models = list(),
    qmatrix = list(covariates_list = NULL, mc.cores = n.cores),
    pmats = list(t_values = time_vec, covariates_list = NULL, mc.cores = n.cores),
    sojourns = list(covariates_list = NULL),
    prevalence = list(
      time_points = time_vec,
      covariates_list = NULL,
      ci = TRUE,
      ci_method = "normal",
      use_approximation = FALSE,
      mc.cores = n.cores
    ),
    hazards = list(skip = TRUE),
    cv = list(skip = TRUE),
    residuals = list(residual_type = "deviance", debug = FALSE),
    auto_covariates = FALSE
  )
)


## CV configurations -----------------------------------------------------------
cv_configs <- list(
  base = list(
    k_folds = 5,
    stratify_by = "final_state",
    prediction_horizon = 365,
    output_dir = here("data", "cv_results", "base"),
    parallel = TRUE,
    n_cores = n.cores
  ),
  
  markov = list(
    k_folds = 5,
    stratify_by = "final_state",
    prediction_horizon = 365,
    output_dir = here("data", "cv_results", "markov"),
    parallel = TRUE,
    n_cores = n.cores
  ),
  
  univar = list(
    k_folds = 5,
    stratify_by = "final_state",
    prediction_horizon = 365,
    output_dir = here("data", "cv_results", "univar"),
    parallel = TRUE,
    n_cores = n.cores
  ),
  
  univar_const = list(
    k_folds = 5,
    stratify_by = "final_state",
    prediction_horizon = 365,
    output_dir = here("data", "cv_results", "univar_const"),
    parallel = TRUE,
    n_cores = n.cores
  ),
  
  univar_deep_dive = list(
    k_folds = 5,
    stratify_by = "final_state",
    prediction_horizon = 365,
    output_dir = here("data", "cv_results", "univar_deep_dive"),
    parallel = TRUE,
    n_cores = n.cores
  ),
  
  time_inhomog = list(
    k_folds = 5,
    stratify_by = "final_state",
    prediction_horizon = 365,
    output_dir = here("data", "cv_results", "time_inhomog"),
    parallel = TRUE,
    n_cores = n.cores
  )
)


# Data prep --------------------------------------------------------------------

cat("\n=== LOADING PREPARED DATA ===\n")
pt_stg <- readRDS(here("data", "pt_stg.rds"))
crude_rates <- readRDS(here("data", "crude_rates.rds"))

cat("Data loaded.\n")
cat("  N patients:", length(unique(pt_stg$deid_enc_id)), "\n")
cat("  N observations:", nrow(pt_stg), "\n")

# Base model structures (no covariates) ----------------------------------------

if ("base" %in% workflow_config$sections) {
  
  base_models <- run_section(
    section_name = "base",
    section_config = workflow_config$base,
    fit_fn = function() {
      
      base_model_specs <- tibble(
        model_name = names(crude_rates),
        model_structure = names(crude_rates)
      )
      
      fit_msm_models(
        patient_data = pt_stg,
        crude_rates = crude_rates,
        model_specs = base_model_specs,
        mc.cores = n.cores
      )
    }
  )
  
  # Compile results
  if (workflow_config$base$compile && !is.null(base_models)) {
    base_comp <- compile_section("base", base_models, analysis_configs$base)
  }
  
  # Cross-validation
  if (workflow_config$base$cv && !is.null(base_models)) {
    base_cv <- run_cv_section("base", base_models, cv_configs$base)
  }
  
  # Compare base models across structures
  if (workflow_config$comparisons$base_structures$run) {
    
    cat("\n=== COMPARING BASE MODELS ACROSS STRUCTURES ===\n")
    
    if (exists("base_models") && !is.null(base_models)) {
      
      # Create all pairwise comparisons
      model_structures <- names(base_models)
      
      if (length(model_structures) > 1) {
        
        model_pairs <- combn(model_structures, 2, simplify = FALSE)
        
        base_structure_comp <- tryCatch({
          compare_across_structures(
            fitted_models = base_models,
            model_pairs = model_pairs,
            methods = workflow_config$comparisons$base_structures$methods,
            mc.cores = n.cores
          )
        }, error = function(e) {
          cat("ERROR in base structure comparison:", e$message, "\n")
          NULL
        })
        
        if (!is.null(base_structure_comp)) {
          saveRDS(base_structure_comp, 
                  here("data", "base_structure_comp.rds"))
          cat("Saved to: base_structure_comp.rds\n")
        }
      } else {
        cat("Only one base model structure - skipping comparison\n")
      }
    } else {
      cat("Base models not available for comparison\n")
    }
  }
  
  gc()
}

# Markov assumption check ------------------------------------------------------

if ("markov" %in% workflow_config$sections) {
  
  markov_models <- run_section(
    section_name = "markov",
    section_config = workflow_config$markov,
    fit_fn = function() {
      
      markov_model_specs <- tibble(
        model_name = c(
          # "base_no_cov",
          # "markov_state",
          "base_time_severe", 
          "base_time_severe_spline",
          "base_time_severe_const",
          "base_time_state",
          "base_time_state_spline",
          "base_time_state_const"
        ),
        model_structure = c(
          # "base_model",
          # "hx_sev",
          "base_model", 
          "base_model",
          "base_model", 
          "base_model", 
          "base_model", 
          "base_model"
        ),
        covariates = list(
          # NULL,
          # NULL,
          c("hx_sev_time"),
          NULL,
          c("hx_sev_time"),
          c("time_in_current_state"),
          NULL,
          c("time_in_current_state")
        ),
        spline_vars = list(
          # NULL,
          # NULL,
          NULL,
          c("hx_sev_time"),
          NULL,
          NULL,
          c("time_in_current_state"),
          NULL
        ),
        spline_df = c(
          # NA, 
          # NA, 
          NA, 
          3, 
          NA,
          NA,
          3,
          NA
          ),
        spline_type = c(
          # NA, 
          # NA, 
          NA, 
          "ns", 
          NA,
          NA,
          "ns",
          NA
          ),
        constraint = c(
          # "transition_specific",
          # "transition_specific",
          "transition_specific",
          "transition_specific",
          "constant_across_transitions",
          "transition_specific",
          "transition_specific",
          "constant_across_transitions"
        )
      )
      
      fit_msm_models(
        patient_data = pt_stg,
        crude_rates = crude_rates,
        model_specs = markov_model_specs,
        mc.cores = n.cores
      )
    }
  )
  
  # Compile results
  if (workflow_config$markov$compile && !is.null(markov_models)) {
    markov_comp <- compile_section("markov", markov_models, analysis_configs$markov)
  }
  
  # Cross-validation
  if (workflow_config$markov$cv && !is.null(markov_models)) {
    markov_cv <- run_cv_section("markov", markov_models, cv_configs$markov)
  }
  
  # Compare Markov models within structure
  if (workflow_config$comparisons$markov_within$run) {
    
    cat("\n=== COMPARING MARKOV MODELS (WITHIN-STRUCTURE) ===\n")
    
    if (exists("markov_models") && !is.null(markov_models)) {
      
      markov_within_comp <- tryCatch({
        compare_within_structure(markov_models)
      }, error = function(e) {
        cat("ERROR in markov within-structure comparison:", e$message, "\n")
        NULL
      })
      
      if (!is.null(markov_within_comp)) {
        saveRDS(markov_within_comp, 
                here("data", "markov_within_comp.rds"))
        cat("Saved to: markov_within_comp.rds\n")
        
        # Print summary
        if (nrow(markov_within_comp) > 0) {
          cat("\nMarkov Model Rankings (by AIC):\n")
          print(markov_within_comp %>% 
                  dplyr::select(model_name, model_structure, AIC, BIC, aic_rank) %>%
                  arrange(aic_rank))
        }
      }
    } else {
      cat("Markov models not available for comparison\n")
    }
  }
  
  gc()
}

# Univariate models ------------------------------------------------------------

## Transition-specific parameters ---------------------------------------------------

if ("univar" %in% workflow_config$sections) {
  
  if (workflow_config$univar$fit) {
    cat("Fitting univariate models (progressive timeout)\n")
    univar_results <- univariate_progressive_timeout(
      patient_data = pt_stg,
      crude_rates = crude_rates,
      covariates = key_covariates,
      spline_vars = continuous_covariates,
      # covariates = c("age"),
      # spline_vars = c("age"),
      model_structures = "base_model",  # ADD THIS LINE
      n.cores = n.cores,
      timeout_vector = c(1, 5, 10),
      constraint = "transition_specific",
      save_prefix = "univar_progressive"
    )
    
    univar_models <- univar_results$models
    
    # Optional unlimited-time retry for failed models
    if (workflow_config$univar$retry_failed && 
        nrow(univar_results$failed_models) > 0) {
      
      cat("\n=== RETRYING FAILED UNIVARIATE MODELS (UNLIMITED TIME) ===\n")
      
      retry_results <- retry_failed_models(
        failed_results = univar_results,
        patient_data = pt_stg,
        crude_rates = crude_rates,
        n.cores = n.cores,
        save_prefix = "univar_retry"
      )
      
      # Simple combination instead of undefined function
      # Merge successful retry results into original results
      for (model_name in names(retry_results$models)) {
        if (is.null(univar_models[[model_name]])) {
          univar_models[[model_name]] <- list()
        }
        for (formula_name in names(retry_results$models[[model_name]])) {
          univar_models[[model_name]][[formula_name]] <- retry_results$models[[model_name]][[formula_name]]
        }
      }
      
      # Create combined results structure
      combined_results <- list(
        models = univar_models,
        original_failed = univar_results$failed_models,
        retry_failed = retry_results$retry_failed,
        settings = univar_results$settings
      )
      saveRDS(combined_results, here("data", "univar_combined.rds"))
    }
    
    # Save final models
    saveRDS(univar_models, here("data", "univar_models.rds"))
    
  } else {
    cat("\n=== LOADING: UNIVARIATE MODELS ===\n")
    univar_models <- readRDS(here("data", "univar_models.rds"))
  }
  
  # Compile results
  if (workflow_config$univar$compile && !is.null(univar_models)) {
    univar_comp <- compile_section("univar", univar_models, analysis_configs$univar)
  }
  
  # Cross-validation
  if (workflow_config$univar$cv && !is.null(univar_models)) {
    univar_cv <- run_cv_section("univar", univar_models, cv_configs$univar)
  }
  
  # Compare univariate models to base
  if (workflow_config$comparisons$univar_to_base$run) {
    
    cat("\n=== COMPARING UNIVARIATE MODELS TO BASE ===\n")
    
    if (exists("univar_models") && !is.null(univar_models) &&
        exists("base_models") && !is.null(base_models)) {
      
      # Get base model for base_model structure
      base_for_comparison <- base_models[["base_model"]]
      
      if (!is.null(base_for_comparison)) {
        
        # Select top univariate models for comparison (by AIC)
        if (exists("univar_comp") && !is.null(univar_comp$model_summary)) {
          
          top_univar <- univar_comp$model_summary %>%
            filter(status == "converged") %>%
            arrange(AIC) %>%
            head(10) %>%
            pull(model_name)
          
          if (length(top_univar) > 0) {
            
            # Create comparison pairs: each top univar vs base
            comparison_pairs <- lapply(top_univar, function(x) {
              c("base_model", x)
            })
            
            # Combine models for comparison
            models_for_comparison <- c(
              list(base_model = base_for_comparison),
              univar_models[top_univar]
            )
            
            univar_base_comp <- tryCatch({
              compare_across_structures(
                fitted_models = models_for_comparison,
                model_pairs = comparison_pairs,
                methods = workflow_config$comparisons$univar_to_base$methods,
                mc.cores = n.cores
              )
            }, error = function(e) {
              cat("ERROR in univar-to-base comparison:", e$message, "\n")
              NULL
            })
            
            if (!is.null(univar_base_comp)) {
              saveRDS(univar_base_comp, 
                      here("data", "univar_base_comp.rds"))
              cat("Saved to: univar_base_comp.rds\n")
              
              # Print summary
              cat("\nTop improvements over base model:\n")
              print(univar_base_comp %>%
                      filter(!is.na(draic)) %>%
                      arrange(draic) %>%
                      select(model2_name, draic, draic_pval) %>%
                      head(5))
            }
            
          } else {
            cat("No converged univariate models for comparison\n")
          }
        } else {
          cat("Univariate model summary not available - run compile first\n")
        }
      } else {
        cat("Base model for base_model structure not found\n")
      }
    } else {
      cat("Required models not available for comparison\n")
    }
  }
  
  
  gc()
}

## Constant transition effects -------------------------------------------------------

if ("univar_const" %in% workflow_config$sections) {
  
  if (workflow_config$univar_const$fit) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("FITTING: UNIVARIATE MODELS - CONSTANT TRANSITIONS (PROGRESSIVE TIMEOUT)\n")
    cat(rep("=", 70), "\n", sep = "")
    
    univar_const_results <- univariate_progressive_timeout(
      patient_data = subset(pt_stg, model == "base_model"),
      crude_rates = list(base_model = crude_rates[["base_model"]]),
      covariates = key_covariates,
      spline_vars = continuous_covariates,
      n.cores = n.cores,
      timeout_vector = univar_timeout_times,
      constraint = "constant_across_transitions",
      save_prefix = "univar_const_progressive"
    )
    
    univar_const_models <- univar_const_results$models
    
    # Optional unlimited-time retry
    if (workflow_config$univar_const$retry_failed && 
        nrow(univar_const_results$failed_models) > 0) {
      
      cat("\n=== RETRYING FAILED CONSTANT-TRANSITION MODELS (UNLIMITED TIME) ===\n")
      
      retry_results <- retry_failed_models(
        failed_results = univar_const_results,
        patient_data = pt_stg,
        crude_rates = crude_rates,
        n.cores = n.cores,
        save_prefix = "univar_const_retry"
      )
      
      combined_results <- combine_univariate_results(univar_const_results, retry_results)
      saveRDS(combined_results, here("data", "univar_const_combined.rds"))
      
      univar_const_models <- combined_results$models
    }
    
    saveRDS(univar_const_models, here("data", "univar_const_models.rds"))
    
  } else {
    cat("\n=== LOADING: UNIVARIATE CONSTANT-TRANSITION MODELS ===\n")
    univar_const_models <- readRDS(here("data", "univar_const_models.rds"))
  }
  
  # Compile results
  if (workflow_config$univar_const$compile && !is.null(univar_const_models)) {
    univar_const_comp <- compile_section(
      "univar_const", 
      univar_const_models, 
      analysis_configs$univar_const
    )
  }
  
  # Cross-validation
  if (workflow_config$univar_const$cv && !is.null(univar_const_models)) {
    univar_const_cv <- run_cv_section(
      "univar_const", 
      univar_const_models, 
      cv_configs$univar_const
    )
  }
  
  gc()
}

## Deep-dive on selected covariates ---------------------------------------------

if ("univar_deep_dive" %in% workflow_config$sections) {
  
  deep_dive_covariates <- c("age", "age_cat")
  
  univar_models <- readRDS(here("data", "univar_models.rds"))
  deep_dive_models <- list(
    "base_model_age_spline" = univar_models[["base_model_age_spline"]],
    "base_model_age_cat" = univar_models[["base_model_age_cat_linear"]],
    "base_model_age_linear" = univar_models[["base_model_age_linear"]]
  )

  # Compile results
  if (workflow_config$univar_deep_dive$compile && !is.null(deep_dive_models)) {
    deep_dive_comp <- compile_section(
      "univar_deep_dive",
      deep_dive_models,
      analysis_configs$univar_deep_dive
    )
  }
  
  # Cross-validation
  if (workflow_config$univar_deep_dive$cv && !is.null(deep_dive_models)) {
    deep_dive_cv <- run_cv_section(
      "univar_deep_dive",
      deep_dive_models,
      cv_configs$univar_deep_dive
    )
  }
  
  gc()
  
}
  
# Multivariable model selection -------------------------------------------------

if ("multivar" %in% workflow_config$sections) {
  cat("RUNNING: MULTIVARIABLE MODEL SELECTION\n")

  forward_results <- forward_selection_msm(
    patient_data = pt_stg,
    crude_rates = crude_rates,
    univar_model_summary = univar_comp$model_summary,  # Now uses 'formula' column
    candidate_vars = key_covariates,
    spline_vars = continuous_covariates,
    model_structures = "base_model",
    mc.cores = n.cores,
    timeout_minutes = 30,
    save_prefix = "forward_selection"
  )

  gc()
}

# Time-homogeneity assumption ------------------------------------------------------

if ("time_inhomog" %in% workflow_config$sections) {
  
  time_inhomog_models <- run_section(
    section_name = "time_inhomog",
    section_config = workflow_config$time_inhomog,
    fit_fn = function() {
      
      run_all_time_models(
        patient_data = subset(pt_stg, model == "base_model"),
        crude_rates = list(base_model = crude_rates[["base_model"]]),
        constraint = "transition_specific",
        mc.cores = n.cores
      )
    }
  )
  
  # Compile results
  if (workflow_config$time_inhomog$compile && !is.null(time_inhomog_models)) {
    time_inhomog_comp <- compile_section(
      "time_inhomog", 
      time_inhomog_models, 
      analysis_configs$time_inhomog
    )
  }
  
  # Cross-validation
  if (workflow_config$time_inhomog$cv && !is.null(time_inhomog_models)) {
    time_inhomog_cv <- run_cv_section(
      "time_inhomog", 
      time_inhomog_models, 
      cv_configs$time_inhomog
    )
  }
  
  gc()
}

# Long-stay sensitivity analysis ---------------------------------------------------

if ("long_stay" %in% workflow_config$sections) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("RUNNING: LONG-STAY SENSITIVITY ANALYSIS\n")
  cat(rep("=", 70), "\n", sep = "")
  
  source(here("R", "06-long-stay.R"))
  
  gc()
}

# Summary -----------------------------------------------------------------------

cat("Workflow summary\n")

cat("\nSections Run:\n")
for (section in workflow_config$sections) {
  cat("  -", section, "\n")
  
  section_config <- workflow_config[[section]]
  if (!is.null(section_config)) {
    steps <- c()
    if (isTRUE(section_config$fit)) steps <- c(steps, "fit")
    if (isTRUE(section_config$compile)) steps <- c(steps, "compile")
    if (isTRUE(section_config$compare)) steps <- c(steps, "compare")
    if (isTRUE(section_config$cv)) steps <- c(steps, "cv")
    if (isTRUE(section_config$retry_failed)) steps <- c(steps, "retry")
    
    if (length(steps) > 0) {
      cat("    Steps:", paste(steps, collapse = ", "), "\n")
    }
  }
}

cat("\nComparisons Run:\n")
comparison_run <- FALSE
for (comp_name in names(workflow_config$comparisons)) {
  if (isTRUE(workflow_config$comparisons[[comp_name]]$run)) {
    cat("  -", comp_name, "\n")
    comparison_run <- TRUE
  }
}
if (!comparison_run) {
  cat("  (none)\n")
}

cat("\nFiles Saved:\n")
temp_files <- list.files(here("data"), pattern = "\\.rds$", full.names = FALSE)
if (length(temp_files) > 0) {
  for (f in temp_files) {
    size <- file.size(here("data", f))
    cat("  -", f, "-", format(size, units = "auto"), "\n")
  }
} else {
  cat("  (no .rds files found)\n")
}

cat("\nCV Results:\n")
cv_dirs <- list.dirs(here("data", "cv_results"), 
                     recursive = FALSE, full.names = FALSE)
if (length(cv_dirs) > 0) {
  for (d in cv_dirs) {
    cv_files <- list.files(here("data", "cv_results", d), 
                           pattern = "\\.rds$")
    cat("  -", d, ":", length(cv_files), "files\n")
  }
} else {
  cat("  (no CV results)\n")
}

