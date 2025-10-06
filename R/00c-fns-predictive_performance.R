
# Cross-validation to assess predictive performance of MSM models -----

## Main CV function -----

#' Main CV function - orchestrates cross-validation across models
#' 
#' @param fitted_models Nested list: fitted_models[[model_structure]][[formula]]
#' @param patient_data Patient-level data with all necessary columns
#' @param crude_rates List of crude rates by model structure
#' @param k_folds Number of CV folds
#' @param stratify_by Variable to stratify folds (default: "final_state")
#' @param prediction_horizon Time horizon for death/severe probabilities (days)
#' @param output_dir Directory to save results
#' @param parallel Use parallel processing
#' @param n_cores Number of cores for parallel processing
#' @param verbose Print progress messages
#' 
#' @return Path to saved summary results file
cv_msm_models <- function(fitted_models,
                          patient_data,
                          crude_rates,
                          k_folds = 5,
                          stratify_by = "final_state",
                          prediction_horizon = 365,
                          output_dir = here::here("data", "temp", "cv_results"),
                          parallel = TRUE,
                          n_cores = min(8, parallel::detectCores() - 1),
                          verbose = TRUE) {
  
  ### Input validation ----

  if (verbose) cat("\n=== Validating Inputs ===\n")
  
  # Check fitted_models structure
  if (!is.list(fitted_models) || length(fitted_models) == 0) {
    stop("fitted_models must be a non-empty list")
  }
  
  validate_model_structure(fitted_models)
  
  # Check patient_data
  required_cols <- c("deid_enc_id", "model", "state", "DaysSinceEntry")
  missing_cols <- setdiff(required_cols, names(patient_data))
  if (length(missing_cols) > 0) {
    stop("patient_data missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.data.frame(patient_data) || nrow(patient_data) == 0) {
    stop("patient_data must be a non-empty data.frame")
  }
  
  # Check crude_rates
  if (!is.list(crude_rates) || length(crude_rates) == 0) {
    stop("crude_rates must be a non-empty list")
  }
  
  # Check stratify_by variable exists
  if (!stratify_by %in% names(patient_data)) {
    # Try to create it
    if (verbose) cat("Creating", stratify_by, "variable from data...\n")
    patient_data <- patient_data %>%
      group_by(deid_enc_id) %>%
      mutate(final_state = last(state)) %>%
      ungroup()
    
    if (!stratify_by %in% names(patient_data)) {
      stop("Could not create stratify_by variable: ", stratify_by)
    }
  }
  
  # Check output directory
  if (!dir.exists(output_dir)) {
    if (verbose) cat("Creating output directory:", output_dir, "\n")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Check parameters
  if (k_folds < 2 || k_folds > 20) {
    stop("k_folds must be between 2 and 20")
  }
  
  if (prediction_horizon <= 0 || prediction_horizon > 1000) {
    stop("prediction_horizon must be between 0 and 1000 days")
  }
  
  ### Setup ----

  start_time <- Sys.time()
  
  # Extract model specifications
  model_specs <- extract_model_specifications(fitted_models)
  n_models <- nrow(model_specs)
  
  if (n_models == 0) {
    stop("No valid models found in fitted_models")
  }
  
  # Create stratified folds
  fold_assignments <- create_stratified_folds(
    patient_data, 
    k_folds, 
    stratify_by,
    verbose
  )
  
  n_patients <- length(unique(fold_assignments$deid_enc_id))
  
  if (verbose) {
    cat("\n=== Cross-Validation Setup ===\n")
    cat("Models to evaluate:", n_models, "\n")
    cat("Folds:", k_folds, "\n")
    cat("Patients:", n_patients, "\n")
    cat("Prediction horizon:", prediction_horizon, "days\n")
    cat("Stratified by:", stratify_by, "\n")
    cat("Parallel:", ifelse(parallel, paste("Yes (", n_cores, "cores)"), "No"), "\n\n")
  }
  
  # Setup parallel processing if requested
  if (parallel && n_cores > 1) {
    if (!requireNamespace("future", quietly = TRUE) || 
        !requireNamespace("future.apply", quietly = TRUE)) {
      warning("Parallel packages not available. Using sequential processing.")
      parallel <- FALSE
    } else {
      library(future)
      library(future.apply)
      plan(multisession, workers = n_cores)
      on.exit(plan(sequential), add = TRUE)
    }
  }

  ### Run CV for each model -----

  
  cv_function <- if (parallel) future.apply::future_lapply else lapply
  
  all_results <- cv_function(1:n_models, function(i) {
    
    spec <- model_specs[i, ]
    
    if (verbose) {
      cat(sprintf("[%d/%d] %s | %s\n", i, n_models, spec$model, spec$formula))
    }
    
    tryCatch({
      # Run CV for this model
      result <- cv_single_model(
        model_structure = spec$model,
        formula_name = spec$formula,
        fitted_models = fitted_models,
        patient_data = patient_data,
        crude_rates = crude_rates,
        fold_assignments = fold_assignments,
        prediction_horizon = prediction_horizon,
        k_folds = k_folds
      )
      
      # Save to file immediately
      filename <- sprintf("cv_%s_%s.rds", 
                          make.names(spec$model), 
                          make.names(spec$formula))
      filepath <- file.path(output_dir, filename)
      
      saveRDS(result, filepath)
      
      # Return just the summary (not full predictions) to save memory
      list(
        model = spec$model,
        formula = spec$formula,
        summary = result$summary,
        file = filepath,
        status = "success"
      )
      
    }, error = function(e) {
      warning(sprintf("Model %s | %s failed: %s", spec$model, spec$formula, e$message))
      list(
        model = spec$model,
        formula = spec$formula,
        summary = NULL,
        file = NA,
        status = "failed",
        error = e$message
      )
    })
  })
  

  ### Compile and save results -----

  
  if (verbose) cat("\n=== Compiling Results ===\n")
  
  # Extract summaries
  all_summaries <- map_dfr(all_results, function(res) {
    if (!is.null(res$summary)) {
      res$summary %>%
        mutate(model = res$model, formula = res$formula, .before = 1)
    } else {
      NULL
    }
  })
  
  # Count successes/failures
  n_success <- sum(sapply(all_results, function(x) x$status == "success"))
  n_failed <- n_models - n_success
  
  # Save combined summary
  summary_file <- file.path(output_dir, "cv_summary_all_models.rds")
  saveRDS(all_summaries, summary_file)
  
  # Create and save metadata
  metadata <- list(
    n_models = n_models,
    n_success = n_success,
    n_failed = n_failed,
    n_patients = n_patients,
    k_folds = k_folds,
    prediction_horizon = prediction_horizon,
    stratify_by = stratify_by,
    runtime = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
    timestamp = Sys.time(),
    individual_files = sapply(all_results, function(x) x$file),
    failed_models = all_results[sapply(all_results, function(x) x$status == "failed")]
  )
  
  saveRDS(metadata, file.path(output_dir, "cv_metadata.rds"))
  
  if (verbose) {
    cat("\nCompleted:", n_success, "/", n_models, "models\n")
    cat("Failed:", n_failed, "models\n")
    cat("Total time:", round(metadata$runtime, 1), "minutes\n")
    cat("Results saved to:", output_dir, "\n")
  }
  
  return(invisible(summary_file))
}


# Helper functions -----

state_to_num <- function(state_name, fitted_model) {
  all_states <- rownames(fitted_model$qmodel$qmatrix)
  state_num <- which(all_states == state_name)
  
  if (length(state_num) == 0) {
    stop("State '", state_name, "' not found in model states: ", 
         paste(all_states, collapse = ", "))
  }
  
  return(state_num)
}

states_to_nums <- function(state_names, fitted_model) {
  sapply(state_names, function(s) state_to_num(s, fitted_model))
}

#' Extract model specifications from nested fitted_models structure
extract_model_specifications <- function(fitted_models) {
  
  specs <- list()
  
  for (model_structure in names(fitted_models)) {
    for (formula_name in names(fitted_models[[model_structure]])) {
      
      model_entry <- fitted_models[[model_structure]][[formula_name]]
      
      # Skip failed models
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && 
           model_entry$status == "failed")) {
        next
      }
      
      # Check for fitted model
      has_model <- if (is.list(model_entry) && "fitted_model" %in% names(model_entry)) {
        !is.null(model_entry$fitted_model)
      } else {
        inherits(model_entry, "msm")
      }
      
      if (has_model) {
        specs[[length(specs) + 1]] <- data.frame(
          model = model_structure,
          formula = formula_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(specs) == 0) {
    return(data.frame(model = character(), formula = character()))
  }
  
  do.call(rbind, specs)
}

#' Create stratified folds ensuring balanced outcomes
create_stratified_folds <- function(patient_data, k_folds, stratify_by, verbose = FALSE) {
  
  # Get one row per patient with stratification variable
  unique_patients <- patient_data %>%
    group_by(deid_enc_id) %>%
    summarise(
      strata = first(!!sym(stratify_by)),
      model = first(model),
      .groups = "drop"
    )
  
  # Check for missing stratification values
  if (any(is.na(unique_patients$strata))) {
    n_missing <- sum(is.na(unique_patients$strata))
    warning(sprintf("%d patients missing %s variable - assigning to random folds", 
                    n_missing, stratify_by))
    unique_patients$strata[is.na(unique_patients$strata)] <- "missing"
  }
  
  # Create stratified folds
  set.seed(123)
  fold_assignments <- unique_patients %>%
    group_by(strata) %>%
    mutate(fold = sample(rep(1:k_folds, length.out = n()))) %>%
    ungroup() %>%
    select(deid_enc_id, fold)
  
  # Validate fold sizes
  fold_sizes <- table(fold_assignments$fold)
  min_fold <- min(fold_sizes)
  max_fold <- max(fold_sizes)
  
  if (max_fold / min_fold > 1.5) {
    warning(sprintf("Unbalanced folds: min=%d, max=%d patients", min_fold, max_fold))
  }
  
  if (verbose) {
    cat("Fold sizes:", paste(fold_sizes, collapse = ", "), "\n")
  }
  
  return(fold_assignments)
}

#' Run CV for a single model across all folds
cv_single_model <- function(model_structure,
                            formula_name,
                            fitted_models,
                            patient_data,
                            crude_rates,
                            fold_assignments,
                            prediction_horizon,
                            k_folds) {
  
  # Extract original fitted model and metadata
  model_entry <- fitted_models[[model_structure]][[formula_name]]
  
  if (is.null(model_entry)) {
    stop("Model entry is NULL")
  }
  
  # Extract metadata
  metadata <- if (is.list(model_entry) && "metadata" %in% names(model_entry)) {
    model_entry$metadata
  } else {
    list()
  }
  
  # Filter data for this model structure
  model_data <- patient_data %>% 
    filter(model == model_structure)
  
  if (nrow(model_data) == 0) {
    stop("No data for model structure: ", model_structure)
  }
  
  # Get fold assignments for patients in this model
  model_fold_assignments <- fold_assignments %>%
    filter(deid_enc_id %in% unique(model_data$deid_enc_id))
  
  # Run each fold
  fold_results <- list()
  fold_errors <- list()
  
  for (fold_num in 1:k_folds) {
    
    fold_result <- tryCatch({
      cv_single_fold(
        fold_num = fold_num,
        model_data = model_data,
        fold_assignments = model_fold_assignments,
        crude_rates = crude_rates[[model_structure]],
        metadata = metadata,
        prediction_horizon = prediction_horizon,
        verbose = FALSE
      )
    }, error = function(e) {
      # ADD MORE DETAIL HERE:
      cat(sprintf("    Fold %d FAILED: %s\n", fold_num, e$message))
      
      fold_errors[[fold_num]] <<- list(
        fold = fold_num,
        error = e$message
      )
      NULL
    })
  }
  
  # Check if any folds succeeded
  n_succeeded <- length(fold_results)
  
  if (n_succeeded == 0) {
    stop("All folds failed for this model")
  }
  
  # Combine fold results
  all_predictions <- map_dfr(fold_results, "predictions")
  all_metrics <- map_dfr(fold_results, "metrics")
  
  # Calculate summary statistics across folds
  summary_stats <- all_metrics %>%
    summarise(
      n_folds_converged = sum(model_converged, na.rm = TRUE),
      n_folds_total = k_folds,
      
      # MAE metrics
      total_los_mae = mean(total_los_mae, na.rm = TRUE),
      total_los_mae_se = sd(total_los_mae, na.rm = TRUE) / sqrt(n()),
      
      days_severe_mae = mean(days_severe_mae, na.rm = TRUE),
      days_severe_mae_se = sd(days_severe_mae, na.rm = TRUE) / sqrt(n()),
      
      # AUC metrics
      death_auc = mean(death_auc, na.rm = TRUE),
      death_auc_se = sd(death_auc, na.rm = TRUE) / sqrt(n()),
      
      severe_auc = mean(severe_auc, na.rm = TRUE),
      severe_auc_se = sd(severe_auc, na.rm = TRUE) / sqrt(n()),
      
      # Calibration metrics
      calibration_slope_death = mean(calibration_slope_death, na.rm = TRUE),
      calibration_slope_severe = mean(calibration_slope_severe, na.rm = TRUE),
      
      # Brier scores
      brier_death = mean(brier_death, na.rm = TRUE),
      brier_severe = mean(brier_severe, na.rm = TRUE)
    )
  
  # Return comprehensive results
  list(
    summary = summary_stats,
    fold_metrics = all_metrics,
    predictions = all_predictions,
    errors = if (length(fold_errors) > 0) do.call(rbind, fold_errors) else NULL,
    metadata = list(
      model = model_structure,
      formula = formula_name,
      k_folds = k_folds,
      n_folds_succeeded = n_succeeded,
      prediction_horizon = prediction_horizon
    )
  )
}

#' Process a single CV fold
cv_single_fold <- function(fold_num,
                           model_data,
                           fold_assignments,
                           crude_rates,
                           metadata,
                           prediction_horizon,
                           verbose = FALSE) {
  
  # Split into train/test
  train_patients <- fold_assignments$deid_enc_id[fold_assignments$fold != fold_num]
  test_patients <- fold_assignments$deid_enc_id[fold_assignments$fold == fold_num]
  
  train_data <- model_data %>% filter(deid_enc_id %in% train_patients)
  test_data <- model_data %>% filter(deid_enc_id %in% test_patients)
  
  # Validate split
  if (nrow(train_data) == 0) {
    stop("Fold ", fold_num, ": No training data")
  }
  
  if (nrow(test_data) == 0) {
    stop("Fold ", fold_num, ": No test data")
  }
  
  # Recalculate crude rates on training data
  train_crude_rates <- tryCatch({
    model_name <- unique(train_data$model)[1]
    temp_qmat_list <- list()
    temp_qmat_list[[model_name]] <- crude_rates$qmat
    calc_crude_init_rates(train_data, temp_qmat_list)
  }, error = function(e) {
    stop("Failed to calculate crude rates for fold ", fold_num, ": ", e$message)
  })
  
  # Refit model on training data
  fitted_model_fold <- tryCatch({
    refit_model_from_metadata(train_data, train_crude_rates, metadata)
  }, error = function(e) {
    stop("Failed to refit model for fold ", fold_num, ": ", e$message)
  })
  
  # Extract the fitted model
  model_structure_name <- names(fitted_model_fold)[1]
  formula_key <- names(fitted_model_fold[[model_structure_name]])[1]
  
  fitted_model <- extract_fitted_model(
    fitted_model_fold, 
    model_structure_name, 
    formula_key,
    paste("in fold", fold_num)
  )
  
  if (is.null(fitted_model)) {
    stop("Fold ", fold_num, ": Failed to extract fitted model")
  }
  
  # Check convergence
  model_converged <- is.null(fitted_model$opt$convergence) || 
    fitted_model$opt$convergence == 0
  
  if (!model_converged) {
    warning("Fold ", fold_num, ": Model did not converge")
  }
  
  if (verbose) {
    cat(sprintf("  Fold %d: Model converged = %s\n", fold_num, model_converged))
    cat(sprintf("  Fold %d: Train=%d, Test=%d patients\n", 
                fold_num, length(train_patients), length(test_patients)))
    cat(sprintf("  Fold %d: Model states = %s\n", fold_num, 
                paste(rownames(fitted_model$qmodel$qmatrix), collapse = ", ")))
  }
  
  # Prepare test data outcomes
  test_outcomes <- test_data %>%
    group_by(deid_enc_id) %>%
    summarise(
      initial_state_name = first(state),
      initial_state_num = state_to_num(first(state), fitted_model),  # ADD
      total_los = max(DaysSinceEntry, na.rm = TRUE),
      days_severe = sum(grepl("S", state), na.rm = TRUE),
      died = any(state == "D"),
      ever_severe = any(grepl("S", state)),
      .groups = "drop"
    )
  
  test_outcomes$initial_state_num <- sapply(test_outcomes$initial_state_name, 
                                            function(s) state_to_num(s, fitted_model))
  
  # Get baseline covariates
  baseline_covariates <- extract_baseline_covariates(test_data, metadata)
  
  # Generate predictions for each patient
  predictions <- generate_predictions_for_fold(
    test_outcomes = test_outcomes,
    baseline_covariates = baseline_covariates,
    fitted_model = fitted_model,
    prediction_horizon = prediction_horizon,
    verbose = verbose
  )
  
  # Calculate metrics
  metrics <- calculate_fold_metrics(
    fold_num = fold_num,
    n_train = length(train_patients),
    n_test = length(test_patients),
    model_converged = model_converged,
    predictions = predictions
  )
  
  # Add fold number to predictions
  predictions$fold <- fold_num
  
  list(
    predictions = predictions,
    metrics = metrics
  )
}

#' Refit model from metadata specifications
refit_model_from_metadata <- function(train_data, train_crude_rates, metadata) {
  
  # Extract parameters with safe defaults
  covariates <- metadata$input_covariates
  spline_vars <- metadata$input_spline_vars
  constraint <- metadata$input_constraint %||% "transition_specific"
  
  spline_df <- if (!is.null(metadata$spline_config$default_df)) {
    metadata$spline_config$default_df
  } else {
    3
  }
  
  spline_type <- if (!is.null(metadata$spline_config$default_type)) {
    metadata$spline_config$default_type
  } else {
    "ns"
  }
  
  time_varying <- if (!is.null(metadata$time_config)) {
    metadata$time_config$type
  } else {
    NULL
  }
  
  time_variable <- if (!is.null(metadata$time_config)) {
    metadata$time_config$variable
  } else {
    "DaysSinceEntry"
  }
  
  time_breakpoints <- if (!is.null(metadata$time_config)) {
    metadata$time_config$breakpoints
  } else {
    NULL
  }
  
  # Call fit_msm_models
  fit_msm_models(
    patient_data = train_data,
    crude_rates = train_crude_rates,
    covariates = covariates,
    spline_vars = spline_vars,
    spline_df = spline_df,
    spline_type = spline_type,
    time_varying = time_varying,
    time_variable = time_variable,
    time_breakpoints = time_breakpoints,
    constraint = constraint,
    mc.cores = 1  # Already parallelizing at model level
  )
}

#' Extract baseline covariates for test patients
extract_baseline_covariates <- function(test_data, metadata) {
  
  # Get list of covariate names from metadata
  base_vars <- extract_base_variables_from_metadata(metadata)
  
  if (length(base_vars) == 0) {
    # No covariates model
    return(test_data %>%
             group_by(deid_enc_id) %>%
             slice_min(DaysSinceEntry, with_ties = FALSE) %>%
             ungroup() %>%
             select(deid_enc_id))
  }
  
  # Extract baseline (first observation) covariates
  baseline <- test_data %>%
    group_by(deid_enc_id) %>%
    slice_min(DaysSinceEntry, with_ties = FALSE) %>%
    ungroup()
  
  # Check which covariates are available
  available_vars <- intersect(base_vars, names(baseline))
  missing_vars <- setdiff(base_vars, names(baseline))
  
  if (length(missing_vars) > 0) {
    warning("Missing baseline covariates: ", paste(missing_vars, collapse = ", "))
  }
  
  if (length(available_vars) == 0) {
    stop("No baseline covariates found in data")
  }
  
  baseline %>%
    select(deid_enc_id, all_of(available_vars))
}

#' Generate predictions for all test patients in a fold
generate_predictions_for_fold <- function(test_outcomes,
                                          baseline_covariates,
                                          fitted_model,
                                          prediction_horizon,
                                          verbose = FALSE) {
  
  # Merge outcomes with baseline covariates
  test_data_combined <- test_outcomes %>%
    left_join(baseline_covariates, by = "deid_enc_id")
  
  if (verbose && nrow(test_data_combined) > 0) {
    cat(sprintf("  Generating predictions for %d patients\n", 
                nrow(test_data_combined)))
    cat(sprintf("  Sample patient - Initial state: %s (num=%d)\n",
                test_data_combined$initial_state_name[1],
                test_data_combined$initial_state_num[1]))
  }
  
  # Generate predictions for each patient
  predictions_list <- lapply(1:nrow(test_data_combined), function(i) {
    
    patient <- test_data_combined[i, ]
    
    tryCatch({
      predict_patient_outcomes(
        fitted_model = fitted_model,
        patient_id = patient$deid_enc_id,
        initial_state_name = patient$initial_state_name,  # CHANGE
        initial_state_num = patient$initial_state_num,     # ADD
        covariates = patient[, !names(patient) %in% c("deid_enc_id", 
                                                      "initial_state_name",   # CHANGE
                                                      "initial_state_num",    # ADD
                                                      "total_los", "days_severe", 
                                                      "died", "ever_severe")],
        prediction_horizon = prediction_horizon
      )
    }, error = function(e) {
      # Return NA predictions if patient prediction fails
      data.frame(
        deid_enc_id = patient$deid_enc_id,
        pred_total_los = NA_real_,
        pred_days_severe = NA_real_,
        pred_prob_death = NA_real_,
        pred_prob_severe = NA_real_,
        error = e$message
      )
    })
  })
  
  predictions <- do.call(rbind, predictions_list)
  
  # Merge predictions with outcomes
  test_outcomes %>%
    left_join(predictions, by = "deid_enc_id")
}

#' Predict outcomes for a single patient using analytical MSM functions
predict_patient_outcomes <- function(fitted_model,
                                     patient_id,
                                     initial_state_name,
                                     initial_state_num,
                                     covariates,
                                     prediction_horizon) {
  
  # Validate inputs
  if (!inherits(fitted_model, "msm")) {
    stop("fitted_model must be an msm object")
  }
  
  all_states <- rownames(fitted_model$qmodel$qmatrix)
  
  if (initial_state_num < 1 || initial_state_num > length(all_states)) {
    stop("initial_state_num out of range: ", initial_state_num)
  }
  
  # Detect severe states
  severe_states <- all_states[grepl("S", all_states)]
  severe_state_nums <- which(all_states %in% severe_states)
  
  # Prepare covariates
  has_covariates <- ncol(covariates) > 0 && any(!is.na(covariates))
  cov_list <- if (has_covariates) as.list(covariates) else NULL
  
  # =============================================================================
  # PREDICTION 1: Total Length of Stay
  # =============================================================================
  
  total_los <- tryCatch({
    if (has_covariates) {
      result <- totlos.msm(fitted_model, start = initial_state_num, covariates = cov_list)
    } else {
      result <- totlos.msm(fitted_model, start = initial_state_num)
    }
    
    # totlos.msm returns a named vector - extract the value for starting state
    unname(result[initial_state_num])
    
  }, error = function(e) {
    warning("Patient ", patient_id, " - totlos.msm failed: ", e$message)
    NA_real_
  })
  
  # =============================================================================
  # PREDICTION 2: Days in Severe States
  # =============================================================================
  
  days_severe <- if (length(severe_state_nums) > 0) {
    tryCatch({
      if (has_covariates) {
        visits <- envisits.msm(fitted_model, start = initial_state_num, 
                               covariates = cov_list)
        sojourns <- sojourn.msm(fitted_model, covariates = cov_list)
      } else {
        visits <- envisits.msm(fitted_model, start = initial_state_num)
        sojourns <- sojourn.msm(fitted_model)
      }
      
      # envisits returns a vector, sojourn returns df with $estimates
      # Extract just the values we need
      visit_values <- unname(visits[severe_state_nums])
      sojourn_values <- sojourns$estimates[severe_state_nums]
      
      sum(visit_values * sojourn_values, na.rm = TRUE)
      
    }, error = function(e) {
      warning("Patient ", patient_id, " - days severe failed: ", e$message)
      NA_real_
    })
  } else {
    0  # No severe states in this model
  }
  
  # =============================================================================
  # PREDICTION 3 & 4: Death and Severe Probabilities (use same pmatrix)
  # =============================================================================
  
  pmat_result <- tryCatch({
    if (has_covariates) {
      pmatrix.msm(fitted_model, t = prediction_horizon, covariates = cov_list)
    } else {
      pmatrix.msm(fitted_model, t = prediction_horizon)
    }
  }, error = function(e) {
    warning("Patient ", patient_id, " - pmatrix.msm failed: ", e$message)
    NULL
  })
  
  # Extract death probability from pmatrix
  prob_death <- if (!is.null(pmat_result) && "D" %in% colnames(pmat_result)) {
    pmat_result[initial_state_num, "D"]
  } else {
    NA_real_
  }
  
  # Extract severe probability from pmatrix
  prob_severe <- if (!is.null(pmat_result) && length(severe_states) > 0) {
    sum(pmat_result[initial_state_num, severe_states], na.rm = TRUE)
  } else if (length(severe_states) == 0) {
    0  # No severe states in model
  } else {
    NA_real_
  }
  
  # Return single-row data frame
  data.frame(
    deid_enc_id = patient_id,
    pred_total_los = as.numeric(total_los),
    pred_days_severe = as.numeric(days_severe),
    pred_prob_death = as.numeric(prob_death),
    pred_prob_severe = as.numeric(prob_severe),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}


#' Calculate performance metrics for a fold
calculate_fold_metrics <- function(fold_num, n_train, n_test, model_converged, predictions) {
  
  # Remove rows with missing predictions
  valid_data <- predictions %>%
    filter(!is.na(pred_total_los) | !is.na(pred_prob_death))
  
  n_valid <- nrow(valid_data)
  
  if (n_valid < 10) {
    warning("Fold ", fold_num, ": Only ", n_valid, " valid predictions")
  }
  
  # Calculate MAE metrics
  total_los_mae <- mean(abs(valid_data$total_los - valid_data$pred_total_los), na.rm = TRUE)
  days_severe_mae <- mean(abs(valid_data$days_severe - valid_data$pred_days_severe), na.rm = TRUE)
  
  # Calculate AUC metrics
  death_auc <- calculate_auc_safe(valid_data$died, valid_data$pred_prob_death)
  severe_auc <- calculate_auc_safe(valid_data$ever_severe, valid_data$pred_prob_severe)
  
  # Calculate calibration metrics
  cal_death <- calculate_calibration_logistic(valid_data$died, valid_data$pred_prob_death)
  cal_severe <- calculate_calibration_logistic(valid_data$ever_severe, valid_data$pred_prob_severe)
  
  # Calculate Brier scores
  brier_death <- mean((valid_data$died - valid_data$pred_prob_death)^2, na.rm = TRUE)
  brier_severe <- mean((valid_data$ever_severe - valid_data$pred_prob_severe)^2, na.rm = TRUE)
  
  # Return metrics for this fold
  data.frame(
    fold = fold_num,
    n_train = n_train,
    n_test = n_test,
    n_valid_predictions = n_valid,
    model_converged = model_converged,
    
    # MAE metrics
    total_los_mae = total_los_mae,
    days_severe_mae = days_severe_mae,
    
    # AUC metrics
    death_auc = as.numeric(death_auc),
    severe_auc = as.numeric(severe_auc),
    
    # Calibration slopes
    calibration_slope_death = cal_death$slope,
    calibration_slope_severe = cal_severe$slope,
    
    # Brier scores
    brier_death = brier_death,
    brier_severe = brier_severe
  )
}

#' Safe AUC calculation with error handling
calculate_auc_safe <- function(observed, predicted) {
  
  # Check for valid inputs
  if (length(observed) == 0 || length(predicted) == 0) {
    return(NA_real_)
  }
  
  # Remove NAs
  valid <- !is.na(observed) & !is.na(predicted)
  
  if (sum(valid) < 10) {
    return(NA_real_)
  }
  
  obs <- observed[valid]
  pred <- predicted[valid]
  
  # Check for variance in outcome
  if (var(obs) == 0) {
    return(NA_real_)
  }
  
  # Check for valid probability range
  if (any(pred < 0 | pred > 1)) {
    warning("Predicted probabilities outside [0,1] range")
    pred <- pmax(0, pmin(1, pred))
  }
  
  # Calculate AUC
  tryCatch({
    if (!requireNamespace("pROC", quietly = TRUE)) {
      stop("pROC package required for AUC calculation")
    }
    pROC::auc(obs, pred, quiet = TRUE)
  }, error = function(e) {
    warning("AUC calculation failed: ", e$message)
    NA_real_
  })
}

#' Calculate logistic calibration slope for binary outcomes
calculate_calibration_logistic <- function(observed, predicted) {
  
  # Check for valid inputs
  if (length(observed) < 10 || length(predicted) < 10) {
    return(list(slope = NA_real_, intercept = NA_real_))
  }
  
  # Remove NAs
  valid <- !is.na(observed) & !is.na(predicted) & 
    predicted > 0 & predicted < 1
  
  if (sum(valid) < 10) {
    return(list(slope = NA_real_, intercept = NA_real_))
  }
  
  obs <- observed[valid]
  pred <- predicted[valid]
  
  # Check for variance
  if (var(obs) == 0 || var(pred) == 0) {
    return(list(slope = NA_real_, intercept = NA_real_))
  }
  
  # Calculate linear predictor (logit)
  linear_predictor <- log(pred / (1 - pred))
  
  # Check for infinite values
  if (any(!is.finite(linear_predictor))) {
    warning("Non-finite linear predictor values detected")
    return(list(slope = NA_real_, intercept = NA_real_))
  }
  
  # Fit logistic regression: observed ~ linear_predictor
  fit <- tryCatch({
    glm(obs ~ linear_predictor, family = binomial())
  }, error = function(e) {
    warning("Calibration GLM failed: ", e$message)
    NULL
  })
  
  if (is.null(fit)) {
    return(list(slope = NA_real_, intercept = NA_real_))
  }
  
  # Extract coefficients
  coefs <- coef(fit)
  
  list(
    slope = coefs[2],
    intercept = coefs[1]
  )
}


# UTILITY FUNCTIONS


#' Load CV results from saved files
#' 
#' @param output_dir Directory where CV results are saved
#' @return List with summary and metadata
load_cv_results <- function(output_dir = here::here("data", "temp")) {
  
  # Check if directory exists
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir)
  }
  
  # Load summary
  summary_file <- file.path(output_dir, "cv_summary_all_models.rds")
  if (!file.exists(summary_file)) {
    stop("Summary file not found: ", summary_file)
  }
  
  summary <- readRDS(summary_file)
  
  # Load metadata
  metadata_file <- file.path(output_dir, "cv_metadata.rds")
  if (file.exists(metadata_file)) {
    metadata <- readRDS(metadata_file)
  } else {
    metadata <- NULL
  }
  
  list(
    summary = summary,
    metadata = metadata
  )
}

#' Load individual model CV results
#' 
#' @param model Model structure name
#' @param formula Formula name
#' @param output_dir Directory where CV results are saved
#' @return Full CV results for specified model
load_model_cv_results <- function(model, formula, 
                                  output_dir = here::here("data", "temp")) {
  
  filename <- sprintf("cv_%s_%s.rds", make.names(model), make.names(formula))
  filepath <- file.path(output_dir, filename)
  
  if (!file.exists(filepath)) {
    stop("Results file not found: ", filepath)
  }
  
  readRDS(filepath)
}

#' Print summary of CV results
#' 
#' @param cv_results Output from load_cv_results()
print_cv_summary <- function(cv_results) {
  
  cat("\n=== Cross-Validation Summary ===\n\n")
  
  if (!is.null(cv_results$metadata)) {
    cat("Models evaluated:", cv_results$metadata$n_models, "\n")
    cat("Successful:", cv_results$metadata$n_success, "\n")
    cat("Failed:", cv_results$metadata$n_failed, "\n")
    cat("Total runtime:", round(cv_results$metadata$runtime, 1), "minutes\n")
    cat("Folds:", cv_results$metadata$k_folds, "\n")
    cat("Patients:", cv_results$metadata$n_patients, "\n\n")
  }
  
  if (!is.null(cv_results$summary) && nrow(cv_results$summary) > 0) {
    cat("Top 5 models by Death AUC:\n")
    top_models <- cv_results$summary %>%
      arrange(desc(death_auc)) %>%
      head(5) %>%
      select(model, formula, death_auc, death_auc_se, 
             total_los_mae, n_folds_converged)
    
    print(top_models, row.names = FALSE)
    
    cat("\n")
  }
}

#' Create comparison plot of model performance
#' 
#' @param cv_summary Summary data frame from CV results
#' @param metric Metric to plot (default: "death_auc")
#' @param top_n Number of top models to show (default: 10)
plot_cv_comparison <- function(cv_summary, metric = "death_auc", top_n = 10) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  
  library(ggplot2)
  
  # Get SE column name
  se_col <- paste0(metric, "_se")
  
  # Check if columns exist
  if (!metric %in% names(cv_summary)) {
    stop("Metric not found: ", metric)
  }
  
  # Prepare data
  plot_data <- cv_summary %>%
    arrange(desc(!!sym(metric))) %>%
    head(top_n) %>%
    mutate(
      model_label = paste(model, formula, sep = "\n"),
      model_label = factor(model_label, levels = rev(model_label))
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = model_label, y = !!sym(metric))) +
    geom_point(size = 3) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste("Top", top_n, "Models by", metric),
      x = "Model",
      y = metric
    )
  
  # Add error bars if SE column exists
  if (se_col %in% names(plot_data)) {
    p <- p + geom_errorbar(
      aes(ymin = !!sym(metric) - 1.96 * !!sym(se_col),
          ymax = !!sym(metric) + 1.96 * !!sym(se_col)),
      width = 0.2
    )
  }
  
  p
}

#' Export CV results to CSV for external analysis
#' 
#' @param cv_results Output from load_cv_results()
#' @param output_file Path to save CSV file
export_cv_results <- function(cv_results, 
                              output_file = here::here("results", "cv_summary.csv")) {
  
  if (is.null(cv_results$summary) || nrow(cv_results$summary) == 0) {
    stop("No summary results to export")
  }
  
  # Ensure directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write to CSV
  write.csv(cv_results$summary, output_file, row.names = FALSE)
  
  cat("Results exported to:", output_file, "\n")
  
  invisible(output_file)
}


