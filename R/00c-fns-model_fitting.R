# Model fitting --------------------------------------------------------------------

## Core fitting function -----------------------------------------------------------

fit_msm_models <- function(patient_data, 
                           crude_rates, 
                           model_specs,
                           
                           # Default specifications (used when columns missing from model_specs)
                           covariates = NULL,
                           spline_vars = NULL,
                           spline_df = 3,
                           spline_type = "ns",
                           time_varying = NULL,
                           time_variable = "DaysSinceEntry",
                           time_breakpoints = NULL,
                           time_spline_df = 3,
                           constraint = "transition_specific",
                           
                           require_se = FALSE,
                           
                           # Processing options
                           mc.cores = min(8, parallel::detectCores() - 1)) {
  
  # Load required packages
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("parallel package is required")
  }
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("splines package is required for spline functionality")  
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required")
  }
  
  # Validate model_specs
  if (!is.data.frame(model_specs)) {
    stop("model_specs must be a data frame or tibble")
  }
  
  required_cols <- c("model_name", "model_structure")
  missing_cols <- setdiff(required_cols, names(model_specs))
  if (length(missing_cols) > 0) {
    stop(paste("model_specs missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check for unique model_name
  if (any(duplicated(model_specs$model_name))) {
    stop("model_name must be unique in model_specs")
  }
  
  # Validate that all model_structure values exist in crude_rates and patient_data
  available_models <- unique(patient_data$model)
  requested_models <- unique(model_specs$model_structure)
  crude_models <- names(crude_rates)
  
  missing_in_data <- setdiff(requested_models, available_models)
  missing_in_crude <- setdiff(requested_models, crude_models)
  
  if (length(missing_in_data) > 0) {
    stop(paste("model_structure values not found in patient_data$model:", 
               paste(missing_in_data, collapse = ", ")))
  }
  
  if (length(missing_in_crude) > 0) {
    stop(paste("model_structure values not found in crude_rates:", 
               paste(missing_in_crude, collapse = ", ")))
  }
  
  # Add defaults for missing columns
  default_cols <- list(
    covariates = covariates,
    spline_vars = spline_vars,
    spline_df = spline_df,
    spline_type = spline_type,
    time_varying = time_varying,
    time_variable = time_variable,
    time_breakpoints = time_breakpoints,
    time_spline_df = time_spline_df,
    constraint = constraint
  )
  
  for (col_name in names(default_cols)) {
    if (!col_name %in% names(model_specs)) {
      model_specs[[col_name]] <- default_cols[[col_name]]
    }
  }
  
  cat("Fitting", nrow(model_specs), "model specifications\n")
  
  # Define function to fit a single row from model_specs
  fit_single_spec <- function(spec_row_idx) {
    # Ensure required packages in worker
    if (!require("msm", quietly = TRUE)) return(list(error = "msm package not available"))
    if (!require("splines", quietly = TRUE)) return(list(error = "splines package not available"))
    
    spec_row <- model_specs[spec_row_idx, ]
    model_name <- spec_row$model_name
    model_structure <- spec_row$model_structure
    
    start_time <- Sys.time()
    
    # Extract spec parameters (handling list columns)
    spec_covariates <- if (is.list(spec_row$covariates)) spec_row$covariates[[1]] else spec_row$covariates
    spec_spline_vars <- if (is.list(spec_row$spline_vars)) spec_row$spline_vars[[1]] else spec_row$spline_vars
    spec_spline_df <- spec_row$spline_df
    spec_spline_type <- spec_row$spline_type
    spec_time_varying <- spec_row$time_varying
    spec_time_variable <- spec_row$time_variable
    spec_time_breakpoints <- if (is.list(spec_row$time_breakpoints)) spec_row$time_breakpoints[[1]] else spec_row$time_breakpoints
    spec_time_spline_df <- spec_row$time_spline_df
    spec_constraint <- spec_row$constraint
    
    # Handle NA values (convert to NULL for our purposes)
    if (length(spec_time_varying) == 1 && is.na(spec_time_varying)) spec_time_varying <- NULL
    
    # Get model data
    model_data <- patient_data[patient_data$model == model_structure, ]
    
    if (nrow(model_data) == 0) {
      return(list(
        model_name = model_name,
        result = NULL,
        error = paste("No data found for model_structure:", model_structure),
        runtime = 0
      ))
    }
    
    # Get crude rates
    crude_result <- crude_rates[[model_structure]]
    if (is.null(crude_result)) {
      return(list(
        model_name = model_name,
        result = NULL,
        error = paste("Crude rates missing for model_structure:", model_structure),
        runtime = 0
      ))
    }
    
    # Preprocess data for time-varying models
    processed_data <- preprocess_time_varying_data(
      model_data, 
      spec_time_varying, 
      spec_time_variable, 
      spec_time_breakpoints
    )[[model_structure]]
    
    # Build model specifications
    spec_result <- build_model_specifications(
      spec_covariates, spec_spline_vars, spec_spline_df, spec_spline_type,
      spec_time_varying, spec_time_variable, spec_time_spline_df,
      setNames(list(processed_data), model_structure)
    )
    
    model_spec <- spec_result$specs[["main"]]
    processed_data <- spec_result$processed_data_list[[model_structure]]
    
    # Create formula
    formula_name <- create_formula_name(model_spec$final_covariates)
    covariate_formula <- if (length(model_spec$final_covariates) > 0) {
      as.formula(formula_name)
    } else {
      NULL
    }
    
    # Handle constraints
    qmat <- crude_result$qmat
    n_transitions <- sum(qmat != 0 & row(qmat) != col(qmat))
    
    constraint_for_msm <- create_constraint_specification(
      model_spec$final_covariates, spec_constraint, n_transitions, processed_data
    )
    
    # Try multiple optimization methods
    fitted_model <- try_multiple_optimizations(
      processed_data, crude_result, covariate_formula, constraint_for_msm, require_se = require_se
    )
    
    # Create comprehensive metadata
    metadata <- create_model_metadata(
      model_spec, spec_covariates, spec_spline_vars, spec_spline_df, spec_spline_type,
      spec_time_varying, spec_time_variable, spec_time_breakpoints, spec_time_spline_df,
      spec_constraint, fitted_model
    )
    
    # Add model_structure to metadata
    metadata$model_structure <- model_structure
    
    # Store result
    result <- if (!is.null(fitted_model$model)) {
      list(
        fitted_model = fitted_model$model,
        status = "converged",
        optimization_method = fitted_model$method,
        metadata = metadata,
        crude_rates = crude_result,
        error_message = NULL,
        runtime = as.numeric(difftime(Sys.time(), start_time, units = "auto"))
      )
    } else {
      list(
        fitted_model = NULL,
        status = "failed",
        optimization_method = fitted_model$method,
        metadata = metadata,
        crude_rates = crude_result,
        error_message = fitted_model$error,
        runtime = as.numeric(difftime(Sys.time(), start_time, units = "auto"))
      )
    }
    
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(list(
      model_name = model_name,
      formula_name = formula_name,
      result = result,
      error = NULL,
      runtime = runtime
    ))
  }
  
  # Parallelize across rows of model_specs
  fitted_results <- parallel::mclapply(
    seq_len(nrow(model_specs)),
    fit_single_spec,
    mc.cores = mc.cores
  )
  
  # Reorganize into nested structure: fitted_models[[model_name]][[formula_name]]
  fitted_msm_models <- list()
  runtimes <- list()
  
  for (result in fitted_results) {
    if (is.null(result$error)) {
      model_name <- result$model_name
      formula_name <- result$formula_name
      
      if (!model_name %in% names(fitted_msm_models)) {
        fitted_msm_models[[model_name]] <- list()
      }
      
      fitted_msm_models[[model_name]][[formula_name]] <- result$result
      runtimes[[model_name]] <- result$runtime
    } else {
      warning(paste("Error fitting", result$model_name, ":", result$error))
    }
  }
  
  # Print runtimes
  cat("\n=== Model Fitting Runtimes ===\n")
  for (model_name in names(runtimes)) {
    runtime_mins <- round(runtimes[[model_name]] / 60, 2)
    runtime_secs <- round(runtimes[[model_name]], 1)
    
    if (runtime_mins >= 1) {
      cat(sprintf("%s: %.2f minutes\n", model_name, runtime_mins))
    } else {
      cat(sprintf("%s: %.1f seconds\n", model_name, runtime_secs))
    }
  }
  
  return(fitted_msm_models)
}

## Univariate models with progressive timeouts -----------------------------

# Helper function to build model_specs tibble for univariate analysis
build_univariate_model_specs <- function(covariate_specs, model_structures, 
                                         spline_df = 3, spline_type = "ns",
                                         time_varying = NULL, time_variable = "DaysSinceEntry",
                                         time_breakpoints = NULL, time_spline_df = 3,
                                         constraint = "transition_specific") {
  
  # Build combinations of covariate specs and model structures
  combinations <- expand.grid(
    model_structure = model_structures,
    spec_name = names(covariate_specs),
    stringsAsFactors = FALSE
  )
  
  # Build model_specs tibble
  model_specs <- tibble(
    model_name = paste(combinations$model_structure, combinations$spec_name, sep = "_"),
    model_structure = combinations$model_structure,
    covariates = list(NULL),
    spline_vars = list(NULL),
    spline_df = spline_df,
    spline_type = spline_type,
    time_varying = if(is.null(time_varying)) NA_character_ else time_varying,
    time_variable = time_variable,
    time_breakpoints = list(time_breakpoints),
    time_spline_df = time_spline_df,
    constraint = constraint
  )
  
  # Fill in covariates and spline_vars based on spec type
  # Pre-allocate lists of the correct length
  covariates_list <- vector("list", nrow(model_specs))
  spline_vars_list <- vector("list", nrow(model_specs))
  
  for (i in seq_len(nrow(model_specs))) {
    spec_name <- combinations$spec_name[i]
    spec <- covariate_specs[[spec_name]]
    
    if (spec$type == "linear") {
      covariates_list[[i]] <- spec$variable
      spline_vars_list[[i]] <- NA
    } else if (spec$type == "spline") {
      covariates_list[[i]] <- NA
      spline_vars_list[[i]] <- spec$variable
    }
  }
  
  # Assign the complete lists
  model_specs$covariates <- covariates_list
  model_specs$spline_vars <- spline_vars_list
  
  for (i in seq_len(nrow(model_specs))) {
    # Handle covariates
    if (length(model_specs$covariates[[i]]) == 1 && is.na(model_specs$covariates[[i]])) {
      model_specs$covariates[[i]] <- as.character(NULL)
    }
    
    # Handle spline_vars
    if (length(model_specs$spline_vars[[i]]) == 1 && is.na(model_specs$spline_vars[[i]])) {
      model_specs$spline_vars[[i]] <- as.character(NULL)
    }
  }
  
  return(model_specs)
}

univariate_progressive_timeout <- function(patient_data, 
                                           crude_rates, 
                                           covariates, 
                                           n.cores,
                                           # Model structure specification
                                           model_structures = NULL,  # NEW: specify which to use
                                           # New parameters with sensible defaults
                                           spline_vars = NULL,
                                           spline_df = 3,
                                           spline_type = "ns",
                                           time_varying = NULL,
                                           time_variable = "DaysSinceEntry",
                                           time_breakpoints = NULL,
                                           time_spline_df = 3,
                                           constraint = "transition_specific",
                                           require_se = FALSE,
                                           # Timeout parameters
                                           timeout_vector = c(5, 15, 30, 60), 
                                           save_prefix = "msm_univar_progress") {
  
  # Determine which model structures to use
  if (is.null(model_structures)) {
    model_structures <- names(crude_rates)
    cat("Using all available model structures:", paste(model_structures, collapse = ", "), "\n")
  } else {
    # Validate specified model structures
    missing_structures <- setdiff(model_structures, names(crude_rates))
    if (length(missing_structures) > 0) {
      stop("Specified model_structures not found in crude_rates: ", 
           paste(missing_structures, collapse = ", "))
    }
    cat("Using specified model structures:", paste(model_structures, collapse = ", "), "\n")
  }
  
  # Build list of all covariate specifications to test
  covariate_specs <- build_univariate_specs(covariates, spline_vars)
  
  # Build model_specs tibble for new fit_msm_models interface
  all_model_specs <- build_univariate_model_specs(
    covariate_specs, model_structures, spline_df, spline_type,
    time_varying, time_variable, time_breakpoints, time_spline_df, constraint
  )
  
  # Initialize tracking
  remaining_model_specs <- all_model_specs
  
  # Initialize results structure
  final_results <- list()
  
  # Track failures with enhanced information
  failed_models <- data.frame(
    model_name = character(),
    model_structure = character(),
    spec_name = character(),
    covariate = character(),
    model_type = character(),
    timeout_minutes = numeric(),
    status = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
  
  cat("Starting progressive timeout analysis\n")
  cat("Timeout sequence:", paste(timeout_vector, "mins", collapse = " -> "), "\n")
  cat("Total model combinations:", nrow(all_model_specs), "\n")
  cat("  Linear models:", sum(sapply(covariate_specs, function(x) x$type == "linear")) * length(model_structures), "\n")
  cat("  Spline models:", sum(sapply(covariate_specs, function(x) x$type == "spline")) * length(model_structures), "\n\n")
  
  # Process each timeout level
  for (timeout_idx in seq_along(timeout_vector)) {
    current_timeout <- timeout_vector[timeout_idx]
    
    if (nrow(remaining_model_specs) == 0) {
      cat("All combinations completed! Stopping early.\n")
      break
    }
    
    cat("=== TIMEOUT ROUND", timeout_idx, "===\n")
    cat("Current timeout:", current_timeout, "minutes\n")
    cat("Model specs to attempt:", nrow(remaining_model_specs), "\n\n")
    
    # Track round start time
    round_start_time <- Sys.time()
    
    # Track what times out in this round
    this_round_timeouts <- tibble()
    
    # Track successes for this round
    this_round_successes <- 0
    this_round_errors <- 0
    
    # Process remaining model specs in batches or individually
    for (i in seq_len(nrow(remaining_model_specs))) {
      current_spec <- remaining_model_specs[i, ]
      model_name <- current_spec$model_name
      model_structure <- current_spec$model_structure
      
      # Extract covariate info for display
      covs <- current_spec$covariates[[1]]
      splines <- current_spec$spline_vars[[1]]
      
      if (!is.null(covs) && length(covs) > 0) {
        display_var <- paste(covs, collapse = "+")
        model_type <- "linear"
      } else if (!is.null(splines) && length(splines) > 0) {
        display_var <- paste(splines, collapse = "+")
        model_type <- "spline"
      } else {
        display_var <- "intercept_only"
        model_type <- "intercept"
      }
      
      cat(sprintf("  [%d/%d] %s: %s (%s, timeout: %d mins)...", 
                  i, nrow(remaining_model_specs), model_structure,
                  display_var, model_type, current_timeout))
      
      # Record start time
      start_time <- Sys.time()
      
      # Try to run with current timeout using single-row model_specs
      model_result <- tryCatch({
        withTimeout({
          fit_msm_models(
            patient_data = patient_data,
            crude_rates = crude_rates,
            model_specs = current_spec,
            mc.cores = n.cores,
            require_se = require_se
          )
        }, timeout = current_timeout * 60)
        
      }, TimeoutException = function(e) {
        # Add to timeout list for next round (if there is one)
        if (timeout_idx < length(timeout_vector)) {
          this_round_timeouts <<- bind_rows(this_round_timeouts, current_spec)
          cat(" TIMEOUT (will retry)\n")
        } else {
          # Final timeout - add to failed models
          failed_models <<- rbind(failed_models, data.frame(
            model_name = model_name,
            model_structure = model_structure,
            spec_name = paste(model_structure, display_var, model_type, sep = "_"),
            covariate = display_var,
            model_type = model_type,
            timeout_minutes = current_timeout,
            status = "final_timeout",
            error_message = paste("Timed out after", current_timeout, "minutes (final attempt)"),
            stringsAsFactors = FALSE
          ))
          cat(" FINAL TIMEOUT\n")
        }
        return(NULL)
        
      }, error = function(e) {
        # Log error - don't retry errors
        this_round_errors <<- this_round_errors + 1
        failed_models <<- rbind(failed_models, data.frame(
          model_name = model_name,
          model_structure = model_structure,
          spec_name = paste(model_structure, display_var, model_type, sep = "_"),
          covariate = display_var,
          model_type = model_type,
          timeout_minutes = current_timeout,
          status = "error",
          error_message = as.character(e$message),
          stringsAsFactors = FALSE
        ))
        cat(" ERROR\n")
        return(NULL)
      })
      
      # Calculate runtime
      end_time <- Sys.time()
      runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
      
      # Store successful result
      if (!is.null(model_result)) {
        # Check if any models actually converged in the result
        any_converged <- FALSE
        
        # Process each model_name in the result (should only be one)
        for (result_model_name in names(model_result)) {
          # Initialize model_name in final_results if it doesn't exist
          if (is.null(final_results[[result_model_name]])) {
            final_results[[result_model_name]] <- list()
          }
          
          # Process each formula_name for this model
          for (formula_name in names(model_result[[result_model_name]])) {
            # Get the model data from new structure
            model_data <- model_result[[result_model_name]][[formula_name]]
            
            # Check if this specific model converged
            if (!is.null(model_data$fitted_model) && model_data$status == "converged") {
              any_converged <- TRUE
              
              # Add run metadata
              model_data$run_metadata <- list(
                timeout_round = timeout_idx,
                actual_runtime_minutes = runtime_minutes,
                timeout_limit_minutes = current_timeout,
                model_name = model_name,
                model_structure = model_structure,
                covariate = display_var,
                model_type = model_type
              )
              
              # Store in the structure
              final_results[[result_model_name]][[formula_name]] <- model_data
            } else {
              # Model failed to converge - treat as error
              this_round_errors <- this_round_errors + 1
              failed_models <<- rbind(failed_models, data.frame(
                model_name = model_name,
                model_structure = model_structure,
                spec_name = paste(model_structure, display_var, model_type, sep = "_"),
                covariate = display_var,
                model_type = model_type,
                timeout_minutes = current_timeout,
                status = "convergence_failed",
                error_message = model_data$error_message %||% "Model did not converge",
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        if (any_converged) {
          this_round_successes <- this_round_successes + 1
          cat(" COMPLETED\n")
        } else {
          cat(" FAILED TO CONVERGE\n")
        }
      }
      
      # Periodic cleanup
      if (i %% 5 == 0) {
        gc()
      }
    }
    
    # Calculate round duration
    round_end_time <- Sys.time()
    round_duration <- difftime(round_end_time, round_start_time, units = "mins")
    
    # Update remaining combinations for next round
    remaining_model_specs <- this_round_timeouts
    
    # Save progress after each timeout round
    progress_data <- list(
      models = final_results,
      failed_models = failed_models,
      timeout_round = timeout_idx,
      current_timeout = current_timeout,
      remaining_combinations = nrow(remaining_model_specs),
      round_duration_minutes = as.numeric(round_duration),
      covariate_specs = covariate_specs,
      all_model_specs = all_model_specs,
      settings = list(
        covariates = covariates,
        spline_vars = spline_vars,
        spline_df = spline_df,
        spline_type = spline_type,
        time_varying = time_varying,
        time_variable = time_variable,
        time_breakpoints = time_breakpoints,
        time_spline_df = time_spline_df,
        constraint = constraint
      )
    )
    
    saveRDS(progress_data, file = here("data", "univar_progressive", paste0(save_prefix, "_round_", timeout_idx, "_", current_timeout, "min.rds")))
    
    # Enhanced round summary
    cat("\n", paste(rep("=", 50), collapse = ""), "\n")
    cat("--- Round", timeout_idx, "Summary ---\n")
    cat("Round duration:", round(as.numeric(round_duration), 1), "minutes\n")
    cat("Completed successfully:", this_round_successes, "\n")
    cat("Errors (not retrying):", this_round_errors, "\n")
    cat("Timed out (will retry):", nrow(this_round_timeouts), "\n")
    cat("Total failed so far:", nrow(failed_models), "\n")
    cat(paste(rep("=", 50), collapse = ""), "\n\n")
  }
  
  # Count successful models for final summary
  total_successful <- 0
  n_linear <- 0
  n_spline <- 0
  
  for (model_name in names(final_results)) {
    for (formula_name in names(final_results[[model_name]])) {
      model_data <- final_results[[model_name]][[formula_name]]
      if (!is.null(model_data$fitted_model) && model_data$status == "converged") {
        total_successful <- total_successful + 1
        if (model_data$run_metadata$model_type == "linear") {
          n_linear <- n_linear + 1
        } else if (model_data$run_metadata$model_type == "spline") {
          n_spline <- n_spline + 1
        }
      }
    }
  }
  
  # Final summary
  cat("=== FINAL SUMMARY ===\n")
  cat("Total successful models:", total_successful, "\n")
  cat("  Linear models:", n_linear, "\n")
  cat("  Spline models:", n_spline, "\n")
  cat("Total failed models:", nrow(failed_models), "\n")
  
  if (nrow(failed_models) > 0) {
    cat("\nFailure breakdown:\n")
    print(table(failed_models$status, failed_models$model_type))
  }
  
  # Save final results
  final_data <- list(
    models = final_results,
    failed_models = failed_models,
    timeout_sequence = timeout_vector,
    covariate_specs = covariate_specs,
    all_model_specs = all_model_specs,
    settings = list(
      covariates = covariates,
      spline_vars = spline_vars,
      spline_df = spline_df,
      spline_type = spline_type,
      time_varying = time_varying,
      time_variable = time_variable,
      time_breakpoints = time_breakpoints,
      time_spline_df = time_spline_df,
      constraint = constraint
    )
  )
  
  saveRDS(final_data, file = here("data", "univar_progressive", paste0(save_prefix, "_full_results.rds")))
  cat("Final results saved\n")
  
  return(final_data)
}

#' Retry failed models without timeout constraints
#' @param failed_results Output from univariate_progressive_timeout containing failed_models
#' @param patient_data Patient data
#' @param crude_rates Crude rates list
#' @param covariate_specs Covariate specifications (from failed_results if available)
#' @param n.cores Number of cores
#' @param save_prefix Prefix for saving results
#' @return List with retry results
retry_failed_models <- function(failed_results,
                                patient_data,
                                crude_rates,
                                covariate_specs = NULL,
                                n.cores = min(8, parallel::detectCores() - 1),
                                save_prefix = "univar_retry") {
  
  # Extract failed models and settings
  failed_models <- failed_results$failed_models
  settings <- failed_results$settings
  
  # Get covariate specs if not provided
  if (is.null(covariate_specs)) {
    if (!is.null(failed_results$covariate_specs)) {
      covariate_specs <- failed_results$covariate_specs
    } else {
      # Rebuild specs from settings
      covariate_specs <- build_univariate_specs(
        settings$covariates,
        settings$spline_vars
      )
    }
  }
  
  # Filter to only final_timeout and error cases (exclude temporary timeouts)
  failed_to_retry <- failed_models %>%
    filter(status %in% c("final_timeout", "error"))
  
  if (nrow(failed_to_retry) == 0) {
    cat("No failed models to retry!\n")
    return(list(
      models = list(),
      failed_models = data.frame(),
      settings = settings,
      message = "No models needed retry"
    ))
  }
  
  cat("=== RETRYING FAILED MODELS (NO TIMEOUT) ===\n")
  cat("Total models to retry:", nrow(failed_to_retry), "\n")
  cat("  Final timeouts:", sum(failed_to_retry$status == "final_timeout"), "\n")
  cat("  Errors:", sum(failed_to_retry$status == "error"), "\n\n")
  
  # Build model_specs for retry
  retry_model_specs <- build_univariate_model_specs(
    covariate_specs, names(crude_rates), 
    settings$spline_df %||% 3, settings$spline_type %||% "ns",
    settings$time_varying, settings$time_variable %||% "DaysSinceEntry",
    settings$time_breakpoints, settings$time_spline_df %||% 3,
    settings$constraint %||% "transition_specific"
  )
  
  # Filter to only the failed models
  retry_specs <- retry_model_specs %>%
    filter(model_name %in% failed_to_retry$model_name)
  
  # Initialize results
  retry_results <- list()
  retry_failed <- data.frame(
    model_name = character(),
    model_structure = character(),
    covariate = character(),
    model_type = character(),
    original_status = character(),
    retry_status = character(),
    error_message = character(),
    runtime_minutes = numeric(),
    stringsAsFactors = FALSE
  )
  
  retry_start_time <- Sys.time()
  successes <- 0
  failures <- 0
  
  # Process each failed model
  for (i in seq_len(nrow(retry_specs))) {
    current_spec <- retry_specs[i, ]
    model_name <- current_spec$model_name
    model_structure <- current_spec$model_structure
    
    # Find original failure info
    original_failure <- failed_to_retry %>%
      filter(model_name == !!model_name) %>%
      slice(1)
    
    original_status <- original_failure$status
    covariate_display <- original_failure$covariate
    model_type <- original_failure$model_type
    
    cat(sprintf("  [%d/%d] %s: %s (%s, originally: %s)...", 
                i, nrow(retry_specs), model_structure,
                covariate_display, model_type, original_status))
    
    # Record start time
    start_time <- Sys.time()
    
    # Try to fit WITHOUT timeout
    model_result <- tryCatch({
      fit_msm_models(
        patient_data = patient_data,
        crude_rates = crude_rates,
        model_specs = current_spec,
        mc.cores = n.cores
      )
    }, error = function(e) {
      cat(" ERROR:", e$message, "\n")
      return(NULL)
    })
    
    # Calculate runtime
    end_time <- Sys.time()
    runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
    
    # Process result with enhanced error handling
    if (!is.null(model_result)) {
      # Check if model actually converged
      model_converged <- FALSE
      convergence_error <- NULL
      
      for (result_model_name in names(model_result)) {
        if (is.null(retry_results[[result_model_name]])) {
          retry_results[[result_model_name]] <- list()
        }
        
        for (formula_name in names(model_result[[result_model_name]])) {
          model_data <- model_result[[result_model_name]][[formula_name]]
          
          if (!is.null(model_data$fitted_model) && model_data$status == "converged") {
            model_converged <- TRUE
            
            # Add retry metadata
            model_data$run_metadata <- list(
              retry_attempt = TRUE,
              original_status = original_status,
              actual_runtime_minutes = runtime_minutes,
              model_name = model_name,
              model_structure = model_structure,
              covariate = covariate_display,
              model_type = model_type
            )
            
            retry_results[[result_model_name]][[formula_name]] <- model_data
          } else {
            # Capture convergence failure details
            convergence_error <- model_data$error_message %||% "Model did not converge"
          }
        }
      }
      
      if (model_converged) {
        successes <- successes + 1
        cat(" SUCCESS (", round(runtime_minutes, 1), "min )\n")
        
        retry_failed <- rbind(retry_failed, data.frame(
          model_name = model_name,
          model_structure = model_structure,
          covariate = covariate_display,
          model_type = model_type,
          original_status = original_status,
          retry_status = "success",
          error_message = NA_character_,
          runtime_minutes = runtime_minutes,
          stringsAsFactors = FALSE
        ))
      } else {
        failures <- failures + 1
        cat(" FAILED TO CONVERGE (", round(runtime_minutes, 1), "min )\n")
        
        retry_failed <- rbind(retry_failed, data.frame(
          model_name = model_name,
          model_structure = model_structure,
          covariate = covariate_display,
          model_type = model_type,
          original_status = original_status,
          retry_status = "failed_convergence",
          error_message = convergence_error %||% "Model did not converge",
          runtime_minutes = runtime_minutes,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      failures <- failures + 1
      
      retry_failed <- rbind(retry_failed, data.frame(
        model_name = model_name,
        model_structure = model_structure,
        covariate = covariate_display,
        model_type = model_type,
        original_status = original_status,
        retry_status = "error",
        error_message = "fit_msm_models returned NULL",
        runtime_minutes = runtime_minutes,
        stringsAsFactors = FALSE
      ))
    }
    
    # Periodic cleanup
    if (i %% 5 == 0) {
      gc()
    }
  }
  
  # Calculate total duration
  retry_end_time <- Sys.time()
  total_duration <- difftime(retry_end_time, retry_start_time, units = "mins")
  
  # Final summary
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("=== RETRY SUMMARY ===\n")
  cat("Total retry duration:", round(as.numeric(total_duration), 1), "minutes\n")
  cat("Models attempted:", nrow(retry_specs), "\n")
  cat("Successful:", successes, "\n")
  cat("Still failed:", failures, "\n")
  cat("Success rate:", round(100 * successes / nrow(retry_specs), 1), "%\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  # Show breakdown by original status
  if (nrow(retry_failed) > 0) {
    cat("Breakdown by original failure reason:\n")
    retry_summary <- retry_failed %>%
      group_by(original_status, retry_status) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(original_status, retry_status)
    print(retry_summary)
    cat("\n")
  }
  
  # Save results
  final_data <- list(
    models = retry_results,
    retry_failed = retry_failed,
    retry_summary = if(exists("retry_summary")) retry_summary else NULL,
    total_duration_minutes = as.numeric(total_duration),
    settings = settings,
    timestamp = Sys.time()
  )
  
  saveRDS(final_data, file = here("data", "temp", paste0(save_prefix, "_results.rds")))
  cat("Retry results saved to:", paste0(save_prefix, "_results.rds"), "\n")
  
  return(final_data)
}

### Multivariable model --------------------------------------------------

# Forward selection for MSM models with parallel processing
forward_selection_msm <- function(patient_data,
                                  crude_rates,
                                  univar_model_summary,  # univar_comp$model_summary dataframe
                                  candidate_vars,  # Pool of variables to consider
                                  spline_vars = NULL,
                                  spline_df = 3,
                                  spline_type = "ns",
                                  model_structures = NULL,
                                  constraint = "transition_specific",
                                  mc.cores = min(8, parallel::detectCores() - 1),
                                  timeout_minutes = 30,
                                  save_prefix = "forward_selection",
                                  verbose = TRUE) {
  
  library(dplyr)
  library(parallel)
  library(msm)
  
  # Validate input
  if (!is.data.frame(univar_model_summary)) {
    stop("univar_model_summary must be a data frame")
  }
  
  required_cols <- c("model_structure", "formula", "AIC", "status")
  missing_cols <- setdiff(required_cols, names(univar_model_summary))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in univar_model_summary: ", paste(missing_cols, collapse = ", "))
  }
  
  # Determine which model structures to use
  if (is.null(model_structures)) {
    model_structures <- unique(univar_model_summary$model_structure)
    if (verbose) cat("Using all available model structures:", paste(model_structures, collapse = ", "), "\n")
  }
  
  # Initialize results storage
  selection_history <- list()
  current_models <- list()
  all_fitted_models <- list()
  
  # Extract best univariate models as starting point
  if (verbose) cat("\n=== INITIALIZING FROM BEST UNIVARIATE MODELS ===\n")
  
  for (structure in model_structures) {
    # Filter to converged models for this structure
    structure_summary <- univar_model_summary %>%
      filter(model_structure == !!structure, 
             status == "converged",
             !is.na(AIC)) %>%
      arrange(AIC)
    
    if (nrow(structure_summary) > 0) {
      # Get best model info
      best_model_info <- structure_summary[1, ]
      best_model_formula <- best_model_info$formula
      best_aic <- best_model_info$AIC
      
      # Extract covariates from formula
      current_vars <- extract_covariates_from_formula(best_model_formula)
      
      current_models[[structure]] <- list(
        variables = current_vars,
        aic = best_aic,
        model_name = best_model_formula,
        formula = best_model_formula
      )
      
      if (verbose) {
        cat("Structure", structure, "- Starting with:", best_model_formula, 
            "| Variables:", paste(current_vars, collapse = ", "),
            "| AIC:", round(best_aic, 2), "\n")
      }
    } else {
      if (verbose) {
        cat("Structure", structure, "- No converged models found\n")
      }
    }
  }
  
  if (length(current_models) == 0) {
    stop("No valid converged univariate models found to start forward selection")
  }
  
  # Store initial step
  step_0 <- data.frame(
    step = 0,
    model_structure = names(current_models),
    variables_added = sapply(current_models, function(x) paste(x$variables, collapse = "+")),
    aic = sapply(current_models, function(x) x$aic),
    model_name = sapply(current_models, function(x) x$model_name),
    stringsAsFactors = FALSE
  )
  selection_history[["step_0"]] <- step_0
  
  # Main forward selection loop
  step <- 1
  continue_selection <- TRUE
  
  while (continue_selection) {
    if (verbose) cat("\n=== FORWARD SELECTION STEP", step, "===\n")
    
    step_improvements <- list()
    step_fitted_models <- list()
    
    # For each current model structure
    for (structure in names(current_models)) {
      current_model_info <- current_models[[structure]]
      current_vars <- current_model_info$variables
      current_aic <- current_model_info$aic
      
      # Find remaining candidate variables
      remaining_vars <- setdiff(candidate_vars, current_vars)
      
      if (length(remaining_vars) == 0) {
        if (verbose) cat("Structure", structure, "- No remaining variables to test\n")
        next
      }
      
      if (verbose) {
        cat("Structure", structure, "- Testing", length(remaining_vars), "additional variables\n")
        cat("  Current variables:", paste(current_vars, collapse = ", "), "\n")
        cat("  Current AIC:", round(current_aic, 2), "\n")
      }
      
      # Build model specs for this step
      step_specs <- build_forward_step_specs(
        current_vars = current_vars,
        candidate_vars = remaining_vars,
        spline_vars = spline_vars,
        spline_df = spline_df,
        spline_type = spline_type,
        model_structure = structure,
        constraint = constraint
      )
      
      # Fit models in parallel for this structure
      if (verbose) cat("  Fitting", nrow(step_specs), "candidate models...\n")
      
      step_results <- fit_forward_step_parallel(
        model_specs = step_specs,
        patient_data = patient_data,
        crude_rates = crude_rates,
        mc.cores = mc.cores,
        timeout_minutes = timeout_minutes,
        verbose = verbose
      )
      
      # Store fitted models
      structure_key <- paste0("step_", step, "_", structure)
      step_fitted_models[[structure_key]] <- step_results$models
      
      # Find best improvement for this structure
      if (length(step_results$models) > 0) {
        aics <- sapply(step_results$models, function(x) {
          if (inherits(x, "msm") && !is.null(x$minus2loglik)) {
            return(x$minus2loglik + 2 * length(x$estimates))
          } else {
            return(Inf)
          }
        })
        
        best_idx <- which.min(aics)
        best_aic <- aics[best_idx]
        best_model_name <- names(step_results$models)[best_idx]
        best_model <- step_results$models[[best_model_name]]
        
        # Check if this is an improvement
        aic_improvement <- current_aic - best_aic
        
        if (best_aic < current_aic) {
          # Get the spec for this model to extract the variables correctly
          model_spec_idx <- which(step_specs$model_name == best_model_name)
          if (length(model_spec_idx) > 0) {
            best_spec <- step_specs[model_spec_idx[1], ]
            new_vars <- best_spec$covariates[[1]]  # Get variables from the spec, NOT from model name
            
            added_var <- setdiff(new_vars, current_vars)
            
            # Handle case where multiple variables differ (shouldn't happen in forward selection)
            if (length(added_var) == 1) {
              added_variable <- added_var[1]
            } else if (length(added_var) > 1) {
              added_variable <- paste(added_var, collapse = "+")
              warning("Multiple variables added in one step: ", added_variable)
            } else {
              added_variable <- "unknown"
              warning("Could not determine added variable")
            }
            
            step_improvements[[structure]] <- list(
              variables = new_vars,
              added_variable = added_variable,
              aic = best_aic,
              aic_improvement = aic_improvement,
              model_name = best_model_name,
              formula = paste("~", paste(new_vars, collapse = " + "))
            )
            
            if (verbose) {
              cat("  Best addition:", added_variable, "| New AIC:", round(best_aic, 2), 
                  "| Improvement:", round(aic_improvement, 2), "\n")
            }
          } else {
            if (verbose) cat("  Could not find model spec for best model\n")
          }
        } else {
          if (verbose) cat("  No improvement found (best AIC:", round(best_aic, 2), ")\n")
        }
      } else {
        if (verbose) cat("  No models converged\n")
      }
    }
    
    # Store all fitted models from this step
    all_fitted_models <- c(all_fitted_models, step_fitted_models)
    
    # Update current models with improvements
    any_improvements <- FALSE
    step_summary <- data.frame(
      step = integer(),
      model_structure = character(),
      variables_added = character(),
      added_variable = character(),
      aic = numeric(),
      aic_improvement = numeric(),
      model_name = character(),
      stringsAsFactors = FALSE
    )
    
    for (structure in names(current_models)) {
      if (structure %in% names(step_improvements)) {
        # Update with improved model
        improvement <- step_improvements[[structure]]
        current_models[[structure]] <- improvement
        any_improvements <- TRUE
        
        step_summary <- rbind(step_summary, data.frame(
          step = step,
          model_structure = structure,
          variables_added = paste(improvement$variables, collapse = "+"),
          added_variable = improvement$added_variable,
          aic = improvement$aic,
          aic_improvement = improvement$aic_improvement,
          model_name = improvement$model_name,
          stringsAsFactors = FALSE
        ))
      } else {
        # No improvement, keep current model
        current_info <- current_models[[structure]]
        step_summary <- rbind(step_summary, data.frame(
          step = step,
          model_structure = structure,
          variables_added = paste(current_info$variables, collapse = "+"),
          added_variable = NA_character_,
          aic = current_info$aic,
          aic_improvement = 0,
          model_name = current_info$model_name,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    selection_history[[paste0("step_", step)]] <- step_summary
    
    # Check stopping criteria
    if (!any_improvements) {
      if (verbose) cat("\nNo improvements found. Stopping forward selection.\n")
      continue_selection <- FALSE
    } else {
      if (verbose) {
        cat("Step", step, "summary:\n")
        improvement_summary <- step_summary %>% 
          filter(!is.na(added_variable)) %>%
          select(model_structure, added_variable, aic, aic_improvement)
        if (nrow(improvement_summary) > 0) {
          print(improvement_summary)
        } else {
          cat("  No improvements found\n")
        }
        
        # Debug: Print current state
        for (structure in names(current_models)) {
          model_info <- current_models[[structure]]
          cat("  Current state for", structure, ":\n")
          cat("    Variables:", paste(model_info$variables, collapse = ", "), "\n")
          cat("    AIC:", round(model_info$aic, 2), "\n")
        }
      }
      step <- step + 1
    }
    
    # Safety check to prevent infinite loops
    if (step > length(candidate_vars)) {
      if (verbose) cat("\nReached maximum possible steps. Stopping.\n")
      continue_selection <- FALSE
    }
  }
  
  # Combine selection history
  final_history <- do.call(rbind, selection_history)
  rownames(final_history) <- NULL
  
  # Extract final models - note: these are just the specifications, not fitted models
  final_models_specs <- list()
  for (structure in names(current_models)) {
    final_models_specs[[structure]] <- current_models[[structure]]
  }
  
  # Save results
  if (!is.null(save_prefix)) {
    saveRDS(list(
      final_models_specs = final_models_specs,
      selection_history = final_history,
      all_fitted_models = all_fitted_models,
      metadata = list(
        candidate_vars = candidate_vars,
        spline_vars = spline_vars,
        model_structures = model_structures,
        constraint = constraint,
        mc.cores = mc.cores,
        timeout_minutes = timeout_minutes,
        completed_steps = step - 1,
        run_time = Sys.time()
      )
    ), file = paste0(save_prefix, "_results.rds"))
    
    if (verbose) cat("\nSaved results to:", paste0(save_prefix, "_results.rds"), "\n")
  }
  
  if (verbose) {
    cat("\n=== FORWARD SELECTION COMPLETE ===\n")
    cat("Total steps completed:", step - 1, "\n")
    cat("Final models:\n")
    for (structure in names(current_models)) {
      model_info <- current_models[[structure]]
      cat("  ", structure, ":", paste(model_info$variables, collapse = " + "), 
          "| AIC:", round(model_info$aic, 2), "\n")
    }
  }
  
  return(list(
    final_models_specs = final_models_specs,
    selection_history = final_history,
    all_fitted_models = all_fitted_models,
    metadata = list(
      candidate_vars = candidate_vars,
      spline_vars = spline_vars, 
      model_structures = model_structures,
      constraint = constraint,
      mc.cores = mc.cores,
      timeout_minutes = timeout_minutes,
      completed_steps = step - 1
    )
  ))
}

# Helper function to extract covariates from formula
extract_covariates_from_formula <- function(formula_string) {
  # Parse formula to extract variables
  # Expecting format like "~ age", "~ age_ns1 + age_ns2 + age_ns3", "~ age + gender"
  
  # Handle special cases first
  if (formula_string == "~ 1" || formula_string == "" || is.na(formula_string)) {
    return(character(0))
  }
  
  # Parse the formula to get all terms
  all_terms <- tryCatch({
    all.vars(as.formula(formula_string))
  }, error = function(e) {
    # Fallback: manually parse the string
    terms_part <- sub("^~\\s*", "", formula_string)
    terms_split <- strsplit(terms_part, "\\s*\\+\\s*")[[1]]
    return(trimws(terms_split))
  })
  
  if (length(all_terms) == 0) {
    return(character(0))
  }
  
  # Extract base variable names by removing spline suffixes
  base_vars <- unique(sapply(all_terms, function(term) {
    # Remove spline term suffixes like "_ns1", "_ns2", "_bs1", etc.
    base_var <- gsub("_(ns|bs)\\d+$", "", term)
    # Remove other common suffixes
    base_var <- gsub("_spline\\d*$", "", base_var)
    base_var <- gsub("_quad$", "", base_var)
    return(base_var)
  }))
  
  return(base_vars)
}

# Helper function to build model specs for forward step
build_forward_step_specs <- function(current_vars, candidate_vars, spline_vars, 
                                     spline_df, spline_type, model_structure,
                                     constraint) {
  
  specs <- data.frame(
    model_name = character(),
    model_structure = character(),
    covariates = I(list()),
    spline_vars = I(list()),
    spline_df = numeric(),
    spline_type = character(),
    constraint = character(),
    stringsAsFactors = FALSE
  )
  
  for (new_var in candidate_vars) {
    # Create new variable set
    new_vars <- c(current_vars, new_var)
    
    # Determine if new variable should be splined
    if (new_var %in% spline_vars) {
      new_spline_vars <- intersect(new_vars, spline_vars)
      model_name <- paste(c(model_structure, new_vars, "spline"), collapse = "_")
    } else {
      new_spline_vars <- intersect(current_vars, spline_vars)  # Keep existing splines
      model_name <- paste(c(model_structure, new_vars), collapse = "_")
    }
    
    specs <- rbind(specs, data.frame(
      model_name = model_name,
      model_structure = model_structure,
      covariates = I(list(new_vars)),
      spline_vars = I(list(new_spline_vars)),
      spline_df = spline_df,
      spline_type = spline_type,
      constraint = constraint,
      stringsAsFactors = FALSE
    ))
  }
  
  return(specs)
}

# Helper function to fit models in parallel for one step
fit_forward_step_parallel <- function(model_specs, patient_data, crude_rates,
                                      mc.cores, timeout_minutes, verbose = FALSE) {
  
  # Fit function for single model spec (similar to your existing pattern)
  fit_single_forward_spec <- function(spec_idx) {
    spec <- model_specs[spec_idx, ]
    
    tryCatch({
      # Set timeout
      setTimeLimit(cpu = timeout_minutes * 60, elapsed = timeout_minutes * 60, 
                   transient = TRUE)
      
      start_time <- Sys.time()
      
      # Get model data and crude rates
      model_data <- subset(patient_data, model == spec$model_structure)
      crude_rate <- crude_rates[[spec$model_structure]]
      
      # Extract variables for this spec
      new_vars <- spec$covariates[[1]]
      new_spline_vars <- spec$spline_vars[[1]]
      
      # Build formula using your existing approach
      spec_result <- build_model_specifications(
        new_vars, new_spline_vars, spec$spline_df, spec$spline_type,
        NULL, NULL, NULL,  # Remove time_varying parameters
        setNames(list(model_data), spec$model_structure)
      )
      
      model_spec <- spec_result$specs[["main"]]
      processed_data <- spec_result$processed_data_list[[spec$model_structure]]
      
      # Create formula
      formula_name <- create_formula_name(model_spec$final_covariates)
      covariate_formula <- if (length(model_spec$final_covariates) > 0) {
        as.formula(formula_name)
      } else {
        NULL
      }
      
      # Handle constraints
      qmat <- crude_rate$qmat
      n_transitions <- sum(qmat != 0 & row(qmat) != col(qmat))
      
      constraint_for_msm <- create_constraint_specification(
        model_spec$final_covariates, spec$constraint, n_transitions, processed_data
      )
      
      # Fit MSM model using your existing approach
      fitted_model_result <- try_multiple_optimizations(
        processed_data, crude_rate, covariate_formula, constraint_for_msm
      )
      
      fitted_model <- fitted_model_result$model
      
      runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      
      return(list(
        model_name = spec$model_name,
        result = fitted_model,
        runtime = runtime,
        error = NULL
      ))
      
    }, error = function(e) {
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      return(list(
        model_name = spec$model_name,
        result = NULL,
        runtime = timeout_minutes * 60,
        error = as.character(e)
      ))
    })
  }
  
  # Run in parallel
  fitted_results <- mclapply(
    seq_len(nrow(model_specs)),
    fit_single_forward_spec,
    mc.cores = mc.cores
  )
  
  # Organize results
  fitted_models <- list()
  for (result in fitted_results) {
    if (is.null(result$error) && inherits(result$result, "msm")) {
      fitted_models[[result$model_name]] <- result$result
    } else if (verbose && !is.null(result$error)) {
      cat("    Failed:", result$model_name, "-", result$error, "\n")
    }
  }
  
  return(list(
    models = fitted_models,
    all_results = fitted_results
  ))
}

## Time-varying models -------------------------------------------------------------

# run_all_time_models <- function(patient_data,
#                                 crude_rates,
#                                 covariates = NULL,
#                                 spline_vars = NULL,
#                                 spline_df = 3,
#                                 spline_type = "ns",
#                                 time_breakpoints = NULL,
#                                 time_spline_df = 3,
#                                 constraint = "transition_specific",
#                                 mc.cores = min(8, parallel::detectCores() - 1)) {
#   
#   # Define the time-varying types and time variables to loop over
#   time_varying_options <- c(
#     "linear", 
#     # "piecewise",
#     "spline"
#     )
#   time_variables <- c(
#     "DaysSinceEntry", 
#     "CalendarTime"
#     )
#   
#   # Initialize a list to store results
#   results <- list()
#   
#   # Loop through each time variable
#   for (time_var in time_variables) {
#     
#     # Create a sublist for each time variable
#     results[[time_var]] <- list()
#     
#     # Loop through each time-varying specification
#     for (tv_type in time_varying_options) {
#       
#       # Optional: message for progress
#       message(paste("Fitting model with time variable =", time_var, 
#                     "and time-varying type =", tv_type))
#       
#       # Fit the model
#       model_result <- fit_msm_models(
#         patient_data = patient_data,
#         crude_rates = crude_rates,
#         
#         # Covariate specifications
#         covariates = covariates,
#         spline_vars = spline_vars,
#         spline_df = spline_df,
#         spline_type = spline_type,
#         
#         # Time-varying specs
#         time_varying = tv_type,
#         time_variable = time_var,
#         time_breakpoints = time_breakpoints,
#         time_spline_df = time_spline_df,
#         
#         # Constraint options
#         constraint = constraint,
#         
#         # Processing
#         mc.cores = mc.cores
#       )
#       
#       # Store in the results list
#       results[[time_var]][[tv_type]] <- model_result
#     }
#   }
#   
#   return(results)
# }
# 

run_all_time_models <- function(patient_data,
                                crude_rates,
                                covariates = NULL,
                                spline_vars = NULL,
                                spline_df = 3,
                                spline_type = "ns",
                                time_breakpoints = NULL,
                                time_spline_df = 3,
                                constraint = "transition_specific",
                                mc.cores = min(8, parallel::detectCores() - 1)) {
  
  # Get unique model structures from patient_data
  model_structures <- unique(patient_data$model)
  
  # Define time-varying options
  time_varying_options <- c("linear", "spline")
  time_variables <- c("DaysSinceEntry", "CalendarTime")
  
  # Create model_specs data frame with all combinations
  model_specs <- expand.grid(
    model_structure = model_structures,
    time_variable = time_variables,
    time_varying = time_varying_options,
    stringsAsFactors = FALSE
  )
  
  # Create unique model names
  model_specs$model_name <- paste(
    model_specs$model_structure,
    model_specs$time_variable,
    model_specs$time_varying,
    sep = "_"
  )
  
  # Add other specification columns
  model_specs$covariates <- list(covariates)
  model_specs$spline_vars <- list(spline_vars)
  model_specs$spline_df <- spline_df
  model_specs$spline_type <- spline_type
  model_specs$time_breakpoints <- list(time_breakpoints)
  model_specs$time_spline_df <- time_spline_df
  model_specs$constraint <- constraint
  
  # Fit all models
  fitted_models <- fit_msm_models(
    patient_data = patient_data,
    crude_rates = crude_rates,
    model_specs = model_specs,
    mc.cores = mc.cores
  )
  
  return(fitted_models)
}


# Helper functions -----------------------------------------------------------------

## Helper functions for model run script ---------------------------------------

#' Run a model fitting section with error handling and timing
run_section <- function(section_name, section_config, fit_fn, ...) {
  
  result_file <- here("data", paste0(section_name, "_models.rds"))
  
  if (section_config$fit) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("FITTING:", toupper(section_name), "MODELS\n")
    cat(rep("=", 70), "\n", sep = "")
    
    start_time <- Sys.time()
    
    models <- tryCatch({
      fit_fn(...)
    }, error = function(e) {
      cat("ERROR in", section_name, "fitting:", e$message, "\n")
      return(NULL)
    })
    
    end_time <- Sys.time()
    runtime <- difftime(end_time, start_time, units = "mins")
    
    if (!is.null(models)) {
      saveRDS(models, result_file)
      cat("Runtime:", round(runtime, 2), "minutes\n")
      cat("Saved to:", result_file, "\n")
    }
    
    return(models)
    
  } else {
    cat("\n=== LOADING:", toupper(section_name), "MODELS ===\n")
    
    if (file.exists(result_file)) {
      models <- readRDS(result_file)
      cat("Loaded from:", result_file, "\n")
      return(models)
    } else {
      cat("WARNING: File not found:", result_file, "\n")
      return(NULL)
    }
  }
}

#' Compile model results with error handling
compile_section <- function(section_name, models, config) {
  
  if (is.null(models)) {
    cat("WARNING: No models to compile for", section_name, "\n")
    return(NULL)
  }
  
  result_file <- here("data", paste0(section_name, "_comp.rds"))
  
  cat("\n=== COMPILING:", toupper(section_name), "RESULTS ===\n")
  
  start_time <- Sys.time()
  
  compiled <- tryCatch({
    run_comprehensive_msm_analysis(
      fitted_msm_models = models,
      patient_data = pt_stg,
      time_vec = time_vec,
      analysis_config = config,
      parallel = TRUE,
      mc.cores = n.cores,
      verbose = TRUE
    )
  }, error = function(e) {
    cat("ERROR in", section_name, "compilation:", e$message, "\n")
    return(NULL)
  })
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  if (!is.null(compiled)) {
    saveRDS(compiled, result_file)
    cat("Runtime:", round(runtime, 2), "minutes\n")
    cat("Saved to:", result_file, "\n")
  }
  
  return(compiled)
}

#' Run cross-validation for a section
run_cv_section <- function(section_name, models, cv_config) {
  
  if (is.null(models)) {
    cat("WARNING: No models for CV in", section_name, "\n")
    return(NULL)
  }
  
  cat("\n=== CROSS-VALIDATION:", toupper(section_name), "===\n")
  
  start_time <- Sys.time()
  
  cv_result <- tryCatch({
    do.call(cv_msm_models, c(
      list(
        fitted_models = models,
        patient_data = pt_stg
      ),
      cv_config
    ))
  }, error = function(e) {
    cat("ERROR in", section_name, "CV:", e$message, "\n")
    return(NULL)
  })
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  if (!is.null(cv_result)) {
    cat("CV Runtime:", round(runtime, 2), "minutes\n")
    cat("Results saved to:", cv_config$output_dir, "\n")
  }
  
  return(cv_result)
}

## Helper functions for model fitting ----------------------------------------------

calc_crude_init_rates <- function(patient_data, qmat_list) {
  
  crude_results <- list()
  for (modelname in unique(patient_data$model)) {
    
    model_data <- patient_data[which(patient_data$model == modelname),]
    qmat <- qmat_list[[modelname]]
    
    crude_result <- tryCatch({
      crudeinits.msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, qmatrix = qmat)
    }, error = function(e) {
      warning(paste("Error calculating crude rates for", modelname, ":", e$message))
      return(NULL)
    })
    
    # Save the crude rates and the model name
    if (!is.null(crude_result)) {
      crude_results[[modelname]] <- list(
        qmat = crude_result,
        modelname = modelname
      )
    }
  }
  
  return(crude_results)
}


# Helper function to preprocess data for time-varying models
preprocess_time_varying_data <- function(patient_data, time_varying, time_variable, time_breakpoints) {
  
  model_names <- intersect(unique(patient_data$model), names(crude_rates))
  processed_list <- list()
  
  for (modelname in model_names) {
    model_data <- patient_data[patient_data$model == modelname, ]
    
    if (!is.null(time_varying)) {
      if (time_varying == "piecewise") {
        # Create piecewise periods
        if (is.null(time_breakpoints)) {
          if (time_variable == "DaysSinceEntry") {
            time_breakpoints <- c(7, 14)
          } else {
            time_breakpoints <- quantile(model_data[[time_variable]], c(0.33, 0.67), na.rm = TRUE)
          }
        }
        
        model_data$time_period <- cut(model_data[[time_variable]], 
                                      breaks = c(-Inf, time_breakpoints, Inf),
                                      labels = paste0("period_", seq_len(length(time_breakpoints) + 1)),
                                      include.lowest = TRUE)
        model_data$time_period <- as.character(model_data$time_period)
        
      } else if (time_varying == "quadratic") {
        # Create quadratic term
        model_data[[paste0(time_variable, "_quad")]] <- model_data[[time_variable]]^2
      }
      # For linear and spline, no preprocessing needed
    }
    
    processed_list[[modelname]] <- model_data
  }
  
  return(processed_list)
}

# Helper function to build model specifications
build_model_specifications <- function(covariates, spline_vars, spline_df, spline_type,
                                       time_varying, time_variable, time_spline_df,
                                       processed_data_list) {
  
  specs <- list()
  
  # Get sample data to create spline bases
  sample_data <- processed_data_list[[1]]
  
  # Build base specification
  spec <- list(
    linear_covariates = covariates,
    spline_covariates = spline_vars,
    time_varying = time_varying,
    time_variable = time_variable,
    spline_metadata = list(),
    final_covariates = character(0)
  )
  
  # Add linear covariates
  if (!is.null(covariates)) {
    spec$final_covariates <- c(spec$final_covariates, covariates)
  }
  
  # Add spline terms
  if (!is.null(spline_vars)) {
    for (spline_var in spline_vars) {
      if (spline_var %in% names(sample_data)) {
        # Create spline basis to get metadata
        spline_basis <- tryCatch({
          create_spline_basis(sample_data[[spline_var]], spline_df, spline_type)
        }, error = function(e) {
          warning(paste("Failed to create spline basis for", spline_var, ":", e$message))
          return(NULL)
        })
        
        if (is.null(spline_basis)) {
          next  # Skip this spline variable
        }
        
        # Store metadata
        spline_terms <- paste0(spline_var, "_", spline_type, 1:spline_df)
        spec$spline_metadata[[spline_var]] <- list(
          original_variable = spline_var,
          type = spline_type,
          df = spline_df,
          knots = attr(spline_basis, "knots"),
          boundary_knots = attr(spline_basis, "Boundary.knots"),
          var_range = range(sample_data[[spline_var]], na.rm = TRUE),
          spline_terms = spline_terms
        )
        
        # Add spline terms to final covariates
        spec$final_covariates <- c(spec$final_covariates, spline_terms)
        
        # Add spline basis to all processed data
        for (modelname in names(processed_data_list)) {
          model_data <- processed_data_list[[modelname]]
          if (spline_var %in% names(model_data)) {
            basis <- create_spline_basis(model_data[[spline_var]], spline_df, spline_type,
                                         knots = attr(spline_basis, "knots"),
                                         boundary_knots = attr(spline_basis, "Boundary.knots"))
            for (i in 1:spline_df) {
              model_data[[spline_terms[i]]] <- basis[, i]
            }
            processed_data_list[[modelname]] <- model_data
            
          }
        }
      }
    }
  }
  
  # Add time-varying terms
  if (!is.null(time_varying)) {
    if (time_varying == "linear") {
      spec$final_covariates <- c(spec$final_covariates, time_variable)
    } else if (time_varying == "spline") {
      # Create time spline
      time_spline_terms <- paste0(time_variable, "_", spline_type, 1:time_spline_df)
      spec$final_covariates <- c(spec$final_covariates, time_spline_terms)
      
      # Add to spline metadata
      sample_time_data <- sample_data[[time_variable]]
      time_basis <- create_spline_basis(sample_time_data, time_spline_df, spline_type)
      
      spec$spline_metadata[[time_variable]] <- list(
        original_variable = time_variable,
        type = spline_type,
        df = time_spline_df,
        knots = attr(time_basis, "knots"),
        boundary_knots = attr(time_basis, "Boundary.knots"),
        var_range = range(sample_time_data, na.rm = TRUE),
        spline_terms = time_spline_terms
      )
      
      # Add to all processed data
      for (modelname in names(processed_data_list)) {
        model_data <- processed_data_list[[modelname]]
        if (time_variable %in% names(model_data)) {
          basis <- create_spline_basis(model_data[[time_variable]], time_spline_df, spline_type,
                                       knots = attr(time_basis, "knots"),
                                       boundary_knots = attr(time_basis, "Boundary.knots"))
          for (i in 1:time_spline_df) {
            model_data[[time_spline_terms[i]]] <- basis[, i]
          }
          processed_data_list[[modelname]] <- model_data
        }
      }
    } else if (time_varying == "piecewise") {
      spec$final_covariates <- c(spec$final_covariates, "time_period")
    } else if (time_varying == "quadratic") {
      spec$final_covariates <- c(spec$final_covariates, time_variable, paste0(time_variable, "_quad"))
    }
  }
  
  specs[["main"]] <- spec
  
  return(list(
    specs = specs,
    processed_data_list = processed_data_list
  ))
}

# Helper function to create spline basis
create_spline_basis <- function(x, df, type, knots = NULL, boundary_knots = NULL) {
  if (type == "ns") {
    if (!is.null(knots) && !is.null(boundary_knots)) {
      splines::ns(x, df = df, knots = knots, Boundary.knots = boundary_knots)
    } else {
      splines::ns(x, df = df)
    }
  } else if (type == "bs") {
    if (!is.null(knots) && !is.null(boundary_knots)) {
      splines::bs(x, df = df, knots = knots, Boundary.knots = boundary_knots)
    } else {
      splines::bs(x, df = df)
    }
  } else {
    stop(paste("Spline type", type, "not supported"))
  }
}

# Helper function to create formula name
create_formula_name <- function(covariates) {
  if (length(covariates) == 0) {
    return("~ 1")
  } else {
    return(paste("~", paste(covariates, collapse = " + ")))
  }
}

# Helper function to create constraint specification
create_constraint_specification <- function(covariates, constraint, n_transitions, model_data) {
  if (constraint == "constant_across_transitions" && length(covariates) > 0) {
    constraint_list <- list()
    for (cov_name in covariates) {
      if (cov_name %in% names(model_data)) {
        if (is.factor(model_data[[cov_name]]) || is.character(model_data[[cov_name]])) {
          factor_levels <- levels(as.factor(model_data[[cov_name]]))
          if (length(factor_levels) > 1) {
            for (i in 2:length(factor_levels)) {
              level_name <- paste0(cov_name, factor_levels[i])
              constraint_list[[level_name]] <- rep(1, n_transitions)
            }
          }
        } else {
          constraint_list[[cov_name]] <- rep(1, n_transitions)
        }
      }
    }
    return(if (length(constraint_list) == 0) NULL else constraint_list)
  } else {
    return(NULL)  # Transition-specific (unconstrained)
  }
}

# Helper function to try multiple optimization methods
try_multiple_optimizations <- function(model_data, crude_result, covariate_formula, constraint_for_msm, require_se = TRUE) {
  
  # Define optimization methods with different strategies
  optimization_methods <- list(
    # First round - standard settings
    list(opt_method = "optim", name = "BFGS_standard", maxit = 20000, fnscale = 10000),
    list(opt_method = "bobyqa", name = "BOBYQA_standard", maxit = 20000, fnscale = 10000),
    
    # Second round - higher fnscale (more aggressive scaling)
    list(opt_method = "optim", name = "BFGS_highscale", maxit = 50000, fnscale = 50000),
    list(opt_method = "bobyqa", name = "BOBYQA_highscale", maxit = 50000, fnscale = 50000),
    
    # Third round - lower fnscale (more conservative)
    list(opt_method = "optim", name = "BFGS_lowscale", maxit = 50000, fnscale = 1000),
    list(opt_method = "bobyqa", name = "BOBYQA_lowscale", maxit = 50000, fnscale = 1000),
    
    # Fourth round - very high iterations
    list(opt_method = "optim", name = "BFGS_100k", maxit = 100000, fnscale = 10000),
    list(opt_method = "bobyqa", name = "BOBYQA_100k", maxit = 100000, fnscale = 10000),
    
    # Fifth round - extreme settings
    list(opt_method = "optim", name = "BFGS_200k", maxit = 200000, fnscale = 50000),
    list(opt_method = "bobyqa", name = "BOBYQA_200k", maxit = 200000, fnscale = 50000)
  )
  
  for (opt_config in optimization_methods) {
    fitted_model <- tryCatch({
      if (opt_config$opt_method == "nlm") {
        # NLM uses different control structure
        control_list <- list(
          stepmax = 10,
          ndigit = 10
        )
      } else {
        control_list <- list(
          fnscale = opt_config$fnscale,
          maxit = opt_config$maxit,
          reltol = 1e-10
        )
      }
      
      msm_args <- list(
        formula = state_num ~ DaysSinceEntry,
        subject = quote(deid_enc_id),
        data = model_data,
        qmatrix = crude_result$qmat,
        opt.method = opt_config$opt_method,
        control = control_list
      )
      
      if (!is.null(covariate_formula)) {
        msm_args$covariates <- covariate_formula
      }
      
      if (!is.null(constraint_for_msm)) {
        msm_args$constraint <- constraint_for_msm
      }
      
      result <- do.call(msm::msm, msm_args)
      
      # Check convergence
      converged <- is.null(result$opt$convergence) || result$opt$convergence == 0
      
      # Check for standard errors
      has_se <- !is.null(result$foundse) && result$foundse == TRUE
      
      # Check if Hessian is positive definite
      hessian_ok <- TRUE
      if (!is.null(result$opt$hessian)) {
        hessian_check <- tryCatch({
          eigen_vals <- eigen(result$opt$hessian, symmetric = TRUE, only.values = TRUE)$values
          all(eigen_vals > 1e-8)  # Allow small negative eigenvalues due to numerical precision
        }, error = function(e) FALSE)
        hessian_ok <- hessian_check
      }
      
      # More lenient acceptance criteria - accept if has SEs even if hessian check fails
      if (converged && has_se) {
        message("Method ", opt_config$name, " succeeded! (maxit=", opt_config$maxit, 
                ", fnscale=", opt_config$fnscale, ", hessian_ok=", hessian_ok, ")")
        return(list(
          model = result, 
          method = opt_config$name, 
          error = NULL,
          has_se = has_se,
          maxit_used = opt_config$maxit,
          fnscale_used = opt_config$fnscale,
          hessian_positive_definite = hessian_ok
        ))
      } else if (converged && !has_se && !require_se) {
        # Accept without SEs if not required
        message("Method ", opt_config$name, " converged without SEs (maxit=", opt_config$maxit, ")")
        return(list(
          model = result, 
          method = opt_config$name, 
          error = NULL,
          has_se = has_se,
          maxit_used = opt_config$maxit,
          fnscale_used = opt_config$fnscale,
          hessian_positive_definite = hessian_ok
        ))
      } else {
        # Log why it failed
        fail_reason <- if (!converged) {
          "not converged"
        } else if (!has_se && require_se) {
          "no SEs"
        } else {
          "hessian not positive definite"
        }
        message("Method ", opt_config$name, " failed: ", fail_reason, 
                " (maxit=", opt_config$maxit, ", fnscale=", opt_config$fnscale, ")")
        NULL
      }
      
    }, error = function(e) {
      message("Method ", opt_config$name, " threw error: ", e$message)
      NULL
    })
    
    if (!is.null(fitted_model)) {
      return(fitted_model)
    }
  }
  
  # If all methods failed
  return(list(
    model = NULL, 
    method = "all_failed", 
    error = if (require_se) {
      "All optimization methods failed or did not produce standard errors"
    } else {
      "All optimization methods failed"
    },
    has_se = FALSE
  ))
}

# Helper function to create comprehensive metadata
create_model_metadata <- function(spec, covariates, spline_vars, spline_df, spline_type,
                                  time_varying, time_variable, time_breakpoints, time_spline_df,
                                  constraint, fitted_model) {
  
  metadata <- list(
    # Input specifications
    input_covariates = covariates,
    input_spline_vars = spline_vars,
    input_constraint = constraint,
    
    # Covariate details
    linear_covariates = spec$linear_covariates,
    spline_covariates = spec$spline_covariates,
    final_covariates = spec$final_covariates,
    
    # Spline specifications
    spline_config = if (length(spec$spline_metadata) > 0) {
      list(
        spline_vars = names(spec$spline_metadata),
        spline_metadata = spec$spline_metadata,
        default_type = spline_type,
        default_df = spline_df
      )
    } else NULL,
    
    # Time-varying specifications
    time_config = if (!is.null(time_varying)) {
      config <- list(
        type = time_varying,
        variable = time_variable
      )
      if (time_varying == "piecewise") {
        config$breakpoints <- time_breakpoints
      } else if (time_varying == "spline") {
        config$spline_df <- time_spline_df
        config$spline_type <- spline_type
        
        if (!is.null(spec$spline_metadata[[time_variable]])) {
          config$knots <- spec$spline_metadata[[time_variable]]$knots
          config$boundary_knots <- spec$spline_metadata[[time_variable]]$boundary_knots
          config$var_range <- spec$spline_metadata[[time_variable]]$var_range
        } else {
          config$knots <- NA
          config$boundary_knots <- NA
          config$var_range <- NA
        }
        
      }
      config
    } else NULL,
    
    # Constraint information
    constraint_config = list(
      type = constraint,
      applied = !is.null(fitted_model$model) && !is.null(fitted_model$model$constraint)
    ),
    
    # Fitting metadata
    fit_timestamp = Sys.time(),
    n_covariates = length(spec$final_covariates),
    has_splines = length(spec$spline_metadata) > 0,
    has_time_varying = !is.null(time_varying)
  )
  
  return(metadata)
}

combine_model_lists <- function(list1, 
                                list2, 
                                model_names, 
                                source = NULL) {
  
  # Validate inputs
  if (!is.list(list1) || !is.list(list2)) {
    stop("list1 and list2 must be lists")
  }
  
  if (!is.character(model_names)) {
    stop("model_names must be a character vector")
  }
  
  # If source is not provided, create it based on where each model is found
  if (is.null(source)) {
    source <- sapply(model_names, function(name) {
      in_list1 <- name %in% names(list1)
      in_list2 <- name %in% names(list2)
      
      if (in_list1 && in_list2) {
        warning(paste("Model", name, "found in both lists. Using list1. ",
                      "Specify 'source' parameter to control selection."))
        return(1)
      } else if (in_list1) {
        return(1)
      } else if (in_list2) {
        return(2)
      } else {
        stop(paste("Model", name, "not found in either list"))
      }
    })
  }
  
  # Validate source vector
  if (length(source) != length(model_names)) {
    stop("source must be the same length as model_names")
  }
  
  if (!all(source %in% c(1, 2))) {
    stop("source must contain only 1 (for list1) or 2 (for list2)")
  }
  
  # Build combined list
  combined_list <- list()
  
  for (i in seq_along(model_names)) {
    model_name <- model_names[i]
    source_list <- if (source[i] == 1) list1 else list2
    
    if (!model_name %in% names(source_list)) {
      stop(paste("Model", model_name, "not found in list", source[i]))
    }
    
    combined_list[[model_name]] <- source_list[[model_name]]
  }
  
  return(combined_list)
}

## Helper functions for univariate models -----------------------------------------

# Helper function to build univariate specifications
build_univariate_specs <- function(covariates, spline_vars) {
  specs <- list()
  
  # Add linear specifications for all covariates
  for (covar in covariates) {
    spec_name <- paste0(covar, "_linear")
    specs[[spec_name]] <- list(
      variable = covar,
      type = "linear"
    )
  }
  
  # Add spline specifications for spline_vars (even if they overlap with covariates)
  if (!is.null(spline_vars)) {
    for (spline_var in spline_vars) {
      spec_name <- paste0(spline_var, "_spline")
      specs[[spec_name]] <- list(
        variable = spline_var,
        type = "spline"
      )
    }
  }
  
  return(specs)
}

#' Combine results from univariate_progressive_timeout and retry_failed_models
#' @param univar_results Results from univariate_progressive_timeout
#' @param retry_results Results from retry_failed_models (optional)
#' @return Combined results list with unified structure
combine_univariate_results <- function(univar_results, retry_results = NULL) {
  
  cat("=== COMBINING UNIVARIATE RESULTS ===\n")
  
  # Start with univar_results models
  combined_models <- univar_results$models
  
  # Count initial models
  initial_count <- 0
  for (model_name in names(combined_models)) {
    initial_count <- initial_count + length(combined_models[[model_name]])
  }
  cat("Initial models (from progressive timeout):", initial_count, "\n")
  
  # Add retry results if provided
  retry_count <- 0
  if (!is.null(retry_results) && !is.null(retry_results$models)) {
    
    for (model_name in names(retry_results$models)) {
      # Initialize model structure if it doesn't exist
      if (is.null(combined_models[[model_name]])) {
        combined_models[[model_name]] <- list()
      }
      
      # Add each formula from retry results
      for (formula_name in names(retry_results$models[[model_name]])) {
        
        # Check if this formula already exists (it shouldn't, but just in case)
        if (!is.null(combined_models[[model_name]][[formula_name]])) {
          warning(paste("Formula", formula_name, "already exists in model", model_name, 
                        "- keeping retry version"))
        }
        
        # Add the retry result
        combined_models[[model_name]][[formula_name]] <- 
          retry_results$models[[model_name]][[formula_name]]
        retry_count <- retry_count + 1
      }
    }
    
    cat("Added from retry:", retry_count, "\n")
  }
  
  # Total count
  total_count <- 0
  for (model_name in names(combined_models)) {
    total_count <- total_count + length(combined_models[[model_name]])
  }
  cat("Total combined models:", total_count, "\n\n")
  
  # Combine failed models
  combined_failed <- univar_results$failed_models
  
  if (!is.null(retry_results) && !is.null(retry_results$retry_failed)) {
    # Only include models that failed in retry (not successes)
    retry_still_failed <- retry_results$retry_failed %>%
      filter(retry_status != "success")
    
    if (nrow(retry_still_failed) > 0) {
      # Remove the original failures that were retried and still failed
      # Keep only the retry version
      retry_specs <- paste(retry_still_failed$crude_rate_name, 
                           retry_still_failed$spec_name, sep = "___")
      original_specs <- paste(combined_failed$crude_rate_name, 
                              combined_failed$spec_name, sep = "___")
      
      # Remove original versions of retried failures
      combined_failed <- combined_failed %>%
        filter(!paste(crude_rate_name, spec_name, sep = "___") %in% retry_specs)
      
      # Add the retry failures with updated information
      retry_failed_formatted <- retry_still_failed %>%
        mutate(
          timeout_minutes = NA,  # Not applicable for retry
          status = retry_status,
          error_message = ifelse(is.na(error_message), 
                                 paste("Failed on retry after", round(runtime_minutes, 1), "minutes"),
                                 error_message)
        ) %>%
        select(crude_rate_name, spec_name, covariate, model_type, 
               timeout_minutes, status, error_message)
      
      combined_failed <- bind_rows(combined_failed, retry_failed_formatted)
    }
    
    # Remove from failed list any that succeeded in retry
    retry_succeeded <- retry_results$retry_failed %>%
      filter(retry_status == "success")
    
    if (nrow(retry_succeeded) > 0) {
      success_specs <- paste(retry_succeeded$crude_rate_name, 
                             retry_succeeded$spec_name, sep = "___")
      combined_failed <- combined_failed %>%
        filter(!paste(crude_rate_name, spec_name, sep = "___") %in% success_specs)
    }
  }
  
  cat("Failed models summary:\n")
  cat("  Total still failed:", nrow(combined_failed), "\n")
  if (nrow(combined_failed) > 0) {
    cat("  By status:\n")
    status_table <- table(combined_failed$status)
    for (status in names(status_table)) {
      cat("    ", status, ":", status_table[status], "\n")
    }
  }
  cat("\n")
  
  # Combine settings (prefer retry settings if available, otherwise use original)
  combined_settings <- univar_results$settings
  if (!is.null(retry_results) && !is.null(retry_results$settings)) {
    # Retry settings should be the same, but just in case
    combined_settings <- retry_results$settings
  }
  
  # Create covariate specs if not present
  if (is.null(univar_results$covariate_specs) && !is.null(retry_results$covariate_specs)) {
    covariate_specs <- retry_results$covariate_specs
  } else {
    covariate_specs <- univar_results$covariate_specs
  }
  
  # Build comprehensive summary
  summary_stats <- create_combined_summary(combined_models, combined_failed, 
                                           univar_results, retry_results)
  
  # Create combined result
  combined_results <- list(
    models = combined_models,
    failed_models = combined_failed,
    covariate_specs = covariate_specs,
    settings = combined_settings,
    summary = summary_stats,
    source_info = list(
      univar_timeout_sequence = univar_results$timeout_sequence,
      retry_performed = !is.null(retry_results),
      combined_timestamp = Sys.time()
    )
  )
  
  # Print summary
  cat("=== COMBINED RESULTS SUMMARY ===\n")
  print_combined_summary(summary_stats)
  
  return(combined_results)
}

#' Create summary statistics for combined results
#' @param combined_models Combined models list
#' @param combined_failed Combined failed models data frame
#' @param univar_results Original univar results
#' @param retry_results Retry results (can be NULL)
#' @return Summary statistics list
create_combined_summary <- function(combined_models, combined_failed, 
                                    univar_results, retry_results) {
  
  # Count converged models by type
  model_counts <- data.frame(
    model_structure = character(),
    formula = character(),
    covariate = character(),
    model_type = character(),
    status = character(),
    source = character(),
    stringsAsFactors = FALSE
  )
  
  for (model_name in names(combined_models)) {
    for (formula_name in names(combined_models[[model_name]])) {
      model_entry <- combined_models[[model_name]][[formula_name]]
      
      # Determine source
      source <- if (!is.null(model_entry$run_metadata$retry_attempt) && 
                    model_entry$run_metadata$retry_attempt) {
        "retry"
      } else {
        "progressive_timeout"
      }
      
      model_counts <- rbind(model_counts, data.frame(
        model_structure = model_name,
        formula = formula_name,
        covariate = model_entry$run_metadata$covariate %||% NA_character_,
        model_type = model_entry$run_metadata$model_type %||% NA_character_,
        status = model_entry$status,
        source = source,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Calculate statistics
  total_models <- nrow(model_counts)
  converged_models <- sum(model_counts$status == "converged")
  from_progressive <- sum(model_counts$source == "progressive_timeout")
  from_retry <- sum(model_counts$source == "retry")
  
  linear_models <- sum(model_counts$model_type == "linear" & model_counts$status == "converged")
  spline_models <- sum(model_counts$model_type == "spline" & model_counts$status == "converged")
  
  # By model structure
  by_structure <- model_counts %>%
    filter(status == "converged") %>%
    group_by(model_structure) %>%
    summarise(
      n_converged = n(),
      n_linear = sum(model_type == "linear"),
      n_spline = sum(model_type == "spline"),
      .groups = "drop"
    )
  
  list(
    total_models = total_models,
    converged_models = converged_models,
    failed_models = nrow(combined_failed),
    from_progressive_timeout = from_progressive,
    from_retry = from_retry,
    linear_models = linear_models,
    spline_models = spline_models,
    by_structure = by_structure,
    model_counts = model_counts
  )
}

#' Print combined summary
#' @param summary_stats Summary statistics from create_combined_summary
print_combined_summary <- function(summary_stats) {
  cat("Total models:", summary_stats$total_models, "\n")
  cat("  Converged:", summary_stats$converged_models, "\n")
  cat("  Failed:", summary_stats$failed_models, "\n")
  cat("\nBy source:\n")
  cat("  From progressive timeout:", summary_stats$from_progressive_timeout, "\n")
  cat("  From retry:", summary_stats$from_retry, "\n")
  cat("\nBy type:\n")
  cat("  Linear models:", summary_stats$linear_models, "\n")
  cat("  Spline models:", summary_stats$spline_models, "\n")
  cat("\nBy model structure:\n")
  print(summary_stats$by_structure)
}


