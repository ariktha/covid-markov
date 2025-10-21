# Data wrangling -------------------------------------------------------------------

add_calendar_time <- function(patient_data, date_column = "Date") {
  # Check if the specified date column exists
  if (!date_column %in% names(patient_data)) {
    stop(paste("Column", date_column, "not found in the data"))
  }
  
  # Convert the specified date column to actual dates if it's not already
  if (!inherits(patient_data[[date_column]], "Date")) {
    patient_data[[date_column]] <- as.Date(patient_data[[date_column]])
  }
  
  # Find the first date in the dataset
  first_date <- min(patient_data[[date_column]], na.rm = TRUE)
  
  # Add CalendarTime column (sequential days from first date)
  patient_data$CalendarTime <- as.numeric(patient_data[[date_column]] - first_date) + 1
  
  return(patient_data)
}

# Model setup ----------------------------------------------------------------------

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
      processed_data, crude_result, covariate_formula, constraint_for_msm
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
        error_message = NULL
      )
    } else {
      list(
        fitted_model = NULL,
        status = "failed",
        optimization_method = fitted_model$method,
        metadata = metadata,
        crude_rates = crude_result,
        error_message = fitted_model$error
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
            mc.cores = n.cores
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

## Time-varying models -------------------------------------------------------------

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
  
  # Define the time-varying types and time variables to loop over
  time_varying_options <- c("linear", "spline", "piecewise")
  time_variables <- c("DaysSinceEntry", "CalendarTime")
  
  # Initialize a list to store results
  results <- list()
  
  # Loop through each time variable
  for (time_var in time_variables) {
    
    # Create a sublist for each time variable
    results[[time_var]] <- list()
    
    # Loop through each time-varying specification
    for (tv_type in time_varying_options) {
      
      # Optional: message for progress
      message(paste("Fitting model with time variable =", time_var, 
                    "and time-varying type =", tv_type))
      
      # Fit the model
      model_result <- fit_msm_models(
        patient_data = patient_data,
        crude_rates = crude_rates,
        
        # Covariate specifications
        covariates = covariates,
        spline_vars = spline_vars,
        spline_df = spline_df,
        spline_type = spline_type,
        
        # Time-varying specs
        time_varying = tv_type,
        time_variable = time_var,
        time_breakpoints = time_breakpoints,
        time_spline_df = time_spline_df,
        
        # Constraint options
        constraint = constraint,
        
        # Processing
        mc.cores = mc.cores
      )
      
      # Store in the results list
      results[[time_var]][[tv_type]] <- model_result
    }
  }
  
  return(results)
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
try_multiple_optimizations <- function(model_data, crude_result, covariate_formula, constraint_for_msm) {
  
  optimization_methods <- list(
    list(opt_method = "optim", method = "BFGS", name = "BFGS"),
    list(opt_method = "bobyqa", name = "BOBYQA"),
    list(opt_method = "optim", method = "Nelder-Mead", name = "Nelder-Mead")
  )
  
  for (opt_config in optimization_methods) {
    fitted_model <- tryCatch({
      control_list <- list(fnscale = 10000, maxit = 20000, reltol = 1e-10)
      
      if (opt_config$opt_method == "optim" && 
          "method" %in% names(opt_config) && 
          !is.null(opt_config$method)) {
        control_list$method <- opt_config$method
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
      if (converged) {
        return(list(model = result, method = opt_config$name, error = NULL))
      } else {
        NULL  # Try next method
      }
      
    }, error = function(e) {
      NULL  # Try next method
    })
    
    if (!is.null(fitted_model)) {
      return(fitted_model)
    }
  }
  
  # If all methods failed
  return(list(model = NULL, method = "all_failed", error = "All optimization methods failed"))
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


