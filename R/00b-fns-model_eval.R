# Extract model results ----------------------------------------------

tidy_msm_models <- function(fitted_msm_models) {
  
  # Check for required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required but not installed")
  }
  
  # Initialize output with tibble for consistency
  tidied_models <- tibble::tibble()
  
  # Input validation
  if (!is.list(fitted_msm_models)) {
    stop("fitted_msm_models must be a list")
  }
  
  validate_model_structure(fitted_msm_models)
  
  # Loop through model structures (outer level)
  for (model_structure in names(fitted_msm_models)) {
    # Loop through covariate formulas (inner level)
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      
      # Check if this is a failed model
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
        
        # Extract metadata and optimization method if available
        opt_method <- if (is.list(model_entry) && !is.null(model_entry$optimization_method)) {
          model_entry$optimization_method
        } else {
          NA_character_
        }
        
        metadata <- if (is.list(model_entry) && !is.null(model_entry$metadata)) {
          model_entry$metadata
        } else {
          list()
        }
        
        # Include failed models in output with NAs
        failed_tidy <- tibble::tibble(
          model = model_structure,
          formula = formula_name,
          covariate = extract_covariate_label_from_metadata(metadata, formula_name),
          base_variables = list(extract_base_variables_from_metadata(metadata)),
          loglik = NA_real_,
          AIC = NA_real_,
          BIC = NA_real_,
          n_params = NA_integer_,
          n_obs = NA_integer_,
          n_states = NA_integer_,
          n_abs_states = NA_integer_,
          n_trans_states = NA_integer_,
          has_covariates = formula_name != "~ 1",
          n_covariates = extract_n_covariates_from_metadata(metadata),
          has_splines = extract_has_splines_from_metadata(metadata),
          spline_types = extract_spline_types_from_metadata(metadata),
          has_transformations = extract_has_transformations_from_metadata(metadata),
          has_time_varying = extract_has_time_varying_from_metadata(metadata),
          time_varying_type = extract_time_varying_type_from_metadata(metadata),
          constraint_type = extract_constraint_type_from_metadata(metadata),
          optimization_method = opt_method,
          status = if (is.list(model_entry) && !is.null(model_entry$status)) model_entry$status else "failed",
          error_message = if (is.list(model_entry) && !is.null(model_entry$error_message)) model_entry$error_message else "Model is NULL"
        )
        tidied_models <- dplyr::bind_rows(tidied_models, failed_tidy)
        next
      }
      
      # Extract the fitted model object and metadata
      fitted_model <- if (is.list(model_entry) && !is.null(model_entry$fitted_model)) {
        model_entry$fitted_model
      } else {
        model_entry  # Backwards compatibility
      }
      
      metadata <- if (is.list(model_entry) && !is.null(model_entry$metadata)) {
        model_entry$metadata
      } else {
        list()
      }
      
      opt_method <- if (is.list(model_entry) && !is.null(model_entry$optimization_method)) {
        model_entry$optimization_method
      } else {
        NA_character_
      }
      
      if (is.null(fitted_model)) {
        warning(paste("Fitted model is NULL for", model_structure, formula_name))
        next
      }
      
      model_tidy <- tryCatch({
        # Calculate model fit statistics
        loglik <- fitted_model$minus2loglik / -2
        aic <- fitted_model$minus2loglik + 2 * length(fitted_model$estimates)
        bic <- calc_bic_msm(fitted_model)
        n_params <- fitted_model$paramdata$npars
        n_obs <- length(unique(fitted_model[["data"]][["mf"]][["(subject)"]]))
        n_abs_states <- length(absorbing.msm(fitted_model))
        n_trans_states <- length(transient.msm(fitted_model))
        n_states <- n_abs_states + n_trans_states
        
        # Check convergence
        converged <- is.null(fitted_model$opt$convergence) || fitted_model$opt$convergence == 0
        
        tidy_result <- tibble::tibble(
          model = model_structure,
          formula = formula_name,
          covariate = extract_covariate_label_from_metadata(metadata, formula_name),
          base_variables = list(extract_base_variables_from_metadata(metadata)),
          loglik = loglik,
          AIC = aic,
          BIC = bic,
          n_params = n_params,
          n_obs = n_obs,
          n_states = n_states,
          n_abs_states = n_abs_states,
          n_trans_states = n_trans_states,
          has_covariates = formula_name != "~ 1",
          n_covariates = extract_n_covariates_from_metadata(metadata),
          has_splines = extract_has_splines_from_metadata(metadata),
          spline_types = extract_spline_types_from_metadata(metadata),
          has_transformations = extract_has_transformations_from_metadata(metadata),
          has_time_varying = extract_has_time_varying_from_metadata(metadata),
          time_varying_type = extract_time_varying_type_from_metadata(metadata),
          constraint_type = extract_constraint_type_from_metadata(metadata),
          optimization_method = opt_method,
          status = if (converged) "converged" else "converged_with_warning",
          error_message = if (!converged) "Model converged with warnings" else NA_character_
        )
      }, error = function(e) {
        warning(paste("Error tidying model for", model_structure, formula_name, ":", e$message))
        
        # Return failed entry with error info
        tibble::tibble(
          model = model_structure,
          formula = formula_name,
          covariate = extract_covariate_label_from_metadata(metadata, formula_name),
          base_variables = list(extract_base_variables_from_metadata(metadata)),
          loglik = NA_real_,
          AIC = NA_real_,
          BIC = NA_real_,
          n_params = NA_integer_,
          n_obs = NA_integer_,
          n_states = NA_integer_,
          n_abs_states = NA_integer_,
          n_trans_states = NA_integer_,
          has_covariates = formula_name != "~ 1",
          n_covariates = extract_n_covariates_from_metadata(metadata),
          has_splines = extract_has_splines_from_metadata(metadata),
          spline_types = extract_spline_types_from_metadata(metadata),
          has_transformations = extract_has_transformations_from_metadata(metadata),
          has_time_varying = extract_has_time_varying_from_metadata(metadata),
          time_varying_type = extract_time_varying_type_from_metadata(metadata),
          constraint_type = extract_constraint_type_from_metadata(metadata),
          optimization_method = opt_method,
          status = "tidy_failed",
          error_message = e$message
        )
      })
      
      if (!is.null(model_tidy)) {
        tidied_models <- dplyr::bind_rows(tidied_models, model_tidy)
      }
    }
  }
  
  return(tidied_models)
}


### Transition intensities --------------------------------------------

tidy_msm_qmatrix <- function(fitted_msm_models, 
                             covariates_list = NULL,
                             ci_type = "normal",
                             mc.cores = min(6, parallel::detectCores() - 1)) {
  
  validate_model_structure(fitted_msm_models)
  
  tidied_qmats <- data.frame()
  
  # Loop through model structures (outer level)
  for (model_structure in names(fitted_msm_models)) {
    # Loop through covariate formulas (inner level)
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      
      # Check if this is a failed model
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
        warning(paste("Skipping tidying for", model_structure, formula_name, "- model failed to fit"))
        next
      }
      
      # Extract the fitted model object
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in tidy_msm_qmatrix")
      if (is.null(fitted_model)) next
      
      # Extract metadata for better covariate info
      metadata <- if (is.list(model_entry) && !is.null(model_entry$metadata)) {
        model_entry$metadata
      } else {
        list()
      }
      
      model_tidy <- tryCatch({
        qmat_result <- qmatrix.msm(fitted_model, ci = ci_type, cores = mc.cores)
        
        # Use metadata for covariate info, fallback to formula parsing
        covariate_vars <- extract_covariate_label_from_metadata(metadata, formula_name)
        
        # Use tidy() directly and add metadata
        tidy(qmat_result) %>%
          mutate(
            model = model_structure,
            formula = formula_name,
            covariate = covariate_vars,
            constraint_type = extract_constraint_type_from_metadata(metadata),
            has_splines = extract_has_splines_from_metadata(metadata),
            has_time_varying = extract_has_time_varying_from_metadata(metadata),
            cov_name = NA,
            cov_value = NA
          )
      }, error = function(e) {
        warning(paste("Error tidying qmatrix for", model_structure, formula_name, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(model_tidy)) {
        tidied_qmats <- bind_rows(tidied_qmats, model_tidy)
      }
    }
  }
  
  return(tidied_qmats)
}

### Transition probabilities ------------------------------------------

tidy_msm_pmats <- function(fitted_msm_models, 
                           t_values = c(1, 5, 14, 30), 
                           covariates_list = NULL,
                           ci_type = "normal",
                           mc.cores = min(6, parallel::detectCores() - 1)) {
  
  validate_model_structure(fitted_msm_models)
  
  tidied_pmats <- data.frame()
  
  # Loop through model structures (outer level)
  for (model_structure in names(fitted_msm_models)) {
    # Loop through covariate formulas (inner level)
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      
      # Check if this is a failed model
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
        warning(paste("Skipping tidying for", model_structure, formula_name, "- model failed to fit"))
        next
      }
      
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in tidy_msm_pmats")
      if (is.null(fitted_model)) next
      
      # Extract metadata for better covariate info
      metadata <- if (is.list(model_entry) && !is.null(model_entry$metadata)) {
        model_entry$metadata
      } else {
        list()
      }
      
      covariate_vars <- extract_covariate_label_from_metadata(metadata, formula_name)
      
      # Handle multiple time points
      for (t_val in t_values) {
        # Handle multiple covariate combinations
        if (is.null(covariates_list)) {
          model_tidy <- tryCatch({
            pmat_result <- pmatrix.msm(fitted_model, t = t_val, ci = ci_type, cores = mc.cores)
            
            # Use tidy() directly and add metadata
            tidy(pmat_result) %>%
              mutate(
                model = model_structure,
                formula = formula_name,
                covariate = covariate_vars,
                constraint_type = extract_constraint_type_from_metadata(metadata),
                has_splines = extract_has_splines_from_metadata(metadata),
                has_time_varying = extract_has_time_varying_from_metadata(metadata),
                cov_name = NA,
                cov_value = NA,
                t_value = t_val
              )
          }, error = function(e) {
            warning(paste("Error tidying pmat for", model_structure, formula_name, "at t =", t_val, ":", e$message))
            return(NULL)
          })
          
          if (!is.null(model_tidy)) {
            tidied_pmats <- bind_rows(tidied_pmats, model_tidy)
          }
        } else {
          # Handle covariate combinations
          for (cov_combo in covariates_list) {
            model_tidy <- tryCatch({
              pmat_result <- pmatrix.msm(fitted_model, t = t_val, 
                                         covariates = cov_combo, ci = ci_type, cores = mc.cores)
              
              # Use tidy() directly and add metadata
              tidy(pmat_result) %>%
                mutate(
                  model = model_structure,
                  formula = formula_name,
                  covariate = covariate_vars,
                  constraint_type = extract_constraint_type_from_metadata(metadata),
                  has_splines = extract_has_splines_from_metadata(metadata),
                  has_time_varying = extract_has_time_varying_from_metadata(metadata),
                  cov_name = names(cov_combo)[1],
                  cov_value = cov_combo[[1]],
                  t_value = t_val
                )
            }, error = function(e) {
              warning(paste("Error tidying pmat for", model_structure, formula_name, ":", e$message))
              return(NULL)
            })
            
            if (!is.null(model_tidy)) {
              tidied_pmats <- bind_rows(tidied_pmats, model_tidy)
            }
          }
        }
      }
    }
  }
  
  return(tidied_pmats)
}

### Sojourn times -----------------------------------------------------

tidy_msm_sojourns <- function(fitted_msm_models, 
                              covariates_list = NULL) {
  
  validate_model_structure(fitted_msm_models)
  tidied_sojourns <- data.frame()
  
  # Loop through model structures (outer level)
  for (model_structure in names(fitted_msm_models)) {
    # Loop through covariate formulas (inner level)
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      
      # Check if this is a failed model
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
        warning(paste("Skipping tidying for", model_structure, formula_name, "- model failed to fit"))
        next
      }
      
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in tidy_msm_sojourns")
      if (is.null(fitted_model)) next
      
      # Extract metadata for better covariate info
      metadata <- if (is.list(model_entry) && !is.null(model_entry$metadata)) {
        model_entry$metadata
      } else {
        list()
      }
      
      covariate_vars <- extract_covariate_label_from_metadata(metadata, formula_name)
      
      # Handle multiple covariate combinations
      if (is.null(covariates_list)) {
        model_tidy <- tryCatch({
          sojourn_result <- sojourn.msm(fitted_model)
          
          states <- rownames(sojourn_result)
          
          # Use tidy() directly and add metadata
          sojourn_result %>%
            mutate(
              state = states,
              model = model_structure,
              formula = formula_name,
              covariate = covariate_vars,
              constraint_type = extract_constraint_type_from_metadata(metadata),
              has_splines = extract_has_splines_from_metadata(metadata),
              has_time_varying = extract_has_time_varying_from_metadata(metadata),
              cov_name = NA,
              cov_value = NA
            )
        }, error = function(e) {
          warning(paste("Error tidying sojourns for", model_structure, formula_name, ":", e$message))
          return(NULL)
        })
        
        if (!is.null(model_tidy)) {
          tidied_sojourns <- bind_rows(tidied_sojourns, model_tidy)
        }
      } else {
        # Handle covariate combinations
        for (cov_combo in covariates_list) {
          model_tidy <- tryCatch({
            sojourn_result <- sojourn.msm(fitted_model, covariates = cov_combo)
            
            # Use tidy() directly and add metadata
            tidy(sojourn_result) %>%
              mutate(
                model = model_structure,
                formula = formula_name,
                covariate = covariate_vars,
                constraint_type = extract_constraint_type_from_metadata(metadata),
                has_splines = extract_has_splines_from_metadata(metadata),
                has_time_varying = extract_has_time_varying_from_metadata(metadata),
                cov_name = names(cov_combo)[1],
                cov_value = cov_combo[[1]]
              )
          }, error = function(e) {
            warning(paste("Error tidying sojourns for", model_structure, formula_name, ":", e$message))
            return(NULL)
          })
          
          if (!is.null(model_tidy)) {
            tidied_sojourns <- bind_rows(tidied_sojourns, model_tidy)
          }
        }
      }
    }
  }
  
  return(tidied_sojourns)
}

### Hazard ratios -----------------------------------------------------

tidy_msm_hazards <- function(fitted_msm_models, hazard_scale = 1) {
  
  validate_model_structure(fitted_msm_models)
  tidied_hrs <- data.frame()
  
  # Loop through model structures (outer level)
  for (model_structure in names(fitted_msm_models)) {
    # Loop through covariate formulas (inner level)
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      
      # Check if this is a failed model
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
        warning(paste("Skipping tidying for", model_structure, formula_name, "- model failed to fit"))
        next
      }
      
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in tidy_msm_hazards")
      
      if (is.null(fitted_model)) {
        warning(paste("Fitted model is NULL for", model_structure, formula_name))
        next
      }
      
      # Extract metadata for better covariate info and checking
      metadata <- if (is.list(model_entry) && !is.null(model_entry$metadata)) {
        model_entry$metadata
      } else {
        list()
      }
      
      # Check if model has covariates using metadata first, then fallback
      has_covariates <- FALSE
      if (!is.null(metadata$has_covariates)) {
        has_covariates <- metadata$has_covariates
      } else if (!is.null(metadata$n_covariates)) {
        has_covariates <- metadata$n_covariates > 0
      } else {
        # Fallback to checking fitted model and formula
        has_covariates <- !is.null(fitted_model$covariates) && formula_name != "~ 1"
      }
      
      if (!has_covariates) {
        warning(paste("Skipping hazards for", model_structure, formula_name, "- model has no covariates"))
        next
      }
      
      # Extract covariate information from metadata
      covariate_vars <- extract_covariate_label_from_metadata(metadata, formula_name)
      base_variables <- extract_base_variables_from_metadata(metadata)
      
      model_tidy <- tryCatch({
        hr_result <- hazard.msm(fitted_model, hazard.scale = hazard_scale)
        
        # Manual tidying of hazard ratios (list of matrices)
        hr_tidy_list <- map_dfr(names(hr_result), function(covariate_name) {
          hr_matrix <- hr_result[[covariate_name]]
          
          # Extract transitions (row names)
          transitions <- rownames(hr_matrix)
          
          # Determine if this is a spline term
          is_spline_term <- is_spline_covariate_term(covariate_name, metadata)
          original_variable <- extract_original_variable_name(covariate_name, metadata)
          
          # Create tidy data frame
          data.frame(
            transition = transitions,
            covariate_term = covariate_name,
            original_variable = original_variable,
            is_spline_term = is_spline_term,
            hazard_ratio = hr_matrix[, "HR"],
            hr_lower = hr_matrix[, "L"], 
            hr_upper = hr_matrix[, "U"],
            log_hr = log(hr_matrix[, "HR"]),
            stringsAsFactors = FALSE
          )
        })
        
        # Add comprehensive metadata
        hr_tidy <- hr_tidy_list %>%
          mutate(
            model = model_structure,
            formula = formula_name,
            covariate = covariate_vars,
            constraint_type = extract_constraint_type_from_metadata(metadata),
            has_splines = extract_has_splines_from_metadata(metadata),
            has_time_varying = extract_has_time_varying_from_metadata(metadata),
            time_varying_type = extract_time_varying_type_from_metadata(metadata),
            hazard_scale = hazard_scale,
            n_base_variables = length(base_variables),
            base_variables_list = list(base_variables)
          )
      }, error = function(e) {
        warning(paste("Error tidying hazards for", model_structure, formula_name, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(model_tidy)) {
        tidied_hrs <- bind_rows(tidied_hrs, model_tidy)
      }
    }
  }
  
  return(tidied_hrs)
}

#' Calculate hazard ratios for spline covariates relative to reference value
#' @param fitted_msm_models Nested list from fit_msm_models
#' @param patient_data Original patient data (to get variable ranges)
#' @param n_points Number of evaluation points along each spline (default: 20)
#' @param reference_level Named list of reference values, or NULL for median (default: NULL)
#' @param ci_method Method for confidence intervals: "delta" or "normal" (default: "delta")
#' @return Data frame with variable values, HRs, CIs, by transition and model
calculate_spline_hazard_ratios <- function(fitted_msm_models,
                                           patient_data,
                                           n_points = 20,
                                           reference_level = NULL,
                                           ci_method = "delta") {
  
  # Validate inputs
  validate_model_structure(fitted_msm_models)
  
  if (!ci_method %in% c("delta", "normal")) {
    stop("ci_method must be 'delta' or 'normal'")
  }
  
  cat("=== CALCULATING SPLINE HAZARD RATIOS ===\n")
  
  all_results <- list()
  
  # Loop through models
  for (model_structure in names(fitted_msm_models)) {
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      
      # Skip failed models
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && 
           model_entry$status == "failed")) {
        next
      }
      
      fitted_model <- model_entry$fitted_model
      metadata <- model_entry$metadata
      
      if (is.null(fitted_model) || is.null(metadata)) next
      
      # Check if model has splines
      if (!extract_has_splines_from_metadata(metadata)) {
        next
      }
      
      cat("Processing:", model_structure, "|", formula_name, "\n")
      
      # Get spline variables and metadata
      spline_config <- metadata$spline_config
      if (is.null(spline_config) || is.null(spline_config$spline_metadata)) {
        next
      }
      
      # Filter patient data for this model
      model_data <- patient_data %>% filter(model == model_structure)
      
      # Process each spline variable
      for (spline_var in names(spline_config$spline_metadata)) {
        
        var_meta <- spline_config$spline_metadata[[spline_var]]
        
        cat("  Variable:", spline_var, "\n")
        
        # Get variable range from data
        var_values <- model_data[[spline_var]]
        if (is.null(var_values) || all(is.na(var_values))) {
          warning(paste("Variable", spline_var, "not found in data"))
          next
        }
        
        var_range <- range(var_values, na.rm = TRUE)
        
        # Set reference level
        ref_value <- if (!is.null(reference_level) && spline_var %in% names(reference_level)) {
          reference_level[[spline_var]]
        } else {
          median(var_values, na.rm = TRUE)
        }
        
        cat("    Range: [", var_range[1], ",", var_range[2], "], Reference:", ref_value, "\n")
        
        # Create evaluation grid
        eval_grid <- seq(var_range[1], var_range[2], length.out = n_points)
        
        # Get base variables (non-spline covariates to hold constant)
        base_variables <- extract_base_variables_from_metadata(metadata)
        other_vars <- setdiff(base_variables, spline_var)
        
        # Create covariate list for reference value
        ref_covariates <- list()
        ref_covariates[[spline_var]] <- ref_value
        
        # Set other covariates to their medians/modes
        for (other_var in other_vars) {
          if (other_var %in% names(model_data)) {
            if (is.numeric(model_data[[other_var]])) {
              ref_covariates[[other_var]] <- median(model_data[[other_var]], na.rm = TRUE)
            } else {
              # For factors, use most common level
              ref_covariates[[other_var]] <- names(sort(table(model_data[[other_var]]), 
                                                        decreasing = TRUE))[1]
            }
          }
        }
        
        # Calculate hazard ratios at reference (should be ~1)
        ref_hazards <- tryCatch({
          hazard.msm(fitted_model, covariates = ref_covariates, 
                     ci = ci_method, hazard.scale = 1)
        }, error = function(e) {
          warning(paste("Failed to calculate reference hazards:", e$message))
          return(NULL)
        })
        
        if (is.null(ref_hazards)) next
        
        # Calculate hazard ratios at each grid point
        grid_results <- list()
        
        for (i in seq_along(eval_grid)) {
          eval_value <- eval_grid[i]
          
          # Create covariate list for this evaluation point
          eval_covariates <- ref_covariates
          eval_covariates[[spline_var]] <- eval_value
          
          # Get hazards at this point
          eval_hazards <- tryCatch({
            hazard.msm(fitted_model, covariates = eval_covariates,
                       ci = ci_method, hazard.scale = 1)
          }, error = function(e) {
            NULL
          })
          
          if (is.null(eval_hazards)) next
          
          # Calculate HR relative to reference for each transition
          for (cov_name in names(eval_hazards)) {
            
            # Get reference hazard for this covariate/transition
            ref_hr <- ref_hazards[[cov_name]]
            eval_hr <- eval_hazards[[cov_name]]
            
            # Extract transition info from covariate name
            # Format is typically "variablename.transition" or just shows transition
            transitions <- rownames(eval_hr)
            
            for (trans in transitions) {
              
              # Calculate HR relative to reference
              # HR = exp(log(eval) - log(ref))
              hr <- eval_hr[trans, "HR"] / ref_hr[trans, "HR"]
              
              # For CIs, use ratio of the HR bounds
              hr_lower <- eval_hr[trans, "L"] / ref_hr[trans, "U"]
              hr_upper <- eval_hr[trans, "U"] / ref_hr[trans, "L"]
              
              grid_results[[length(grid_results) + 1]] <- data.frame(
                model = model_structure,
                formula = formula_name,
                spline_variable = spline_var,
                variable_value = eval_value,
                reference_value = ref_value,
                transition = trans,
                hazard_ratio = hr,
                hr_lower = hr_lower,
                hr_upper = hr_upper,
                log_hr = log(hr),
                log_hr_lower = log(hr_lower),
                log_hr_upper = log(hr_upper),
                stringsAsFactors = FALSE
              )
            }
          }
        }
        
        if (length(grid_results) > 0) {
          var_results <- do.call(rbind, grid_results)
          
          # Add metadata
          var_results$constraint_type <- extract_constraint_type_from_metadata(metadata)
          var_results$spline_df <- var_meta$df
          var_results$spline_type <- var_meta$type
          
          combo_key <- paste(model_structure, formula_name, spline_var, sep = "___")
          all_results[[combo_key]] <- var_results
        }
      }
    }
  }
  
  if (length(all_results) == 0) {
    warning("No spline hazard ratios calculated")
    return(data.frame())
  }
  
  # Combine all results
  final_results <- do.call(rbind, all_results)
  rownames(final_results) <- NULL
  
  cat("=== COMPLETED ===\n")
  cat("Calculated HRs for", length(unique(final_results$spline_variable)), 
      "spline variables\n")
  cat("Across", length(unique(final_results$model)), "model structures\n")
  cat("Total rows:", nrow(final_results), "\n")
  
  return(final_results)
}

### State prevalence --------------------------------------------------

tidy_msm_prevalences <- function(fitted_msm_models, 
                                 time_points = seq(1, 30, by = 1),
                                 covariates_list = NULL, 
                                 mc.cores = min(8, parallel::detectCores() - 1),
                                 ci = TRUE,
                                 ci_method = "normal",
                                 use_approximation = FALSE,
                                 approx_method = "midpoint",
                                 show_progress = TRUE) {
  
  # Remove all future dependencies
  library(dplyr)
  
  validate_model_structure(fitted_msm_models)
  
  cat("=== MSM PREVALENCE PROCESSING ===\n")
  cat("CI:", ifelse(ci, paste("enabled (", ci_method, ")", sep = ""), "disabled"), "\n")
  cat("Approximation:", ifelse(use_approximation, paste("enabled (", approx_method, ")", sep = ""), "disabled"), "\n")
  
  # Validate approximation method
  if (use_approximation && !approx_method %in% c("start", "midpoint")) {
    cat("Invalid approximation method '", approx_method, "'. Using 'midpoint' instead.\n", sep = "")
    approx_method <- "midpoint"
  }
  
  start_time <- Sys.time()
  
  # Build model combinations - more memory efficient
  model_combinations <- list()
  combination_count <- 0
  
  for (model_structure in names(fitted_msm_models)) {
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      if (!is.null(model_entry) && 
          is.list(model_entry) && 
          model_entry$status == "converged") {
        combination_count <- combination_count + 1
        # Store only the essential info, not the full model objects
        model_combinations[[combination_count]] <- list(
          model_structure = model_structure, 
          formula_name = formula_name,
          index = combination_count
        )
      }
    }
  }
  
  cat("Processing", length(model_combinations), "model combinations with", 
      length(time_points), "time points\n")
  
  # Define the worker function - keep it lightweight
  process_model_combination <- function(combo_info) {
    # Load required packages in worker
    if (!require("msm", quietly = TRUE)) return(NULL)
    if (!require("dplyr", quietly = TRUE)) return(NULL)
    
    model_structure <- combo_info$model_structure
    formula_name <- combo_info$formula_name
    
    # Extract model object and metadata inside worker to minimize data transfer
    model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
    fitted_model <- model_entry$fitted_model
    metadata <- model_entry$metadata
    
    if (is.null(fitted_model)) return(NULL)
    
    # Get covariate information from metadata
    covariate_vars <- extract_covariate_label_from_metadata(metadata, formula_name)
    
    # Process prevalences
    tryCatch({
      if (is.null(covariates_list)) {
        # Build arguments for prevalence.msm
        prev_args <- list(x = fitted_model, times = time_points)
        if (ci) prev_args$ci <- ci_method
        
        if (use_approximation) {
          if ("interp" %in% names(formals(prevalence.msm))) {
            valid_methods <- c("start", "midpoint")
            if (approx_method %in% valid_methods) {
              prev_args$interp <- approx_method
            } else {
              prev_args$interp <- "midpoint"
            }
          }
        }
        
        prev_result <- do.call(prevalence.msm, prev_args)
        
        prevalence_to_tib(prev_result, fitted_model, has_ci = ci) %>%
          mutate(
            model = model_structure,
            formula = formula_name,
            covariate = covariate_vars,
            constraint_type = extract_constraint_type_from_metadata(metadata),
            has_splines = extract_has_splines_from_metadata(metadata),
            has_time_varying = extract_has_time_varying_from_metadata(metadata),
            time_varying_type = extract_time_varying_type_from_metadata(metadata),
            cov_name = NA,
            cov_value = NA,
            .after = time
          )
      } else {
        # Handle multiple covariate combinations
        cov_results <- list()
        for (i in seq_along(covariates_list)) {
          cov_combo <- covariates_list[[i]]
          
          prev_args <- list(x = fitted_model, times = time_points, covariates = cov_combo)
          if (ci) prev_args$ci <- ci_method
          
          if (use_approximation) {
            if ("interp" %in% names(formals(prevalence.msm))) {
              valid_methods <- c("start", "midpoint")
              if (approx_method %in% valid_methods) {
                prev_args$interp <- approx_method
              } else {
                prev_args$interp <- "midpoint"
              }
            }
          }
          
          prev_result <- do.call(prevalence.msm, prev_args)
          
          cov_results[[i]] <- prevalence_to_tib(prev_result, fitted_model, has_ci = ci) %>%
            mutate(
              model = model_structure,
              formula = formula_name,
              covariate = covariate_vars,
              constraint_type = extract_constraint_type_from_metadata(metadata),
              has_splines = extract_has_splines_from_metadata(metadata),
              has_time_varying = extract_has_time_varying_from_metadata(metadata),
              time_varying_type = extract_time_varying_type_from_metadata(metadata),
              cov_name = names(cov_combo)[1],
              cov_value = cov_combo[[1]],
              .after = time
            )
        }
        do.call(bind_rows, cov_results)
      }
    }, error = function(e) {
      # Return minimal error info to reduce memory usage
      return(list(error = e$message, model = model_structure, formula = formula_name))
    })
  }
  
  # Choose processing method based on mc.cores
  if (mc.cores == 1) {
    cat("Using sequential processing\n")
    
    # True sequential processing - most memory efficient
    prevalence_list <- vector("list", length(model_combinations))
    for (i in seq_along(model_combinations)) {
      prevalence_list[[i]] <- process_model_combination(model_combinations[[i]])
      
      # Optional progress reporting
      if (i %% 5 == 0) {
        cat("Completed", i, "of", length(model_combinations), "combinations\n")
      }
    }
    
  } else {
    cat("Using", mc.cores, "parallel workers (mclapply)\n")
    
    # Use mclapply instead of future - more memory efficient
    prevalence_list <- parallel::mclapply(
      model_combinations, 
      process_model_combination,
      mc.cores = mc.cores,
      mc.set.seed = TRUE
    )
  }
  
  # Check for errors and combine results - memory efficient
  error_count <- 0
  valid_results <- vector("list", length(prevalence_list))
  valid_count <- 0
  
  for (i in seq_along(prevalence_list)) {
    result <- prevalence_list[[i]]
    if (is.list(result) && "error" %in% names(result)) {
      error_count <- error_count + 1
      cat("Error in combination", i, ":", result$model, result$formula, ":", result$error, "\n")
    } else if (is.data.frame(result)) {
      valid_count <- valid_count + 1
      valid_results[[valid_count]] <- result
    }
  }
  
  if (error_count > 0) {
    cat("Found", error_count, "errors in prevalence calculations\n")
  }
  
  if (valid_count == 0) {
    warning("No valid prevalence results calculated")
    return(data.frame())
  }
  
  # Combine results efficiently
  valid_results <- valid_results[1:valid_count]  # Remove empty slots
  final_result <- do.call(bind_rows, valid_results)
  
  total_time <- Sys.time() - start_time
  cat("=== COMPLETED in", as.numeric(total_time, units = "secs"), "seconds ===\n")
  cat("Final dimensions:", nrow(final_result), "x", ncol(final_result), "\n")
  
  return(final_result)
}

prevalence_to_tib <- function(prevalence_result, fitted_model, has_ci = TRUE) {
  # Get basic info with error checking
  total_at_risk <- prevalence_result$Observed[1, "Total"]
  time_points <- as.numeric(rownames(prevalence_result$Observed))
  
  # Debug info - check dimensions
  cat("DEBUG prevalence_to_tib:\n")
  cat("  Observed dims:", dim(prevalence_result$Observed), "\n")
  cat("  Time points:", length(time_points), "\n")
  
  # Handle different structures based on CI computation
  observed_cols <- colnames(prevalence_result$Observed)
  cat("  Observed columns:", paste(observed_cols, collapse = ", "), "\n")
  
  # Expected data structure varies with CI
  if (is.list(prevalence_result$Expected) && "estimates" %in% names(prevalence_result$Expected)) {
    expected_cols <- colnames(prevalence_result$Expected$estimates)
    expected_data <- prevalence_result$Expected$estimates
    has_ci_data <- has_ci && !is.null(prevalence_result$Expected$ci)
  } else {
    expected_cols <- colnames(prevalence_result$Expected)
    expected_data <- prevalence_result$Expected
    has_ci_data <- FALSE
  }
  
  cat("  Expected columns:", paste(expected_cols, collapse = ", "), "\n")
  cat("  Expected dims:", dim(expected_data), "\n")
  
  # Get state columns (exclude "Total") with bounds checking
  observed_states <- observed_cols[observed_cols != "Total"]
  expected_states <- expected_cols[expected_cols != "Total"]
  
  cat("  Observed states:", paste(observed_states, collapse = ", "), "\n")
  cat("  Expected states:", paste(expected_states, collapse = ", "), "\n")
  
  # Ensure we have matching states
  if (length(observed_states) != length(expected_states)) {
    cat("  WARNING: Mismatch in state counts - observed:", length(observed_states), "expected:", length(expected_states), "\n")
    # Use the shorter list to avoid out of bounds
    n_states <- min(length(observed_states), length(expected_states))
    observed_states <- observed_states[1:n_states]
    expected_states <- expected_states[1:n_states]
  }
  
  # Get state names from fitted model for consistency
  if (!is.null(fitted_model) && !is.null(fitted_model$qmodel)) {
    model_states <- rownames(fitted_model$qmodel$qmatrix)
    cat("  Model states from qmatrix:", paste(model_states, collapse = ", "), "\n")
  } else {
    model_states <- expected_states
    cat("  Using expected states as model states\n")
  }
  
  # Build result data frame with bounds checking
  result_list <- list()
  
  for (t_idx in seq_along(time_points)) {
    for (s_idx in seq_along(observed_states)) {
      # Safe matrix access with error handling
      tryCatch({
        # Use character names directly for matrix subsetting
        obs_state_name <- observed_states[s_idx]
        exp_state_name <- expected_states[s_idx]
        
        obs_count <- prevalence_result$Observed[t_idx, obs_state_name]
        obs_pct <- prevalence_result$`Observed percentages`[t_idx, obs_state_name]
        exp_count <- expected_data[t_idx, exp_state_name]
        
        # Handle expected percentages structure
        if (is.list(prevalence_result$`Expected percentages`)) {
          exp_pct <- prevalence_result$`Expected percentages`$estimates[t_idx, exp_state_name]
        } else {
          exp_pct <- prevalence_result$`Expected percentages`[t_idx, exp_state_name]
        }
        
        # Use expected state as state name
        state_name <- exp_state_name

        # Handle confidence intervals with bounds checking
        if (has_ci_data && !is.null(prevalence_result$Expected$ci)) {
          if (dim(prevalence_result$Expected$ci)[3] >= 2) {
            exp_count_ll <- unname(prevalence_result$Expected$ci[t_idx, s_idx, 1])
            exp_count_ul <- unname(prevalence_result$Expected$ci[t_idx, s_idx, 2])
          } else {
            exp_count_ll <- exp_count_ul <- NA_real_
          }
          
          if (is.list(prevalence_result$`Expected percentages`) && 
              !is.null(prevalence_result$`Expected percentages`$ci) &&
              dim(prevalence_result$`Expected percentages`$ci)[3] >= 2) {
            exp_pct_ll <- unname(prevalence_result$`Expected percentages`$ci[t_idx, s_idx, 1])
            exp_pct_ul <- unname(prevalence_result$`Expected percentages`$ci[t_idx, s_idx, 2])
          } else {
            exp_pct_ll <- exp_pct_ul <- NA_real_
          }
        } else {
          exp_count_ll <- exp_count_ul <- exp_pct_ll <- exp_pct_ul <- NA_real_
        }
        
        result_list[[length(result_list) + 1]] <- data.frame(
          time = time_points[t_idx],
          state_num = obs_state_name,  # Changed from observed_states[s_idx]
          state = state_name,
          observed_count = obs_count,
          observed_percentage = obs_pct,
          expected_count = exp_count,
          expected_percentage = exp_pct,
          expected_count_ll = exp_count_ll,
          expected_count_ul = exp_count_ul,
          expected_percentage_ll = exp_pct_ll,
          expected_percentage_ul = exp_pct_ul,
          stringsAsFactors = FALSE
        )
        
      }, error = function(e) {
        cat("  ERROR accessing matrix element t_idx:", t_idx, "s_idx:", s_idx, 
            "obs_state:", observed_states[s_idx], "exp_state:", expected_states[s_idx],
            "error:", e$message, "\n")
      })
    }
  }
  
  if (length(result_list) == 0) {
    cat("  ERROR: No valid results created\n")
    return(data.frame())
  }
  
  # Combine results
  combined_prev <- do.call(bind_rows, result_list)
  
  # Calculate N_at_risk with error handling
  death_recovery_states <- c("D", "R")
  
  N_at_risk_df <- tryCatch({
    combined_prev %>%
      filter(state %in% death_recovery_states) %>%
      group_by(time) %>%
      summarise(N_at_risk = total_at_risk - sum(observed_count, na.rm = TRUE), .groups = "drop")
  }, error = function(e) {
    cat("  ERROR calculating N_at_risk:", e$message, "\n")
    data.frame(time = time_points, N_at_risk = total_at_risk)
  })
  
  # Add final calculations with error handling
  final_result <- tryCatch({
    combined_prev %>%
      left_join(N_at_risk_df, by = "time") %>%
      mutate(
        count_residual = observed_count - expected_count,
        percentage_residual = observed_percentage - expected_percentage,
        binomial_variance = N_at_risk * (expected_percentage/100) * (1 - expected_percentage/100),
        standardized_residual = ifelse(!is.na(binomial_variance) & binomial_variance > 0, 
                                       count_residual / sqrt(binomial_variance), NA_real_),
        residual_category = case_when(
          is.na(standardized_residual) ~ "Unknown",
          abs(standardized_residual) > 2 ~ "Large",
          abs(standardized_residual) > 1 ~ "Moderate", 
          TRUE ~ "Small"
        )
      ) %>%
      arrange(time, state)
  }, error = function(e) {
    cat("  ERROR in final calculations:", e$message, "\n")
    combined_prev
  })
  
  return(final_result)
}


### Predictive performance --------------------------------------------

#' Cross-validation with covariate-specific calibration assessment
#' @param fitted_models Nested list of fitted models: fitted_models[[model_structure]][[formula]]
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param k_folds Number of CV folds
#' @param prediction_times Vector of prediction horizons
#' @param stratify_by Variable to stratify folds
#' @param calibration_covariates List of covariate values to evaluate calibration at
#' @param calibration_subgroups List of variables to stratify calibration by
#' @param parallel Logical, whether to use parallel processing (auto-detect if NULL)
#' @param n_cores Number of cores for parallel processing
#' @return Data frame with CV results including covariate-specific calibration
calculate_predictive_performance <- function(patient_data, 
                                             fitted_models, crude_rates,
                                             k_folds = 5, 
                                             prediction_times = c(14, 30),
                                             stratify_by = "final_state",
                                             calibration_covariates = NULL,
                                             calibration_subgroups = NULL,
                                             parallel = NULL, 
                                             n_cores = min(6, parallel::detectCores() - 1)) {
  
  validate_model_structure(fitted_models)
  
  # Auto-detect parallel processing based on total number of model/formula combinations
  total_combinations <- sum(sapply(fitted_models, length))
  if (is.null(parallel)) {
    parallel <- total_combinations > 2 && parallel::detectCores() > 2
  }
  
  # Setup parallel backend if requested
  if (parallel) {
    if (!requireNamespace("furrr", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      warning("furrr/future packages not available. Falling back to sequential processing.")
      parallel <- FALSE
    } else {
      library(furrr)
      library(future)
      plan(multisession, workers = n_cores)
      on.exit(plan(sequential), add = TRUE)
      cat("Using parallel processing with", n_cores, "cores\n")
    }
  }
  
  if (!parallel) cat("Using sequential processing\n")
  
  # Create stratified fold assignments
  unique_patients <- patient_data %>%
    group_by(deid_enc_id) %>%
    summarise(
      final_state = last(state),
      model = first(model),
      .groups = "drop"
    )
  
  set.seed(123)
  fold_assignments <- unique_patients %>%
    group_by(!!sym(stratify_by)) %>%
    mutate(fold = sample(rep(1:k_folds, length.out = n()))) %>%
    ungroup() %>%
    select(deid_enc_id, fold)
  
  cv_results <- list()
  
  # Loop through model structures
  for (model_structure in names(fitted_models)) {
    # Loop through formulas within each model structure
    for (formula_name in names(fitted_models[[model_structure]])) {
      model_entry <- fitted_models[[model_structure]][[formula_name]]
      
      # Check if this is a failed model
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
        warning(paste("Skipping CV for", model_structure, formula_name, "- model failed to fit"))
        next
      }
      
      # Extract the fitted model object and metadata
      fitted_model <- if (is.list(model_entry) && !is.null(model_entry$fitted_model)) {
        model_entry$fitted_model
      } else {
        model_entry  # Backwards compatibility
      }
      
      metadata <- if (is.list(model_entry) && !is.null(model_entry$metadata)) {
        model_entry$metadata
      } else {
        list()
      }
      
      if (is.null(fitted_model)) {
        warning(paste("Fitted model is NULL for", model_structure, formula_name))
        next
      }
      
      cat("Running CV for model:", model_structure, "formula:", formula_name, "\n")
      
      model_data <- patient_data %>% filter(model == model_structure)
      model_fold_assignments <- fold_assignments %>%
        filter(deid_enc_id %in% unique(model_data$deid_enc_id))
      
      # Extract covariate formula from the formula name
      covariate_formula <- if (formula_name == "~ 1") {
        NULL
      } else {
        as.formula(formula_name)
      }
      
      combination_key <- paste(model_structure, formula_name, sep = "_")
      
      # Run CV folds with enhanced calibration
      fold_results <- if (parallel) {
        # First attempt with current settings
        tryCatch({
          future_map_dfr(1:k_folds, function(fold) {
            cv_fold_core(fold, model_fold_assignments, model_data, crude_rates[[model_structure]], 
                         covariate_formula, prediction_times, calibration_covariates, calibration_subgroups, model_entry)
          }, .options = furrr_options(seed = TRUE))
        }, error = function(e) {
          if (grepl("future.globals.maxSize", e$message)) {
            cat("  Globals size exceeded for", model_structure, formula_name, "- increasing to 750 MiB and retrying...\n")
            
            # Store current setting
            current_maxsize <- getOption("future.globals.maxSize")
            
            # Increase to 750 MiB temporarily
            options(future.globals.maxSize = 750 * 1024^2)
            
            tryCatch({
              result <- future_map_dfr(1:k_folds, function(fold) {
                cv_fold_core(fold, model_fold_assignments, model_data, crude_rates[[model_structure]], 
                             covariate_formula, prediction_times, calibration_covariates, calibration_subgroups, model_entry)
              }, .options = furrr_options(seed = TRUE))
              
              cat("  Successfully completed", model_structure, formula_name, "with increased maxSize\n")
              return(result)
              
            }, error = function(e2) {
              cat("  Failed even with increased maxSize for", model_structure, formula_name, ":", e2$message, "\n")
              stop("Model ", model_structure, " ", formula_name, " failed: ", e2$message)
            }, finally = {
              # Restore original setting
              options(future.globals.maxSize = current_maxsize)
            })
          } else {
            # Re-throw if it's a different error
            stop("Model ", model_structure, " ", formula_name, " failed: ", e$message)
          }
        })
      } else {
        map_dfr(1:k_folds, function(fold) {
          cv_fold_core(fold, model_fold_assignments, model_data, crude_rates[[model_structure]], 
                       covariate_formula, prediction_times, calibration_covariates, calibration_subgroups, model_entry)
        })
      }
      
      # Extract covariate information from metadata
      covariate_vars <- extract_covariate_label_from_metadata(metadata, formula_name)
      
      # Aggregate results with enhanced metadata
      base_results <- fold_results %>%
        filter(is.na(prediction_time) & is.na(subgroup_var)) %>%
        summarise(
          model = model_structure,
          formula = formula_name,
          covariate = covariate_vars,
          constraint_type = extract_constraint_type_from_metadata(metadata),
          has_splines = extract_has_splines_from_metadata(metadata),
          has_time_varying = extract_has_time_varying_from_metadata(metadata),
          time_varying_type = extract_time_varying_type_from_metadata(metadata),
          n_base_variables = length(extract_base_variables_from_metadata(metadata)),
          n_folds_converged = sum(model_converged, na.rm = TRUE),
          total_los_mae_cv = mean(total_los_mae, na.rm = TRUE),
          total_los_mae_se = sd(total_los_mae, na.rm = TRUE) / sqrt(n()),
          total_los_mae_soj_cv = mean(total_los_mae_soj, na.rm = TRUE),
          total_los_mae_soj_se = sd(total_los_mae_soj, na.rm = TRUE) / sqrt(n()),
          total_los_mae_vis_soj_cv = mean(total_los_mae_vis_soj, na.rm = TRUE),
          total_los_mae_vis_soj_se = sd(total_los_mae_vis_soj, na.rm = TRUE) / sqrt(n()),
          days_severe_mae_cv = mean(days_severe_mae, na.rm = TRUE),
          days_severe_mae_se = sd(days_severe_mae, na.rm = TRUE) / sqrt(n()),
          days_severe_mae_soj_cv = mean(days_severe_mae_soj, na.rm = TRUE),
          days_severe_mae_soj_se = sd(days_severe_mae_soj, na.rm = TRUE) / sqrt(n()),
          days_severe_mae_vis_soj_cv = mean(days_severe_mae_vis_soj, na.rm = TRUE),
          days_severe_mae_vis_soj_se = sd(days_severe_mae_vis_soj, na.rm = TRUE) / sqrt(n()),
          death_auc_cv = mean(death_auc, na.rm = TRUE),
          death_auc_se = sd(death_auc, na.rm = TRUE) / sqrt(n()),
          severe_auc_cv = mean(severe_auc, na.rm = TRUE),
          severe_auc_se = sd(severe_auc, na.rm = TRUE) / sqrt(n()),
          calibration_slope_death_cv = mean(calibration_slope_death, na.rm = TRUE),
          calibration_intercept_death_cv = mean(calibration_intercept_death, na.rm = TRUE),
          brier_death_cv = mean(brier_death, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Time-specific calibration results
      time_specific_results <- fold_results %>%
        filter(!is.na(prediction_time) & is.na(subgroup_var)) %>%
        group_by(prediction_time) %>%
        summarise(
          model = model_structure,
          formula = formula_name,
          covariate = covariate_vars,
          constraint_type = extract_constraint_type_from_metadata(metadata),
          has_splines = extract_has_splines_from_metadata(metadata),
          has_time_varying = extract_has_time_varying_from_metadata(metadata),
          time_varying_type = extract_time_varying_type_from_metadata(metadata),
          prediction_time = first(prediction_time),
          calibration_slope_death_time = mean(calibration_slope_death, na.rm = TRUE),
          calibration_intercept_death_time = mean(calibration_intercept_death, na.rm = TRUE),
          brier_death_time = mean(brier_death, na.rm = TRUE),
          death_auc_time = mean(death_auc, na.rm = TRUE),
          severe_auc_time = mean(severe_auc, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Subgroup-specific calibration results
      subgroup_specific_results <- fold_results %>%
        filter(is.na(prediction_time) & !is.na(subgroup_var)) %>%
        group_by(subgroup_var, subgroup_value) %>%
        summarise(
          model = model_structure,
          formula = formula_name,
          covariate = covariate_vars,
          constraint_type = extract_constraint_type_from_metadata(metadata),
          has_splines = extract_has_splines_from_metadata(metadata),
          has_time_varying = extract_has_time_varying_from_metadata(metadata),
          time_varying_type = extract_time_varying_type_from_metadata(metadata),
          subgroup_var = first(subgroup_var),
          subgroup_value = first(subgroup_value),
          calibration_slope_death_subgroup = mean(calibration_slope_death, na.rm = TRUE),
          calibration_intercept_death_subgroup = mean(calibration_intercept_death, na.rm = TRUE),
          brier_death_subgroup = mean(brier_death, na.rm = TRUE),
          death_auc_subgroup = mean(death_auc, na.rm = TRUE),
          n_subgroup = sum(!is.na(calibration_slope_death)),
          .groups = "drop"
        )
      
      cv_results[[combination_key]] <- list(
        base = base_results,
        time_specific = time_specific_results,
        subgroup_specific = subgroup_specific_results
      )
    }
  }
  
  return(cv_results)
}

#' Enhanced core function for single fold cross-validation with calibration assessment
cv_fold_core <- function(fold_num, fold_assignments, 
                         model_data, crude_rate, 
                         covariate_formula = NULL, 
                         prediction_times = c(14, 30),
                         calibration_covariates = NULL,
                         calibration_subgroups = NULL,
                         original_model_entry = NULL) {
  
  # Split data
  train_patients <- fold_assignments$deid_enc_id[fold_assignments$fold != fold_num]
  test_patients <- fold_assignments$deid_enc_id[fold_assignments$fold == fold_num]
  
  train_data <- model_data %>% filter(deid_enc_id %in% train_patients)
  test_data <- model_data %>% filter(deid_enc_id %in% test_patients)
  
  if (nrow(train_data) == 0 || nrow(test_data) == 0) {
    return(create_empty_fold_result(fold_num, nrow(train_data), nrow(test_data)))
  }
  
  # Recompute crude rates for training data 
  train_crude_rates <- tryCatch({
    temp_qmat_list <- list()
    temp_qmat_list[[unique(train_data$model)[1]]] <- crude_rate$qmat
    calc_crude_init_rates(train_data, temp_qmat_list)
  }, error = function(e) {
    warning(paste("Crude rate calculation failed in fold", fold_num, ":", e$message))
    temp_list <- list()
    temp_list[[unique(train_data$model)[1]]] <- crude_rate
    return(temp_list)
  })
  
  # Extract comprehensive metadata from the original fitted model
  metadata <- if (!is.null(original_model_entry) && is.list(original_model_entry) && !is.null(original_model_entry$metadata)) {
    original_model_entry$metadata
  } else {
    list()
  }
  
  # Refit model on training data using the unified fit_msm_models function
  fitted_models_fold <- refit_model_for_fold_unified(train_data, train_crude_rates, metadata)
  
  # Extract the fitted model from nested structure
  model_structure_name <- names(fitted_models_fold)[1]
  if (is.null(model_structure_name)) {
    return(create_empty_fold_result(fold_num, nrow(train_data), nrow(test_data)))
  }
  
  formula_key <- names(fitted_models_fold[[model_structure_name]])[1]
  model_entry <- fitted_models_fold[[model_structure_name]][[formula_key]]
  
  fitted_model <- extract_fitted_model(fitted_models_fold, model_structure_name, formula_key, "in cv_fold_core")
  if (is.null(fitted_model)) {
    return(create_empty_fold_result(fold_num, nrow(train_data), nrow(test_data)))
  }
  
  if (is.null(fitted_model) || 
      (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
    return(create_empty_fold_result(fold_num, nrow(train_data), nrow(test_data)))
  }
  
  model_converged <- is.null(fitted_model$opt$convergence) || fitted_model$opt$convergence == 0
  
  # Extract baseline covariates and outcomes
  test_baseline <- test_data %>%
    group_by(deid_enc_id) %>%
    slice_min(DaysSinceEntry, with_ties = FALSE) %>%
    ungroup()
  
  test_outcomes <- test_data %>%
    group_by(deid_enc_id) %>%
    summarise(
      total_los = max(DaysSinceEntry, na.rm = TRUE),
      days_severe = sum(state %in% c("S", "S1", "S2"), na.rm = TRUE),
      died = any(state == "D"),
      ever_severe = any(state %in% c("S", "S1", "S2")),
      final_state = last(state),
      initial_state = first(state),
      .groups = "drop"
    )
  
  # Generate predictions using metadata-aware approach
  all_predictions <- generate_enhanced_predictions_unified(
    test_outcomes, test_baseline, fitted_model, 
    metadata, prediction_times, calibration_covariates
  )
  
  # Calculate base performance metrics
  base_results <- calculate_base_performance_metrics(
    fold_num, nrow(train_data), nrow(test_data), model_converged,
    test_outcomes, all_predictions$base_predictions
  )
  
  return(base_results)
}

# New helper function to refit model using unified structure
refit_model_for_fold_unified <- function(train_data, train_crude_rates, metadata) {
  
  # Extract fitting parameters from metadata
  covariates <- metadata$input_covariates
  spline_vars <- metadata$input_spline_vars
  constraint <- metadata$input_constraint
  
  # Extract spline specifications
  spline_df <- if (!is.null(metadata$spline_config) && !is.null(metadata$spline_config$default_df)) {
    metadata$spline_config$default_df
  } else {
    3
  }
  
  spline_type <- if (!is.null(metadata$spline_config) && !is.null(metadata$spline_config$default_type)) {
    metadata$spline_config$default_type
  } else {
    "ns"
  }
  
  # Extract time-varying specifications
  time_varying <- if (!is.null(metadata$time_config)) metadata$time_config$type else NULL
  time_variable <- if (!is.null(metadata$time_config)) metadata$time_config$variable else "DaysSinceEntry"
  time_breakpoints <- if (!is.null(metadata$time_config)) metadata$time_config$breakpoints else NULL
  
  # Call the unified fit_msm_models function
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
    constraint = constraint %||% "transition_specific"
  )
}

# New helper function for metadata-aware predictions
generate_enhanced_predictions_unified <- function(test_outcomes, test_baseline, fitted_model, 
                                                  metadata, prediction_times, calibration_covariates) {
  
  # Determine which variables to extract from baseline data
  base_variables <- extract_base_variables_from_metadata(metadata)
  
  # Base predictions using metadata-aware approach
  base_predictions <- map_dfr(test_outcomes$deid_enc_id, function(patient_id) {
    baseline_covs <- test_baseline %>%
      filter(deid_enc_id == patient_id) %>%
      select(any_of(base_variables)) %>%
      as.list()
    
    initial_state <- test_outcomes$initial_state[test_outcomes$deid_enc_id == patient_id]
    max_time <- max(prediction_times)
    
    generate_patient_predictions_unified(patient_id, baseline_covs, initial_state, fitted_model, 
                                         metadata, max_time)
  })
  
  return(list(
    base_predictions = base_predictions,
    time_specific_predictions = NULL,  # Can be extended later
    covariate_specific_predictions = NULL  # Can be extended later
  ))
}

# Enhanced patient prediction function using metadata
generate_patient_predictions_unified <- function(patient_id, baseline_covs, initial_state, 
                                                 fitted_model, metadata, prediction_time) {
  
  state_names <- rownames(fitted_model$qmodel$qmatrix)
  state_mapping <- state_name_to_num(initial_state, fitted_model)
  initial_state_num <- state_mapping[initial_state]
  if (is.na(initial_state_num)) {
    warning(paste("Initial state", initial_state, "not found in model, using state 1"))
    initial_state_num <- 1
  }
  
  tryCatch({
    # Generate predictions based on whether model has covariates
    has_covariates <- extract_has_splines_from_metadata(metadata) || 
      extract_has_time_varying_from_metadata(metadata) ||
      length(extract_base_variables_from_metadata(metadata)) > 0
    
    pmat <- if (has_covariates && length(baseline_covs) > 0) {
      pmatrix.msm(fitted_model, t = prediction_time, covariates = baseline_covs)
    } else {
      pmatrix.msm(fitted_model, t = prediction_time)
    }
    
    probs <- pmat[initial_state_num, ]
    
    # Enhanced LOS prediction using metadata about model type
    sojourn_result <- if (has_covariates && length(baseline_covs) > 0) {
      sojourn.msm(fitted_model, covariates = baseline_covs)
    } else {
      sojourn.msm(fitted_model)
    }
    
    # LOS calculations (keeping existing complex logic but enhanced)
    transient_states <- state_names[!state_names %in% c("D", "R")]
    Q_matrix <- fitted_model$qmodel$qmatrix
    transient_indices <- which(rownames(Q_matrix) %in% transient_states)
    
    if (length(transient_indices) > 1) {
      Q_transient <- Q_matrix[transient_indices, transient_indices]
      I_matrix <- diag(nrow(Q_transient))
      
      tryCatch({
        N_matrix <- solve(I_matrix - Q_transient)
        initial_transient_idx <- which(transient_states == initial_state)
        if (length(initial_transient_idx) == 0) initial_transient_idx <- 1
        
        expected_visits <- N_matrix[initial_transient_idx, ]
        pred_total_los_vis_soj <- sum(expected_visits * sojourn_result$estimates[transient_indices])
        pred_total_los_soj <- sum(sojourn_result$estimates[transient_indices])
        pred_total_los <- NA_real_
        
      }, error = function(e) {
        pred_total_los <- sum(sojourn_result$estimates[transient_indices])
        pred_total_los_vis_soj <- NA_real_
        pred_total_los_soj <- NA_real_
      })
    } else {
      pred_total_los <- sojourn_result$estimates[transient_indices[1]] %||% NA_real_
      pred_total_los_soj <- pred_total_los
      pred_total_los_vis_soj <- pred_total_los
    }
    
    # Severe state predictions
    severe_states <- state_names[grepl("S", state_names)]
    if (length(severe_states) > 0) {
      severe_indices <- which(transient_states %in% severe_states)
      if (exists("expected_visits") && length(severe_indices) > 0) {
        pred_days_severe_vis_soj <- sum(expected_visits[severe_indices] * 
                                          sojourn_result$estimates[severe_indices])
        pred_days_severe_soj <- sum(probs[severe_indices]) * 
          mean(sojourn_result$estimates[severe_indices], na.rm = TRUE)
        pred_days_severe <- NA_real_  
      } else {
        pred_days_severe_soj <- NA_real_
        pred_days_severe_vis_soj <- NA_real_
        pred_days_severe <- sum(probs[severe_indices]) * 
          mean(sojourn_result$estimates[severe_indices], na.rm = TRUE)
      }
    } else {
      pred_days_severe <- 0
      pred_days_severe_soj <- 0
      pred_days_severe_vis_soj <- 0
    }
    
    data.frame(
      deid_enc_id = patient_id,
      pred_prob_death = probs["D"] %||% 0,
      pred_prob_recovery = probs["R"] %||% 0,
      pred_prob_severe = sum(probs[severe_states]) %||% 0,
      pred_total_los = pred_total_los,
      pred_total_los_soj = pred_total_los_soj,
      pred_total_los_vis_soj = pred_total_los_vis_soj,
      pred_days_severe = pred_days_severe,
      pred_days_severe_soj = pred_days_severe_soj,
      pred_days_severe_vis_soj = pred_days_severe_vis_soj
    )
    
  }, error = function(e) {
    data.frame(
      deid_enc_id = patient_id,
      pred_prob_death = NA,
      pred_prob_recovery = NA, 
      pred_prob_severe = NA,
      pred_total_los = NA,
      pred_total_los_soj = NA,
      pred_total_los_vis_soj = NA,
      pred_days_severe = NA,
      pred_days_severe_soj = NA,
      pred_days_severe_vis_soj = NA
    )
  })
}

# Model fit: ------------------------------------------------------------

#' Calculate transition residuals using pre-computed transition probabilities
#' @param fitted_msm_models Nested list of fitted models: fitted_models[[model_structure]][[formula]]
#' @param patient_data Patient data used to fit the models
#' @param pmats_results Output from tidy_msm_pmats with pre-computed transition probabilities
#' @param residual_type Type of residuals ("deviance", "pearson", "raw")
#' @param debug Logical, whether to print debug information
#' @return Data frame with transition residuals for all model/formula combinations
calculate_transition_residuals <- function(fitted_msm_models, patient_data, 
                                           pmats_results,
                                           residual_type = "deviance",
                                           debug = FALSE) {
  
  validate_model_structure(fitted_msm_models)
  
  # Check if pmats_results is empty
  if (is.null(pmats_results) || nrow(pmats_results) == 0) {
    warning("pmats_results is empty - cannot calculate transition residuals")
    return(data.frame())
  }
  
  if (debug) cat("=== DEBUG: Starting transition residuals calculation ===\n")
  
  all_residuals <- list()
  
  # Process one model at a time to minimize memory usage
  for (model_structure in names(fitted_msm_models)) {
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      
      if (debug) cat("Processing:", model_structure, "with formula:", formula_name, "\n")
      
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      
      # Check if this is a failed model
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
        if (debug) cat("Skipping failed model:", model_structure, formula_name, "\n")
        next
      }
      
      # Extract the fitted model object
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in calculate_transition_residuals")
      if (is.null(fitted_model)) next
      
      if (!inherits(fitted_model, "msm")) {
        if (debug) cat("Invalid model object for:", model_structure, formula_name, "\n")
        next
      }
      
      # Filter patient data for this model structure - do this early to reduce memory
      model_data <- patient_data[patient_data$model == model_structure, ]
      
      if (debug) {
        cat("Model data dimensions:", nrow(model_data), "x", ncol(model_data), "\n")
        cat("Unique patients:", length(unique(model_data$deid_enc_id)), "\n")
        cat("Model convergence:", fitted_model$opt$convergence == 0, "\n")
        cat("Model has covariates:", !is.null(fitted_model$covariates), "\n")
      }
      
      # Get observed transitions - memory efficient approach
      observed_transitions <- model_data %>%
        arrange(deid_enc_id, DaysSinceEntry) %>%
        group_by(deid_enc_id) %>%
        mutate(
          from_state = state,
          to_state = lead(state),
          time_diff = lead(DaysSinceEntry) - DaysSinceEntry,
          transition_time = DaysSinceEntry,
          calendar_time = if("CalendarTime" %in% names(.)) CalendarTime else NA
        ) %>%
        filter(!is.na(to_state), time_diff > 0) %>%
        ungroup()
      
      # Standardize the state columns after observed_transitions is created
      observed_transitions <- standardize_state_columns(observed_transitions, fitted_model)
      
      if (debug) {
        cat("Observed transitions created. Dimensions:", nrow(observed_transitions), "x", ncol(observed_transitions), "\n")
        cat("Unique transitions:", length(unique(paste(observed_transitions$from_state, "->", observed_transitions$to_state))), "\n")
        cat("Unique time_diffs:", paste(unique(observed_transitions$time_diff), collapse = ", "), "\n")
      }
      
      if (nrow(observed_transitions) == 0) {
        if (debug) cat("No transitions found for:", model_structure, formula_name, "\n")
        next
      }
      
      # Filter pmats for this model/formula combination
      pmats_filtered <- pmats_results %>%
        filter(model == model_structure, formula == formula_name)
      
      if (nrow(pmats_filtered) == 0) {
        if (debug) cat("No pmats results found for:", model_structure, formula_name, "\n")
        next
      }
      
      if (debug) {
        cat("Pmats results filtered. Dimensions:", nrow(pmats_filtered), "x", ncol(pmats_filtered), "\n")
        cat("Available t_values:", paste(unique(pmats_filtered$t_value), collapse = ", "), "\n")
      }
      
      # Join observed transitions with expected probabilities
      # Match on: from_state, to_state, time_diff (= t_value), model, formula
      transition_residuals <- observed_transitions %>%
        left_join(
          pmats_filtered %>%
            select(statename, tostatename, estimate, t_value, model, formula) %>%
            rename(from_state = statename, 
                   to_state = tostatename, 
                   expected_prob = estimate, 
                   time_diff = t_value),
          by = c("from_state", "to_state", "time_diff", "model")
        ) %>%
        # Filter out transitions where we don't have matching pmats
        filter(!is.na(expected_prob), expected_prob > 0) %>%
        mutate(
          # Calculate residuals
          raw_residual = 1 - expected_prob,
          pearson_residual = raw_residual / sqrt(expected_prob * (1 - expected_prob)),
          deviance_residual = sign(raw_residual) * sqrt(pmax(-2 * log(pmax(expected_prob, 1e-10)), 0)),
          transition = paste(from_state, "->", to_state),
          formula = formula_name,
          covariate = if (formula_name == "~ 1") {
            "None"
          } else {
            tryCatch({
              paste(all.vars(as.formula(formula_name)), collapse = ", ")
            }, error = function(e) "Error parsing formula")
          }
        )
      
      if (nrow(transition_residuals) == 0) {
        if (debug) cat("No valid transitions after filtering for:", model_structure, formula_name, "\n")
        next
      }
      
      # Select residual type and standardize
      residual_col <- switch(residual_type,
                             "deviance" = "deviance_residual",
                             "pearson" = "pearson_residual", 
                             "raw" = "raw_residual")
      
      transition_residuals <- transition_residuals %>%
        mutate(
          residual = .data[[residual_col]],
          residual_std = residual / sd(residual, na.rm = TRUE),
          residual_magnitude = abs(residual_std),
          residual_category = case_when(
            residual_magnitude > 2.5 ~ "Very Large",
            residual_magnitude > 2 ~ "Large", 
            residual_magnitude > 1.5 ~ "Moderate",
            TRUE ~ "Small"
          )
        )
      
      if (debug) {
        cat("Calculated residuals for:", model_structure, formula_name, "\n")
        cat("Final dimensions:", nrow(transition_residuals), "x", ncol(transition_residuals), "\n")
        cat("Transitions without matching pmats:", 
            nrow(observed_transitions) - nrow(transition_residuals), "\n")
      }
      
      # Store results
      combo_key <- paste(model_structure, formula_name, sep = "___")
      all_residuals[[combo_key]] <- transition_residuals
      
      # Clean up to free memory
      rm(model_data, observed_transitions, pmats_filtered, transition_residuals)
      gc(verbose = FALSE)
    }
  }
  
  # Combine all results efficiently
  if (length(all_residuals) == 0) {
    warning("No valid residuals calculated for any model/formula combination")
    return(data.frame())
  }
  
  final_residuals <- do.call(bind_rows, all_residuals)
  
  if (debug) {
    cat("=== DEBUG: Final Combined Results ===\n")
    cat("Total transitions across all models:", nrow(final_residuals), "\n")
    if (nrow(final_residuals) > 0) {
      cat("Model structures:", paste(unique(final_residuals$model), collapse = ", "), "\n")
      cat("Formulas:", paste(unique(final_residuals$formula), collapse = ", "), "\n")
      cat("Time horizons represented:", paste(sort(unique(final_residuals$time_diff)), collapse = ", "), "\n")
    }
  }
  
  return(final_residuals)
}

# Compare models --------------------------------------------

#' Comprehensive model comparison combining within and across structure comparisons
#' @param fitted_models Nested list of fitted MSM models
#' @param include_within Include within-structure comparisons (default: TRUE)
#' @param include_across Include cross-structure comparisons (default: TRUE)  
#' @param cross_structure_methods Methods for cross-structure comparison
#' @param cores Number of cores for DRLCV
#' @return List with within_structure and across_structure comparison results
comprehensive_model_comparison <- function(fitted_models, include_within = TRUE, 
                                           include_across = FALSE,
                                           cross_structure_methods = c("draic"),
                                           mc.cores = min(6, detectCores() - 1)) {
  
  results <- list()
  
  if (include_within) {
    cat("=== Within-Structure Comparisons ===\n")
    results$within_structure <- compare_within_structure(fitted_models)
    
    # Summary of within-structure results
    within_summary <- results$within_structure %>%
      group_by(model) %>%
      summarise(
        n_models = n(),
        n_converged = sum(status == "converged", na.rm = TRUE),
        best_aic_formula = formula[which.min(AIC)],
        best_aic_value = min(AIC, na.rm = TRUE),
        .groups = "drop"
      )
    
    cat("Within-structure summary:\n")
    print(within_summary)
    
    best_models <- results$within_structure %>%
      filter(status == "converged") %>%
      group_by(model) %>%
      slice_min(AIC, n = 1) %>%
      ungroup() %>%
      arrange(AIC)
    
    cat("\n=== Best Model per Structure (by AIC) ===\n")
    print(best_models %>% select(model, formula, covariate, AIC, BIC, n_params))
    
  }
  
  if (include_across) {
    cat("\n=== Cross-Structure Comparisons ===\n")
    results$across_structure <- compare_across_structures(
      fitted_models, 
      methods = cross_structure_methods, 
      cores = mc.cores
    )
    
    # Summary of cross-structure results
    if ("draic" %in% cross_structure_methods) {
      draic_summary <- results$across_structure %>%
        filter(!is.na(draic)) %>%
        arrange(draic) %>%
        head(5)
      
      cat("Top 5 model pairs by DRAIC (lower is better for model1):\n")
      print(draic_summary %>% select(model1_name, model2_name, draic, draic_pval))
    }
  }
  
  return(results)
}

#' Compare models within the same structure using AIC, BIC, and LRT with nesting direction detection
#' @param fitted_models Nested list of fitted MSM models
#' @param model_structure Which model structure to compare (if NULL, compares all)
#' @param reference_formula Reference formula for LRT (default: "~ 1")
#' @param use_tidy_models Whether to use tidy_msm_models for extraction (default: TRUE)
#' @return Data frame with within-structure comparison metrics
compare_within_structure <- function(fitted_models, model_structure = NULL, 
                                     reference_formula = "~ 1", use_tidy_models = TRUE) {
  
  validate_model_structure(fitted_models)
  
  # Detect nesting direction
  nesting_info <- detect_nesting_direction(fitted_models)
  
  if (nesting_info$direction == "unknown") {
    warning("Unable to determine nesting direction. No valid models found.")
    return(data.frame())
  }
  
  cat("Detected nesting direction:", nesting_info$direction, "\n")
  cat("Structure levels:", paste(nesting_info$structure_levels, collapse = ", "), "\n")
  cat("Formula levels found:", paste(head(nesting_info$formula_levels, 3), collapse = ", "), 
      if(length(nesting_info$formula_levels) > 3) "..." else "", "\n")
  
  # If no specific structure provided, do this for all structures
  if (is.null(model_structure)) {
    all_comparisons <- map_dfr(nesting_info$structure_levels, function(struct) {
      compare_within_structure(fitted_models, struct, reference_formula, use_tidy_models)
    })
    return(all_comparisons)
  }
  
  # Extract models for the specified structure based on nesting direction
  if (nesting_info$direction == "structure_formula") {
    # fitted_models[[structure]][[formula]]
    if (!model_structure %in% names(fitted_models)) {
      stop(paste("Model structure", model_structure, "not found in fitted_models"))
    }
    structure_models <- fitted_models[[model_structure]]
    
  } else if (nesting_info$direction == "formula_structure") {
    # fitted_models[[formula]][[structure]] - need to reorganize
    structure_models <- list()
    for (formula_name in names(fitted_models)) {
      if (model_structure %in% names(fitted_models[[formula_name]])) {
        structure_models[[formula_name]] <- fitted_models[[formula_name]][[model_structure]]
      }
    }
    
    if (length(structure_models) == 0) {
      warning(paste("No models found for structure", model_structure))
      return(data.frame())
    }
  }
  
  if (use_tidy_models) {
    # Create properly nested structure for tidy_msm_models
    temp_nested <- if (nesting_info$direction == "structure_formula") {
      setNames(list(structure_models), model_structure)
    } else {
      # Reorganize for tidy_msm_models which expects structure_formula nesting
      temp_structure <- list()
      temp_structure[[model_structure]] <- structure_models
      temp_structure
    }
    
    tidied <- tidy_msm_models(temp_nested) %>%
      filter(model == model_structure, status == "converged")
    
    if (nrow(tidied) == 0) {
      warning(paste("No converged models found for structure", model_structure))
      return(data.frame())
    }
    
    # Get reference model stats
    ref_stats <- tidied %>% filter(formula == reference_formula)
    if (nrow(ref_stats) == 0) {
      warning(paste("Reference formula", reference_formula, "not found. Using first model as reference."))
      ref_stats <- tidied[1, ]
      reference_formula <- ref_stats$formula
    }
    
    # Calculate deltas from reference
    comparison_results <- tidied %>%
      mutate(
        delta_AIC = AIC - ref_stats$AIC,
        delta_BIC = BIC - ref_stats$BIC,
        delta_loglik = loglik - ref_stats$loglik,
        aic_rank = rank(AIC),
        bic_rank = rank(BIC)
      )
    
    # Add LRT results for non-reference models
    lrt_results <- map_dfr(setdiff(tidied$formula, reference_formula), function(formula) {
      
      # Extract model objects based on nesting direction
      if (nesting_info$direction == "structure_formula") {
        ref_model_entry <- structure_models[[reference_formula]]
        test_model_entry <- structure_models[[formula]]
      } else {
        ref_model_entry <- structure_models[[reference_formula]]
        test_model_entry <- structure_models[[formula]]
      }
      
      ref_model <- if (is.list(ref_model_entry) && !is.null(ref_model_entry$fitted_model)) {
        ref_model_entry$fitted_model
      } else {
        ref_model_entry
      }
      
      test_model <- if (is.list(test_model_entry) && !is.null(test_model_entry$fitted_model)) {
        test_model_entry$fitted_model
      } else {
        test_model_entry
      }
      
      # Replace the LRT calculation section in compare_within_structure function:
      
      lrt_result <- tryCatch({
        if (!is.null(ref_model) && !is.null(test_model) && 
            inherits(ref_model, "msm") && inherits(test_model, "msm")) {
          
          lrt <- lrtest.msm(ref_model, test_model)
          
          # Handle different return structures from lrtest.msm
          if (is.matrix(lrt)) {
            # Matrix format: extract from first row
            c(lrt[1, 1], lrt[1, 2], lrt[1, 3])
          } else if (is.list(lrt)) {
            # List format: extract named elements
            c(lrt$statistic, lrt$df, lrt$p.value)
          } else if (is.vector(lrt) && length(lrt) >= 3) {
            # Vector format: use first three elements
            c(lrt[1], lrt[2], lrt[3])
          } else {
            # Unknown format - return NAs
            c(NA, NA, NA)
          }
        } else {
          c(NA, NA, NA)
        }
      }, error = function(e) {
        warning(paste("LRT failed for", formula, "vs", reference_formula, ":", e$message))
        c(NA, NA, NA)
      })      
      data.frame(
        formula = formula,
        lrt_stat = lrt_result[1],
        lrt_df = lrt_result[2],
        lrt_pval = lrt_result[3]
      )
    })
    
    # Add reference model LRT (all zeros/NAs)
    ref_lrt <- data.frame(
      formula = reference_formula,
      lrt_stat = 0,
      lrt_df = 0,
      lrt_pval = 1
    )
    
    lrt_results <- bind_rows(lrt_results, ref_lrt)
    
    # Merge results
    final_results <- comparison_results %>%
      left_join(lrt_results, by = "formula") %>%
      mutate(reference_model = reference_formula) %>%
      arrange(AIC)
    
  } else {
    # Manual extraction (fallback)
    warning("Manual extraction not fully implemented - using tidy_msm_models")
    return(compare_within_structure(fitted_models, model_structure, reference_formula, use_tidy_models = TRUE))
  }
  
  return(final_results)
}

#' Compare models across different structures using DRAIC and DRLCV
#' @param fitted_models Nested list of fitted MSM models
#' @param model_pairs Data frame with columns model1_struct, formula1, model2_struct, formula2
#' @param methods Vector of comparison methods ("draic", "drlcv", or both)
#' @param cores Number of cores for parallel computation
#' @return Data frame with cross-structure comparison metrics
compare_across_structures <- function(fitted_models, model_pairs = NULL, 
                                      methods = c("draic", "drlcv"), mc.cores = min(6, detectCores() - 1)) {
  
  validate_model_structure(fitted_models)
  
  # If no model pairs specified, create all pairwise comparisons
  if (is.null(model_pairs)) {
    # Use tidy_msm_models to get all converged models
    tidied <- tidy_msm_models(fitted_models) %>%
      filter(status == "converged")
    
    if (nrow(tidied) < 2) {
      stop("Need at least 2 converged models for comparison")
    }
    
    # Create all pairwise combinations
    model_pairs <- expand_grid(
      model1_struct = tidied$model,
      formula1 = tidied$formula,
      model2_struct = tidied$model,
      formula2 = tidied$formula
    ) %>%
      filter(!(model1_struct == model2_struct & formula1 == formula2)) %>%  # Remove self-comparisons
      filter(model1_struct != model2_struct) %>%  # Only cross-structure comparisons
      distinct()
  }
  
  # Validate methods
  valid_methods <- c("draic", "drlcv")
  methods <- intersect(methods, valid_methods)
  if (length(methods) == 0) {
    stop("At least one valid method must be specified: 'draic' or 'drlcv'")
  }
  
  # Safe wrappers
  safe_draic <- safely(draic.msm)
  safe_drlcv <- safely(drlcv.msm)
  
  cat("Comparing", nrow(model_pairs), "model pairs across structures\n")
  
  comparison_results <- model_pairs %>%
    mutate(
      
      model1_obj = map2(model1_struct, formula1, function(struct, form) {
        extract_fitted_model(fitted_models, struct, form, "in compare_across_structures model1")
      }),
      
      model2_obj = map2(model2_struct, formula1, function(struct, form) {
        extract_fitted_model(fitted_models, struct, form, "in compare_across_structures model2")
      })
    )
  
  # DRAIC comparisons
  if ("draic" %in% methods) {
    cat("Computing DRAIC comparisons...\n")
    comparison_results <- comparison_results %>%
      mutate(
        draic_result = map2(model1_obj, model2_obj, function(m1, m2) {
          if (is.null(m1) || is.null(m2)) return(NULL)
          
          result <- safe_draic(m1, m2)
          if (!is.null(result)) {
            list(
              draic = result$result$draic[1],
              draic_ll = result$result$ti["2.5%"],
              draic_ul = result$result$ti["97.5%"],
              draic_pval = result$result$ti["Prob<0"]
            )
          } else {
            warning(paste("DRAIC failed:", result$error$message))
            NULL
          }
        }),
        
        draic      = map_dbl(draic_result, pluck, "draic", .default = NA_real_),
        draic_ll   = map_dbl(draic_result, pluck, "draic_ll", .default = NA_real_),
        draic_ul   = map_dbl(draic_result, pluck, "draic_ul", .default = NA_real_),
        draic_pval = map_dbl(draic_result, pluck, "draic_pval", .default = NA_real_)
      ) %>%
      select(-draic_result)
  }
  
  # DRLCV comparisons
  if ("drlcv" %in% methods) {
    cat("Computing DRLCV comparisons (this may take a while)...\n")
    comparison_results <- comparison_results %>%
      mutate(
        drlcv_result = map2(model1_obj, model2_obj, function(m1, m2) {
          if (is.null(m1) || is.null(m2)) return(NULL)
          
          result <- safe_drlcv(m1, m2, cores = mc.cores, verbose = FALSE)
          if (!is.null(result$result)) {
            list(
              drlcv = result$result$drlcv[1],
              drlcv_ll = result$result$ti["2.5%"],
              drlcv_ul = result$result$ti["97.5%"],
              drlcv_pval = result$result$ti["Prob<0"]
            )
          } else {
            warning(paste("DRLCV failed:", result$error$message))
            NULL
          }
        }),
        
        drlcv      = map_dbl(drlcv_result, pluck, "drlcv",     .default = NA_real_),
        drlcv_ll   = map_dbl(drlcv_result, pluck, "drlcv_ll",  .default = NA_real_),
        drlcv_ul   = map_dbl(drlcv_result, pluck, "drlcv_ul",  .default = NA_real_),
        drlcv_pval = map_dbl(drlcv_result, pluck, "drlcv_pval",.default = NA_real_)
      ) %>%
      select(-drlcv_result)
  }
  
  # Clean up and add interpretive columns
  final_results <- comparison_results %>%
    select(-model1_obj, -model2_obj) %>%
    mutate(
      model1_name = paste(model1_struct, formula1, sep = " | "),
      model2_name = paste(model2_struct, formula2, sep = " | "),
      comparison_type = "cross_structure"
    )
  
  return(final_results)
}

#' Detect the nesting direction of fitted models by examining structure
#' @param fitted_models Nested list structure to analyze
#' @return List with direction, structure_levels, and formula_levels
detect_nesting_direction <- function(fitted_models) {
  
  if (!is.list(fitted_models) || length(fitted_models) == 0) {
    return(list(
      direction = "unknown", 
      structure_levels = character(), 
      formula_levels = character(),
      confidence = 0
    ))
  }
  
  # Get the first level names
  first_level_names <- names(fitted_models)
  if (is.null(first_level_names) || length(first_level_names) == 0) {
    return(list(
      direction = "unknown", 
      structure_levels = character(), 
      formula_levels = character(),
      confidence = 0
    ))
  }
  
  # Examine several elements to get a robust assessment
  sample_size <- min(3, length(fitted_models))
  sample_elements <- fitted_models[1:sample_size]
  
  structure_evidence <- list()
  
  for (i in seq_along(sample_elements)) {
    element_name <- names(sample_elements)[i]
    element <- sample_elements[[i]]
    
    if (!is.list(element)) {
      structure_evidence[[i]] <- list(
        first_level = element_name,
        is_list = FALSE,
        direction = "unknown"
      )
      next
    }
    
    second_level_names <- names(element)
    if (is.null(second_level_names) || length(second_level_names) == 0) {
      structure_evidence[[i]] <- list(
        first_level = element_name,
        second_level = character(),
        is_nested = FALSE,
        direction = "unknown"
      )
      next
    }
    
    # Examine the deepest accessible level
    deep_samples <- list()
    deep_sample_count <- min(2, length(element))
    
    for (j in 1:deep_sample_count) {
      deep_element <- element[[j]]
      deep_name <- names(element)[j]
      
      deep_samples[[j]] <- list(
        name = deep_name,
        is_msm = inherits(deep_element, "msm"),
        has_fitted_model = is.list(deep_element) && "fitted_model" %in% names(deep_element),
        has_status = is.list(deep_element) && "status" %in% names(deep_element),
        is_list = is.list(deep_element),
        class = class(deep_element)[1]
      )
    }
    
    # Determine what this element suggests about nesting direction
    has_msm_objects <- any(sapply(deep_samples, function(x) x$is_msm))
    has_fitted_model_structure <- any(sapply(deep_samples, function(x) x$has_fitted_model))
    
    # Check if second level names look like formulas
    second_level_is_formula <- any(grepl("^~", second_level_names)) ||
      any(second_level_names %in% c("~ 1", "~1"))
    
    # Check if first level name looks like a model structure identifier
    first_level_looks_structural <- grepl("model|mod_|base|sev|trans|hx_", element_name, ignore.case = TRUE)
    
    structure_evidence[[i]] <- list(
      first_level = element_name,
      second_level = second_level_names,
      has_msm_objects = has_msm_objects,
      has_fitted_model_structure = has_fitted_model_structure,
      second_level_is_formula = second_level_is_formula,
      first_level_looks_structural = first_level_looks_structural,
      deep_samples = deep_samples
    )
  }
  
  # Aggregate evidence to determine direction
  votes_structure_formula <- 0
  votes_formula_structure <- 0
  confidence_score <- 0
  
  for (evidence in structure_evidence) {
    # Vote based on multiple criteria
    
    # Strong evidence for structure -> formula nesting
    if (evidence$has_fitted_model_structure && evidence$second_level_is_formula) {
      votes_structure_formula <- votes_structure_formula + 3
    }
    
    if (evidence$first_level_looks_structural && evidence$second_level_is_formula) {
      votes_structure_formula <- votes_structure_formula + 2  
    }
    
    if (evidence$has_msm_objects) {
      votes_structure_formula <- votes_structure_formula + 1
    }
    
    # Evidence for formula -> structure nesting would need different checks
    # (This is less common in the MSM context, but we can add checks if needed)
    
    # Add to confidence based on quality of evidence
    if (evidence$has_fitted_model_structure) confidence_score <- confidence_score + 0.4
    if (evidence$has_msm_objects) confidence_score <- confidence_score + 0.3
    if (evidence$second_level_is_formula) confidence_score <- confidence_score + 0.2
    if (evidence$first_level_looks_structural) confidence_score <- confidence_score + 0.1
  }
  
  # Determine final direction
  if (votes_structure_formula > votes_formula_structure) {
    direction <- "structure_formula"
    structure_levels <- first_level_names
    formula_levels <- unique(unlist(lapply(fitted_models, names)))
  } else if (votes_formula_structure > votes_structure_formula) {
    direction <- "formula_structure"
    formula_levels <- first_level_names  
    structure_levels <- unique(unlist(lapply(fitted_models, names)))
  } else {
    # Default assumption based on MSM patterns
    direction <- "structure_formula"
    structure_levels <- first_level_names
    formula_levels <- unique(unlist(lapply(fitted_models, names)))
    confidence_score <- confidence_score * 0.5  # Lower confidence for default
  }
  
  return(list(
    direction = direction,
    structure_levels = structure_levels,
    formula_levels = formula_levels,
    confidence = min(confidence_score, 1.0),
    evidence_summary = list(
      votes_structure_formula = votes_structure_formula,
      votes_formula_structure = votes_formula_structure,
      sample_count = length(structure_evidence)
    ),
    debug_info = structure_evidence
  ))
}


# Comprehensive wrapper for base performance metrics calculation ----

run_comprehensive_msm_analysis <- function(fitted_msm_models, 
                                           patient_data, 
                                           crude_rates,
                                           analysis_config = list(),
                                           parallel = TRUE,
                                           mc.cores = NULL,
                                           verbose = TRUE) {
  
  # Validate inputs
  if (!is.list(fitted_msm_models) || length(fitted_msm_models) == 0) {
    stop("fitted_msm_models must be a non-empty list")
  }
  
  validate_model_structure(fitted_msm_models)
  
  if (!is.data.frame(patient_data)) {
    stop("patient_data must be a data frame")
  }
  
  if (!is.list(crude_rates)) {
    stop("crude_rates must be a list")
  }
  
  # Setup parallel processing
  if (is.null(mc.cores)) {
    mc.cores <- min(8, parallel::detectCores() - 1)
  }
  
  # Default configuration
  default_config <- list(
    # Model summary parameters
    tidy_models = list(),
    
    # Q-matrix parameters
    qmatrix = list(
      covariates_list = NULL,
      mc.cores = mc.cores
    ),
    
    # P-matrix parameters  
    pmats = list(
      t_values = time_vec,
      covariates_list = NULL,
      mc.cores = mc.cores
    ),
    
    # Sojourn time parameters
    sojourns = list(
      covariates_list = NULL
    ),
    
    # Prevalence parameters
    prevalence = list(
      time_points = time_vec,
      covariates_list = NULL,
      ci = TRUE,
      ci_method = "normal",
      use_approximation = FALSE,
      mc.cores = mc.cores
    ),
    
    # Hazard ratio parameters
    hazards = list(
      hazard_scale = 1
    ),
    
    # Cross-validation parameters
    cv = list(
      k_folds = 3,
      prediction_times = time_vec,
      stratify_by = "final_state",
      calibration_covariates = NULL,
      calibration_subgroups = NULL,
      parallel = parallel,
      n_cores = mc.cores
    ),
    
    # Residual analysis parameters
    residuals = list(
      residual_type = "deviance",
      debug = verbose
    ),
    
    # Model comparison parameters
    comparison = list(
      include_within = TRUE,
      include_across = TRUE,
      cross_structure_methods = c("draic", "drlcv"),
      mc.cores = mc.cores
    )
  )
  
  # Merge user config with defaults
  config <- modifyList(default_config, analysis_config)
  
  # Initialize results list
  results <- list(
    metadata = list(
      run_time = Sys.time(),
      n_structures = length(fitted_msm_models),
      total_models = sum(sapply(fitted_msm_models, length)),
      config = config,
      parallel = parallel,
      mc.cores = mc.cores
    )
  )
  
  if (verbose) {
    cat("=== COMPREHENSIVE MSM ANALYSIS ===\n")
    cat("Structures:", length(fitted_msm_models), "\n")
    cat("Total models:", sum(sapply(fitted_msm_models, length)), "\n")
    cat("Parallel processing:", ifelse(parallel, paste("enabled (", mc.cores, "cores)"), "disabled"), "\n")
    cat("\n")
  }
  
  # 1. Model Summary
  if (verbose) cat("1. Tidying model summaries...\n")
  results$model_summary <- tryCatch({
    tidy_msm_models(fitted_msm_models)
  }, error = function(e) {
    warning(paste("tidy_msm_models failed:", e$message))
    data.frame()
  })
  
  if (verbose && nrow(results$model_summary) > 0) {
    converged <- sum(results$model_summary$status == "converged", na.rm = TRUE)
    cat("   -", converged, "of", nrow(results$model_summary), "models converged\n")
  }
  
  # 2. Transition Intensities (Q-matrix)
  if (verbose) cat("2. Extracting transition intensities...\n")
  results$qmatrix <- tryCatch({
    tidy_msm_qmatrix(fitted_msm_models, 
                     covariates_list = config$qmatrix$covariates_list,
                     ci_type = "normal",
                     mc.cores = config$qmatrix$mc.cores)
  }, error = function(e) {
    warning(paste("tidy_msm_qmatrix failed:", e$message))
    data.frame()
  })
  
  # 3. Transition Probabilities (P-matrix)
  if (verbose) cat("3. Calculating transition probabilities...\n")
  results$pmats <- tryCatch({
    tidy_msm_pmats(fitted_msm_models,
                   t_values = config$pmats$t_values,
                   covariates_list = config$pmats$covariates_list,
                   ci_type = "normal",
                   mc.cores = config$pmats$mc.cores)
  }, error = function(e) {
    warning(paste("tidy_msm_pmats failed:", e$message))
    data.frame()
  })
  
  # 4. Sojourn Times
  if (verbose) cat("4. Extracting sojourn times...\n")
  results$sojourns <- tryCatch({
    tidy_msm_sojourns(fitted_msm_models, 
                      covariates_list = config$sojourns$covariates_list)
  }, error = function(e) {
    warning(paste("tidy_msm_sojourns failed:", e$message))
    data.frame()
  })
  
  # 5. State Prevalences (most computationally intensive)
  if (!isTRUE(config$prevalence$skip)) {
    if (verbose) cat("5. Computing state prevalences...\n")
    prevalence_start_time <- Sys.time()
    results$prevalence <- tryCatch({
      do.call(tidy_msm_prevalences, c(list(fitted_msm_models), config$prevalence))
    }, error = function(e) {
      warning(paste("tidy_msm_prevalences failed:", e$message))
      data.frame()
    })
    prevalence_duration <- Sys.time() - prevalence_start_time
    if (verbose) cat("   - Completed in", round(as.numeric(prevalence_duration, units = "secs"), 1), "seconds\n")
  } else {
    if (verbose) cat("5. Skipping state prevalences (skip = TRUE)\n")
    results$prevalence <- data.frame()
  }
  
  # 6. Hazard Ratios (only for models with covariates)
  if (!isTRUE(config$hazards$skip)) {
    if (verbose) cat("6. Calculating hazard ratios...\n")
    results$hazards <- tryCatch({
      tidy_msm_hazards(fitted_msm_models, hazard_scale = config$hazards$hazard_scale)
    }, error = function(e) {
      warning(paste("tidy_msm_hazards failed:", e$message))
      data.frame()
    })
  } else {
    if (verbose) cat("6. Skipping hazard ratios (skip = TRUE)\n")
    results$hazards <- data.frame()
  }
  
  # 7. Cross-Validation and Predictive Performance
  if (!isTRUE(config$cv$skip)) {
    if (verbose) cat("7. Running cross-validation...\n")
    cv_start_time <- Sys.time()
    results$predictive_performance <- tryCatch({
      calculate_predictive_performance(
        patient_data = patient_data, 
        fitted_models = fitted_msm_models, 
        crude_rates = crude_rates,
        k_folds = config$cv$k_folds,
        prediction_times = config$cv$prediction_times,
        stratify_by = config$cv$stratify_by,
        calibration_covariates = config$cv$calibration_covariates,
        calibration_subgroups = config$cv$calibration_subgroups,
        parallel = config$cv$parallel,
        n_cores = config$cv$n_cores
      )
    }, error = function(e) {
      warning(paste("calculate_predictive_performance failed:", e$message))
      list()
    })
    cv_duration <- Sys.time() - cv_start_time
    if (verbose) cat("   - Completed in", round(as.numeric(cv_duration, units = "mins"), 1), "minutes\n")
  } else {
    if (verbose) cat("7. Skipping cross-validation (skip = TRUE)\n")
    results$predictive_performance <- list()
  }
  
  # 8. Transition Residuals
  if (!isTRUE(config$residuals$skip)) {
    if (verbose) cat("8. Computing transition residuals...\n")
    results$residuals <- tryCatch({
      do.call(calculate_transition_residuals, 
              c(list(fitted_msm_models = fitted_msm_models,
                     patient_data = patient_data,
                     pmats_results = results$pmats),  # Pass pre-computed pmats
                config$residuals))
    }, error = function(e) {
      warning(paste("calculate_transition_residuals failed:", e$message))
      data.frame()
    })
  } else {
    if (verbose) cat("8. Skipping transition residuals (skip = TRUE)\n")
    results$residuals <- data.frame()
  }
  
  # 9. Comprehensive Model Comparison
  if (verbose) cat("9. Running model comparisons...\n")
  comparison_start_time <- Sys.time()
  results$model_comparison <- tryCatch({
    comprehensive_model_comparison(
      fitted_models = fitted_msm_models,
      include_within = config$comparison$include_within,
      include_across = config$comparison$include_across,
      cross_structure_methods = config$comparison$cross_structure_methods,
      mc.cores = config$comparison$mc.cores
    )
  }, error = function(e) {
    warning(paste("comprehensive_model_comparison failed:", e$message))
    list()
  })
  comparison_duration <- Sys.time() - comparison_start_time
  if (verbose) cat("   - Completed in", round(as.numeric(comparison_duration, units = "secs"), 1), "seconds\n")
  
  # Final metadata
  total_duration <- Sys.time() - results$metadata$run_time
  results$metadata$total_runtime <- total_duration
  results$metadata$components_completed <- sum(!sapply(results[2:length(results)], function(x) {
    is.data.frame(x) && nrow(x) == 0 || is.list(x) && length(x) == 0
  }))
  
  if (verbose) {
    cat("\n=== ANALYSIS COMPLETE ===\n")
    cat("Total runtime:", round(as.numeric(total_duration, units = "mins"), 1), "minutes\n")
    cat("Components completed:", results$metadata$components_completed, "of 9\n")
    cat("Result object size:", format(object.size(results), units = "Mb"), "\n")
  }
  
  # Add helper attributes for easier access
  class(results) <- c("msm_analysis_results", "list")
  attr(results, "summary") <- create_analysis_summary(results)
  
  return(results)
}


## Print analysis ----------------------------------------------------------

#' Create a summary of analysis results
#' @param results Output from run_comprehensive_msm_analysis
#' @return Data frame summarizing key findings
create_analysis_summary <- function(results) {
  
  summary_list <- list()
  
  # Model convergence summary
  if (!is.null(results$model_summary) && nrow(results$model_summary) > 0) {
    model_stats <- results$model_summary %>%
      group_by(model) %>%
      summarise(
        total_models = n(),
        converged = sum(status == "converged", na.rm = TRUE),
        best_aic = min(AIC[status == "converged"], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(convergence_rate = converged / total_models)
    
    summary_list$model_convergence <- model_stats
  }
  
  # Cross-validation performance summary
  if (!is.null(results$predictive_performance) && length(results$predictive_performance) > 0) {
    cv_summary <- map_dfr(names(results$predictive_performance), function(model_key) {
      cv_data <- results$predictive_performance[[model_key]]
      if (!is.null(cv_data$base)) {
        cv_data$base %>%
          select(model, formula, death_auc_cv, calibration_slope_death_cv, n_folds_converged) %>%
          mutate(model_key = model_key)
      }
    })
    
    if (nrow(cv_summary) > 0) {
      summary_list$predictive_performance <- cv_summary %>%
        arrange(desc(death_auc_cv)) %>%
        head(10)  # Top 10 by AUC
    }
  }
  
  # Model comparison summary
  if (!is.null(results$model_comparison) && !is.null(results$model_comparison$within_structure)) {
    comparison_summary <- results$model_comparison$within_structure %>%
      filter(status == "converged") %>%
      group_by(model) %>%
      slice_min(AIC, n = 1) %>%
      ungroup() %>%
      arrange(AIC) %>%
      select(model, formula, covariate, AIC, BIC, n_params)
    
    summary_list$best_models <- comparison_summary
  }
  
  return(summary_list)
}

#' Print method for MSM analysis results
#' @param x MSM analysis results object
#' @param ... Additional arguments
print.msm_analysis_results <- function(x, ...) {
  cat("MSM Analysis Results\n")
  cat("====================\n")
  cat("Run time:", format(x$metadata$run_time), "\n")
  cat("Total runtime:", round(as.numeric(x$metadata$total_runtime, units = "mins"), 1), "minutes\n")
  cat("Components completed:", x$metadata$components_completed, "of 9\n\n")
  
  # Print component sizes
  component_sizes <- sapply(x[2:length(x)], function(comp) {
    if (is.data.frame(comp)) {
      paste0(nrow(comp), " rows")
    } else if (is.list(comp)) {
      paste0(length(comp), " elements")
    } else {
      "unknown"
    }
  })
  
  cat("Components:\n")
  for (i in seq_along(component_sizes)) {
    cat("  ", names(component_sizes)[i], ":", component_sizes[i], "\n")
  }
  
  # Print summary if available
  if (!is.null(attr(x, "summary"))) {
    cat("\nKey Findings:\n")
    summary_obj <- attr(x, "summary")
    
    if (!is.null(summary_obj$model_convergence)) {
      cat("- Model convergence rates:\n")
      print(summary_obj$model_convergence)
    }
    
    if (!is.null(summary_obj$best_models)) {
      cat("\n- Best models by AIC:\n")
      print(head(summary_obj$best_models, 3))
    }
  }
}

#' Extract specific analysis component with error handling
#' @param results MSM analysis results object
#' @param component Component name to extract
#' @return Requested component or informative error
extract_component <- function(results, component) {
  
  valid_components <- c("model_summary", "qmatrix", "pmats", "sojourns", 
                        "prevalence", "hazards", "predictive_performance", 
                        "residuals", "model_comparison")
  
  if (!component %in% valid_components) {
    stop(paste("Invalid component. Must be one of:", 
               paste(valid_components, collapse = ", ")))
  }
  
  if (!component %in% names(results)) {
    stop(paste("Component", component, "not found in results"))
  }
  
  comp_data <- results[[component]]
  
  # Add informative messages for empty results
  if (is.data.frame(comp_data) && nrow(comp_data) == 0) {
    message(paste("Component", component, "is empty - check for analysis errors"))
  } else if (is.list(comp_data) && length(comp_data) == 0) {
    message(paste("Component", component, "is empty - check for analysis errors"))
  }
  
  return(comp_data)
}


# Helper functions ------------------------------------------------------

calc_bic_msm <- function(fitted_model) {
  if (is.null(fitted_model)) {
    return(NA_real_)
  }
  
  # Check if minus2loglik is valid
  if (!is.finite(fitted_model$minus2loglik)) {
    return(NA_real_)
  }
  
  # Get number of parameters
  n_params <- length(fitted_model$estimates)
  if (n_params <= 0) {
    return(NA_real_)
  }
  
  # Get sample size - MSM models store this in different places
  n_obs <- length(unique(fitted_model[["data"]][["mf"]][["(subject)"]]))
  
  bic <- fitted_model$minus2loglik + n_params * log(n_obs)
  
  return(bic)
}

calculate_auc_safe <- function(outcome, prediction) {
  if (var(outcome, na.rm = TRUE) > 0 && sum(!is.na(prediction)) > 0) {
    tryCatch({
      pROC::auc(outcome, prediction, quiet = TRUE)
    }, error = function(e) NA)
  } else {
    NA
  }
}

extract_fitted_model <- function(fitted_models, structure, formula, context = "") {
  if (!structure %in% names(fitted_models)) {
    warning(paste("Structure", structure, "not found in fitted_models", context))
    return(NULL)
  }
  
  if (!formula %in% names(fitted_models[[structure]])) {
    warning(paste("Formula", formula, "not found in structure", structure, context))
    return(NULL)
  }
  
  model_entry <- fitted_models[[structure]][[formula]]
  
  # Handle the new unified structure
  if (!is.list(model_entry)) {
    warning(paste("Invalid model structure for", structure, formula, 
                  "- expected list with $fitted_model component", context))
    return(NULL)
  }
  
  if (!"fitted_model" %in% names(model_entry)) {
    warning(paste("Missing fitted_model component for", structure, formula, context))
    return(NULL)
  }
  
  # Check if model failed
  if (!is.null(model_entry$status) && model_entry$status == "failed") {
    warning(paste("Model failed for", structure, formula, context))
    return(NULL)
  }
  
  # Additional check that the fitted_model is actually an msm object
  fitted_model <- model_entry$fitted_model
  if (is.null(fitted_model)) {
    warning(paste("fitted_model is NULL for", structure, formula, context))
    return(NULL)
  }
  
  if (!inherits(fitted_model, "msm")) {
    warning(paste("fitted_model is not an msm object for", structure, formula, context))
    return(NULL)
  }
  
  return(fitted_model)
}

validate_model_structure <- function(fitted_models) {
  if (!is.list(fitted_models) || length(fitted_models) == 0) {
    stop("fitted_models must be a non-empty list")
  }
  
  for (structure in names(fitted_models)) {
    structure_models <- fitted_models[[structure]]
    
    if (!is.list(structure_models)) {
      stop(paste("Invalid structure: fitted_models[[", structure, 
                 "]] must be a list of model entries"))
    }
    
    if (length(structure_models) == 0) {
      warning(paste("Empty structure:", structure, "- no models found"))
      next
    }
    
    for (formula in names(structure_models)) {
      entry <- structure_models[[formula]]
      
      if (!is.list(entry)) {
        stop(paste("Invalid entry: fitted_models[[", structure, "]][[", formula, 
                   "]] must be a list with required components"))
      }
      
      # Check for required components in the new unified structure
      required_components <- c("fitted_model", "status", "optimization_method", "metadata")
      missing_components <- setdiff(required_components, names(entry))
      
      if (length(missing_components) > 0) {
        stop(paste("Missing components in fitted_models[[", structure, "]][[", formula, 
                   "]]: ", paste(missing_components, collapse = ", ")))
      }
      
      # Validate status values
      valid_statuses <- c("converged", "failed", "converged_with_warning")
      if (!is.null(entry$status) && !entry$status %in% valid_statuses) {
        warning(paste("Invalid status '", entry$status, "' for", structure, formula, 
                      ". Expected one of:", paste(valid_statuses, collapse = ", ")))
      }
      
      # Check that successful models have valid msm objects
      if (entry$status == "converged") {
        if (is.null(entry$fitted_model)) {
          stop(paste("Converged model has NULL fitted_model:", structure, formula))
        }
        
        if (!inherits(entry$fitted_model, "msm")) {
          stop(paste("fitted_model is not an msm object:", structure, formula))
        }
      }
      
      # Validate metadata structure
      if (!is.null(entry$metadata) && !is.list(entry$metadata)) {
        warning(paste("Metadata should be a list for", structure, formula))
      }
    }
  }
  
  return(TRUE)
}

# Additional helper function to get model summary statistics
summarize_model_structure <- function(fitted_models) {
  validate_model_structure(fitted_models)
  
  structure_summary <- list()
  
  for (structure in names(fitted_models)) {
    structure_models <- fitted_models[[structure]]
    
    total_models <- length(structure_models)
    converged_models <- sum(sapply(structure_models, function(x) x$status == "converged"))
    failed_models <- sum(sapply(structure_models, function(x) x$status == "failed"))
    
    # Count model types
    has_splines <- sum(sapply(structure_models, function(x) {
      if (!is.null(x$metadata) && !is.null(x$metadata$has_splines)) {
        return(x$metadata$has_splines)
      }
      return(FALSE)
    }))
    
    has_time_varying <- sum(sapply(structure_models, function(x) {
      if (!is.null(x$metadata) && !is.null(x$metadata$has_time_varying)) {
        return(x$metadata$has_time_varying)
      }
      return(FALSE)
    }))
    
    structure_summary[[structure]] <- list(
      total_models = total_models,
      converged = converged_models,
      failed = failed_models,
      convergence_rate = converged_models / total_models,
      models_with_splines = has_splines,
      models_with_time_varying = has_time_varying
    )
  }
  
  return(structure_summary)
}

# Helper function to extract all successfully fitted models
extract_successful_models <- function(fitted_models) {
  validate_model_structure(fitted_models)
  
  successful_models <- list()
  
  for (structure in names(fitted_models)) {
    successful_models[[structure]] <- list()
    
    for (formula in names(fitted_models[[structure]])) {
      model_entry <- fitted_models[[structure]][[formula]]
      
      if (model_entry$status == "converged" && !is.null(model_entry$fitted_model)) {
        successful_models[[structure]][[formula]] <- model_entry
      }
    }
    
    # Remove empty structures
    if (length(successful_models[[structure]]) == 0) {
      successful_models[[structure]] <- NULL
    }
  }
  
  return(successful_models)
}

state_name_to_num <- function(state_names, fitted_model) {
  if (is.null(fitted_model) || !inherits(fitted_model, "msm")) {
    warning("Invalid fitted_model provided")
    return(setNames(seq_along(unique(state_names)), unique(state_names)))
  }
  
  # Get state names from Q matrix
  model_states <- rownames(fitted_model$qmodel$qmatrix)
  state_mapping <- setNames(seq_along(model_states), model_states)
  
  # Map input states
  result <- state_mapping[state_names]
  names(result) <- state_names
  
  # Warn about unmapped states
  unmapped <- state_names[is.na(result)]
  if (length(unmapped) > 0) {
    warning(paste("States not found in model:", paste(unmapped, collapse = ", ")))
  }
  
  return(result)
}

state_num_to_name <- function(state_nums, fitted_model) {
  if (is.null(fitted_model) || !inherits(fitted_model, "msm")) {
    warning("Invalid fitted_model provided")
    return(setNames(as.character(state_nums), state_nums))
  }
  
  model_states <- rownames(fitted_model$qmodel$qmatrix)
  
  # Validate state numbers
  invalid_nums <- state_nums[state_nums < 1 | state_nums > length(model_states)]
  if (length(invalid_nums) > 0) {
    warning(paste("Invalid state numbers:", paste(invalid_nums, collapse = ", ")))
  }
  
  result <- model_states[state_nums]
  names(result) <- state_nums
  
  return(result)
}

standardize_state_columns <- function(data, fitted_model, config = NULL) {
  
  # If we have state but not state_num
  if ("state" %in% names(data) && !"state_num" %in% names(data)) {
    state_mapping <- state_name_to_num(data$state, fitted_model)
    data$state_num <- state_mapping[data$state]
  }
  
  # If we have state_num but not state  
  if ("state_num" %in% names(data) && !"state" %in% names(data)) {
    state_mapping <- state_num_to_name(data$state_num, fitted_model)
    data$state <- state_mapping[as.character(data$state_num)]
  }
  
  # Validate consistency if both exist
  if ("state" %in% names(data) && "state_num" %in% names(data)) {
    inconsistent <- which(data$state != state_num_to_name(data$state_num, fitted_model)[as.character(data$state_num)])
    if (length(inconsistent) > 0) {
      warning(paste("Inconsistent state/state_num mapping in", length(inconsistent), "rows"))
    }
  }
  
  return(data)
}


### Extract from metadata ------------------------------------------------

# Helper functions to extract information from metadata
extract_covariate_label_from_metadata <- function(metadata, formula_name) {
  if (!is.null(metadata$input_covariates) || !is.null(metadata$input_spline_vars) || !is.null(metadata$time_config)) {
    # Build label from metadata
    parts <- c()
    if (!is.null(metadata$input_covariates)) {
      parts <- c(parts, paste(metadata$input_covariates, collapse = ", "))
    }
    if (!is.null(metadata$input_spline_vars)) {
      spline_part <- paste0(metadata$input_spline_vars, " (", 
                            ifelse(is.null(metadata$spline_config$default_type), "spline", metadata$spline_config$default_type), 
                            ")")
      parts <- c(parts, paste(spline_part, collapse = ", "))
    }
    if (!is.null(metadata$time_config)) {
      time_part <- paste0(metadata$time_config$variable, " (", metadata$time_config$type, ")")
      parts <- c(parts, time_part)
    }
    
    if (length(parts) > 0) {
      return(paste(parts, collapse = " + "))
    }
  }
  
  # Fallback to formula parsing
  if (formula_name == "~ 1") {
    return("None")
  } else {
    vars <- tryCatch({
      all.vars(as.formula(formula_name))
    }, error = function(e) {
      return(character(0))
    })
    return(if (length(vars) > 0) paste(vars, collapse = ", ") else "Unknown")
  }
}

extract_base_variables_from_metadata <- function(metadata) {
  vars <- c()
  if (!is.null(metadata$input_covariates)) vars <- c(vars, metadata$input_covariates)
  if (!is.null(metadata$input_spline_vars)) vars <- c(vars, metadata$input_spline_vars)
  if (!is.null(metadata$time_config)) vars <- c(vars, metadata$time_config$variable)
  return(unique(vars))
}

extract_n_covariates_from_metadata <- function(metadata) {
  if (!is.null(metadata$n_covariates)) {
    return(metadata$n_covariates)
  }
  return(length(extract_base_variables_from_metadata(metadata)))
}

extract_has_splines_from_metadata <- function(metadata) {
  if (!is.null(metadata$has_splines)) {
    return(metadata$has_splines)
  }
  return(!is.null(metadata$spline_config) && length(metadata$spline_config$spline_vars) > 0)
}

extract_spline_types_from_metadata <- function(metadata) {
  if (!is.null(metadata$spline_config) && !is.null(metadata$spline_config$spline_metadata)) {
    types <- sapply(metadata$spline_config$spline_metadata, function(x) x$type)
    return(paste(unique(types), collapse = ","))
  }
  return(NA_character_)
}

extract_has_transformations_from_metadata <- function(metadata) {
  # For now, return FALSE - could be extended if transformation info is stored
  return(FALSE)
}

extract_has_time_varying_from_metadata <- function(metadata) {
  if (!is.null(metadata$has_time_varying)) {
    return(metadata$has_time_varying)
  }
  return(!is.null(metadata$time_config))
}

extract_time_varying_type_from_metadata <- function(metadata) {
  if (!is.null(metadata$time_config) && !is.null(metadata$time_config$type)) {
    return(metadata$time_config$type)
  }
  return(NA_character_)
}

extract_constraint_type_from_metadata <- function(metadata) {
  if (!is.null(metadata$constraint_config) && !is.null(metadata$constraint_config$type)) {
    return(metadata$constraint_config$type)
  }
  return(NA_character_)
}

### Helper functions to extract hazard ratios from fitted model ----------

# Helper function to check if a covariate term is a spline term
is_spline_covariate_term <- function(covariate_name, metadata) {
  if (is.null(metadata$spline_config) || is.null(metadata$spline_config$spline_metadata)) {
    return(FALSE)
  }
  
  # Check if this covariate name matches any spline terms
  for (spline_var in names(metadata$spline_config$spline_metadata)) {
    spline_meta <- metadata$spline_config$spline_metadata[[spline_var]]
    if (!is.null(spline_meta$spline_terms)) {
      if (covariate_name %in% spline_meta$spline_terms) {
        return(TRUE)
      }
    }
  }
  
  # Also check for common spline patterns in case metadata is incomplete
  spline_patterns <- c("_ns\\d+$", "_bs\\d+$", "_spline\\d+$")
  for (pattern in spline_patterns) {
    if (grepl(pattern, covariate_name)) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}

# Helper function to extract original variable name from a covariate term
extract_original_variable_name <- function(covariate_name, metadata) {
  # First, check spline metadata for exact matches
  if (!is.null(metadata$spline_config) && !is.null(metadata$spline_config$spline_metadata)) {
    for (spline_var in names(metadata$spline_config$spline_metadata)) {
      spline_meta <- metadata$spline_config$spline_metadata[[spline_var]]
      if (!is.null(spline_meta$spline_terms)) {
        if (covariate_name %in% spline_meta$spline_terms) {
          return(spline_var)
        }
      }
    }
  }
  
  # Check for time-varying terms
  if (!is.null(metadata$time_config)) {
    time_var <- metadata$time_config$variable
    # Check for various time-varying patterns
    if (covariate_name == time_var || 
        covariate_name == "time_period" ||
        grepl(paste0("^", time_var, "_"), covariate_name)) {
      return(time_var)
    }
  }
  
  # For regular variables and factor levels, extract base name
  # Handle factor levels (e.g., "sexMale" -> "sex")
  if (grepl("^[a-zA-Z_][a-zA-Z0-9_]*[A-Z]", covariate_name)) {
    # Try to find the base variable from input covariates
    if (!is.null(metadata$input_covariates)) {
      for (base_var in metadata$input_covariates) {
        if (startsWith(covariate_name, base_var)) {
          return(base_var)
        }
      }
    }
  }
  
  # Fallback: remove common suffixes and patterns
  original_name <- covariate_name
  
  # Remove spline suffixes
  original_name <- gsub("_(ns|bs|spline)\\d+$", "", original_name)
  
  # Remove quadratic suffix
  original_name <- gsub("_quad$", "", original_name)
  
  # If no transformation detected, return original
  if (original_name == covariate_name) {
    return(covariate_name)
  }
  
  return(original_name)
}

# Enhanced helper function for grouping hazard ratios by original variable
group_hazard_ratios_by_variable <- function(tidied_hrs) {
  if (nrow(tidied_hrs) == 0) {
    return(tidied_hrs)
  }
  
  # Add grouping information
  tidied_hrs %>%
    group_by(model, formula, transition, original_variable) %>%
    mutate(
      n_terms_for_variable = n(),
      is_multi_term_variable = n_terms_for_variable > 1,
      variable_type = case_when(
        is_spline_term & is_multi_term_variable ~ "spline",
        is_multi_term_variable & !is_spline_term ~ "factor_or_interaction", 
        !is_multi_term_variable ~ "single_term",
        TRUE ~ "unknown"
      )
    ) %>%
    ungroup()
}
