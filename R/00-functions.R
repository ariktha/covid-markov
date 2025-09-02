
# Model fitting functions ---------------------------------------------

## Calculate crude initial rates --------------------------------------

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

## Fit multi-state model ----------------------------------------------

fit_msm_models <- function(patient_data, 
                           crude_rates, 
                           covariates = NULL, 
                           spline_vars = NULL, 
                           spline_df = 3,
                           spline_type = "ns",
                           spline_degree = 3,
                           constraint = "default",
                           nest_name = NULL, 
                           mc.cores = min(8, parallel::detectCores() - 1)) {
  
  # Load required packages
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("parallel package is required for parallelization")
  }
  
  # Warn if non-ns spline types are requested and fall back to ns
  if (!is.null(spline_type)) {
    if (is.list(spline_type)) {
      non_ns_types <- unlist(spline_type)[unlist(spline_type) != "ns"]
    } else {
      non_ns_types <- if (spline_type != "ns") spline_type else character(0)
    }
    
    if (length(non_ns_types) > 0) {
      warning(paste("Only natural splines (ns) are currently implemented. Falling back to ns for requested types:", 
                    paste(unique(non_ns_types), collapse = ", ")))
    }
  }
  
  # Handle the case where covariates is NULL
  if (is.null(covariates)) {
    covariates <- list("no_covariates" = NULL)
  }
  
  # Ensure spline_df is a named list (can be variable-based or combination-based)
  if (!is.null(spline_df) && !is.list(spline_df)) {
    # If spline_df is a single value, convert to list
    spline_df <- list(default = spline_df)
  }
  
  # Handle constraint parameter - ALWAYS define constraint_for_msm
  constraint_for_msm <- NULL  # Default value
  
  if (!is.null(constraint)) {
    if (constraint == "default") {
      constraint_for_msm <- NULL  # Use MSM default behavior
    } else if (constraint == "unconstrained") {
      constraint_for_msm <- NULL
    } else if (constraint == "same") {
      constraint_for_msm <- NULL
    } else if (is.matrix(constraint)) {
      constraint_for_msm <- constraint
    } else {
      constraint_for_msm <- NULL
    }
  }
  
  # Get unique model names
  model_names <- unique(patient_data$model)
  
  # Define function to fit models for a single model structure
  fit_single_model_structure <- function(modelname) {
    # Ensure required packages are loaded in parallel worker
    if (!require("msm", quietly = TRUE)) {
      return(list(error = "msm package not available"))
    }
    if (!require("splines", quietly = TRUE)) {
      return(list(error = "splines package not available"))
    }
    
    # Explicitly attach splines to make ns() available in formula environment
    library(splines)
    
    model_data <- patient_data[which(patient_data$model == modelname), ]
    crude_result <- crude_rates[[modelname]]
    
    if (is.null(crude_result)) {
      warning(paste("Crude rates missing for", modelname, "- skipping model fitting"))
      return(NULL)
    }
    
    # Validate that required variables exist in the data
    all_vars_needed <- unique(unlist(c(covariates, spline_vars)))
    all_vars_needed <- all_vars_needed[!is.null(all_vars_needed)]
    if (length(all_vars_needed) > 0) {
      missing_vars <- setdiff(all_vars_needed, names(model_data))
      if (length(missing_vars) > 0) {
        warning(paste("Variables not found in data for model", modelname, ":", 
                      paste(missing_vars, collapse = ", ")))
        # Continue but exclude missing variables
      }
    }
    
    model_results <- list()
    
    # Loop through each covariate combination
    for (cov_name in names(covariates)) {
      current_covariates <- covariates[[cov_name]]
      current_splines <- if (!is.null(spline_vars)) {
        if (is.list(spline_vars)) spline_vars[[cov_name]] else spline_vars
      } else {
        NULL
      }
      
      # Filter out missing variables
      if (!is.null(current_covariates)) {
        current_covariates <- intersect(current_covariates, names(model_data))
      }
      if (!is.null(current_splines)) {
        current_splines <- intersect(current_splines, names(model_data))
      }
      
      # Create formula string for this combination
      formula_name <- "~ 1"  # Default intercept-only
      covariate_formula <- NULL
      
      if (!is.null(current_covariates) || !is.null(current_splines)) {
        formula_terms <- c()
        
        # Add regular covariates (those not in spline list)
        if (!is.null(current_covariates)) {
          non_spline_covariates <- setdiff(current_covariates, current_splines)
          if (length(non_spline_covariates) > 0) {
            formula_terms <- c(formula_terms, non_spline_covariates)
          }
        }
        
        # Add spline terms (only ns splines)
        if (!is.null(current_splines) && length(current_splines) > 0) {
          for (spline_var in current_splines) {
            # Get df value for this variable/combination
            df_value <- if (!is.null(spline_df) && !is.null(spline_df[[cov_name]])) {
              spline_df[[cov_name]]
            } else if (!is.null(spline_df) && !is.null(spline_df[[spline_var]])) {
              spline_df[[spline_var]]
            } else if (!is.null(spline_df) && !is.null(spline_df[["default"]])) {
              spline_df[["default"]]
            } else {
              3  # Final fallback
            }
            
            # Only use natural splines
            spline_term <- paste0("ns(", spline_var, ", df = ", df_value, ")")
            formula_terms <- c(formula_terms, spline_term)
          }
        }
        
        if (length(formula_terms) > 0) {
          formula_name <- paste("~", paste(formula_terms, collapse = " + "))
          # Create formula in an environment where splines functions are available
          covariate_formula <- as.formula(formula_name, env = environment())
        }
      }
      
      # Try different optimization methods if convergence fails
      optimization_methods <- list(
        list(opt_method = "optim", method = "BFGS", name = "BFGS"),
        list(opt_method = "bobyqa", name = "BOBYQA"),
        list(opt_method = "optim", method = "Nelder-Mead", name = "Nelder-Mead")
      )
      
      fitted_model <- NULL
      last_error <- NULL
      successful_method <- NULL
      
      for (opt_config in optimization_methods) {
        
        fitted_model <- tryCatch({
          # Base control parameters
          control_list <- list(fnscale = 10000, maxit = 5000, reltol = 1e-10)
          
          # Only add method to control if using optim
          if (opt_config$opt_method == "optim" && 
              "method" %in% names(opt_config) && 
              !is.null(opt_config$method) && 
              opt_config$method != "") {
            control_list$method <- opt_config$method
          }
          
          # Pre-validate spline basis if splines are used
          if (!is.null(covariate_formula)) {
            # Test that the spline basis can be created successfully
            basis_matrix <- tryCatch({
              model.matrix(covariate_formula, data = model_data)
            }, error = function(e) {
              stop(paste("Error creating spline basis matrix:", e$message))
            })
            
            # Debug: Check dimensions
            cat("Spline basis matrix dimensions:", dim(basis_matrix), "for formula:", formula_name, "\n")
            cat("Data dimensions:", nrow(model_data), "rows,", ncol(model_data), "cols\n")
            cat("Q-matrix dimensions:", dim(crude_result$qmat), "\n")
          }
          
          # Build MSM arguments
          msm_args <- list(
            formula = state_num ~ DaysSinceEntry,
            subject = quote(deid_enc_id),
            data = model_data,
            qmatrix = crude_result$qmat,
            opt.method = opt_config$opt_method,
            control = control_list
          )
          
          # Add covariates if they exist
          if (!is.null(covariate_formula)) {
            msm_args$covariates <- covariate_formula
          }
          
          # Add constraint if specified
          if (!is.null(constraint_for_msm)) {
            msm_args$constraint <- constraint_for_msm
          }
          
          # Fit the model
          cat("Attempting to fit MSM model with:", length(msm_args), "arguments\n")
          cat("Formula:", deparse(msm_args$formula), "\n")
          if (!is.null(msm_args$covariates)) {
            cat("Covariates formula:", deparse(msm_args$covariates), "\n")
          }
          
          result <- do.call(msm::msm, msm_args)
          cat("MSM fitting successful\n")
          result
          
        }, error = function(e) {
          last_error <<- e$message
          return(NULL)
        })
        
        # Check if model fitted successfully and converged
        if (!is.null(fitted_model)) {
          converged <- is.null(fitted_model$opt$convergence) || fitted_model$opt$convergence == 0
          
          if (converged) {
            message(paste("Model", modelname, "with formula", formula_name, "converged successfully using", opt_config$name))
            successful_method <- opt_config$name
            break
          } else {
            message(paste("Model", modelname, "with formula", formula_name, "failed to converge with", opt_config$name, "- trying next method"))
            fitted_model <- NULL  # Reset for next iteration
          }
        } else {
          message(paste("Model", modelname, "with formula", formula_name, "failed to fit with", opt_config$name, "- trying next method"))
        }
      }
      
      # Store result or failure information
      if (is.null(fitted_model)) {
        failure_msg <- paste("All optimization methods failed for", modelname, "with formula", formula_name)
        if (!is.null(last_error)) {
          failure_msg <- paste(failure_msg, "- Last error:", last_error)
        }
        warning(failure_msg)
        model_results[[formula_name]] <- list(
          fitted_model = NULL,
          error_message = failure_msg,
          status = "failed",
          optimization_method = NA,
          spline_config = NULL,
          constraint_type = NA
        )
      } else {
        # Create spline configuration metadata (fixed scope issues)
        spline_config <- list(
          spline_vars = current_splines,
          spline_type = if (!is.null(current_splines) && length(current_splines) > 0) "ns" else NULL,
          spline_df = if (!is.null(current_splines) && length(current_splines) > 0) {
            # Get the df value that was actually used
            if (!is.null(spline_df) && !is.null(spline_df[[cov_name]])) {
              spline_df[[cov_name]]
            } else if (!is.null(spline_df) && !is.null(spline_df[["default"]])) {
              spline_df[["default"]]
            } else {
              3
            }
          } else NULL,
          spline_degree = NULL  # Not applicable for ns splines
        )
        
        model_results[[formula_name]] <- list(
          fitted_model = fitted_model,
          error_message = NULL,
          status = "converged",
          optimization_method = successful_method,
          spline_config = spline_config,
          constraint_type = constraint
        )        
      }
    }
    
    return(model_results)
  }
  
  # Parallelize across model structures
  fitted_models_list <- parallel::mclapply(
    model_names,
    fit_single_model_structure,
    mc.cores = mc.cores
  )
  
  # Convert back to named list and remove NULL entries (but keep failed models)
  fitted_msm_models <- setNames(fitted_models_list, model_names)
  fitted_msm_models <- fitted_msm_models[!sapply(fitted_msm_models, is.null)]
  
  if (!is.null(nest_name)) {
    fitted_models_nest <- list()
    fitted_models_nest[[nest_name]] <- fitted_msm_models
    return(fitted_models_nest)
  }
  
  return(fitted_msm_models)
}

# Tidy model output ---------------------------------------------------

## Wrappers -----------------------------------------------------------

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
  
  # Helper function to detect spline usage and type
  detect_splines <- function(formula_str) {
    spline_info <- list(
      has_splines = FALSE,
      spline_types = c(),
      has_transformations = FALSE,
      transformation_types = c()
    )
    
    # Check for different spline types (with or without package prefix)
    if (grepl("\\bns\\(", formula_str)) {
      spline_info$has_splines <- TRUE
      spline_info$spline_types <- c(spline_info$spline_types, "ns")
    }
    if (grepl("\\bbs\\(", formula_str)) {
      spline_info$has_splines <- TRUE
      spline_info$spline_types <- c(spline_info$spline_types, "bs")
    }
    if (grepl("\\bpspline\\(", formula_str)) {
      spline_info$has_splines <- TRUE
      spline_info$spline_types <- c(spline_info$spline_types, "pspline")
    }
    
    # Check for other transformations
    if (grepl("\\bpoly\\(", formula_str)) {
      spline_info$has_transformations <- TRUE
      spline_info$transformation_types <- c(spline_info$transformation_types, "poly")
    }
    if (grepl("\\bI\\(", formula_str)) {
      spline_info$has_transformations <- TRUE
      spline_info$transformation_types <- c(spline_info$transformation_types, "arithmetic")
    }
    if (grepl("\\b(log|sqrt|exp)\\(", formula_str)) {
      spline_info$has_transformations <- TRUE
      spline_info$transformation_types <- c(spline_info$transformation_types, "math")
    }
    
    return(spline_info)
  }
  
  # Helper function to extract covariate info
  extract_covariate_info <- function(formula_str) {
    if (formula_str == "~ 1") {
      list(
        covariate_label = "None",
        has_splines = FALSE,
        spline_types = NA_character_,
        has_transformations = FALSE,
        n_covariates = 0,
        base_variables = character(0)
      )
    } else {
      # Extract variable names
      vars <- tryCatch({
        all.vars(as.formula(formula_str))
      }, error = function(e) {
        # If formula parsing fails, try to extract variables manually
        # This handles cases where the formula might be malformed
        vars_pattern <- "\\b[a-zA-Z_][a-zA-Z0-9_]*\\b"
        potential_vars <- regmatches(formula_str, gregexpr(vars_pattern, formula_str))[[1]]
        # Filter out function names and common R keywords
        function_names <- c("ns", "bs", "pspline", "poly", "I", "log", "sqrt", "exp", 
                            "splines", "survival", "df", "degree", "nterm")
        potential_vars[!potential_vars %in% function_names]
      })
      
      # Get spline information
      spline_info <- detect_splines(formula_str)
      
      # Create covariate label
      covariate_label <- if (length(vars) > 0) {
        base_label <- paste(vars, collapse = ", ")
        
        # Add spline type information if present
        if (spline_info$has_splines) {
          spline_types_str <- paste(unique(spline_info$spline_types), collapse = "/")
          base_label <- paste0(base_label, " (", spline_types_str, " splines)")
        } else if (spline_info$has_transformations) {
          base_label <- paste0(base_label, " (transformed)")
        }
        base_label
      } else {
        "Unknown"
      }
      
      list(
        covariate_label = covariate_label,
        has_splines = spline_info$has_splines,
        spline_types = if (spline_info$has_splines) paste(spline_info$spline_types, collapse = ",") else NA_character_,
        has_transformations = spline_info$has_transformations,
        n_covariates = length(vars),
        base_variables = vars
      )
    }
  }
  
  # Loop through model structures (outer level)
  for (model_structure in names(fitted_msm_models)) {
    # Loop through covariate formulas (inner level)
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      
      # Extract covariate information
      cov_info <- extract_covariate_info(formula_name)
      
      # Check if this is a failed model
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && model_entry$status == "failed")) {
        
        # Extract optimization method if available
        opt_method <- if (is.list(model_entry) && !is.null(model_entry$optimization_method)) {
          model_entry$optimization_method
        } else {
          NA_character_
        }
        
        # Include failed models in output with NAs
        failed_tidy <- tibble::tibble(
          model = model_structure,
          formula = formula_name,
          covariate = cov_info$covariate_label,
          base_variables = list(cov_info$base_variables),
          loglik = NA_real_,
          AIC = NA_real_,
          BIC = NA_real_,
          n_params = NA_integer_,
          n_obs = NA_integer_,
          has_covariates = formula_name != "~ 1",
          n_covariates = cov_info$n_covariates,
          has_splines = cov_info$has_splines,
          spline_types = cov_info$spline_types,
          has_transformations = cov_info$has_transformations,
          optimization_method = opt_method,
          status = if (is.list(model_entry) && !is.null(model_entry$status)) model_entry$status else "failed",
          error_message = if (is.list(model_entry) && !is.null(model_entry$error_message)) model_entry$error_message else "Model is NULL"
        )
        tidied_models <- dplyr::bind_rows(tidied_models, failed_tidy)
        next
      }
      
      # Extract the fitted model object
      fitted_model <- if (is.list(model_entry) && !is.null(model_entry$fitted_model)) {
        model_entry$fitted_model
      } else {
        model_entry  # Backwards compatibility if structure is different
      }
      
      # Extract optimization method if available
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
        n_params <- length(fitted_model$estimates)
        n_obs <- length(unique(fitted_model[["data"]][["mf"]][["(subject)"]]))
        
        # Check convergence
        converged <- is.null(fitted_model$opt$convergence) || fitted_model$opt$convergence == 0
        
        tidy_result <- tibble::tibble(
          model = model_structure,
          formula = formula_name,
          covariate = cov_info$covariate_label,
          base_variables = list(cov_info$base_variables),
          loglik = loglik,
          AIC = aic,
          BIC = bic,
          n_params = n_params,
          n_obs = n_obs,
          has_covariates = formula_name != "~ 1",
          n_covariates = cov_info$n_covariates,
          has_splines = cov_info$has_splines,
          spline_types = cov_info$spline_types,
          has_transformations = cov_info$has_transformations,
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
          covariate = cov_info$covariate_label,
          base_variables = list(cov_info$base_variables),
          loglik = NA_real_,
          AIC = NA_real_,
          BIC = NA_real_,
          n_params = NA_integer_,
          n_obs = NA_integer_,
          has_covariates = formula_name != "~ 1",
          n_covariates = cov_info$n_covariates,
          has_splines = cov_info$has_splines,
          spline_types = cov_info$spline_types,
          has_transformations = cov_info$has_transformations,
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

#' Comprehensive MSM Model Analysis Wrapper
#' 
#' Runs all major MSM analysis functions on fitted model list and returns organized results
#' 
#' @param fitted_msm_models Nested list of fitted MSM models (structure[[formula]]$fitted_model)
#' @param patient_data Original patient data used to fit models
#' @param crude_rates List of crude rates used in model fitting
#' @param analysis_config List of analysis configuration parameters
#' @param parallel Logical, whether to use parallel processing where available
#' @param mc.cores Number of cores for parallel processing (auto-detect if NULL)
#' @param verbose Logical, whether to print progress messages
#' 
#' @return Named list with all analysis results
#' 
#' @examples
#' # Basic usage with defaults
#' results <- run_comprehensive_msm_analysis(fitted_models, patient_data, crude_rates)
#' 
#' # Customized analysis
#' config <- list(
#'   prevalence = list(time_points = seq(1, 60, by = 2), ci = TRUE),
#'   cv = list(k_folds = 10, prediction_times = c(7, 14, 30))
#' )
#' results <- run_comprehensive_msm_analysis(fitted_models, patient_data, crude_rates, 
#'                                          analysis_config = config)
#' Comprehensive MSM Model Analysis Wrapper
#' 
#' Runs all major MSM analysis functions on fitted model list and returns organized results
#' 
#' @param fitted_msm_models Nested list of fitted MSM models (structure[[formula]]$fitted_model)
#' @param patient_data Original patient data used to fit models
#' @param crude_rates List of crude rates used in model fitting
#' @param analysis_config List of analysis configuration parameters
#' @param parallel Logical, whether to use parallel processing where available
#' @param mc.cores Number of cores for parallel processing (auto-detect if NULL)
#' @param verbose Logical, whether to print progress messages
#' 
#' @return Named list with all analysis results
#' 
#' @examples
#' # Basic usage with defaults
#' results <- run_comprehensive_msm_analysis(fitted_models, patient_data, crude_rates)
#' 
#' # Customized analysis
#' config <- list(
#'   prevalence = list(time_points = seq(1, 60, by = 2), ci = TRUE),
#'   cv = list(k_folds = 10, prediction_times = c(7, 14, 30))
#' )
#' results <- run_comprehensive_msm_analysis(fitted_models, patient_data, crude_rates, 
#'                                          analysis_config = config)
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
      covariates_list = NULL
    ),
    
    # P-matrix parameters  
    pmats = list(
      t_values = c(1, 7, 14, 30),
      covariates_list = NULL
    ),
    
    # Sojourn time parameters
    sojourns = list(
      covariates_list = NULL
    ),
    
    # Prevalence parameters
    prevalence = list(
      time_points = seq(1, 30, by = 1),
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
      k_folds = 5,
      prediction_times = c(14, 30),
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
      cores = 1  # DRAIC/DRLCV are computationally intensive
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
    do.call(tidy_msm_models, c(list(fitted_msm_models), config$tidy_models))
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
    do.call(tidy_msm_qmatrix, c(list(fitted_msm_models), config$qmatrix))
  }, error = function(e) {
    warning(paste("tidy_msm_qmatrix failed:", e$message))
    data.frame()
  })
  
  # 3. Transition Probabilities (P-matrix)
  if (verbose) cat("3. Calculating transition probabilities...\n")
  results$pmats <- tryCatch({
    do.call(tidy_msm_pmats, c(list(fitted_msm_models), config$pmats))
  }, error = function(e) {
    warning(paste("tidy_msm_pmats failed:", e$message))
    data.frame()
  })
  
  # 4. Sojourn Times
  if (verbose) cat("4. Extracting sojourn times...\n")
  results$sojourns <- tryCatch({
    do.call(tidy_msm_sojourns, c(list(fitted_msm_models), config$sojourns))
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
      do.call(tidy_msm_hazards, c(list(fitted_msm_models), config$hazards))
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
      do.call(calculate_predictive_performance, 
              c(list(patient_data = patient_data, 
                     fitted_models = fitted_msm_models, 
                     crude_rates = crude_rates), 
                config$cv))
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
                     patient_data = patient_data), 
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
    do.call(comprehensive_model_comparison, c(list(fitted_msm_models), config$comparison))
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

## Model output wrangling ---------------------------------------------

### Transition intensities --------------------------------------------

tidy_msm_qmatrix <- function(fitted_msm_models, 
                             covariates_list = NULL) {
  
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
      
      if (is.null(fitted_model)) {
        warning(paste("Fitted model is NULL for", model_structure, formula_name))
        next
      }
      
      model_tidy <- tryCatch({
        qmat_result <- qmatrix.msm(fitted_model, ci = "normal")
        
        # Extract variable names from formula for covariate info
        covariate_vars <- if (formula_name == "~ 1") {
          "None"
        } else {
          vars_in_formula <- all.vars(as.formula(formula_name))
          paste(vars_in_formula, collapse = ", ")
        }
        
        # Use tidy() directly and add metadata
        tidy(qmat_result) %>%
          mutate(
            model = model_structure,
            formula = formula_name,
            covariate = covariate_vars,
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
                           covariates_list = NULL) {
  
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
      
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in tidy_msm_qmatrix")
      if (is.null(fitted_model)) next
      
      if (is.null(fitted_model)) {
        warning(paste("Fitted model is NULL for", model_structure, formula_name))
        next
      }
      
      # Extract variable names from formula for covariate info
      covariate_vars <- if (formula_name == "~ 1") {
        "None"
      } else {
        vars_in_formula <- all.vars(as.formula(formula_name))
        paste(vars_in_formula, collapse = ", ")
      }
      
      # Handle multiple time points
      for (t_val in t_values) {
        # Handle multiple covariate combinations
        if (is.null(covariates_list)) {
          model_tidy <- tryCatch({
            pmat_result <- pmatrix.msm(fitted_model, t = t_val, ci = "normal")
            
            # Use tidy() directly and add metadata
            tidy(pmat_result) %>%
              mutate(
                model = model_structure,
                formula = formula_name,
                covariate = covariate_vars,
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
                                         covariates = cov_combo, ci = "normal")
              
              # Use tidy() directly and add metadata
              tidy(pmat_result) %>%
                mutate(
                  model = model_structure,
                  formula = formula_name,
                  covariate = covariate_vars,
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
      
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in tidy_msm_qmatrix")
      if (is.null(fitted_model)) next
      
      if (is.null(fitted_model)) {
        warning(paste("Fitted model is NULL for", model_structure, formula_name))
        next
      }
      
      # Extract variable names from formula for covariate info
      covariate_vars <- if (formula_name == "~ 1") {
        "None"
      } else {
        vars_in_formula <- all.vars(as.formula(formula_name))
        paste(vars_in_formula, collapse = ", ")
      }
      
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

### State prevalence --------------------------------------------------

tidy_msm_prevalences_old <- function(fitted_msm_models, 
                                 time_points = seq(1, 30, by = 1),
                                 covariates_list = NULL, 
                                 mc.cores = min(8, parallel::detectCores() - 1),
                                 ci = TRUE,
                                 ci_method = "normal",
                                 use_approximation = FALSE,
                                 approx_method = "midpoint") {
  
  # Load required packages
  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install required packages: install.packages(c('future', 'future.apply'))")
  }
  
  library(future)
  library(future.apply)
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
  
  # Build model combinations
  model_combinations <- list()
  for (model_structure in names(fitted_msm_models)) {
    for (formula_name in names(fitted_msm_models[[model_structure]])) {
      model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
      if (!is.null(model_entry) && 
          ((is.list(model_entry) && is.null(model_entry$status)) || 
           (is.list(model_entry) && model_entry$status == "converged"))) {
        model_combinations[[paste(model_structure, formula_name, sep = "___")]] <- 
          list(model_structure = model_structure, formula_name = formula_name)
      }
    }
  }
  
  cat("Processing", length(model_combinations), "model combinations with", 
      length(time_points), "time points\n")
  
  # Set up parallel processing
  plan(multisession, workers = mc.cores)
  cat("Using", mc.cores, "parallel workers\n")
  
  # Process combinations in parallel
  prevalence_list <- future_lapply(model_combinations, function(combo_info) {
    library(msm)
    library(dplyr)
    
    model_structure <- combo_info$model_structure
    formula_name <- combo_info$formula_name
    
    # Get the model
    model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
    fitted_model <- if (is.list(model_entry) && !is.null(model_entry$fitted_model)) {
      model_entry$fitted_model
    } else {
      model_entry
    }
    
    if (is.null(fitted_model)) return(NULL)
    
    # Get covariate variable names
    covariate_vars <- if (formula_name == "~ 1") {
      "None"
    } else {
      paste(all.vars(as.formula(formula_name)), collapse = ", ")
    }
    
    # Process prevalences
    tryCatch({
      if (is.null(covariates_list)) {
        # Build arguments for prevalence.msm
        prev_args <- list(x = fitted_model, times = time_points)
        if (ci) prev_args$ci <- ci_method
        
        # Check if approximation is supported and use valid method
        if (use_approximation) {
          # Check valid options for interp parameter
          if ("interp" %in% names(formals(prevalence.msm))) {
            valid_methods <- c("start", "midpoint")
            if (approx_method %in% valid_methods) {
              prev_args$interp <- approx_method
            } else {
              prev_args$interp <- "midpoint"  # Default safe option
            }
          }
        }
        
        prev_result <- do.call(prevalence.msm, prev_args)
        
        prevalence_to_tib(prev_result, has_ci = ci) %>%
          mutate(
            model = model_structure,
            formula = formula_name,
            covariate = covariate_vars,
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
          
          cov_results[[i]] <- prevalence_to_tib(prev_result, has_ci = ci) %>%
            mutate(
              model = model_structure,
              formula = formula_name,
              covariate = covariate_vars,
              cov_name = names(cov_combo)[1],
              cov_value = cov_combo[[1]],
              .after = time
            )
        }
        do.call(bind_rows, cov_results)
      }
    }, error = function(e) {
      # Return error info for debugging
      return(list(error = e$message, model = model_structure, formula = formula_name))
    })
  }, future.seed = TRUE)
  
  # Clean up parallel processing
  plan(sequential)
  
  # Check for errors and combine results
  error_results <- prevalence_list[sapply(prevalence_list, function(x) is.list(x) && "error" %in% names(x))]
  if (length(error_results) > 0) {
    cat("Errors found:\n")
    for (err in error_results) {
      cat("  -", err$model, err$formula, ":", err$error, "\n")
    }
  }
  
  # Filter out errors and NULLs
  valid_results <- prevalence_list[sapply(prevalence_list, function(x) is.data.frame(x))]
  
  if (length(valid_results) == 0) {
    warning("No valid prevalence results calculated")
    return(data.frame())
  }
  
  final_result <- do.call(bind_rows, valid_results)
  
  total_time <- Sys.time() - start_time
  cat("=== COMPLETED in", as.numeric(total_time, units = "secs"), "seconds ===\n")
  cat("Final dimensions:", nrow(final_result), "x", ncol(final_result), "\n")
  
  return(final_result)
}

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
          ((is.list(model_entry) && is.null(model_entry$status)) || 
           (is.list(model_entry) && model_entry$status == "converged"))) {
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
    
    # Extract model object inside worker to minimize data transfer
    model_entry <- fitted_msm_models[[model_structure]][[formula_name]]
    fitted_model <- if (is.list(model_entry) && !is.null(model_entry$fitted_model)) {
      model_entry$fitted_model
    } else {
      model_entry
    }
    
    if (is.null(fitted_model)) return(NULL)
    
    # Get covariate variable names
    covariate_vars <- if (formula_name == "~ 1") {
      "None"
    } else {
      paste(all.vars(as.formula(formula_name)), collapse = ", ")
    }
    
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

# FIXED: Added fitted_model parameter
prevalence_to_tib <- function(prevalence_result, fitted_model, has_ci = TRUE) {
  # Get basic info
  total_at_risk <- prevalence_result$Observed[1, "Total"]
  time_points <- as.numeric(rownames(prevalence_result$Observed))
  
  # Handle different structures based on CI computation
  observed_cols <- colnames(prevalence_result$Observed)
  
  # Expected data structure varies with CI:
  if (is.list(prevalence_result$Expected) && "estimates" %in% names(prevalence_result$Expected)) {
    expected_cols <- colnames(prevalence_result$Expected$estimates)
    expected_data <- prevalence_result$Expected$estimates
    has_ci_data <- has_ci && !is.null(prevalence_result$Expected$ci)
  } else {
    expected_cols <- colnames(prevalence_result$Expected)
    expected_data <- prevalence_result$Expected
    has_ci_data <- FALSE
  }
  
  # Same for Expected percentages
  if (is.list(prevalence_result$`Expected percentages`) && "estimates" %in% names(prevalence_result$`Expected percentages`)) {
    expected_pct_data <- prevalence_result$`Expected percentages`$estimates
  } else {
    expected_pct_data <- prevalence_result$`Expected percentages`
  }
  
  # Get state columns (exclude "Total")
  observed_states <- observed_cols[observed_cols != "Total"]
  expected_states <- expected_cols[expected_cols != "Total"]
  
  # Build result data frame
  result_list <- list()
  
  for (t_idx in seq_along(time_points)) {
    for (s_idx in seq_along(observed_states)) {
      # Extract values
      obs_count <- prevalence_result$Observed[t_idx, observed_states[s_idx]]
      obs_pct <- prevalence_result$`Observed percentages`[t_idx, observed_states[s_idx]]
      exp_count <- expected_data[t_idx, expected_states[s_idx]]
      exp_pct <- expected_pct_data[t_idx, expected_states[s_idx]]
      
      # FIXED: Now fitted_model is available as a parameter
      # Ensure consistent state naming
      state_name <- if (is.numeric(observed_states[s_idx])) {
        # If we have state numbers, convert to names
        tryCatch({
          rownames(fitted_model$qmodel$qmatrix)[as.numeric(observed_states[s_idx])]
        }, error = function(e) {
          # Fallback to the original state identifier
          observed_states[s_idx]
        })
      } else {
        observed_states[s_idx]
      }
      
      # Handle confidence intervals
      if (has_ci_data) {
        exp_count_ll <- prevalence_result$Expected$ci[t_idx, expected_states[s_idx], 1]
        exp_count_ul <- prevalence_result$Expected$ci[t_idx, expected_states[s_idx], 2]
        exp_pct_ll <- prevalence_result$`Expected percentages`$ci[t_idx, expected_states[s_idx], 1]
        exp_pct_ul <- prevalence_result$`Expected percentages`$ci[t_idx, expected_states[s_idx], 2]
      } else {
        exp_count_ll <- exp_count_ul <- exp_pct_ll <- exp_pct_ul <- NA_real_
      }
      
      result_list[[length(result_list) + 1]] <- data.frame(
        time = time_points[t_idx],
        state_num = observed_states[s_idx],
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
    }
  }
  
  # Combine and add derived variables
  combined_prev <- do.call(bind_rows, result_list)
  
  # Calculate N_at_risk
  death_recovery_states <- c("D", "R")
  N_at_risk_df <- combined_prev %>%
    filter(state %in% death_recovery_states) %>%
    group_by(time) %>%
    summarise(N_at_risk = total_at_risk - sum(observed_count, na.rm = TRUE), .groups = "drop")
  
  # Add final calculations
  final_result <- combined_prev %>%
    left_join(N_at_risk_df, by = "time") %>%
    mutate(
      count_residual = observed_count - expected_count,
      percentage_residual = observed_percentage - expected_percentage,
      binomial_variance = N_at_risk * (expected_percentage/100) * (1 - expected_percentage/100),
      standardized_residual = ifelse(binomial_variance > 0 & !is.na(binomial_variance), 
                                     count_residual / sqrt(binomial_variance), NA_real_),
      residual_category = case_when(
        is.na(standardized_residual) ~ "Unknown",
        abs(standardized_residual) > 2 ~ "Large",
        abs(standardized_residual) > 1 ~ "Moderate", 
        TRUE ~ "Small"
      )
    ) %>%
    arrange(time, state)
  
  return(final_result)
}

# Cached version for repeated analyses
tidy_msm_prevalences_cached <- function(fitted_msm_models, 
                                        time_points = seq(1, 30, by = 1),
                                        cache_dir = "msm_cache",
                                        ...) {
  
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package required for caching. Install with: install.packages('digest')")
  }
  
  # Create cache directory
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Generate cache key
  cache_key <- digest::digest(list(
    models = lapply(fitted_msm_models, names),
    time_points = time_points,
    args = list(...)
  ))
  
  cache_file <- file.path(cache_dir, paste0("prevalence_", cache_key, ".rds"))
  
  # Check cache
  if (file.exists(cache_file)) {
    cat("Loading cached result...\n")
    return(readRDS(cache_file))
  }
  
  # Compute and cache
  cat("Computing prevalences...\n")
  result <- tidy_msm_prevalences(fitted_msm_models, time_points, ...)
  
  saveRDS(result, cache_file)
  cat("Result cached\n")
  
  return(result)
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
      
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in tidy_msm_qmatrix")

      if (is.null(fitted_model)) {
        warning(paste("Fitted model is NULL for", model_structure, formula_name))
        next
      }
      
      # Check if model has covariates (hazards only meaningful for models with covariates)
      if (is.null(fitted_model$covariates) || formula_name == "~ 1") {
        warning(paste("Skipping hazards for", model_structure, formula_name, "- model has no covariates"))
        next
      }
      
      # Extract variable names from formula for covariate info
      covariate_vars <- if (formula_name == "~ 1") {
        "None"
      } else {
        vars_in_formula <- all.vars(as.formula(formula_name))
        paste(vars_in_formula, collapse = ", ")
      }
      
      model_tidy <- tryCatch({
        hr_result <- hazard.msm(fitted_model, hazard.scale = hazard_scale)
        
        # Manual tidying of hazard ratios (list of matrices)
        hr_tidy_list <- map_dfr(names(hr_result), function(covariate_name) {
          hr_matrix <- hr_result[[covariate_name]]
          
          # Extract transitions (row names)
          transitions <- rownames(hr_matrix)
          
          # Create tidy data frame
          data.frame(
            transition = transitions,
            covariate_term = covariate_name,
            hazard_ratio = hr_matrix[, "HR"],
            hr_lower = hr_matrix[, "L"], 
            hr_upper = hr_matrix[, "U"],
            log_hr = log(hr_matrix[, "HR"]),
            stringsAsFactors = FALSE
          )
        })
        
        # Add metadata
        hr_tidy <- hr_tidy_list %>%
          mutate(
            model = model_structure,
            formula = formula_name,
            covariate = covariate_vars,
            hazard_scale = hazard_scale
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

## Predictive performance evaluation ----------------------------------

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
                                             n_cores = min(8, parallel::detectCores() - 1)) {
  
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
      
      # Extract the fitted model object
      fitted_model <- if (is.list(model_entry) && !is.null(model_entry$fitted_model)) {
        model_entry$fitted_model
      } else {
        model_entry  # Backwards compatibility
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
        future_map_dfr(1:k_folds, function(fold) {
          cv_fold_core(fold, model_fold_assignments, model_data, crude_rates[[model_structure]], 
                       covariate_formula, prediction_times, calibration_covariates, calibration_subgroups, model_entry)
        }, .options = furrr_options(seed = TRUE))
      } else {
        map_dfr(1:k_folds, function(fold) {
          cv_fold_core(fold, model_fold_assignments, model_data, crude_rates[[model_structure]], 
                       covariate_formula, prediction_times, calibration_covariates, calibration_subgroups, model_entry)
        })
      }
      
      # Extract variable names from formula for covariate info
      covariate_vars <- if (formula_name == "~ 1") {
        "None"
      } else {
        vars_in_formula <- all.vars(as.formula(formula_name))
        paste(vars_in_formula, collapse = ", ")
      }
      
      # Aggregate results with time-specific and covariate-specific calibration
      base_results <- fold_results %>%
        filter(is.na(prediction_time) & is.na(subgroup_var)) %>%
        summarise(
          model = model_structure,
          formula = formula_name,
          covariate = covariate_vars,
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
#' @param fold_num Fold number
#' @param fold_assignments Data frame with patient fold assignments
#' @param model_data Data for specific model
#' @param crude_rate Crude rates for model initialization
#' @param covariate_formula Formula for covariates (NULL for none)
#' @param prediction_times Vector of prediction horizons
#' @param calibration_covariates List of covariate values to evaluate calibration at
#' @param calibration_subgroups List of variables to stratify calibration by
#' @param original_model_entry Original model entry containing spline metadata
#' @return Data frame with fold results including time and subgroup-specific calibration
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
  
  # Extract spline metadata from the original fitted model
  spline_metadata <- if (!is.null(original_model_entry) && is.list(original_model_entry) && !is.null(original_model_entry$spline_config)) {
    original_model_entry$spline_config
  } else {
    NULL
  }
  
  # Refit model on training data with spline metadata
  fitted_models_fold <- refit_model_for_fold(train_data, train_crude_rates, covariate_formula, spline_metadata)
  
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
  
  # Generate predictions for multiple time points and covariate values
  all_predictions <- generate_enhanced_predictions(
    test_outcomes, test_baseline, fitted_model, 
    covariate_formula, prediction_times, calibration_covariates
  )
  
  # Calculate base performance metrics (original functionality)
  base_results <- calculate_base_performance_metrics(
    fold_num, nrow(train_data), nrow(test_data), model_converged,
    test_outcomes, all_predictions$base_predictions
  )
  
  # Calculate time-specific calibration metrics
  time_specific_results <- calculate_time_specific_calibration(
    fold_num, nrow(train_data), nrow(test_data), model_converged,
    all_predictions$time_specific_predictions, test_outcomes, test_baseline
  )
  
  # Calculate subgroup-specific calibration metrics
  subgroup_specific_results <- calculate_subgroup_calibration(
    fold_num, nrow(train_data), nrow(test_data), model_converged,
    all_predictions$base_predictions, test_outcomes, test_baseline, calibration_subgroups
  )
  
  # Combine all results
  combined_results <- bind_rows(
    base_results,
    time_specific_results,
    subgroup_specific_results
  )
  
  return(combined_results)
}

### Helper functions for calibration/predictive performance -----------------

# Helper function to create empty fold result
create_empty_fold_result <- function(fold_num, n_train, n_test) {
  data.frame(
    fold = fold_num,
    n_train = n_train,
    n_test = n_test,
    model_converged = FALSE,
    total_los_mae = NA,
    total_los_mae_soj = NA,
    total_los_mae_vis_soj = NA,
    days_severe_mae = NA,
    days_severe_mae_soj = NA,
    days_severe_mae_vis_soj = NA,
    death_auc = NA,
    severe_auc = NA,
    calibration_slope_death = NA,
    calibration_intercept_death = NA,
    brier_death = NA,
    prediction_time = NA,
    subgroup_var = NA,
    subgroup_value = NA
  )
}

# Helper function to refit model for CV fold using stored spline metadata
refit_model_for_fold <- function(train_data, train_crude_rates, covariate_formula, spline_metadata = NULL) {
  if (!is.null(covariate_formula)) {
    all_vars <- all.vars(covariate_formula)
    
    # Use stored spline metadata if available
    if (!is.null(spline_metadata)) {
      covariates_list <- list("formula_combination" = all_vars)
      
      # Extract spline configuration from metadata
      spline_vars_list <- if (!is.null(spline_metadata$spline_vars)) {
        list("formula_combination" = spline_metadata$spline_vars)
      } else {
        NULL
      }
      
      spline_type_list <- if (!is.null(spline_metadata$spline_type)) {
        list("formula_combination" = spline_metadata$spline_type)
      } else {
        NULL
      }
      
      spline_df_list <- if (!is.null(spline_metadata$spline_df)) {
        setNames(list(spline_metadata$spline_df), "formula_combination")
      } else {
        NULL
      }
      
      spline_degree_list <- if (!is.null(spline_metadata$spline_degree)) {
        setNames(list(spline_metadata$spline_degree), "formula_combination")
      } else {
        NULL
      }
      
      # Call fit_msm_models with stored spline configuration
      fit_msm_models(train_data, train_crude_rates, 
                     covariates = covariates_list,
                     spline_vars = spline_vars_list,
                     spline_type = spline_type_list,
                     spline_df = spline_df_list,
                     spline_degree = spline_degree_list)
    } else {
      # Fallback to basic parsing for backwards compatibility
      formula_str <- deparse(covariate_formula)
      
      # Enhanced pattern matching for different spline types
      ns_pattern <- "ns\\(([^,]+),\\s*df\\s*=\\s*(\\d+)\\)"
      bs_pattern <- "bs\\(([^,]+),\\s*df\\s*=\\s*(\\d+)(?:,\\s*degree\\s*=\\s*(\\d+))?\\)"
      pspline_pattern <- "pspline\\(([^,]+),\\s*nterm\\s*=\\s*(\\d+)\\)"
      
      # Check for any spline terms
      ns_matches <- regmatches(formula_str, gregexpr(ns_pattern, formula_str))[[1]]
      bs_matches <- regmatches(formula_str, gregexpr(bs_pattern, formula_str))[[1]]
      pspline_matches <- regmatches(formula_str, gregexpr(pspline_pattern, formula_str))[[1]]
      
      if (length(ns_matches) > 0 || length(bs_matches) > 0 || length(pspline_matches) > 0) {
        # Parse spline information
        spline_vars <- character(0)
        spline_df <- list()
        spline_type <- list()
        spline_degree <- list()
        
        # Process ns() splines
        if (length(ns_matches) > 0) {
          ns_info <- str_match_all(formula_str, ns_pattern)[[1]]
          if (nrow(ns_info) > 0) {
            spline_vars <- c(spline_vars, ns_info[, 2])
            for (i in seq_len(nrow(ns_info))) {
              var_name <- ns_info[i, 2]
              spline_df[[var_name]] <- as.numeric(ns_info[i, 3])
              spline_type[[var_name]] <- "ns"
            }
          }
        }
        
        # Process bs() splines
        if (length(bs_matches) > 0) {
          bs_info <- str_match_all(formula_str, bs_pattern)[[1]]
          if (nrow(bs_info) > 0) {
            spline_vars <- c(spline_vars, bs_info[, 2])
            for (i in seq_len(nrow(bs_info))) {
              var_name <- bs_info[i, 2]
              spline_df[[var_name]] <- as.numeric(bs_info[i, 3])
              spline_type[[var_name]] <- "bs"
              if (!is.na(bs_info[i, 4])) {
                spline_degree[[var_name]] <- as.numeric(bs_info[i, 4])
              }
            }
          }
        }
        
        # Process pspline() splines
        if (length(pspline_matches) > 0) {
          pspline_info <- str_match_all(formula_str, pspline_pattern)[[1]]
          if (nrow(pspline_info) > 0) {
            spline_vars <- c(spline_vars, pspline_info[, 2])
            for (i in seq_len(nrow(pspline_info))) {
              var_name <- pspline_info[i, 2]
              spline_df[[var_name]] <- as.numeric(pspline_info[i, 3])
              spline_type[[var_name]] <- "pspline"
            }
          }
        }
        
        # Create configuration lists
        regular_vars <- setdiff(all_vars, spline_vars)
        covariates_list <- list("formula_combination" = c(regular_vars, spline_vars))
        spline_vars_list <- list("formula_combination" = spline_vars)
        
        fit_msm_models(train_data, train_crude_rates, 
                       covariates = covariates_list,
                       spline_vars = spline_vars_list,
                       spline_type = spline_type,
                       spline_df = spline_df,
                       spline_degree = spline_degree)
      } else {
        # No splines detected
        covariates_list <- list("formula_combination" = all_vars)
        fit_msm_models(train_data, train_crude_rates, 
                       covariates = covariates_list)
      }
    }
  } else {
    fit_msm_models(train_data, train_crude_rates)
  }
}

# Helper function to generate enhanced predictions
generate_enhanced_predictions <- function(test_outcomes, test_baseline, fitted_model, 
                                          covariate_formula, prediction_times, calibration_covariates) {
  
  # Base predictions (original functionality)
  base_predictions <- map_dfr(test_outcomes$deid_enc_id, function(patient_id) {
    baseline_covs <- test_baseline %>%
      filter(deid_enc_id == patient_id) %>%
      select(any_of(all.vars(covariate_formula %||% ~ 1))) %>%
      as.list()
    
    initial_state <- test_outcomes$initial_state[test_outcomes$deid_enc_id == patient_id]
    max_time <- max(prediction_times)
    
    generate_patient_predictions(patient_id, baseline_covs, initial_state, fitted_model, 
                                 covariate_formula, max_time)
  })
  
  # Time-specific predictions
  time_specific_predictions <- map_dfr(prediction_times, function(time_point) {
    map_dfr(test_outcomes$deid_enc_id, function(patient_id) {
      baseline_covs <- test_baseline %>%
        filter(deid_enc_id == patient_id) %>%
        select(any_of(all.vars(covariate_formula %||% ~ 1))) %>%
        as.list()
      
      initial_state <- test_outcomes$initial_state[test_outcomes$deid_enc_id == patient_id]
      
      pred_result <- generate_patient_predictions(patient_id, baseline_covs, initial_state, 
                                                  fitted_model, covariate_formula, time_point)
      pred_result$prediction_time <- time_point
      return(pred_result)
    })
  })
  
  # Covariate-specific predictions (if specified)
  covariate_specific_predictions <- if (!is.null(calibration_covariates)) {
    map_dfr(names(calibration_covariates), function(cov_set_name) {
      cov_values <- calibration_covariates[[cov_set_name]]
      
      map_dfr(test_outcomes$deid_enc_id, function(patient_id) {
        initial_state <- test_outcomes$initial_state[test_outcomes$deid_enc_id == patient_id]
        max_time <- max(prediction_times)
        
        pred_result <- generate_patient_predictions(patient_id, cov_values, initial_state, 
                                                    fitted_model, covariate_formula, max_time)
        pred_result$covariate_set <- cov_set_name
        return(pred_result)
      })
    })
  } else {
    NULL
  }
  
  return(list(
    base_predictions = base_predictions,
    time_specific_predictions = time_specific_predictions,
    covariate_specific_predictions = covariate_specific_predictions
  ))
}

# Helper function to generate predictions for a single patient
generate_patient_predictions <- function(patient_id, baseline_covs, initial_state, 
                                         fitted_model, covariate_formula, prediction_time) {
  
  state_names <- rownames(fitted_model$qmodel$qmatrix)
  state_mapping <- state_name_to_num(initial_state, fitted_model)
  initial_state_num <- state_mapping[initial_state]
  if (is.na(initial_state_num)) {
    warning(paste("Initial state", initial_state, "not found in model, using state 1"))
    initial_state_num <- 1
  }
  
  tryCatch({
    pmat <- if (length(baseline_covs) > 0 && !is.null(covariate_formula)) {
      pmatrix.msm(fitted_model, t = prediction_time, covariates = baseline_covs)
    } else {
      pmatrix.msm(fitted_model, t = prediction_time)
    }
    
    probs <- pmat[initial_state_num, ]
    
    # Enhanced LOS prediction (keeping existing logic)
    sojourn_result <- if (length(baseline_covs) > 0 && !is.null(covariate_formula)) {
      sojourn.msm(fitted_model, covariates = baseline_covs)
    } else {
      sojourn.msm(fitted_model)
    }
    
    # LOS calculations (keeping existing complex logic)
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

# Helper function to calculate base performance metrics
calculate_base_performance_metrics <- function(fold_num, n_train, n_test, model_converged,
                                               test_outcomes, predictions) {
  
  fold_data <- test_outcomes %>%
    left_join(predictions, by = "deid_enc_id")
  
  tryCatch({
    # MAE calculations (keeping existing logic)
    total_los_mae <- mean(abs(fold_data$total_los - fold_data$pred_total_los), na.rm = TRUE)
    total_los_mae_soj <- mean(abs(fold_data$total_los - fold_data$pred_total_los_soj), na.rm = TRUE)
    total_los_mae_vis_soj <- mean(abs(fold_data$total_los - fold_data$pred_total_los_vis_soj), na.rm = TRUE)
    
    days_severe_mae <- mean(abs(fold_data$days_severe - fold_data$pred_days_severe), na.rm = TRUE)
    days_severe_mae_soj <- mean(abs(fold_data$days_severe - fold_data$pred_days_severe_soj), na.rm = TRUE)
    days_severe_mae_vis_soj <- mean(abs(fold_data$days_severe - fold_data$pred_days_severe_vis_soj), na.rm = TRUE)
    
    # AUC calculations (keeping existing logic)
    death_auc <- calculate_auc_safe(fold_data$died, fold_data$pred_prob_death)
    severe_auc <- calculate_auc_safe(fold_data$ever_severe, fold_data$pred_prob_severe)
    
    # Calibration metrics (keeping existing logic)
    calibration_death <- calculate_calibration_metrics(fold_data$died, fold_data$pred_prob_death)
    
    # Brier score
    brier_death <- if (sum(!is.na(fold_data$pred_prob_death)) > 0) {
      mean((fold_data$died - fold_data$pred_prob_death)^2, na.rm = TRUE)
    } else NA
    
    data.frame(
      fold = fold_num,
      n_train = n_train,
      n_test = n_test,
      model_converged = model_converged,
      total_los_mae = total_los_mae,
      total_los_mae_soj = total_los_mae_soj,
      total_los_mae_vis_soj = total_los_mae_vis_soj,
      days_severe_mae = days_severe_mae,
      days_severe_mae_soj = days_severe_mae_soj,
      days_severe_mae_vis_soj = days_severe_mae_vis_soj,
      death_auc = as.numeric(death_auc),
      severe_auc = as.numeric(severe_auc),
      calibration_slope_death = calibration_death$slope,
      calibration_intercept_death = calibration_death$intercept,
      brier_death = brier_death,
      prediction_time = NA,
      subgroup_var = NA,
      subgroup_value = NA
    )
    
  }, error = function(e) {
    create_empty_fold_result(fold_num, n_train, n_test)
  })
}

# Helper function to calculate time-specific calibration
calculate_time_specific_calibration <- function(fold_num, n_train, n_test, model_converged,
                                                time_predictions, test_outcomes, test_baseline) {
  
  if (is.null(time_predictions) || nrow(time_predictions) == 0) {
    return(data.frame())
  }
  
  # Create observed outcomes at each time point
  time_specific_results <- map_dfr(unique(time_predictions$prediction_time), function(time_point) {
    
    # Filter predictions for this time point
    time_preds <- time_predictions %>%
      filter(prediction_time == time_point)
    
    # Create time-specific outcomes (died/severe by time_point)
    time_outcomes <- test_outcomes %>%
      left_join(time_preds, by = "deid_enc_id") %>%
      mutate(
        # For time-specific calibration, we need outcomes by specific time
        died_by_time = died,  # This would need actual time-to-event data
        severe_by_time = ever_severe  # This would need time-specific severe outcomes
      )
    
    if (nrow(time_outcomes) == 0) return(data.frame())
    
    # Calculate time-specific calibration
    death_auc_time <- calculate_auc_safe(time_outcomes$died_by_time, time_outcomes$pred_prob_death)
    severe_auc_time <- calculate_auc_safe(time_outcomes$severe_by_time, time_outcomes$pred_prob_severe)
    
    calibration_death_time <- calculate_calibration_metrics(time_outcomes$died_by_time, time_outcomes$pred_prob_death)
    
    brier_death_time <- if (sum(!is.na(time_outcomes$pred_prob_death)) > 0) {
      mean((time_outcomes$died_by_time - time_outcomes$pred_prob_death)^2, na.rm = TRUE)
    } else NA
    
    data.frame(
      fold = fold_num,
      n_train = n_train,
      n_test = n_test,
      model_converged = model_converged,
      total_los_mae = NA,
      total_los_mae_soj = NA,
      total_los_mae_vis_soj = NA,
      days_severe_mae = NA,
      days_severe_mae_soj = NA,
      days_severe_mae_vis_soj = NA,
      death_auc = as.numeric(death_auc_time),
      severe_auc = as.numeric(severe_auc_time),
      calibration_slope_death = calibration_death_time$slope,
      calibration_intercept_death = calibration_death_time$intercept,
      brier_death = brier_death_time,
      prediction_time = time_point,
      subgroup_var = NA,
      subgroup_value = NA
    )
  })
  
  return(time_specific_results)
}

# Helper function to calculate subgroup-specific calibration
calculate_subgroup_calibration <- function(fold_num, n_train, n_test, model_converged,
                                           predictions, test_outcomes, test_baseline, calibration_subgroups) {
  
  if (is.null(calibration_subgroups) || length(calibration_subgroups) == 0) {
    return(data.frame())
  }
  
  fold_data <- test_outcomes %>%
    left_join(predictions, by = "deid_enc_id") %>%
    left_join(test_baseline %>% select(deid_enc_id, any_of(calibration_subgroups)), by = "deid_enc_id")
  
  subgroup_results <- map_dfr(calibration_subgroups, function(subgroup_var) {
    
    if (!subgroup_var %in% names(fold_data)) {
      return(data.frame())
    }
    
    # Get unique values of subgroup variable
    subgroup_values <- unique(fold_data[[subgroup_var]])
    subgroup_values <- subgroup_values[!is.na(subgroup_values)]
    
    map_dfr(subgroup_values, function(subgroup_value) {
      
      subgroup_data <- fold_data %>%
        filter(.data[[subgroup_var]] == subgroup_value)
      
      if (nrow(subgroup_data) < 5) return(data.frame())  # Skip small subgroups
      
      # Calculate subgroup-specific metrics
      death_auc_subgroup <- calculate_auc_safe(subgroup_data$died, subgroup_data$pred_prob_death)
      calibration_death_subgroup <- calculate_calibration_metrics(subgroup_data$died, subgroup_data$pred_prob_death)
      
      brier_death_subgroup <- if (sum(!is.na(subgroup_data$pred_prob_death)) > 0) {
        mean((subgroup_data$died - subgroup_data$pred_prob_death)^2, na.rm = TRUE)
      } else NA
      
      data.frame(
        fold = fold_num,
        n_train = n_train,
        n_test = n_test,
        model_converged = model_converged,
        total_los_mae = NA,
        total_los_mae_soj = NA,
        total_los_mae_vis_soj = NA,
        days_severe_mae = NA,
        days_severe_mae_soj = NA,
        days_severe_mae_vis_soj = NA,
        death_auc = as.numeric(death_auc_subgroup),
        severe_auc = NA,
        calibration_slope_death = calibration_death_subgroup$slope,
        calibration_intercept_death = calibration_death_subgroup$intercept,
        brier_death = brier_death_subgroup,
        prediction_time = NA,
        subgroup_var = subgroup_var,
        subgroup_value = as.character(subgroup_value)
      )
    })
  })
  
  return(subgroup_results)
}

# Helper function to calculate calibration metrics
calculate_calibration_metrics <- function(outcome, prediction) {
  if (var(outcome, na.rm = TRUE) > 0 && sum(!is.na(prediction)) > 0) {
    tryCatch({
      # Convert probabilities to logit scale for calibration regression
      pred_logit <- qlogis(pmax(pmin(prediction, 0.999), 0.001))
      cal_model <- glm(outcome ~ pred_logit, family = binomial())
      list(slope = coef(cal_model)[2], intercept = coef(cal_model)[1])
    }, error = function(e) {
      list(slope = NA, intercept = NA)
    })
  } else {
    list(slope = NA, intercept = NA)
  }
}

## Residual diagnostics -----------------------------------------------------

### Transition residuals ----------------------------------------------------

#' Calculate transition residuals for nested model structure
#' @param fitted_msm_models Nested list of fitted models: fitted_models[[model_structure]][[formula]]
#' @param patient_data Patient data used to fit the models
#' @param residual_type Type of residuals ("deviance", "pearson", "raw")
#' @param debug Logical, whether to print debug information
#' @return Data frame with transition residuals for all model/formula combinations
calculate_transition_residuals_old <- function(fitted_msm_models, patient_data, 
                                           residual_type = "deviance",
                                           debug = TRUE) {
  
  validate_model_structure(fitted_msm_models)
  
  if (debug) cat("=== DEBUG: Starting transition residuals calculation for nested models ===\n")
  
  all_residuals <- list()
  
  # Loop through model structures
  for (model_structure in names(fitted_msm_models)) {
    # Loop through formulas within each model structure
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
      fitted_model <- extract_fitted_model(fitted_msm_models, model_structure, formula_name, "in tidy_msm_qmatrix")
      if (is.null(fitted_model)) next
      
      if (is.null(fitted_model) || !inherits(fitted_model, "msm")) {
        if (debug) cat("Invalid model object for:", model_structure, formula_name, "\n")
        next
      }
      
      # Filter patient data for this model structure
      model_data <- patient_data %>% filter(model == model_structure)
      
      if (debug) {
        cat("Model data dimensions:", nrow(model_data), "x", ncol(model_data), "\n")
        cat("Unique patients:", length(unique(model_data$deid_enc_id)), "\n")
        cat("Model convergence:", fitted_model$opt$convergence == 0, "\n")
        cat("Model has covariates:", !is.null(fitted_model$covariates), "\n")
      }
      
      # Get observed transitions
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
      
      observed_transitions <- standardize_state_columns(observed_transitions, fitted_model)
      
      if (debug) {
        cat("Observed transitions created. Dimensions:", nrow(observed_transitions), "x", ncol(observed_transitions), "\n")
        cat("Unique transitions:", length(unique(paste(observed_transitions$from_state, "->", observed_transitions$to_state))), "\n")
      }
      
      if (nrow(observed_transitions) == 0) {
        if (debug) cat("No transitions found for:", model_structure, formula_name, "\n")
        next
      }
      
      # Calculate expected probabilities using tidy_msm_pmats
      unique_time_diffs <- unique(observed_transitions$time_diff)
      
      # Create temporary nested structure for tidy_msm_pmats
      temp_nested_model <- list()
      temp_nested_model[[model_structure]] <- list()
      temp_nested_model[[model_structure]][[formula_name]] <- list(
        fitted_model = fitted_model,
        status = "converged"
      )
      
      expected_probs_df <- tryCatch({
        tidy_msm_pmats(
          fitted_msm_models = temp_nested_model,
          t_values = unique_time_diffs,
          covariates_list = NULL
        )
      }, error = function(e) {
        if (debug) cat("tidy_msm_pmats ERROR:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(expected_probs_df) || nrow(expected_probs_df) == 0) {
        if (debug) cat("No expected probabilities calculated for:", model_structure, formula_name, "\n")
        next
      }
      
      if (debug) {
        cat("Expected probs calculated. Dimensions:", nrow(expected_probs_df), "x", ncol(expected_probs_df), "\n")
        cat("Column names in expected_probs_df:", paste(names(expected_probs_df), collapse = ", "), "\n")
        cat("Sample of expected_probs_df:\n")
        print(head(expected_probs_df))
      }
      
      # Join expected probabilities with observed transitions
      # Convert state numbers to state names to match observed transitions
      transition_residuals <- observed_transitions %>%
        left_join(
          expected_probs_df %>% 
            select(state, tostate, estimate, t_value, statename, tostatename) %>%
            rename(from_state = statename, to_state = tostatename, expected_prob = estimate, time_diff = t_value) %>%
            select(-state, -tostate),  # Remove the numeric columns
          by = c("from_state", "to_state", "time_diff")
        )
      
      transition_residuals <- standardize_state_columns(transition_residuals, fitted_model)
      
      if (debug) {
        cat("After joining - Non-NA expected probs:", sum(!is.na(transition_residuals$expected_prob)), "/", nrow(transition_residuals), "\n")
      }
      
      # Filter and calculate residuals
      transition_residuals <- transition_residuals %>%
        filter(!is.na(expected_prob), expected_prob > 0) %>%
        mutate(
          # Raw residual (observed = 1 for actual transitions)
          raw_residual = 1 - expected_prob,
          
          # Pearson residual
          pearson_residual = raw_residual / sqrt(expected_prob * (1 - expected_prob)),
          
          # Deviance residual with safety check
          deviance_residual = sign(raw_residual) * sqrt(pmax(-2 * log(pmax(expected_prob, 1e-10)), 0)),
          
          # Transition identifier
          transition = paste(from_state, "->", to_state),
          
          # Add model metadata
          model = model_structure,
          formula = formula_name,
          covariate = if (formula_name == "~ 1") {
            "None"
          } else {
            tryCatch({
              vars_in_formula <- all.vars(as.formula(formula_name))
              paste(vars_in_formula, collapse = ", ")
            }, error = function(e) {
              "Error parsing formula"
            })
          }
        )
      
      if (nrow(transition_residuals) == 0) {
        if (debug) cat("No valid transitions after filtering for:", model_structure, formula_name, "\n")
        next
      }
      
      # Select and standardize residuals
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
      }
      
      # Store results
      combo_key <- paste(model_structure, formula_name, sep = "___")
      all_residuals[[combo_key]] <- transition_residuals
    }
  }
  
  # Combine all results
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
    }
  }
  
  return(final_residuals)
}

# Optimized calculate_transition_residuals function
calculate_transition_residuals <- function(fitted_msm_models, patient_data, 
                                           residual_type = "deviance",
                                           debug = FALSE) {
  
  validate_model_structure(fitted_msm_models)
  
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
      
      # Fix: NOW standardize the state columns after observed_transitions is created
      observed_transitions <- standardize_state_columns(observed_transitions, fitted_model)
      
      if (debug) {
        cat("Observed transitions created. Dimensions:", nrow(observed_transitions), "x", ncol(observed_transitions), "\n")
        cat("Unique transitions:", length(unique(paste(observed_transitions$from_state, "->", observed_transitions$to_state))), "\n")
      }
      
      if (nrow(observed_transitions) == 0) {
        if (debug) cat("No transitions found for:", model_structure, formula_name, "\n")
        next
      }
      
      # Calculate expected probabilities more efficiently
      unique_time_diffs <- unique(observed_transitions$time_diff)
      
      # Create minimal temporary structure for tidy_msm_pmats
      temp_nested_model <- list()
      temp_nested_model[[model_structure]] <- list()
      temp_nested_model[[model_structure]][[formula_name]] <- list(
        fitted_model = fitted_model,
        status = "converged"
      )
      
      expected_probs_df <- tryCatch({
        # Calculate probabilities for unique time differences only
        tidy_msm_pmats(
          fitted_msm_models = temp_nested_model,
          t_values = unique_time_diffs,
          covariates_list = NULL
        )
      }, error = function(e) {
        if (debug) cat("tidy_msm_pmats ERROR:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(expected_probs_df) || nrow(expected_probs_df) == 0) {
        if (debug) cat("No expected probabilities calculated for:", model_structure, formula_name, "\n")
        next
      }
      
      if (debug) {
        cat("Expected probs calculated. Dimensions:", nrow(expected_probs_df), "x", ncol(expected_probs_df), "\n")
      }
      
      # Join and calculate residuals more efficiently
      transition_residuals <- observed_transitions %>%
        left_join(
          expected_probs_df %>% 
            select(state, tostate, estimate, t_value, statename, tostatename) %>%
            rename(from_state = statename, to_state = tostatename, expected_prob = estimate, time_diff = t_value) %>%
            select(-state, -tostate),
          by = c("from_state", "to_state", "time_diff")
        ) %>%
        filter(!is.na(expected_prob), expected_prob > 0) %>%
        mutate(
          # Calculate residuals
          raw_residual = 1 - expected_prob,
          pearson_residual = raw_residual / sqrt(expected_prob * (1 - expected_prob)),
          deviance_residual = sign(raw_residual) * sqrt(pmax(-2 * log(pmax(expected_prob, 1e-10)), 0)),
          transition = paste(from_state, "->", to_state),
          model = model_structure,
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
      }
      
      # Store results
      combo_key <- paste(model_structure, formula_name, sep = "___")
      all_residuals[[combo_key]] <- transition_residuals
      
      # Clean up to free memory
      rm(model_data, observed_transitions, expected_probs_df, transition_residuals, temp_nested_model)
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
    }
  }
  
  return(final_residuals)
}

# Covariate effects -------------------------------------------------------

## Transition-specific effects --------------------------------------------

#' Fit models with transition-specific covariate effects
#' @param patient_data Patient data with transition trend categories
#' @param crude_rates List of crude rates
#' @param covariates List structure compatible with fit_msm_models (or vector for backwards compatibility)
#' @param constraint_configs List specifying transition-specific constraints
#' @return Nested list of fitted models compatible with enhanced structure
fit_transition_specific_models <- function(patient_data, crude_rates, covariates,
                                           constraint_configs = list(
                                             "constant" = "default",
                                             "by_trend" = "by_trend",
                                             "by_transition" = NULL
                                           )) {
  
  transition_models <- list()
  
  # Convert covariates to list format if needed
  if (!is.list(covariates)) {
    covariates_list <- setNames(list(covariates), "main_covariates")
  } else {
    covariates_list <- covariates
  }
  
  for (cov_set_name in names(covariates_list)) {
    cov_set <- covariates_list[[cov_set_name]]
    
    cat("Fitting transition-specific models for covariate set:", cov_set_name, "\n")
    
    # Fit models by constraint type
    for (constraint_name in names(constraint_configs)) {
      constraint_type <- constraint_configs[[constraint_name]]
      
      cat("  - Fitting constraint type:", constraint_name, "\n")
      
      if (constraint_name == "constant") {
        # Standard model with constant effects across transitions
        model_key <- paste(cov_set_name, "constant", sep = "_")
        transition_models[[model_key]] <- 
          fit_msm_models(
            patient_data = patient_data, 
            crude_rates = crude_rates, 
            covariates = setNames(list(cov_set), cov_set_name),
            nest_name = model_key
          )
        
      } else if (constraint_name == "by_trend") {
        # Effects vary by trend type (requires custom constraint matrix)
        # This is a placeholder - actual implementation would need
        # careful constraint matrix specification based on your transition structure
        model_key <- paste(cov_set_name, "by_trend", sep = "_")
        
        # For now, fit as unconstrained and add metadata about intended constraint
        models_result <- fit_msm_models(
          patient_data = patient_data, 
          crude_rates = crude_rates, 
          covariates = setNames(list(cov_set), cov_set_name),
          nest_name = model_key
        )
        
        # Add constraint metadata
        if (!is.null(models_result) && length(models_result) > 0) {
          for (nested_model_name in names(models_result)) {
            if (!is.null(models_result[[nested_model_name]])) {
              for (formula_name in names(models_result[[nested_model_name]])) {
                if (is.list(models_result[[nested_model_name]][[formula_name]])) {
                  models_result[[nested_model_name]][[formula_name]]$constraint_type <- "by_trend"
                }
              }
            }
          }
        }
        
        transition_models[[model_key]] <- models_result
        
      } else if (constraint_name == "by_transition") {
        # Fully flexible - different effect for each transition
        model_key <- paste(cov_set_name, "by_transition", sep = "_")
        
        # Fit unconstrained model (constraint = NULL in msm)
        # Note: This requires modifying fit_msm_models to accept constraint parameter
        models_result <- fit_msm_models(
          patient_data = patient_data, 
          crude_rates = crude_rates, 
          covariates = setNames(list(cov_set), cov_set_name),
          nest_name = model_key
        )
        
        # Add constraint metadata
        if (!is.null(models_result) && length(models_result) > 0) {
          for (nested_model_name in names(models_result)) {
            if (!is.null(models_result[[nested_model_name]])) {
              for (formula_name in names(models_result[[nested_model_name]])) {
                if (is.list(models_result[[nested_model_name]][[formula_name]])) {
                  models_result[[nested_model_name]][[formula_name]]$constraint_type <- "by_transition"
                }
              }
            }
          }
        }
        
        transition_models[[model_key]] <- models_result
      }
    }
  }
  
  return(transition_models)
}

## Non-linear effects -----------------------------------------------------

#' Test for non-linear relationships using splines with enhanced fit_msm_models
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param covariate Name of continuous covariate
#' @param spline_df Degrees of freedom for spline (default: 3)
#' @param spline_type Type of spline ("ns", "bs", "pspline")
#' @param comparison_method Method for model comparison ("AIC", "LRT", "both")
#' @return List with linear and spline models plus comparison
test_nonlinear_relationship <- function(patient_data, crude_rates, covariate, 
                                        spline_df = 3, spline_type = "ns",
                                        comparison_method = "both") {
  
  cat("Testing nonlinear relationship for:", covariate, "\n")
  cat("Using spline type:", spline_type, "with df =", spline_df, "\n")
  
  # Validate that covariate exists in data
  if (!covariate %in% names(patient_data)) {
    stop(paste("Covariate", covariate, "not found in patient_data"))
  }
  
  # Check for sufficient variation in covariate
  covariate_summary <- patient_data %>%
    group_by(model) %>%
    summarise(
      n_unique = length(unique(.data[[covariate]])),
      n_missing = sum(is.na(.data[[covariate]])),
      .groups = "drop"
    )
  
  models_to_fit <- covariate_summary %>%
    filter(n_unique >= 5, n_missing < 0.8 * length(unique(patient_data$deid_enc_id))) %>%
    pull(model)
  
  if (length(models_to_fit) == 0) {
    warning(paste("Insufficient variation in", covariate, "across all models"))
    return(NULL)
  }
  
  # Fit linear model
  cat("Fitting linear models...\n")
  linear_models <- fit_msm_models(
    patient_data = patient_data, 
    crude_rates = crude_rates, 
    covariates = list("linear" = covariate),
    nest_name = paste0(covariate, "_linear")
  )
  
  # Fit spline models
  cat("Fitting spline models...\n")
  spline_models <- fit_msm_models(
    patient_data = patient_data, 
    crude_rates = crude_rates, 
    covariates = list("spline" = covariate),
    spline_vars = list("spline" = covariate),
    spline_df = setNames(list(spline_df), covariate),
    spline_type = setNames(list(spline_type), covariate),
    nest_name = paste0(covariate, "_spline")
  )
  
  # Extract model objects for comparison
  comparison_results <- list()
  
  for (model_structure in models_to_fit) {
    # Extract linear model
    linear_model_obj <- NULL
    if (!is.null(linear_models) && !is.null(linear_models[[paste0(covariate, "_linear")]])) {
      linear_entry <- linear_models[[paste0(covariate, "_linear")]][[model_structure]]
      if (!is.null(linear_entry)) {
        linear_formula_key <- names(linear_entry)[1]
        if (!is.null(linear_formula_key)) {
          linear_model_entry <- linear_entry[[linear_formula_key]]
          linear_model_obj <- if (is.list(linear_model_entry) && !is.null(linear_model_entry$fitted_model)) {
            linear_model_entry$fitted_model
          } else {
            linear_model_entry
          }
        }
      }
    }
    
    # Extract spline model
    spline_model_obj <- NULL
    if (!is.null(spline_models) && !is.null(spline_models[[paste0(covariate, "_spline")]])) {
      spline_entry <- spline_models[[paste0(covariate, "_spline")]][[model_structure]]
      if (!is.null(spline_entry)) {
        spline_formula_key <- names(spline_entry)[1]
        if (!is.null(spline_formula_key)) {
          spline_model_entry <- spline_entry[[spline_formula_key]]
          spline_model_obj <- if (is.list(spline_model_entry) && !is.null(spline_model_entry$fitted_model)) {
            spline_model_entry$fitted_model
          } else {
            spline_model_entry
          }
        }
      }
    }
    
    # Compare models if both fitted successfully
    if (!is.null(linear_model_obj) && !is.null(spline_model_obj)) {
      comparison_result <- compare_linear_vs_spline(
        linear_model_obj, spline_model_obj, model_structure, covariate, 
        spline_df, comparison_method
      )
      comparison_results[[model_structure]] <- comparison_result
    } else {
      warning(paste("Could not extract models for comparison in", model_structure))
      comparison_results[[model_structure]] <- data.frame(
        model = model_structure,
        covariate = covariate,
        linear_converged = !is.null(linear_model_obj),
        spline_converged = !is.null(spline_model_obj),
        comparison_possible = FALSE
      )
    }
  }
  
  # Combine comparison results
  comparison_df <- bind_rows(comparison_results)
  
  return(list(
    linear_models = linear_models,
    spline_models = spline_models,
    comparison = comparison_df,
    covariate = covariate,
    spline_df = spline_df,
    spline_type = spline_type,
    summary = list(
      models_tested = length(models_to_fit),
      models_compared = sum(comparison_df$comparison_possible, na.rm = TRUE),
      prefer_spline = sum(comparison_df$prefer_spline, na.rm = TRUE)
    )
  ))
}

#' Compare linear vs spline models
compare_linear_vs_spline <- function(linear_model, spline_model, model_name, covariate, 
                                     spline_df, comparison_method) {
  
  # Check convergence
  linear_converged <- is.null(linear_model$opt$convergence) || linear_model$opt$convergence == 0
  spline_converged <- is.null(spline_model$opt$convergence) || spline_model$opt$convergence == 0
  
  if (!linear_converged || !spline_converged) {
    return(data.frame(
      model = model_name,
      covariate = covariate,
      linear_converged = linear_converged,
      spline_converged = spline_converged,
      comparison_possible = FALSE,
      linear_AIC = ifelse(linear_converged, AIC(linear_model), NA),
      spline_AIC = ifelse(spline_converged, AIC(spline_model), NA),
      delta_AIC = NA,
      lrt_stat = NA,
      lrt_df = spline_df - 1,
      lrt_pval = NA,
      prefer_spline = NA
    ))
  }
  
  # Calculate AICs
  linear_aic <- AIC(linear_model)
  spline_aic <- AIC(spline_model)
  delta_aic <- spline_aic - linear_aic
  
  # Likelihood ratio test
  lrt_result <- tryCatch({
    test_result <- lrtest.msm(linear_model, spline_model)
    c(test_result$statistic, test_result$df, test_result$p.value)
  }, error = function(e) {
    warning(paste("LRT failed for", model_name, ":", e$message))
    c(NA, NA, NA)
  })
  
  # Decision criteria
  prefer_spline_aic <- delta_aic < -2
  prefer_spline_lrt <- !is.na(lrt_result[3]) && lrt_result[3] < 0.05
  
  prefer_spline <- if (comparison_method == "AIC") {
    prefer_spline_aic
  } else if (comparison_method == "LRT") {
    prefer_spline_lrt
  } else {  # both
    prefer_spline_aic && prefer_spline_lrt
  }
  
  return(data.frame(
    model = model_name,
    covariate = covariate,
    linear_converged = linear_converged,
    spline_converged = spline_converged,
    comparison_possible = TRUE,
    linear_AIC = linear_aic,
    spline_AIC = spline_aic,
    delta_AIC = delta_aic,
    lrt_stat = lrt_result[1],
    lrt_df = lrt_result[2],
    lrt_pval = lrt_result[3],
    prefer_spline_aic = prefer_spline_aic,
    prefer_spline_lrt = prefer_spline_lrt,
    prefer_spline = prefer_spline
  ))
}

## Multivariate model selection ------------------------------------------------

#' Fit multivariate models with forward selection using enhanced structure
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param candidate_covariates Vector of candidate covariates
#' @param method Selection method ("forward", "backward", "stepwise")
#' @param alpha_enter P-value threshold for entering (default: 0.05)
#' @param alpha_remove P-value threshold for removal (default: 0.10)
#' @param max_variables Maximum number of variables to include (default: 5)
#' @return Nested list of selected models with metadata
multivariate_selection <- function(patient_data, crude_rates, candidate_covariates, 
                                   method = "forward", alpha_enter = 0.05, alpha_remove = 0.10,
                                   max_variables = 5) {
  
  # Validate inputs
  missing_covariates <- setdiff(candidate_covariates, names(patient_data))
  if (length(missing_covariates) > 0) {
    warning(paste("Covariates not found in data:", paste(missing_covariates, collapse = ", ")))
    candidate_covariates <- intersect(candidate_covariates, names(patient_data))
  }
  
  if (length(candidate_covariates) == 0) {
    stop("No valid candidate covariates provided")
  }
  
  selected_models <- list()
  model_names <- names(crude_rates)
  
  for (model_name in model_names) {
    cat("Running", method, "selection for", model_name, "\n")
    
    model_data <- patient_data %>% filter(model == model_name)
    if (nrow(model_data) == 0) {
      warning(paste("No data for model", model_name))
      next
    }
    
    # Check covariate availability for this model
    available_covariates <- candidate_covariates[
      sapply(candidate_covariates, function(cov) {
        sum(!is.na(model_data[[cov]])) > 10  # At least 10 non-missing values
      })
    ]
    
    if (length(available_covariates) == 0) {
      warning(paste("No covariates with sufficient data for model", model_name))
      next
    }
    
    if (method == "forward") {
      selection_result <- forward_selection(
        model_data, crude_rates[[model_name]], available_covariates, 
        alpha_enter, max_variables
      )
    } else if (method == "backward") {
      selection_result <- enhanced_backward_elimination(
        model_data, crude_rates[[model_name]], available_covariates, alpha_remove
      )
    } else if (method == "stepwise") {
      selection_result <- enhanced_stepwise_selection(
        model_data, crude_rates[[model_name]], available_covariates, 
        alpha_enter, alpha_remove, max_variables
      )
    } else {
      stop(paste("Unknown selection method:", method))
    }
    
    # Store results in nested structure
    if (!is.null(selection_result)) {
      selected_models[[model_name]] <- selection_result
    }
  }
  
  return(selected_models)
}

#' Forward selection compatible with fit_msm_models structure
forward_selection <- function(model_data, crude_rate, candidate_covariates, 
                              alpha_enter, max_variables) {
  
  selected_covariates <- character(0)
  remaining_covariates <- candidate_covariates
  selection_history <- list()
  step <- 0
  
  # Fit base model (no covariates)
  base_models <- fit_msm_models(
    patient_data = model_data,
    crude_rates = setNames(list(crude_rate), unique(model_data$model)[1])
  )
  
  # Extract base model object
  model_name <- unique(model_data$model)[1]
  formula_key <- names(base_models[[model_name]])[1]
  base_model_obj <- extract_fitted_model(base_models, model_name, formula_key, "in forward_selection base")  
  if (is.null(base_model_obj)) {
    warning("Base model failed to fit")
    return(NULL)
  }
  
  current_aic <- AIC(base_model_obj)
  current_model <- base_model_obj
  
  cat("Base model AIC:", round(current_aic, 2), "\n")
  
  while (length(remaining_covariates) > 0 && length(selected_covariates) < max_variables) {
    step <- step + 1
    best_covariate <- NULL
    best_pvalue <- 1
    best_aic <- current_aic
    best_model <- NULL
    
    cat("Step", step, ": Testing", length(remaining_covariates), "covariates\n")
    
    for (cov in remaining_covariates) {
      test_covariates <- c(selected_covariates, cov)
      
      cat("  Testing covariate:", cov, "\n")
      
      # Fit model with test covariates
      test_models <- fit_msm_models(
        patient_data = model_data,
        crude_rates = setNames(list(crude_rate), unique(model_data$model)[1]),
        covariates = list("test" = test_covariates)
      )
      
      test_formula_key <- names(test_models[[model_name]])[1]  # Get the actual formula key for this test model
      test_model_obj <- extract_fitted_model(test_models, model_name, test_formula_key, "in forward_selection test")
      
      if (is.null(test_model_obj)) {
        cat("    Model failed to fit\n")
        next
      }
      
      test_aic <- AIC(test_model_obj)
      cat("    Test AIC:", round(test_aic, 2), "\n")
      
      # Likelihood ratio test
      lrt_result <- tryCatch({
        lrt <- lrtest.msm(current_model, test_model_obj)
        
        # Extract from matrix structure
        statistic <- lrt[1, 1]  # First row, first column ("-2 log LR")
        df <- lrt[1, 2]         # First row, second column ("df") 
        pvalue <- lrt[1, 3]     # First row, third column ("p")
        
        cat("    LRT statistic:", round(statistic, 4), "df:", df, "p-value:", round(pvalue, 6), "\n")
        c(statistic, df, pvalue)
      }, error = function(e) {
        cat("    LRT failed:", e$message, "\n")
        c(NA, NA, 1)
      })
      
      if (!is.na(lrt_result[3]) && lrt_result[3] < best_pvalue && test_aic < best_aic) {
        best_covariate <- cov
        best_pvalue <- lrt_result[3]
        best_aic <- test_aic
        best_model <- test_model_obj
      }
    }
    
    if (!is.null(best_covariate) && best_pvalue < alpha_enter) {
      selected_covariates <- c(selected_covariates, best_covariate)
      remaining_covariates <- setdiff(remaining_covariates, best_covariate)
      current_aic <- best_aic
      current_model <- best_model
      
      selection_history[[step]] <- list(
        action = "add",
        variable = best_covariate,
        pvalue = best_pvalue,
        aic = best_aic,
        selected_vars = selected_covariates
      )
      
      cat("  Added:", best_covariate, "(p =", round(best_pvalue, 4), ", AIC =", round(best_aic, 2), ")\n")
    } else {
      cat("  No variables meet entry criterion (alpha =", alpha_enter, ")\n")
      break
    }
  }
  
  # Fit final model with selected covariates
  if (length(selected_covariates) > 0) {
    final_models <- fit_msm_models(
      patient_data = model_data,
      crude_rates = setNames(list(crude_rate), unique(model_data$model)[1]),
      covariates = list("final" = selected_covariates)
    )
    final_formula_key <- names(final_models[[model_name]])[1]
    final_model_obj <- extract_fitted_model(final_models, model_name, final_formula_key, "final model")
  } else {
    final_models <- base_models
    final_model_obj <- base_model_obj
  }
  
  return(list(
    final_models = final_models,
    selected_covariates = selected_covariates,
    final_aic = current_aic,
    selection_history = selection_history,
    method = "forward",
    alpha_enter = alpha_enter,
    n_steps = step
  ))
}

# Time-varying models -----------------------------------------------------

## Model drifts over time -------------------------------------------------

#' Fit models with time-dependent transition rates
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param time_covariates Vector of time covariate specifications
#' @param time_term Which time variable to use ("days_since_entry" or "calendar_time")
#' @param spline_df Degrees of freedom for spline models (default: 3)
#' @param spline_type Type of spline for time trends ("ns", "bs", "pspline")
#' @param piecewise_breakpoints Vector of breakpoints for piecewise models
#' @param constraint Constraint specification for time effects across transitions
#' @return Nested list of fitted time-varying models
fit_time_varying_models <- function(patient_data, crude_rates, 
                                    time_covariates = c("linear", "spline", "piecewise"),
                                    time_term = "days_since_entry",
                                    spline_df = 3,
                                    spline_type = "ns",
                                    piecewise_breakpoints = NULL,
                                    constraint = "default") {
  
  # Validate time_term parameter
  if (!time_term %in% c("days_since_entry", "calendar_time")) {
    stop("time_term must be either 'days_since_entry' or 'calendar_time'")
  }
  
  # Select the appropriate time variable
  if (time_term == "days_since_entry") {
    time_var <- "DaysSinceEntry"
    if (!time_var %in% names(patient_data)) {
      stop("DaysSinceEntry column not found in patient_data")
    }
  } else {
    time_var <- "CalendarTime"
    if (!time_var %in% names(patient_data)) {
      stop("CalendarTime column not found in patient_data. Please run add_calendar_time() first.")
    }
  }
  
  # Validate that there's sufficient time variation
  time_summary <- patient_data %>%
    group_by(model) %>%
    summarise(
      min_time = min(.data[[time_var]], na.rm = TRUE),
      max_time = max(.data[[time_var]], na.rm = TRUE),
      n_unique_times = length(unique(.data[[time_var]])),
      .groups = "drop"
    )
  
  insufficient_models <- time_summary %>%
    filter(n_unique_times < 5 | (max_time - min_time) < 2) %>%
    pull(model)
  
  if (length(insufficient_models) > 0) {
    warning(paste("Insufficient time variation for models:", 
                  paste(insufficient_models, collapse = ", ")))
  }
  
  time_models <- list()
  
  for (time_type in time_covariates) {
    cat("Fitting", time_type, "time-varying model using", time_term, "\n")
    
    if (time_type == "linear") {
      # Linear time trend
      model_key <- paste0("time_linear_", time_term)
      
      time_models[[model_key]] <- fit_msm_models(
        patient_data = patient_data, 
        crude_rates = crude_rates, 
        covariates = list("linear_time" = time_var),
        constraint = constraint,
        nest_name = model_key
      )
      
    } else if (time_type == "spline") {
      # Spline time trend using enhanced spline functionality
      model_key <- paste0("time_spline_", time_term)
      
      time_models[[model_key]] <- fit_msm_models(
        patient_data = patient_data,
        crude_rates = crude_rates,
        covariates = list("spline_time" = time_var),
        spline_vars = list("spline_time" = time_var),
        spline_df = setNames(list(spline_df), time_var),
        spline_type = setNames(list(spline_type), time_var),
        constraint = constraint,
        nest_name = model_key
      )
      
    } else if (time_type == "piecewise") {
      # Piecewise constant models
      model_key <- paste0("time_piecewise_", time_term)
      
      # Determine breakpoints
      if (is.null(piecewise_breakpoints)) {
        if (time_term == "days_since_entry") {
          breakpoints <- c(7, 14)  # Early, middle, late periods
        } else {
          # For calendar time, use quantiles or predefined dates
          breakpoints <- quantile(patient_data[[time_var]], c(0.33, 0.67), na.rm = TRUE)
        }
      } else {
        breakpoints <- piecewise_breakpoints
      }
      
      # Create time periods
      patient_data_piece <- patient_data %>%
        mutate(
          time_period = cut(.data[[time_var]], 
                            breaks = c(-Inf, breakpoints, Inf),
                            labels = paste0("period_", seq_len(length(breakpoints) + 1)),  # Fixed this line
                            include.lowest = TRUE)
        ) %>%
        # Convert to character to avoid factor issues in msm
        mutate(time_period = as.character(time_period))
      
      time_models[[model_key]] <- fit_msm_models(
        patient_data = patient_data_piece,
        crude_rates = crude_rates,
        covariates = list("piecewise_time" = "time_period"),
        constraint = constraint,
        nest_name = model_key
      )
      
    } else if (time_type == "quadratic") {
      # Quadratic time trend (additional option)
      model_key <- paste0("time_quadratic_", time_term)
      
      # Create quadratic term
      patient_data_quad <- patient_data %>%
        mutate(
          !!paste0(time_var, "_quad") := (.data[[time_var]])^2
        )
      
      quad_covariates <- c(time_var, paste0(time_var, "_quad"))
      
      time_models[[model_key]] <- fit_msm_models(
        patient_data = patient_data_quad,
        crude_rates = crude_rates,
        covariates = list("quadratic_time" = quad_covariates),
        constraint = constraint,
        nest_name = model_key
      )
      
    } else if (time_type == "interaction") {
      # Time interactions with other covariates (if available)
      # This requires specifying which covariate to interact with
      warning("Interaction models require additional covariate specification - skipping")
      next
      
    } else {
      warning(paste("Unknown time covariate type:", time_type, "- skipping"))
      next
    }
  }
  
  # # Add model comparison if multiple types fitted
  # if (length(time_models) > 1) {
  #   time_models$comparison <- compare_time_models(time_models, time_term)
  # }
  
  return(time_models)
}


#' Fit models with calendar time effects (COVID waves, seasonal patterns)
#' @param patient_data Patient data with CalendarTime variable
#' @param crude_rates List of crude rates
#' @param wave_periods Named list of COVID wave periods
#' @param seasonal Include seasonal effects (default: FALSE)
#' @param constraint Constraint specification
#' @return Nested list of calendar time models
fit_calendar_time_models <- function(patient_data, crude_rates, 
                                     wave_periods = list(
                                       "wave1" = c(0, 100),
                                       "wave2" = c(200, 350),
                                       "wave3" = c(450, 600)
                                     ),
                                     seasonal = FALSE,
                                     constraint = "default") {
  
  # Validate CalendarTime exists
  if (!"CalendarTime" %in% names(patient_data)) {
    stop("CalendarTime column not found. Please run add_calendar_time() first.")
  }
  
  calendar_models <- list()
  
  # Wave-specific models
  if (!is.null(wave_periods)) {
    cat("Fitting COVID wave models\n")
    
    # Create wave indicators
    patient_data_waves <- patient_data
    for (wave_name in names(wave_periods)) {
      wave_range <- wave_periods[[wave_name]]
      patient_data_waves[[paste0("wave_", wave_name)]] <- 
        as.numeric(patient_data_waves$CalendarTime >= wave_range[1] & 
                     patient_data_waves$CalendarTime <= wave_range[2])
    }
    
    # Create wave factor
    patient_data_waves$covid_wave <- "baseline"
    for (wave_name in names(wave_periods)) {
      wave_range <- wave_periods[[wave_name]]
      patient_data_waves$covid_wave[
        patient_data_waves$CalendarTime >= wave_range[1] & 
          patient_data_waves$CalendarTime <= wave_range[2]
      ] <- wave_name
    }
    
    calendar_models[["covid_waves"]] <- fit_msm_models(
      patient_data = patient_data_waves,
      crude_rates = crude_rates,
      covariates = list("waves" = "covid_wave"),
      constraint = constraint,
      nest_name = "covid_waves"
    )
  }
  
  # Seasonal models
  if (seasonal) {
    cat("Fitting seasonal models\n")
    
    # Create seasonal terms (assuming CalendarTime is days since some reference)
    patient_data_seasonal <- patient_data %>%
      mutate(
        # Convert to approximately monthly cycles
        month_cycle = (CalendarTime %% 365.25) / 30.44,
        season_sin = sin(2 * pi * month_cycle / 12),
        season_cos = cos(2 * pi * month_cycle / 12)
      )
    
    calendar_models[["seasonal"]] <- fit_msm_models(
      patient_data = patient_data_seasonal,
      crude_rates = crude_rates,
      covariates = list("seasonal" = c("season_sin", "season_cos")),
      constraint = constraint,
      nest_name = "seasonal"
    )
  }
  
  return(calendar_models)
}

## Time-varying covariate effects ------------------------------------------

#' Test for time-varying effects of specific covariates
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param covariate Covariate to test for time-varying effects
#' @param time_var Time variable ("DaysSinceEntry" or "CalendarTime")
#' @param interaction_types Types of interactions to test
#' @return List with constant and time-varying models plus comparison
test_time_varying_covariate_effects <- function(patient_data, crude_rates, covariate,
                                                time_var = "DaysSinceEntry",
                                                interaction_types = c("linear", "spline")) {
  
  # Validate inputs
  if (!covariate %in% names(patient_data)) {
    stop(paste("Covariate", covariate, "not found in patient_data"))
  }
  
  if (!time_var %in% names(patient_data)) {
    stop(paste("Time variable", time_var, "not found in patient_data"))
  }
  
  cat("Testing time-varying effects for:", covariate, "\n")
  
  time_varying_models <- list()
  
  # Constant effect model (baseline)
  time_varying_models[["constant"]] <- fit_msm_models(
    patient_data = patient_data,
    crude_rates = crude_rates,
    covariates = list("constant" = covariate),
    nest_name = paste0(covariate, "_constant")
  )
  
  for (interaction_type in interaction_types) {
    
    if (interaction_type == "linear") {
      # Linear interaction with time
      interaction_data <- patient_data %>%
        mutate(
          !!paste0(covariate, "_time_interaction") := .data[[covariate]] * .data[[time_var]]
        )
      
      interaction_covs <- c(covariate, time_var, paste0(covariate, "_time_interaction"))
      
      time_varying_models[[paste0("linear_", time_var)]] <- fit_msm_models(
        patient_data = interaction_data,
        crude_rates = crude_rates,
        covariates = list("linear_interaction" = interaction_covs),
        nest_name = paste0(covariate, "_linear_time")
      )
      
    } else if (interaction_type == "spline") {
      # Spline interaction - covariate effect varies smoothly over time
      # This is complex and would require custom spline basis construction
      warning("Spline interactions not yet implemented")
      next
    }
  }
  
  # Compare models
  comparison <- compare_time_varying_effects(time_varying_models, covariate)
  
  return(list(
    models = time_varying_models,
    comparison = comparison,
    covariate = covariate,
    time_var = time_var
  ))
}

# Model comparison functions ----------------------------------------------

## Basic comparison functions ---------------------------------------------

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

#' Compare models across different structures using DRAIC and DRLCV
#' @param fitted_models Nested list of fitted MSM models
#' @param model_pairs Data frame with columns model1_struct, formula1, model2_struct, formula2
#' @param methods Vector of comparison methods ("draic", "drlcv", or both)
#' @param cores Number of cores for parallel computation
#' @return Data frame with cross-structure comparison metrics
compare_across_structures <- function(fitted_models, model_pairs = NULL, 
                                      methods = c("draic", "drlcv"), cores = 1) {
  
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
          
          result <- safe_drlcv(m1, m2, cores = cores, verbose = FALSE)
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


#' Comprehensive model comparison combining within and across structure comparisons
#' @param fitted_models Nested list of fitted MSM models
#' @param include_within Include within-structure comparisons (default: TRUE)
#' @param include_across Include cross-structure comparisons (default: TRUE)  
#' @param cross_structure_methods Methods for cross-structure comparison
#' @param cores Number of cores for DRLCV
#' @return List with within_structure and across_structure comparison results
comprehensive_model_comparison <- function(fitted_models, include_within = TRUE, 
                                           include_across = TRUE,
                                           cross_structure_methods = c("draic"),
                                           cores = 1) {
  
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
  }
  
  if (include_across) {
    cat("\n=== Cross-Structure Comparisons ===\n")
    results$across_structure <- compare_across_structures(
      fitted_models, 
      methods = cross_structure_methods, 
      cores = cores
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
  
  # Overall best models summary
  if (include_within) {
    best_models <- results$within_structure %>%
      filter(status == "converged") %>%
      group_by(model) %>%
      slice_min(AIC, n = 1) %>%
      ungroup() %>%
      arrange(AIC)
    
    cat("\n=== Best Model per Structure (by AIC) ===\n")
    print(best_models %>% select(model, formula, covariate, AIC, BIC, n_params))
  }
  
  return(results)
}

## Bespoke comparison functions ------------------------------------------

#' Compare time-varying models using the same framework
#' @param time_models Output from fit_time_varying_models()
#' @param include_comparison Use built-in comparison if available (default: TRUE)
#' @return Comparison results using consistent framework
compare_time_varying_models <- function(time_models, include_comparison = TRUE) {
  
  # Remove comparison entry if it exists to avoid confusion
  model_entries <- time_models[names(time_models) != "comparison"]
  
  # Use comprehensive comparison framework
  time_comparison <- comprehensive_model_comparison(
    model_entries,
    include_within = TRUE,
    include_across = TRUE,
    cross_structure_methods = c("draic"),
    cores = 1
  )
  
  # Add time-specific interpretation
  if (!is.null(time_comparison$within_structure)) {
    time_comparison$within_structure <- time_comparison$within_structure %>%
      mutate(
        time_model_type = case_when(
          grepl("linear", model) ~ "Linear time trend",
          grepl("spline", model) ~ "Spline time trend", 
          grepl("quadratic", model) ~ "Quadratic time trend",
          grepl("piecewise", model) ~ "Piecewise time periods",
          TRUE ~ "Other time model"
        )
      )
  }
  
  return(time_comparison)
}

#' Extract and compare covariate effects across nested models
#' @param fitted_models Nested list of fitted models
#' @param base_formula Formula to use as baseline (default: "~ 1")
#' @return Data frame with covariate effect comparisons
compare_covariate_effects <- function(fitted_models, base_formula = "~ 1") {
  
  # Use tidy_msm_models for consistent extraction
  tidied <- tidy_msm_models(fitted_models) %>%
    filter(status == "converged")
  
  if (nrow(tidied) == 0) {
    stop("No converged models found")
  }
  
  # Get baseline models for each structure
  base_models <- tidied %>%
    filter(formula == base_formula) %>%
    select(model, AIC_base = AIC, BIC_base = BIC, loglik_base = loglik)
  
  # Compare all models to their respective baselines
  covariate_comparison <- tidied %>%
    left_join(base_models, by = "model") %>%
    mutate(
      delta_AIC_from_base = AIC - AIC_base,
      delta_BIC_from_base = BIC - BIC_base,
      delta_loglik_from_base = loglik - loglik_base,
      improvement_AIC = delta_AIC_from_base < -2,
      improvement_BIC = delta_BIC_from_base < -2,
      substantial_improvement = delta_AIC_from_base < -10
    ) %>%
    arrange(model, AIC)
  
  # Add LRT results comparing to baseline
  lrt_results <- map_dfr(1:nrow(covariate_comparison), function(i) {
    row <- covariate_comparison[i, ]
    
    # Skip if this is the baseline model
    if (row$formula == base_formula) {
      return(data.frame(
        model = row$model,
        formula = row$formula,
        lrt_vs_base_stat = 0,
        lrt_vs_base_df = 0,
        lrt_vs_base_pval = 1
      ))
    }
    
    # Extract model objects using proper structure
    base_model <- extract_fitted_model(fitted_models, row$model, base_formula, 
                                       "in compare_covariate_effects base")
    test_model <- extract_fitted_model(fitted_models, row$model, row$formula, 
                                       "in compare_covariate_effects test")
    
    if (is.null(base_model) || is.null(test_model)) {
      return(data.frame(
        model = row$model,
        formula = row$formula,
        lrt_vs_base_stat = NA,
        lrt_vs_base_df = NA,
        lrt_vs_base_pval = NA
      ))
    }
    
    lrt_result <- tryCatch({
      lrt <- lrtest.msm(base_model, test_model)
      c(lrt$statistic, lrt$df, lrt$p.value)
    }, error = function(e) {
      c(NA, NA, NA)
    })
    
    data.frame(
      model = row$model,
      formula = row$formula,
      lrt_vs_base_stat = lrt_result[1],
      lrt_vs_base_df = lrt_result[2],
      lrt_vs_base_pval = lrt_result[3]
    )
  })

  # Merge with comparison results
  final_comparison <- covariate_comparison %>%
    left_join(lrt_results, by = c("model", "formula")) %>%
    mutate(
      significant_improvement = lrt_vs_base_pval < 0.05 & !is.na(lrt_vs_base_pval),
      effect_strength = case_when(
        delta_AIC_from_base < -10 ~ "Strong",
        delta_AIC_from_base < -2 ~ "Moderate", 
        delta_AIC_from_base < 2 ~ "Weak",
        TRUE ~ "None"
      )
    )
  
  return(final_comparison)
}

# Helper functions --------------------------------------------------------

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

#' Helper function to extract model object from nested structure
extract_model_object_old <- function(nested_models, model_name) {
  if (is.null(nested_models) || length(nested_models) == 0) return(NULL)
  
  # Handle different nesting levels
  if (!is.null(nested_models[[model_name]])) {
    model_entry <- nested_models[[model_name]]
    formula_key <- names(model_entry)[1]
    if (!is.null(formula_key)) {
      model_data <- model_entry[[formula_key]]
      if (is.list(model_data) && !is.null(model_data$fitted_model)) {
        return(model_data$fitted_model)
      } else {
        return(model_data)
      }
    }
  }
  
  return(NULL)
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
  
  # Only accept the correct nested structure
  if (!is.list(model_entry) || !"fitted_model" %in% names(model_entry)) {
    warning(paste("Invalid model structure for", structure, formula, 
                  "- expected fitted_models[[structure]][[formula]]$fitted_model", context))
    return(NULL)
  }
  
  # Check if model failed
  if (!is.null(model_entry$status) && model_entry$status == "failed") {
    warning(paste("Model failed for", structure, formula, context))
    return(NULL)
  }
  
  return(model_entry$fitted_model)
}

validate_model_structure <- function(fitted_models) {
  for (structure in names(fitted_models)) {
    for (formula in names(fitted_models[[structure]])) {
      entry <- fitted_models[[structure]][[formula]]
      
      if (!is.list(entry)) {
        stop(paste("Invalid structure: fitted_models[[", structure, "]][[", formula, 
                   "]] must be a list with $fitted_model component"))
      }
      
      if (!"fitted_model" %in% names(entry)) {
        stop(paste("Missing fitted_model: fitted_models[[", structure, "]][[", formula, 
                   "]]$fitted_model not found"))
      }
    }
  }
  return(TRUE)
}

#' Convert state names to state numbers based on model structure
#' @param state_names Vector of state names
#' @param fitted_model MSM model object to get state mapping from
#' @return Named vector mapping state names to numbers
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

#' Convert state numbers to state names
#' @param state_nums Vector of state numbers  
#' @param fitted_model MSM model object to get state mapping from
#' @return Named vector mapping state numbers to names
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

#' Standardize state columns in data
#' @param data Data frame with state information
#' @param fitted_model MSM model object for mapping
#' @param config MSM configuration (optional)
#' @return Data frame with both state and state_num columns
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

# Function to add calendar time column to patient data
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

#' Safe extraction of list elements
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

#' Format p-values for reporting
format_pvalue <- function(p, digits = 3) {
  case_when(
    is.na(p) ~ "—",
    p < 0.001 ~ "<0.001",
    p < 0.01 ~ paste0("<0.01"),
    TRUE ~ as.character(round(p, digits))
  )
}

#' Format confidence intervals
format_ci <- function(estimate, lower, upper, digits = 2) {
  paste0(
    round(estimate, digits), 
    " (", 
    round(lower, digits), 
    "—", 
    round(upper, digits), 
    ")"
  )
}


