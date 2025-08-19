# Additional utility functions for MSM analysis
# Add to end of 00-functions.R

\# ============================================================================
# ADVANCED DIAGNOSTIC FUNCTIONS
# ============================================================================

#' Calculate Pearson residuals for multi-state models
#' @param fitted_model A fitted msm model
#' @param patient_data Patient data
#' @return Data frame with Pearson residuals
calculate_pearson_residuals <- function(fitted_model, patient_data) {
  
  if (is.null(fitted_model)) return(NULL)
  
  # Get observed and expected transition counts
  obs_expected <- calculate_transition_deviance(fitted_model, patient_data)
  
  if (is.null(obs_expected)) return(NULL)
  
  # Calculate Pearson residuals
  pearson_residuals <- obs_expected %>%
    mutate(
      pearson_residual = (observed - expected_count) / sqrt(expected_count),
      standardized_residual = pearson_residual / sqrt(1 - 1/person_days)  # Approximate
    )
  
  return(pearson_residuals)
}

#' Perform goodness-of-fit test using chi-square
#' @param fitted_model A fitted msm model
#' @param patient_data Patient data
#' @return List with test statistic and p-value
goodness_of_fit_test <- function(fitted_model, patient_data) {
  
  residuals_data <- calculate_pearson_residuals(fitted_model, patient_data)
  
  if (is.null(residuals_data) || nrow(residuals_data) == 0) {
    return(list(chi_square = NA, df = NA, p_value = NA))
  }
  
  # Calculate chi-square statistic
  chi_square <- sum(residuals_data$pearson_residual^2, na.rm = TRUE)
  
  # Degrees of freedom (number of transitions - number of parameters)
  n_transitions <- nrow(residuals_data)
  n_params <- length(fitted_model$estimates)
  df <- max(1, n_transitions - n_params)
  
  # P-value
  p_value <- 1 - pchisq(chi_square, df)
  
  return(list(
    chi_square = chi_square,
    df = df,
    p_value = p_value,
    residuals = residuals_data
  ))
}

# ============================================================================
# BOOTSTRAP CONFIDENCE INTERVALS
# ============================================================================

#' Calculate bootstrap confidence intervals for model parameters
#' @param fitted_model A fitted msm model
#' @param patient_data Patient data  
#' @param n_bootstrap Number of bootstrap samples
#' @param confidence_level Confidence level (default: 0.95)
#' @return Data frame with bootstrap CIs
bootstrap_ci <- function(fitted_model, patient_data, n_bootstrap = 200, confidence_level = 0.95) {
  
  if (is.null(fitted_model)) return(NULL)
  
  alpha <- 1 - confidence_level
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2
  
  # Get unique patients for resampling
  unique_patients <- unique(patient_data$deid_enc_id)
  n_patients <- length(unique_patients)
  
  # Bootstrap resampling
  bootstrap_results <- map_dfr(1:n_bootstrap, function(b) {
    # Resample patients with replacement
    boot_patients <- sample(unique_patients, n_patients, replace = TRUE)
    
    # Create bootstrap dataset
    boot_data <- map_dfr(seq_along(boot_patients), function(i) {
      original_id <- boot_patients[i]
      patient_rows <- patient_data[patient_data$deid_enc_id == original_id, ]
      patient_rows$deid_enc_id <- i  # Assign new patient ID
      return(patient_rows)
    })
    
    # Fit model to bootstrap sample
    boot_model <- tryCatch({
      update(fitted_model, data = boot_data)
    }, error = function(e) NULL)
    
    if (!is.null(boot_model) && boot_model$opt$convergence == 0) {
      # Extract parameter estimates
      params <- boot_model$estimates
      data.frame(
        bootstrap = b,
        parameter = names(params),
        estimate = as.numeric(params)
      )
    } else {
      NULL
    }
  })
  
  if (nrow(bootstrap_results) > 0) {
    # Calculate confidence intervals
    ci_results <- bootstrap_results %>%
      group_by(parameter) %>%
      summarise(
        n_successful = n(),
        mean_estimate = mean(estimate, na.rm = TRUE),
        sd_estimate = sd(estimate, na.rm = TRUE),
        lower_ci = quantile(estimate, lower_q, na.rm = TRUE),
        upper_ci = quantile(estimate, upper_q, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        original_estimate = fitted_model$estimates[parameter],
        bias = mean_estimate - original_estimate,
        bootstrap_success_rate = n_successful / n_bootstrap
      )
    
    return(ci_results)
  } else {
    return(NULL)
  }
}

# ============================================================================
# MODEL VALIDATION FUNCTIONS
# ============================================================================

#' Perform leave-one-out cross-validation for model selection
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param model_specifications List of model specifications to compare
#' @return Data frame with LOO-CV results
leave_one_out_cv <- function(patient_data, crude_rates, model_specifications) {
  
  unique_patients <- unique(patient_data$deid_enc_id)
  n_patients <- length(unique_patients)
  
  cv_results <- list()
  
  for (spec_name in names(model_specifications)) {
    cat("Running LOO-CV for", spec_name, "\n")
    
    spec <- model_specifications[[spec_name]]
    fold_results <- list()
    
    for (i in seq_along(unique_patients)) {
      if (i %% 50 == 0) cat("  Fold", i, "of", n_patients, "\n")
      
      # Leave out patient i
      train_data <- patient_data %>% filter(deid_enc_id != unique_patients[i])
      test_data <- patient_data %>% filter(deid_enc_id == unique_patients[i])
      
      # Fit model on training data
      train_model <- tryCatch({
        if (!is.null(spec$covariates)) {
          fit_msm_models(train_data, crude_rates, covariates = spec$covariates)
        } else {
          fit_msm_models(train_data, crude_rates)
        }
      }, error = function(e) NULL)
      
      if (!is.null(train_model)) {
        # Calculate log-likelihood on test data
        test_loglik <- map_dbl(names(train_model), function(model_name) {
          model_test_data <- test_data %>% filter(model == model_name)
          fitted_model <- train_model[[model_name]]
          
          if (is.null(fitted_model) || nrow(model_test_data) == 0) {
            return(NA)
          }
          
          tryCatch({
            # Calculate log-likelihood for test patient
            logLik(fitted_model, newdata = model_test_data)
          }, error = function(e) NA)
        })
        
        fold_results[[i]] <- data.frame(
          fold = i,
          patient_id = unique_patients[i],
          model_name = names(train_model),
          test_loglik = test_loglik
        )
      }
    }
    
    if (length(fold_results) > 0) {
      cv_results[[spec_name]] <- bind_rows(fold_results) %>%
        mutate(specification = spec_name)
    }
  }
  
  # Combine and summarize results
  if (length(cv_results) > 0) {
    combined_results <- bind_rows(cv_results) %>%
      group_by(specification, model_name) %>%
      summarise(
        n_folds = sum(!is.na(test_loglik)),
        mean_test_loglik = mean(test_loglik, na.rm = TRUE),
        se_test_loglik = sd(test_loglik, na.rm = TRUE) / sqrt(n_folds),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_test_loglik))
    
    return(combined_results)
  } else {
    return(NULL)
  }
}


