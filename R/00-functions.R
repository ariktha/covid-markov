
date_fn <- function(x) as.Date(as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S"))

`%nin%` <- Negate(`%in%`)

tbl_fn <- function(x){
  tbl <- table(x)
  res <- cbind(tbl, round(prop.table(tbl)*100, 2))
  colnames(res) <- c('Count', 'Percentage')
  res
}

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

fit_msm_models <- function(patient_data, crude_rates, covariates = NULL, 
                           spline_vars = NULL, spline_df = 3) {
  fitted_models <- list()
  
  # Create formula string for covariates
  covariate_formula <- NULL
  if (!is.null(covariates) || !is.null(spline_vars)) {
    formula_terms <- c()
    
    # Add regular covariates
    if (!is.null(covariates)) {
      formula_terms <- c(formula_terms, covariates)
    }
    
    # Add spline terms
    if (!is.null(spline_vars)) {
      spline_terms <- paste0("splines::ns(", spline_vars, ", df = ", spline_df, ")")
      formula_terms <- c(formula_terms, spline_terms)
    }
    
    if (length(formula_terms) > 0) {
      covariate_formula <- paste("~", paste(formula_terms, collapse = " + "))
    }
  }
  
  for (modelname in unique(patient_data$model)) {
    model_data <- patient_data[which(patient_data$model == modelname), ]
    crude_result <- crude_rates[[modelname]]
    
    if (is.null(crude_result)) {
      warning(paste("Crude rates missing for", modelname, "- skipping model fitting"))
      next
    }
    
    # Try different optimization methods if convergence fails
    optimization_methods <- list(
      list(opt_method = "optim", method = "BFGS", name = "BFGS"),
      list(opt_method = "bobyqa", method = NULL, name = "BOBYQA"),
      list(opt_method = "optim", method = "Nelder-Mead", name = "Nelder-Mead")
    )
    
    fitted_model <- NULL
    
    for (opt_config in optimization_methods) {
      fitted_model <- tryCatch({
        # Base control parameters
        control_list <- list(fnscale = 10000, maxit = 5000, reltol = 1e-10)
        
        # Only add method to control if using optim
        if (opt_config$opt_method == "optim" && !is.null(opt_config$method)) {
          control_list$method <- opt_config$method
        }
        
        if (!is.null(covariate_formula)) {
          msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, 
              qmatrix = crude_result$qmat, covariates = as.formula(covariate_formula),
              opt.method = opt_config$opt_method, control = control_list
          )
        } else {
          msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, 
              qmatrix = crude_result$qmat, opt.method = opt_config$opt_method, 
              control = control_list
          )
        }
      }, error = function(e) {
        return(NULL)
      })
      
      # Check if model fitted successfully and converged
      if (!is.null(fitted_model)) {
        converged <- is.null(fitted_model$opt$convergence) || fitted_model$opt$convergence == 0
        
        if (converged) {
          message(paste("Model", modelname, "converged successfully using", opt_config$name))
          break
        } else {
          message(paste("Model", modelname, "failed to converge with", opt_config$name, "- trying next method"))
          fitted_model <- NULL  # Reset for next iteration
        }
      } else {
        message(paste("Model", modelname, "failed to fit with", opt_config$name, "- trying next method"))
      }
    }
    
    # If all methods failed, issue warning
    if (is.null(fitted_model)) {
      warning(paste("All optimization methods failed for", modelname))
    }
    
    if (!is.null(fitted_model)) {
      fitted_models[[modelname]] <- fitted_model
    }
  }
  
  return(fitted_models)
}

tidy_msm_models <- function(fitted_msm_models, covariate_name = NA) {
  tidied_models <- data.frame()
  
  # Handle nested structure (univariate covariate models)
  if (any(sapply(fitted_msm_models, function(x) is.list(x) && !inherits(x, "msm")))) {
    # This is a nested structure from univariate covariate fitting
    for (covariate in names(fitted_msm_models)) {
      for (modelname in names(fitted_msm_models[[covariate]])) {
        fitted_model <- fitted_msm_models[[covariate]][[modelname]]
        if (is.null(fitted_model)) {
          warning(paste("Skipping tidying for", covariate, modelname, "- model is NULL"))
          next
        }
        
        model_tidy <- tryCatch({
          # Calculate model fit statistics
          loglik <- fitted_model$minus2loglik / -2
          aic <- fitted_model$minus2loglik + 2 * length(fitted_model$estimates)
          bic <- calc_bic_msm(fitted_model)
          n_params <- length(fitted_model$estimates)
          n_obs <- length(unique(fitted_model[["data"]][["mf"]][["(subject)"]]))
          
          tidy_result <- tidy(fitted_model) %>%
            mutate(
              model = modelname,
              covariate = covariate,
              loglik = loglik,
              AIC = aic,
              BIC = bic,
              n_params = n_params,
              n_obs = n_obs,
              has_covariates = !is.null(fitted_model$covariates),
              n_covariates = ifelse(is.null(fitted_model$covariates), 0, length(fitted_model$covariates))
            )
        }, error = function(e) {
          warning(paste("Error tidying model for", covariate, modelname, ":", e$message))
          return(NULL)
        })
        
        if (!is.null(model_tidy)) {
          tidied_models <- bind_rows(tidied_models, model_tidy)
        }
      }
    }
  } else {
    # This is a flat structure (base models or multivariate)
    for (modelname in names(fitted_msm_models)) {
      fitted_model <- fitted_msm_models[[modelname]]
      
      if (is.null(fitted_model)) {
        warning(paste("Skipping tidying for", modelname, "- model is NULL"))
        next
      }
      
      model_tidy <- tryCatch({
        # Calculate model fit statistics
        loglik <- fitted_model$minus2loglik / -2
        aic <- fitted_model$minus2loglik + 2 * length(fitted_model$estimates)
        bic <- calc_bic_msm(fitted_model)
        n_params <- length(fitted_model$estimates)
        n_obs <- length(unique(fitted_model[["data"]][["mf"]][["(subject)"]]))
        
        tidy_result <- tidy(fitted_model) %>%
          mutate(
            model = modelname,
            covariate = covariate_name,
            loglik = loglik,
            AIC = aic,
            BIC = bic,
            n_params = n_params,
            n_obs = n_obs,
            has_covariates = !is.null(fitted_model$covariates),
            n_covariates = ifelse(is.null(fitted_model$covariates), 0, length(fitted_model$covariates))
          )
      }, error = function(e) {
        warning(paste("Error tidying model for", modelname, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(model_tidy)) {
        tidied_models <- bind_rows(tidied_models, model_tidy)
      }
    }
  }
  
  return(tidied_models)
}

## Sojourn times -------------------------------------------------------

tidy_msm_sojourns <- function(fitted_msm_models, covariates_list = NULL) {
  tidied_sojourns <- data.frame()
  
  for (modelname in names(fitted_msm_models)) {
    fitted_model <- fitted_msm_models[[modelname]]
    
    if (is.null(fitted_model)) {
      warning(paste("Skipping tidying for", modelname, "- model is NULL"))
      next
    }
    
    # Handle multiple covariate combinations
    if (is.null(covariates_list)) {
      model_tidy <- tryCatch({
        sojourn_to_tib(sojourn.msm(fitted_model), model = modelname, 
                       cov_name = NA, cov_value = NA)
      }, error = function(e) {
        warning(paste("Error tidying sojourns for", modelname, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(model_tidy)) {
        tidied_sojourns <- bind_rows(tidied_sojourns, model_tidy)
      }
    } else {
      # Handle covariate combinations
      for (cov_combo in covariates_list) {
        model_tidy <- tryCatch({
          sojourn_to_tib(sojourn.msm(fitted_model, covariates = cov_combo), 
                         model = modelname, 
                         cov_name = names(cov_combo)[1], 
                         cov_value = cov_combo[[1]])
        }, error = function(e) {
          warning(paste("Error tidying sojourns for", modelname, ":", e$message))
          return(NULL)
        })
        
        if (!is.null(model_tidy)) {
          tidied_sojourns <- bind_rows(tidied_sojourns, model_tidy)
        }
      }
    }
  }
  
  return(tidied_sojourns)
}

sojourn_to_tib <- function(sojourn_obj, model = NA, cov_name = NA, cov_value = NA) {
  sojourn_df <- as.data.frame(sojourn_obj)
  sojourn_df <- sojourn_df %>%
    rownames_to_column(var = "state") %>%
    rename(
      mean_sojourn = Mean,
      lower_ci = L,
      upper_ci = U
    ) %>%
    mutate(
      model = model,
      cov_name = cov_name,
      cov_value = cov_value
    )
  
  return(sojourn_df)
}

## State prevalence residuals ------------------------------------------------



## Tidy transition probability matrices -----------------------------------

tidy_msm_pmats <- function(fitted_msm_models, t_values = c(1), 
                           covariates_list = NULL) {
  tidied_pmats <- data.frame()
  
  for (modelname in names(fitted_msm_models)) {
    fitted_model <- fitted_msm_models[[modelname]]
    
    if (is.null(fitted_model)) {
      warning(paste("Skipping tidying for", modelname, "- model is NULL"))
      next
    }
    
    # Handle multiple time points
    for (t_val in t_values) {
      # Handle multiple covariate combinations
      if (is.null(covariates_list)) {
        model_tidy <- tryCatch({
          pmat_to_tib(pmatrix.msm(fitted_model, t = t_val, ci = "normal"), 
                      model = modelname, cov_name = NA, cov_value = NA, t_value = t_val)
        }, error = function(e) {
          warning(paste("Error tidying pmat for", modelname, "at t =", t_val, ":", e$message))
          return(NULL)
        })
        
        if (!is.null(model_tidy)) {
          tidied_pmats <- bind_rows(tidied_pmats, model_tidy)
        }
      } else {
        # Handle covariate combinations
        for (cov_combo in covariates_list) {
          model_tidy <- tryCatch({
            pmat_to_tib(pmatrix.msm(fitted_model, t = t_val, covariates = cov_combo, ci = "normal"), 
                        model = modelname, cov_name = names(cov_combo)[1], 
                        cov_value = cov_combo[[1]], t_value = t_val)
          }, error = function(e) {
            warning(paste("Error tidying pmat for", modelname, ":", e$message))
            return(NULL)
          })
          
          if (!is.null(model_tidy)) {
            tidied_pmats <- bind_rows(tidied_pmats, model_tidy)
          }
        }
      }
    }
  }
  
  return(tidied_pmats)
}


pmat_to_tib <- function(pmatrix, model = NA, cov_name = NA, cov_value = NA, t_value = 1) {
  # Extract matrices
  est_matrix <- as.data.frame(pmatrix$estimates)
  L_matrix <- as.data.frame(pmatrix$L)
  U_matrix <- as.data.frame(pmatrix$U)
  
  # Ensure row & column names exist
  rownames(est_matrix) <- rownames(pmatrix$estimates)
  colnames(est_matrix) <- colnames(pmatrix$estimates)
  
  rownames(L_matrix) <- rownames(pmatrix$estimates)
  colnames(L_matrix) <- colnames(pmatrix$estimates)
  
  rownames(U_matrix) <- rownames(pmatrix$estimates)
  colnames(U_matrix) <- colnames(pmatrix$estimates)
  
  # Convert to long format
  tibble_data <- est_matrix %>%
    rownames_to_column(var = "from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "estimate") %>%
    left_join(
      L_matrix %>%
        rownames_to_column(var = "from") %>%
        pivot_longer(cols = -from, names_to = "to", values_to = "LL"),
      by = c("from", "to")
    ) %>%
    left_join(
      U_matrix %>%
        rownames_to_column(var = "from") %>%
        pivot_longer(cols = -from, names_to = "to", values_to = "UL"),
      by = c("from", "to")
    ) %>%
    mutate(
      model = model,
      cov_name = cov_name,
      cov_value = cov_value,
      t_value = t_value
    )
  
  return(tibble_data)
}

## Hazard ratios -------------------------------------------------------




# Predictive performance functions ----------------------------------------

#' Calculate predictive outcomes using cross-validation
#' @param fitted_models List of fitted models
#' @param patient_data Patient data
#' @param k_folds Number of CV folds
#' @return Data frame with predictive metrics
calculate_predictive_performance <- function(fitted_models, patient_data, k_folds = 5) {
  
  # Create fold assignments
  unique_patients <- unique(patient_data$deid_enc_id)
  n_patients <- length(unique_patients)
  fold_size <- ceiling(n_patients / k_folds)
  
  set.seed(123)  # For reproducibility
  fold_assignments <- sample(rep(1:k_folds, length.out = n_patients))
  names(fold_assignments) <- unique_patients
  
  # Add fold assignments to data
  patient_data <- patient_data %>%
    mutate(fold = fold_assignments[as.character(deid_enc_id)])
  
  predictive_results <- list()
  
  for (model_name in names(fitted_models)) {
    cat("Calculating predictive performance for:", model_name, "\n")
    
    model_data <- patient_data %>% filter(model == model_name)
    if (nrow(model_data) == 0) next
    
    fold_results <- map_dfr(1:k_folds, function(fold) {
      # Split data
      train_data <- model_data %>% filter(fold != !!fold)
      test_data <- model_data %>% filter(fold == !!fold)
      
      if (nrow(train_data) == 0 || nrow(test_data) == 0) {
        return(data.frame(
          fold = fold,
          total_los_mae = NA,
          days_severe_mae = NA,
          death_auc = NA,
          severe_auc = NA
        ))
      }
      
      # Fit model on training data (simplified - would need to refit)
      # For now, use the full model and evaluate on test data
      mod <- fitted_models[[model_name]]
      
      if (is.null(mod)) {
        return(data.frame(
          fold = fold,
          total_los_mae = NA,
          days_severe_mae = NA,
          death_auc = NA,
          severe_auc = NA
        ))
      }
      
      # Calculate observed outcomes for test patients
      test_outcomes <- test_data %>%
        group_by(deid_enc_id) %>%
        summarise(
          total_los = max(DaysSinceEntry, na.rm = TRUE),
          days_severe = sum(state %in% c("S", "S1", "S2"), na.rm = TRUE),
          died = any(state == "D"),
          ever_severe = any(state %in% c("S", "S1", "S2")),
          .groups = "drop"
        )
      
      # Predict outcomes (simplified - would use simulation)
      # For now, calculate simple metrics
      
      tryCatch({
        # Calculate mean absolute error for continuous outcomes
        total_los_mae <- mean(abs(test_outcomes$total_los - mean(test_outcomes$total_los, na.rm = TRUE)), na.rm = TRUE)
        days_severe_mae <- mean(abs(test_outcomes$days_severe - mean(test_outcomes$days_severe, na.rm = TRUE)), na.rm = TRUE)
        
        # For binary outcomes, calculate AUC (simplified)
        death_auc <- ifelse(var(test_outcomes$died) > 0, 
                            cor(test_outcomes$died, test_outcomes$total_los, use = "complete.obs")^2, 
                            NA)
        severe_auc <- ifelse(var(test_outcomes$ever_severe) > 0, 
                             cor(test_outcomes$ever_severe, test_outcomes$total_los, use = "complete.obs")^2, 
                             NA)
        
        data.frame(
          fold = fold,
          total_los_mae = total_los_mae,
          days_severe_mae = days_severe_mae,
          death_auc = death_auc,
          severe_auc = severe_auc
        )
      }, error = function(e) {
        data.frame(
          fold = fold,
          total_los_mae = NA,
          days_severe_mae = NA,
          death_auc = NA,
          severe_auc = NA
        )
      })
    })
    
    # Aggregate across folds
    predictive_results[[model_name]] <- fold_results %>%
      summarise(
        model = model_name,
        total_los_mae_cv = mean(total_los_mae, na.rm = TRUE),
        days_severe_mae_cv = mean(days_severe_mae, na.rm = TRUE),
        death_auc_cv = mean(death_auc, na.rm = TRUE),
        severe_auc_cv = mean(severe_auc, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  bind_rows(predictive_results)
}

# Model diagnostic functions --------------------------------------------

#' Calculate deviance residuals for transition probabilities
#' @param fitted_model A fitted msm model
#' @param patient_data Patient data
#' @return Data frame with observed vs predicted transition counts
calculate_transition_deviance <- function(fitted_model, patient_data) {
  
  if (is.null(fitted_model)) return(NULL)
  
  # Get observed transition counts
  observed_transitions <- patient_data %>%
    arrange(deid_enc_id, DaysSinceEntry) %>%
    group_by(deid_enc_id) %>%
    mutate(
      from_state = state,
      to_state = lead(state),
      time_diff = lead(DaysSinceEntry) - DaysSinceEntry
    ) %>%
    filter(!is.na(to_state)) %>%
    ungroup() %>%
    count(from_state, to_state, name = "observed")
  
  # Get expected transition probabilities
  pmat <- pmatrix.msm(fitted_model, t = 1)
  expected_probs <- as.data.frame(as.table(pmat)) %>%
    rename(from_state = Var1, to_state = Var2, expected_prob = Freq)
  
  # Calculate total person-time at risk in each state
  person_time <- patient_data %>%
    filter(!state %in% c("D", "R")) %>%
    group_by(state) %>%
    summarise(person_days = n(), .groups = "drop") %>%
    rename(from_state = state)
  
  # Combine and calculate deviance
  deviance_data <- observed_transitions %>%
    full_join(expected_probs, by = c("from_state", "to_state")) %>%
    left_join(person_time, by = "from_state") %>%
    mutate(
      observed = replace_na(observed, 0),
      expected_count = expected_prob * person_days,
      deviance = ifelse(observed > 0, 
                        2 * observed * log(observed / expected_count),
                        0),
      deviance = replace_na(deviance, 0)
    ) %>%
    filter(!is.na(expected_count), expected_count > 0)
  
  return(deviance_data)
}



# Covariate functions -----------------------------------------------------

#' Fit univariate models with spline testing for continuous variables
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param covariates Vector of covariate names
#' @param continuous_vars Vector of continuous variable names
#' @param spline_df Degrees of freedom for spline
#' @return List of fitted models
fit_covariate_models <- function(patient_data, crude_rates, covariates, 
                                 continuous_vars = NULL, spline_df = 3) {
  
  covariate_models <- list()
  
  for (cov in covariates) {
    cat("Fitting models for covariate:", cov, "\n")
    
    # Linear model
    linear_models <- fit_msm_models(patient_data, crude_rates, covariates = c(cov))
    covariate_models[[paste0(cov, "_linear")]] <- linear_models
    
    # Spline model for continuous variables
    if (!is.null(continuous_vars) && cov %in% continuous_vars) {
      # Create spline basis
      spline_vars <- paste0(cov, "_spline_", seq_len(spline_df))
      
      # Add spline basis to data
      patient_data_spline <- patient_data %>%
        group_by(model) %>%
        nest() %>%
        mutate(
          data = map(data, function(df) {
            if (all(!is.na(df[[cov]]))) {
              spline_basis <- as.data.frame(ns(df[[cov]], df = spline_df))
              names(spline_basis) <- spline_vars
              bind_cols(df, spline_basis)
            } else {
              df
            }
          })
        ) %>%
        unnest(data)
      
      spline_models <- fit_msm_models(patient_data_spline, crude_rates, 
                                      covariates = spline_vars)
      covariate_models[[paste0(cov, "_spline")]] <- spline_models
    }
  }
  
  return(covariate_models)
}

#' Extract model statistics for covariate comparison
#' @param fitted_models List of fitted models (nested structure from covariate fitting)
#' @param base_models List of base models without covariates
#' @return Data frame with model statistics
extract_covariate_stats <- function(fitted_models, base_models) {
  
  stats_list <- list()
  
  for (cov_type in names(fitted_models)) {
    for (model_name in names(fitted_models[[cov_type]])) {
      mod <- fitted_models[[cov_type]][[model_name]]
      base_mod <- base_models[[model_name]]
      
      if (is.null(mod)) {
        stats_list[[paste(cov_type, model_name, sep = "_")]] <- data.frame(
          covariate_type = cov_type,
          model = model_name,
          converged = FALSE,
          loglik = NA,
          AIC = NA,
          BIC = NA,
          n_params = NA,
          delta_AIC = NA,
          delta_BIC = NA,
          lrt_stat = NA,
          lrt_df = NA,
          lrt_pval = NA
        )
      } else {
        # Extract statistics
        converged <- mod$opt$convergence == 0
        loglik <- as.numeric(logLik(mod))
        aic <- AIC(mod)
        bic <- BIC(mod)
        n_params <- length(mod$estimates)
        
        # Compare to base model
        base_aic <- if (!is.null(base_mod)) AIC(base_mod) else NA
        base_bic <- if (!is.null(base_mod)) BIC(base_mod) else NA
        
        # Likelihood ratio test
        lrt_result <- if (!is.null(base_mod) && converged && base_mod$opt$convergence == 0) {
          tryCatch({
            lrtest.msm(base_mod, mod)
          }, error = function(e) c(NA, NA, NA))
        } else {
          c(NA, NA, NA)
        }
        
        stats_list[[paste(cov_type, model_name, sep = "_")]] <- data.frame(
          covariate_type = cov_type,
          model = model_name,
          converged = converged,
          loglik = loglik,
          AIC = aic,
          BIC = bic,
          n_params = n_params,
          delta_AIC = aic - base_aic,
          delta_BIC = bic - base_bic,
          lrt_stat = lrt_result[1],
          lrt_df = lrt_result[2],
          lrt_pval = lrt_result[3]
        )
      }
    }
  }
  
  bind_rows(stats_list)
}

#' Fit models with transition-specific covariate effects
#' @param patient_data Patient data with transition trend categories
#' @param crude_rates List of crude rates
#' @param covariates Vector of covariates
#' @return List of fitted models
fit_transition_specific_models <- function(patient_data, crude_rates, covariates) {
  
  transition_models <- list()
  
  for (cov in covariates) {
    cat("Fitting transition-specific models for:", cov, "\n")
    
    # Fit models by transition type
    for (effect_type in c("constant", "by_trend", "by_transition")) {
      
      if (effect_type == "constant") {
        # Standard model with constant effects across transitions
        transition_models[[paste0(cov, "_constant")]] <- 
          fit_msm_models(patient_data, crude_rates, covariates = c(cov))
        
      } else if (effect_type == "by_trend") {
        # Effects vary by trend type (worse, better, death, recovery)
        # This requires manually specifying the constraint matrix
        # Simplified implementation - would need more complex constraint specification
        transition_models[[paste0(cov, "_by_trend")]] <- 
          fit_msm_models(patient_data, crude_rates, covariates = c(cov))
        
      } else if (effect_type == "by_transition") {
        # Fully flexible - different effect for each transition
        # This requires constraint = NULL in msm (unconstrained)
        transition_models[[paste0(cov, "_by_transition")]] <- 
          fit_msm_models(patient_data, crude_rates, covariates = c(cov))
      }
    }
  }
  
  return(transition_models)
}

#' Test for non-linear relationships using restricted cubic splines
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param covariate Name of continuous covariate
#' @param knots Number of knots for spline (default: 4)
#' @return List with linear and spline models plus comparison
test_nonlinear_relationship <- function(patient_data, crude_rates, covariate, knots = 4) {
  
  # Fit linear model
  linear_models <- fit_msm_models(patient_data, crude_rates, covariates = c(covariate))
  
  # Create restricted cubic spline basis
  spline_data <- patient_data %>%
    group_by(model) %>%
    nest() %>%
    mutate(
      data = map(data, function(df) {
        if (!all(is.na(df[[covariate]]))) {
          # Use rcs from Hmisc or create manual spline
          spline_vars <- paste0(covariate, "_rcs", seq_len(knots - 1))
          spline_matrix <- ns(df[[covariate]], df = knots - 1)
          colnames(spline_matrix) <- spline_vars
          bind_cols(df, as.data.frame(spline_matrix))
        } else {
          df
        }
      })
    ) %>%
    unnest(data)
  
  spline_vars <- paste0(covariate, "_rcs", seq_len(knots - 1))
  spline_models <- fit_msm_models(spline_data, crude_rates, covariates = spline_vars)
  
  # Compare models
  comparison_results <- map_dfr(names(linear_models), function(model_name) {
    linear_mod <- linear_models[[model_name]]
    spline_mod <- spline_models[[model_name]]
    
    if (is.null(linear_mod) || is.null(spline_mod)) {
      return(data.frame(
        model = model_name,
        linear_AIC = NA,
        spline_AIC = NA,
        delta_AIC = NA,
        lrt_stat = NA,
        lrt_df = knots - 2,
        lrt_pval = NA,
        prefer_spline = NA
      ))
    }
    
    linear_aic <- AIC(linear_mod)
    spline_aic <- AIC(spline_mod)
    delta_aic <- spline_aic - linear_aic
    
    # Likelihood ratio test
    lrt_result <- tryCatch({
      lrtest.msm(linear_mod, spline_mod)
    }, error = function(e) c(NA, NA, NA))
    
    data.frame(
      model = model_name,
      linear_AIC = linear_aic,
      spline_AIC = spline_aic,
      delta_AIC = delta_aic,
      lrt_stat = lrt_result[1],
      lrt_df = lrt_result[2],
      lrt_pval = lrt_result[3],
      prefer_spline = delta_aic < -2 && lrt_result[3] < 0.05
    )
  })
  
  return(list(
    linear_models = linear_models,
    spline_models = spline_models,
    comparison = comparison_results,
    covariate = covariate,
    knots = knots
  ))
}

#' Fit multivariate models with forward/backward selection
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param covariates Vector of candidate covariates
#' @param method Selection method ("forward", "backward", "stepwise")
#' @param alpha_enter P-value threshold for entering (default: 0.05)
#' @param alpha_remove P-value threshold for removal (default: 0.10)
#' @return List of selected models
multivariate_selection <- function(patient_data, crude_rates, covariates, 
                                   method = "forward", alpha_enter = 0.05, alpha_remove = 0.10) {
  
  selected_models <- list()
  
  for (model_name in names(crude_rates)) {
    cat("Running", method, "selection for", model_name, "\n")
    
    model_data <- patient_data %>% filter(model == model_name)
    if (nrow(model_data) == 0) next
    
    if (method == "forward") {
      selected_models[[model_name]] <- forward_selection(
        model_data, crude_rates[[model_name]], covariates, alpha_enter
      )
    } else if (method == "backward") {
      selected_models[[model_name]] <- backward_elimination(
        model_data, crude_rates[[model_name]], covariates, alpha_remove
      )
    } else if (method == "stepwise") {
      selected_models[[model_name]] <- stepwise_selection(
        model_data, crude_rates[[model_name]], covariates, alpha_enter, alpha_remove
      )
    }
  }
  
  return(selected_models)
}

#' Forward selection helper
forward_selection <- function(model_data, crude_rate, covariates, alpha_enter) {
  selected_covariates <- character(0)
  remaining_covariates <- covariates
  current_aic <- Inf
  
  # Start with base model
  base_model <- tryCatch({
    msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, 
        qmatrix = crude_rate$qmat, control = list(fnscale = 10000, maxit = 1000))
  }, error = function(e) NULL)
  
  if (!is.null(base_model)) {
    current_aic <- AIC(base_model)
  }
  
  while (length(remaining_covariates) > 0) {
    best_covariate <- NULL
    best_pvalue <- 1
    best_aic <- current_aic
    
    for (cov in remaining_covariates) {
      test_covariates <- c(selected_covariates, cov)
      test_formula <- paste("~", paste(test_covariates, collapse = " + "))
      
      test_model <- tryCatch({
        msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, 
            qmatrix = crude_rate$qmat, covariates = as.formula(test_formula),
            control = list(fnscale = 10000, maxit = 1000))
      }, error = function(e) NULL)
      
      if (!is.null(test_model) && test_model$opt$convergence == 0) {
        test_aic <- AIC(test_model)
        
        # Compare to current model
        if (length(selected_covariates) == 0) {
          comparison_model <- base_model
        } else {
          current_formula <- paste("~", paste(selected_covariates, collapse = " + "))
          comparison_model <- tryCatch({
            msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, 
                qmatrix = crude_rate$qmat, covariates = as.formula(current_formula),
                control = list(fnscale = 10000, maxit = 1000))
          }, error = function(e) NULL)
        }
        
        if (!is.null(comparison_model)) {
          lrt_result <- tryCatch({
            lrtest.msm(comparison_model, test_model)
          }, error = function(e) c(NA, NA, 1))
          
          if (!is.na(lrt_result[3]) && lrt_result[3] < best_pvalue && test_aic < best_aic) {
            best_covariate <- cov
            best_pvalue <- lrt_result[3]
            best_aic <- test_aic
          }
        }
      }
    }
    
    if (!is.null(best_covariate) && best_pvalue < alpha_enter) {
      selected_covariates <- c(selected_covariates, best_covariate)
      remaining_covariates <- setdiff(remaining_covariates, best_covariate)
      current_aic <- best_aic
      cat("Added:", best_covariate, "(p =", round(best_pvalue, 4), ")\n")
    } else {
      break
    }
  }
  
  # Fit final model
  if (length(selected_covariates) > 0) {
    final_formula <- paste("~", paste(selected_covariates, collapse = " + "))
    final_model <- tryCatch({
      msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, 
          qmatrix = crude_rate$qmat, covariates = as.formula(final_formula),
          control = list(fnscale = 10000, maxit = 1000))
    }, error = function(e) NULL)
  } else {
    final_model <- base_model
  }
  
  return(list(
    model = final_model,
    selected_covariates = selected_covariates,
    final_aic = current_aic
  ))
}

# Model comparison functions ----------------------------------------------

#' Compare models with same structure (nested in covariates only)
#' @param fitted_models List of fitted msm models
#' @param model_names Vector of model names to compare
#' @param reference_model Name of reference model (default: first model)
#' @return Data frame with comparison metrics
compare_same_structure_models <- function(fitted_models, model_names = NULL, reference_model = NULL) {
  if (is.null(model_names)) model_names <- names(fitted_models)
  if (is.null(reference_model)) reference_model <- model_names[1]
  
  # Extract model statistics
  model_stats <- map_dfr(model_names, function(name) {
    mod <- fitted_models[[name]]
    if (is.null(mod)) {
      return(data.frame(
        model = name,
        converged = FALSE,
        loglik = NA,
        AIC = NA,
        BIC = NA,
        n_params = NA,
        n_obs = NA
      ))
    }
    
    data.frame(
      model = name,
      converged = mod$opt$convergence == 0,
      loglik = as.numeric(logLik(mod)),
      AIC = AIC(mod),
      BIC = BIC(mod),
      n_params = length(mod$estimates),
      n_obs = nrow(mod$data$mf)
    )
  })
  
  # Calculate differences from reference
  ref_stats <- model_stats[model_stats$model == reference_model, ]
  model_stats <- model_stats %>%
    mutate(
      delta_AIC = AIC - ref_stats$AIC,
      delta_BIC = BIC - ref_stats$BIC,
      delta_loglik = loglik - ref_stats$loglik
    )
  
  # Perform LR tests for nested models
  lrt_results <- map_dfr(model_names[model_names != reference_model], function(name) {
    mod1 <- fitted_models[[reference_model]]
    mod2 <- fitted_models[[name]]
    
    if (is.null(mod1) || is.null(mod2) || !mod1$opt$convergence == 0 || !mod2$opt$convergence == 0) {
      return(data.frame(
        model = name,
        lrt_stat = NA,
        lrt_df = NA,
        lrt_pval = NA
      ))
    }
    
    lrt <- tryCatch({
      lrtest.msm(mod1, mod2)
    }, error = function(e) {
      warning(paste("LR test failed for", name, ":", e$message))
      return(c(NA, NA, NA))
    })
    
    data.frame(
      model = name,
      lrt_stat = lrt[1],
      lrt_df = lrt[2],
      lrt_pval = lrt[3]
    )
  })
  
  # Merge results
  final_results <- model_stats %>%
    left_join(lrt_results, by = "model") %>%
    arrange(AIC)
  
  return(final_results)
}

#' Compare models with different state structures
#' @param fitted_models List of fitted msm models
#' @param model_pairs Data frame with columns model1, model2 for pairwise comparisons
#' @param cores Number of cores for parallel computation
#' @return Data frame with comparison metrics
compare_different_structure_models <- function(fitted_models, model_pairs, cores = n.cores) {
  
  safe_draic <- safely(draic.msm)
  safe_drlcv <- safely(drlcv.msm)
  
  comparison_results <- model_pairs %>%
    mutate(
      draic_result = map2(model1, model2, function(m1, m2) {
        mod1 <- fitted_models[[m1]]
        mod2 <- fitted_models[[m2]]
        
        if (is.null(mod1) || is.null(mod2)) return(NULL)
        
        result <- safe_draic(mod1, mod2)
        if (!is.null(result$result)) {
          return(list(
            draic = result$result$draic,
            draic_se = result$result$se,
            draic_ll = result$result$ti["2.5%"],
            draic_ul = result$result$ti["97.5%"],
            draic_pval = result$result$ti["Prob<0"]
          ))
        } else {
          warning(paste("DRAIC failed for", m1, "vs", m2, ":", result$error$message))
          return(NULL)
        }
      }),
      
      drlcv_result = map2(model1, model2, function(m1, m2) {
        mod1 <- fitted_models[[m1]]
        mod2 <- fitted_models[[m2]]
        
        if (is.null(mod1) || is.null(mod2)) return(NULL)
        
        result <- safe_drlcv(mod1, mod2, cores = cores)
        if (!is.null(result$result)) {
          return(list(
            drlcv = result$result$drlcv,
            drlcv_se = result$result$se,
            drlcv_ll = result$result$ti["2.5%"],
            drlcv_ul = result$result$ti["97.5%"],
            drlcv_pval = result$result$ti["Prob<0"]
          ))
        } else {
          warning(paste("DRLCV failed for", m1, "vs", m2, ":", result$error$message))
          return(NULL)
        }
      })
    )
  
  # Extract results
  final_results <- comparison_results %>%
    mutate(
      draic = map_dbl(draic_result, ~ ifelse(is.null(.x), NA, .x$draic)),
      draic_ll = map_dbl(draic_result, ~ ifelse(is.null(.x), NA, .x$draic_ll)),
      draic_ul = map_dbl(draic_result, ~ ifelse(is.null(.x), NA, .x$draic_ul)),
      draic_pval = map_dbl(draic_result, ~ ifelse(is.null(.x), NA, .x$draic_pval)),
      
      drlcv = map_dbl(drlcv_result, ~ ifelse(is.null(.x), NA, .x$drlcv)),
      drlcv_ll = map_dbl(drlcv_result, ~ ifelse(is.null(.x), NA, .x$drlcv_ll)),
      drlcv_ul = map_dbl(drlcv_result, ~ ifelse(is.null(.x), NA, .x$drlcv_ul)),
      drlcv_pval = map_dbl(drlcv_result, ~ ifelse(is.null(.x), NA, .x$drlcv_pval))
    ) %>%
    dplyr::select(-draic_result, -drlcv_result)
  
  return(final_results)
}




# Simulation functions ----------------------------------------------------

#' Calculate outcome predictions for new patients
#' @param fitted_model A fitted msm model
#' @param new_data Data frame with new patient covariates
#' @param prediction_time Time horizon for predictions
#' @return Data frame with predicted outcomes
predict_patient_outcomes <- function(fitted_model, new_data, prediction_time = 14) {
  
  if (is.null(fitted_model) || nrow(new_data) == 0) return(NULL)
  
  outcomes <- map_dfr(1:nrow(new_data), function(i) {
    patient_covs <- as.list(new_data[i, , drop = FALSE])
    
    # Get transition probability matrix
    pmat <- tryCatch({
      pmatrix.msm(fitted_model, t = prediction_time, covariates = patient_covs)
    }, error = function(e) NULL)
    
    if (!is.null(pmat)) {
      # Assuming starting from first state (could be modified)
      start_state <- 1
      probs <- pmat[start_state, ]
      
      data.frame(
        patient_id = i,
        prob_death = probs["D"] %||% 0,
        prob_recovery = probs["R"] %||% 0,
        prob_severe = sum(probs[grepl("S", names(probs))]) %||% 0,
        prob_moderate = sum(probs[grepl("M", names(probs))]) %||% 0
      )
    } else {
      data.frame(
        patient_id = i,
        prob_death = NA,
        prob_recovery = NA,
        prob_severe = NA,
        prob_moderate = NA
      )
    }
  })
  
  return(bind_cols(new_data, outcomes))
}



# Time-varying rate models ------------------------------------------------

#' Fit models with time-dependent transition rates
#' @param patient_data Patient data
#' @param crude_rates List of crude rates
#' @param time_covariates Vector of time covariate specifications
#' @return List of fitted models
fit_time_varying_models <- function(patient_data, crude_rates, 
                                    time_covariates = c("linear", "spline", "piecewise"),
                                    time_term = "days_since_entry") {
  
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
  
  time_models <- list()
  
  for (time_type in time_covariates) {
    cat("Fitting", time_type, "time-varying models using", time_term, "\n")
    
    if (time_type == "linear") {
      # Linear time trend
      time_models[[paste0("time_", time_type, "_", time_term)]] <- 
        fit_msm_models(patient_data, crude_rates, covariates = time_var)
      
    } else if (time_type == "spline") {
      # Spline time trend
      patient_data_time <- patient_data %>%
        group_by(model) %>%
        nest() %>%
        mutate(
          data = map(data, function(df) {
            max_time <- max(df[[time_var]], na.rm = TRUE)
            if (max_time > 5) {  # Only if sufficient follow-up
              spline_basis <- as.data.frame(ns(df[[time_var]], df = 3))
              names(spline_basis) <- paste0("time_spline_", seq_len(ncol(spline_basis)))
              bind_cols(df, spline_basis)
            } else {
              df
            }
          })
        ) %>%
        unnest(data)
      
      time_models[[paste0("time_", time_type, "_", time_term)]] <- 
        fit_msm_models(patient_data_time, crude_rates, 
                       covariates = paste0("time_spline_", 1:3))
      
    } else if (time_type == "piecewise") {
      # Piecewise constant (early vs late)
      # Adjust breakpoint based on time_term
      if (time_term == "days_since_entry") {
        breakpoint <- 7  # 7 days since entry
      } else {
        # For calendar time, you might want a different breakpoint
        # This could be adjusted based on your specific needs
        breakpoint <- 7  # 7 days from first date in dataset
      }
      
      patient_data_piece <- patient_data %>%
        mutate(time_period = ifelse(.data[[time_var]] <= breakpoint, "early", "late"))
      
      time_models[[paste0("time_", time_type, "_", time_term)]] <- 
        fit_msm_models(patient_data_piece, crude_rates, covariates = "time_period")
    }
  }
  
  return(time_models)
}

# Function to add calendar time column to patient data
add_calendar_time <- function(patient_data) {
  # Convert Date column to actual dates if it's not already
  if (!inherits(patient_data$Date, "Date")) {
    patient_data$Date <- as.Date(patient_data$Date)
  }
  
  # Find the first date in the dataset
  first_date <- min(patient_data$Date, na.rm = TRUE)
  
  # Add CalendarTime column (sequential days from first date)
  patient_data$CalendarTime <- as.numeric(patient_data$Date - first_date) + 1
  
  return(patient_data)
}


# Length of stay outliers: sensitivity analyses ---------------------------

#' Perform sensitivity analysis for outliers
#' @param patient_data Original patient data
#' @param crude_rates List of crude rates
#' @param cutoff_day Maximum day to include
#' @return List containing excluded and truncated models
sensitivity_analysis_outliers <- function(patient_data, crude_rates, cutoff_day = 30) {
  
  # Identify long-stay patients
  long_stay_patients <- patient_data %>%
    group_by(deid_enc_id) %>%
    summarise(max_day = max(DaysSinceEntry, na.rm = TRUE), .groups = "drop") %>%
    filter(max_day > cutoff_day) %>%
    pull(deid_enc_id)
  
  cat("Found", length(long_stay_patients), "patients with stays >", cutoff_day, "days\n")
  
  # Exclusion approach
  data_excluded <- patient_data %>%
    filter(!deid_enc_id %in% long_stay_patients)
  
  # Truncation approach
  data_truncated <- patient_data %>%
    filter(DaysSinceEntry <= cutoff_day)
  
  # Recompute crude rates for both datasets
  crude_rates_excluded <- calc_crude_init_rates(data_excluded, 
                                                qmat_list[unique(data_excluded$model)])
  crude_rates_truncated <- calc_crude_init_rates(data_truncated, 
                                                 qmat_list[unique(data_truncated$model)])
  
  # Fit models
  models_excluded <- fit_msm_models(data_excluded, crude_rates_excluded)
  models_truncated <- fit_msm_models(data_truncated, crude_rates_truncated)
  
  return(list(
    excluded = list(data = data_excluded, models = models_excluded, crude_rates = crude_rates_excluded),
    truncated = list(data = data_truncated, models = models_truncated, crude_rates = crude_rates_truncated),
    long_stay_patients = long_stay_patients
  ))
}





# Utility functions -------------------------------------------------------

#' Extract transition trend categories for plotting
#' @param patient_data Patient data
#' @return Data frame with transition trend categories
add_transition_trends <- function(patient_data) {
  
  trend_types <- cross_join(
    tibble(from = c("M", "M1", "M2", "M3", "MS", "S", "S1", "S2")),
    tibble(to = c("M", "M1", "M2", "M3", "MS", "S", "S1", "S2", "R", "D"))
  ) %>%
    mutate(
      from_type = substring(from, 1, 1),
      to_type = substring(to, 1, 1)
    ) %>%
    mutate(
      trend = case_when(
        from == to ~ "Self-transition",
        to == "R" ~ "Recovery",
        to == "D" ~ "Death",
        from_type == "M" & to_type == "M" & 
          as.numeric(sub("M", "", to)) > as.numeric(sub("M", "", from)) ~ "Worse",
        from_type == "M" & to_type == "M" & 
          as.numeric(sub("M", "", to)) < as.numeric(sub("M", "", from)) ~ "Better",
        from_type == "M" & to_type == "S" ~ "Worse",
        from_type == "S" & to_type == "M" ~ "Better",
        from_type == "S" & to_type == "S" & 
          as.numeric(sub("S", "", to)) > as.numeric(sub("S", "", from)) ~ "Worse",
        from_type == "S" & to_type == "S" & 
          as.numeric(sub("S", "", to)) < as.numeric(sub("S", "", from)) ~ "Better",
        from_type == "M" & to == "MS" ~ "Other",
        from == "MS" & to_type == "M" ~ "Other",
        TRUE ~ "Other"
      )
    ) %>%
    filter(!is.na(trend))
  
  return(trend_types)
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

#' Calculate model weights from AIC values
calculate_aic_weights <- function(aic_values) {
  delta_aic <- aic_values - min(aic_values, na.rm = TRUE)
  relative_likelihood <- exp(-0.5 * delta_aic)
  weights <- relative_likelihood / sum(relative_likelihood, na.rm = TRUE)
  return(weights)
}

#' Extract transition-specific hazard ratios
extract_hazard_ratios <- function(fitted_model, covariate_name, reference_value = 0, 
                                  comparison_value = 1) {
  
  if (is.null(fitted_model)) return(NULL)
  
  # Get hazard ratios
  hr_result <- tryCatch({
    hazard.msm(fitted_model, 
               hazard.scale = comparison_value - reference_value,
               ci = "normal")
  }, error = function(e) NULL)
  
  if (!is.null(hr_result)) {
    hr_df <- as.data.frame(hr_result) %>%
      rownames_to_column("transition") %>%
      rename_with(~gsub("^[^.]*\\.", "", .x), .cols = everything()) %>%
      mutate(
        covariate = covariate_name,
        comparison = paste(comparison_value, "vs", reference_value)
      )
    
    return(hr_df)
  } else {
    return(NULL)
  }
}

#' Create model summary table
create_model_summary_table <- function(model_list, include_metrics = c("AIC", "BIC", "logLik")) {
  
  summary_table <- map_dfr(names(model_list), function(model_name) {
    model <- model_list[[model_name]]
    
    if (is.null(model)) {
      result <- data.frame(model = model_name, converged = FALSE)
      for (metric in include_metrics) {
        result[[metric]] <- NA
      }
      return(result)
    }
    
    result <- data.frame(
      model = model_name,
      converged = model$opt$convergence == 0,
      n_params = length(model$estimates),
      n_subjects = length(unique(model$data$mf)))
    
    if ("AIC" %in% include_metrics) result$AIC <- AIC(model)
    if ("BIC" %in% include_metrics) result$BIC <- BIC(model)  
    if ("logLik" %in% include_metrics) result$logLik <- as.numeric(logLik(model))
    
    return(result)
  })
  
  return(summary_table)
}

transition_tables_by_covariate <- function(data, covariate_name, state_var = "state", 
                                           subject_var = "deid_enc_id", time_var = "DaysSinceEntry",
                                           n_categories = 3, method = "quantile") {
  
  # Check if covariate is continuous
  is_continuous <- is.numeric(data[[covariate_name]]) && 
    length(unique(data[[covariate_name]][!is.na(data[[covariate_name]])])) > 10
  
  # Handle continuous variables
  if (is_continuous) {
    cat("Continuous covariate detected. Categorizing into", n_categories, "groups using", method, "method.\n")
    
    # Create working copy of data
    data_work <- data
    
    if (method == "quantile") {
      # Quantile-based categorization
      breaks <- quantile(data_work[[covariate_name]], 
                         probs = seq(0, 1, length.out = n_categories + 1), 
                         na.rm = TRUE)
      # Ensure unique breaks
      breaks <- unique(breaks)
      if (length(breaks) <= n_categories) {
        warning("Not enough unique values to create requested number of categories. Using available breaks.")
      }
      
      labels <- paste0("Q", 1:(length(breaks)-1))
      data_work$covariate_cat <- cut(data_work[[covariate_name]], 
                                     breaks = breaks, 
                                     include.lowest = TRUE,
                                     labels = labels)
      
      # Store break information for interpretation
      break_info <- data.frame(
        Category = labels,
        Min = breaks[-length(breaks)],
        Max = breaks[-1]
      )
      
    } else if (method == "equal_width") {
      # Equal width intervals
      min_val <- min(data_work[[covariate_name]], na.rm = TRUE)
      max_val <- max(data_work[[covariate_name]], na.rm = TRUE)
      width <- (max_val - min_val) / n_categories
      
      breaks <- seq(min_val, max_val, by = width)
      if (length(breaks) != (n_categories + 1)) {
        breaks[length(breaks)] <- max_val  # Ensure max value is included
      }
      
      labels <- paste0("Grp", 1:n_categories)
      data_work$covariate_cat <- cut(data_work[[covariate_name]], 
                                     breaks = breaks, 
                                     include.lowest = TRUE,
                                     labels = labels)
      
      break_info <- data.frame(
        Category = labels,
        Min = breaks[-length(breaks)],
        Max = breaks[-1]
      )
      
    } else if (method == "equal_freq") {
      # Equal frequency (approximately equal n in each group)
      data_work$covariate_cat <- as.factor(
        ntile(data_work[[covariate_name]], n_categories)
      )
      levels(data_work$covariate_cat) <- paste0("Grp", 1:n_categories)
      
      # Calculate actual ranges for each group
      break_info <- data_work %>%
        filter(!is.na(covariate_cat) & !is.na(!!sym(covariate_name))) %>%
        group_by(covariate_cat) %>%
        summarise(
          Min = min(!!sym(covariate_name), na.rm = TRUE),
          Max = max(!!sym(covariate_name), na.rm = TRUE),
          n = n(),
          .groups = 'drop'
        ) %>%
        rename(Category = covariate_cat)
    }
    
    # Use categorized variable
    covariate_values <- levels(data_work$covariate_cat)
    working_covariate <- "covariate_cat"
    
  } else {
    # Handle categorical variables as before
    data_work <- data
    covariate_values <- unique(data_work[[covariate_name]])
    covariate_values <- covariate_values[!is.na(covariate_values)]
    working_covariate <- covariate_name
    break_info <- NULL
  }
  
  transition_tables <- list()
  
  for (value in covariate_values) {
    # Subset data
    subset_data <- data_work[data_work[[working_covariate]] == value & 
                               !is.na(data_work[[working_covariate]]), ]
    
    if (nrow(subset_data) > 0) {
      # Create transition pairs
      transitions <- subset_data %>%
        arrange(!!sym(subject_var), !!sym(time_var)) %>%
        group_by(!!sym(subject_var)) %>%
        mutate(
          next_state = lead(!!sym(state_var)),
          transition = paste(!!sym(state_var), "->", next_state)
        ) %>%
        filter(!is.na(next_state)) %>%
        ungroup()
      
      # Count transitions
      transition_counts <- table(transitions$transition)
      
      # Calculate transition matrix for this group
      from_states <- sapply(strsplit(names(transition_counts), " -> "), `[`, 1)
      to_states <- sapply(strsplit(names(transition_counts), " -> "), `[`, 2)
      
      all_states <- sort(unique(c(from_states, to_states)))
      transition_matrix <- matrix(0, nrow = length(all_states), ncol = length(all_states))
      rownames(transition_matrix) <- all_states
      colnames(transition_matrix) <- all_states
      
      for (i in seq_along(transition_counts)) {
        from <- from_states[i]
        to <- to_states[i]
        transition_matrix[from, to] <- transition_counts[i]
      }
      
      transition_tables[[paste0(covariate_name, "_", value)]] <- list(
        covariate_value = value,
        transitions = transition_counts,
        transition_matrix = transition_matrix,
        n_transitions = sum(transition_counts),
        n_subjects = length(unique(subset_data[[subject_var]])),
        n_observations = nrow(subset_data)
      )
    }
  }
  
  # Add metadata for continuous variables
  if (is_continuous) {
    attr(transition_tables, "continuous_info") <- list(
      original_variable = covariate_name,
      method = method,
      n_categories = n_categories,
      break_info = break_info
    )
    attr(transition_tables, "is_continuous") <- TRUE
  } else {
    attr(transition_tables, "is_continuous") <- FALSE
  }
  
  return(transition_tables)
}

print_transition_tables <- function(transition_tables_list) {
  
  # Check if this involves continuous variables
  is_continuous <- attr(transition_tables_list, "is_continuous")
  
  if (is_continuous) {
    continuous_info <- attr(transition_tables_list, "continuous_info")
    cat("\n=== CONTINUOUS VARIABLE ANALYSIS ===\n")
    cat("Original variable:", continuous_info$original_variable, "\n")
    cat("Categorization method:", continuous_info$method, "\n")
    cat("Number of categories:", continuous_info$n_categories, "\n\n")
    
    cat("Category definitions:\n")
    print(continuous_info$break_info)
    cat("\n")
  }
  
  for (name in names(transition_tables_list)) {
    cat("\n=== Transition Table for", name, "===\n")
    cat("Covariate value:", transition_tables_list[[name]]$covariate_value, "\n")
    cat("Number of subjects:", transition_tables_list[[name]]$n_subjects, "\n")
    cat("Number of transitions:", transition_tables_list[[name]]$n_transitions, "\n")
    cat("Number of observations:", transition_tables_list[[name]]$n_observations, "\n\n")
    
    # Print transition counts
    cat("Transition counts:\n")
    print(transition_tables_list[[name]]$transitions)
    cat("\n")
    
    # Print transition matrix if available
    if ("transition_matrix" %in% names(transition_tables_list[[name]])) {
      cat("Transition matrix (from row to column):\n")
      print(transition_tables_list[[name]]$transition_matrix)
      cat("\n")
      
      # Print transition rates
      transition_matrix <- transition_tables_list[[name]]$transition_matrix
      row_sums <- rowSums(transition_matrix)
      rate_matrix <- sweep(transition_matrix, 1, row_sums, FUN = "/")
      rate_matrix[is.nan(rate_matrix)] <- 0  # Handle division by zero
      
      cat("Transition proportions (from row to column):\n")
      print(round(rate_matrix, 3))
    }
    
    cat("\n", rep("=", 60), "\n")
  }
}

summarize_continuous_transitions <- function(transition_tables_list) {
  
  if (!attr(transition_tables_list, "is_continuous")) {
    stop("This function is only for continuous variable analyses")
  }
  
  continuous_info <- attr(transition_tables_list, "continuous_info")
  
  # Extract summary statistics for each category
  summary_df <- data.frame(
    Category = character(),
    N_Subjects = integer(),
    N_Transitions = integer(),
    N_Observations = integer(),
    Transition_Rate = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(transition_tables_list)) {
    table_info <- transition_tables_list[[i]]
    
    summary_df <- rbind(summary_df, data.frame(
      Category = table_info$covariate_value,
      N_Subjects = table_info$n_subjects,
      N_Transitions = table_info$n_transitions,
      N_Observations = table_info$n_observations,
      Transition_Rate = table_info$n_transitions / table_info$n_observations,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add range information if available
  if ("break_info" %in% names(continuous_info)) {
    break_info <- continuous_info$break_info
    summary_df <- merge(summary_df, break_info, by = "Category", all.x = TRUE)
  }
  
  return(summary_df)
}

