# Additional utility functions for MSM analysis
# Add to end of 00-functions.R

# ============================================================================
# ADVANCED COVARIATE TESTING FUNCTIONS
# ============================================================================

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

# ============================================================================
# SIMULATION AND PREDICTION FUNCTIONS
# ============================================================================

#' Simulate patient trajectories from fitted model
#' @param fitted_model A fitted msm model
#' @param n_patients Number of patients to simulate
#' @param max_time Maximum follow-up time
#' @param start_state Starting state (default: first state)
#' @param covariates Named list of covariate values
#' @return Data frame with simulated trajectories
simulate_trajectories <- function(fitted_model, n_patients = 100, max_time = 30, 
                                  start_state = NULL, covariates = NULL) {
  
  if (is.null(fitted_model)) return(NULL)
  
  # Get model states
  states <- fitted_model$qmodel$state.names
  if (is.null(start_state)) start_state <- states[1]
  
  # Simulate trajectories
  sim_results <- list()
  
  for (i in 1:n_patients) {
    sim_traj <- tryCatch({
      sim.msm(fitted_model, 
              qmatrix = qmatrix.msm(fitted_model, covariates = covariates),
              maxtime = max_time,
              start = match(start_state, states))
    }, error = function(e) NULL)
    
    if (!is.null(sim_traj)) {
      sim_results[[i]] <- data.frame(
        patient_id = i,
        time = sim_traj$times,
        state = states[sim_traj$states]
      )
    }
  }
  
  if (length(sim_results) > 0) {
    return(bind_rows(sim_results))
  } else {
    return(NULL)
  }
}

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

# ============================================================================
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

# ============================================================================
# REAL-TIME RISK CALCULATOR FUNCTIONS  
# ============================================================================

#' Create risk calculator function from fitted model
#' @param fitted_model A fitted msm model
#' @param prediction_times Vector of time points for predictions
#' @return Function that calculates risks for new patients
create_risk_calculator <- function(fitted_model, prediction_times = c(1, 3, 7, 14)) {
  
  if (is.null(fitted_model)) {
    return(function(...) stop("Model not available"))
  }
  
  states <- fitted_model$qmodel$state.names
  
  risk_function <- function(patient_covariates, current_state = states[1], 
                            time_in_current_state = 0) {
    
    # Input validation
    if (!current_state %in% states) {
      stop("Invalid current state: ", current_state)
    }
    
    # Calculate risk predictions for each time horizon
    risk_predictions <- map_dfr(prediction_times, function(t) {
      
      # Get transition probability matrix
      pmat <- tryCatch({
        pmatrix.msm(fitted_model, t = t, covariates = patient_covariates)
      }, error = function(e) {
        warning("Failed to calculate probabilities for time ", t)
        return(NULL)
      })
      
      if (!is.null(pmat)) {
        current_state_idx <- match(current_state, states)
        probs <- pmat[current_state_idx, ]
        
        data.frame(
          time_horizon = t,
          prob_death = probs["D"] %||% 0,
          prob_recovery = probs["R"] %||% 0,
          prob_severe = sum(probs[grepl("S", names(probs))], na.rm = TRUE),
          prob_moderate = sum(probs[grepl("M", names(probs))], na.rm = TRUE),
          prob_stay_current = probs[current_state] %||% 0
        )
      } else {
        data.frame(
          time_horizon = t,
          prob_death = NA,
          prob_recovery = NA,
          prob_severe = NA,
          prob_moderate = NA,
          prob_stay_current = NA
        )
      }
    })
    
    # Add risk categories
    risk_predictions <- risk_predictions %>%
      mutate(
        death_risk_category = case_when(
          prob_death < 0.05 ~ "Low",
          prob_death < 0.15 ~ "Moderate", 
          prob_death < 0.30 ~ "High",
          TRUE ~ "Very High"
        ),
        progression_risk = prob_severe + prob_death,
        improvement_chance = prob_recovery + ifelse(current_state != "M", prob_moderate, 0)
      )
    
    # Return results with metadata
    return(list(
      predictions = risk_predictions,
      current_state = current_state,
      time_in_state = time_in_current_state,
      patient_covariates = patient_covariates,
      model_states = states,
      calculation_time = Sys.time()
    ))
  }
  
  return(risk_function)
}

# ============================================================================
# VISUALIZATION HELPER FUNCTIONS
# ============================================================================

#' Create transition flow diagram
#' @param qmatrix Transition intensity matrix
#' @param model_name Name of the model
#' @return ggplot object
plot_transition_diagram <- function(qmatrix, model_name = "") {
  
  # Extract non-zero transitions
  transitions <- as.data.frame(as.table(qmatrix)) %>%
    rename(from = Var1, to = Var2, rate = Freq) %>%
    filter(rate > 0, from != to) %>%
    mutate(
      rate_category = case_when(
        rate < 0.01 ~ "Very Low",
        rate < 0.05 ~ "Low", 
        rate < 0.1 ~ "Moderate",
        TRUE ~ "High"
      )
    )
  
  # Create state positions (simple layout)
  states <- unique(c(as.character(transitions$from), as.character(transitions$to)))
  n_states <- length(states)
  
  # Arrange states in a circle
  angles <- seq(0, 2*pi, length.out = n_states + 1)[1:n_states]
  state_positions <- data.frame(
    state = states,
    x = cos(angles),
    y = sin(angles)
  )
  
  # Add positions to transitions
  transitions <- transitions %>%
    left_join(state_positions, by = c("from" = "state")) %>%
    rename(x_from = x, y_from = y) %>%
    left_join(state_positions, by = c("to" = "state")) %>%
    rename(x_to = x, y_to = y)
  
  # Create plot
  p <- ggplot() +
    # Draw transitions as arrows
    geom_segment(data = transitions,
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to, 
                     color = rate_category, size = rate),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
    # Draw states as circles
    geom_point(data = state_positions, aes(x = x, y = y), 
               size = 15, color = "white", fill = "lightblue", shape = 21, stroke = 2) +
    # Add state labels
    geom_text(data = state_positions, aes(x = x, y = y, label = state), 
              size = 4, fontface = "bold") +
    # Styling
    scale_color_manual(values = c("Very Low" = "gray70", "Low" = "steelblue", 
                                  "Moderate" = "orange", "High" = "red")) +
    scale_size_continuous(range = c(0.5, 2)) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(
      title = paste("Transition Diagram:", model_name),
      color = "Transition Rate",
      size = "Rate Magnitude"
    )
  
  return(p)
}

#' Create survival curves from multi-state model
#' @param fitted_model A fitted msm model
#' @param times Vector of time points
#' @param covariates List of covariate scenarios
#' @return ggplot object with survival curves
plot_survival_curves <- function(fitted_model, times = seq(0, 30, 1), covariates = list()) {
  
  if (is.null(fitted_model)) return(NULL)
  
  states <- fitted_model$qmodel$state.names
  absorbing_states <- c("D", "R")
  
  # If no covariates specified, use mean values
  if (length(covariates) == 0) {
    covariates <- list("Mean values" = NULL)
  }
  
  # Calculate survival probabilities for each covariate scenario
  survival_data <- map_dfr(names(covariates), function(scenario_name) {
    scenario_covs <- covariates[[scenario_name]]
    
    scenario_probs <- map_dfr(times, function(t) {
      if (t == 0) {
        # Starting probabilities (assume start in first non-absorbing state)
        start_state <- states[!states %in% absorbing_states][1]
        probs <- setNames(rep(0, length(states)), states)
        probs[start_state] <- 1
      } else {
        # Get transition probabilities
        pmat <- tryCatch({
          pmatrix.msm(fitted_model, t = t, covariates = scenario_covs)
        }, error = function(e) NULL)
        
        if (!is.null(pmat)) {
          start_state_idx <- which(!states %in% absorbing_states)[1]
          probs <- pmat[start_state_idx, ]
        } else {
          probs <- setNames(rep(NA, length(states)), states)
        }
      }
      
      data.frame(
        time = t,
        scenario = scenario_name,
        prob_alive = 1 - (probs["D"] %||% 0),
        prob_recovered = probs["R"] %||% 0,
        prob_hospitalized = 1 - sum(probs[absorbing_states], na.rm = TRUE)
      )
    })
    
    return(scenario_probs)
  })
  
  # Create survival plot
  survival_long <- survival_data %>%
    pivot_longer(cols = starts_with("prob_"), names_to = "outcome", values_to = "probability") %>%
    mutate(
      outcome = case_when(
        outcome == "prob_alive" ~ "Survival (alive)",
        outcome == "prob_recovered" ~ "Recovery",
        outcome == "prob_hospitalized" ~ "Still hospitalized"
      )
    )
  
  p <- ggplot(survival_long, aes(x = time, y = probability, color = outcome, linetype = scenario)) +
    geom_line(size = 1) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_color_manual(values = c("Survival (alive)" = "darkgreen", 
                                  "Recovery" = "blue", 
                                  "Still hospitalized" = "orange")) +
    labs(
      title = "Survival and Outcome Curves",
      x = "Time (days)",
      y = "Probability",
      color = "Outcome",
      linetype = "Scenario"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# ============================================================================
# HELPER UTILITY FUNCTIONS
# ============================================================================

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