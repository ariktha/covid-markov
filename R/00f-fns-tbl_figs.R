#' Extract transition trend categories for plotting
#' @param patient_data Patient data
#' @return Data frame with transition trend categories
add_transition_trends <- function(patient_data) {
  
  
  
  patient_data <- patient_data %>%
    left_join(trend_types)
  
  return(patient_data)
}

#' Format confidence intervals
format_ci <- function(estimate, lower, upper, digits = 2) {
  paste0(
    round(estimate, digits), 
    " (95% CI: ", 
    round(lower, digits), 
    "—", 
    round(upper, digits), 
    ")"
  )
}

add_trend_type <- function(df, from_name = "from_state", to_name = "to_state") {
  df %>%
    left_join(trend_types, by = join_by(from_name == "from", to_name == "to"))
}

add_covariate_labels <- function(df, cov_col_name = "covariate") {
  cov_sub <- unique(df[[cov_col_name]])
  df %>%
    mutate(covariate = factor())
}

## Hazard ratio plotting functions

#' Plot smooth hazard ratio curves for spline variables
#' 
#' Creates smooth curves showing how hazard ratios change across the range 
#' of a spline variable (e.g., age vs hazard ratio).
#'
#' @param fitted_msm_models Your fitted models
#' @param patient_data Your patient data 
#' @param spline_var Variable name (e.g., "age", "BMI")
#' @param transitions Which transitions to plot (optional)
#' @param reference_value Reference value for HR=1 (default: median)
#' @param n_points Number of points along the curve (default: 50)
#' @return ggplot object
plot_spline_hazard_ratios <- function(fitted_msm_models,
                                      patient_data,
                                      spline_var = "age",
                                      transitions = NULL,
                                      reference_value = NULL,
                                      n_points = 50) {
  
  library(ggplot2)
  library(dplyr)
  
  # Find models with this spline variable
  spline_results <- list()
  
  for (model_name in names(fitted_msm_models)) {
    for (formula_name in names(fitted_msm_models[[model_name]])) {
      model_entry <- fitted_msm_models[[model_name]][[formula_name]]
      
      # Skip failed models
      if (is.null(model_entry) || 
          (is.list(model_entry) && !is.null(model_entry$status) && 
           model_entry$status == "failed")) {
        next
      }
      
      # Extract metadata and check for splines
      metadata <- if (is.list(model_entry) && !is.null(model_entry$metadata)) {
        model_entry$metadata
      } else {
        list()
      }
      
      # Check if this model has our spline variable
      if (!is.null(metadata$spline_config) && 
          !is.null(metadata$spline_config$spline_metadata) &&
          spline_var %in% names(metadata$spline_config$spline_metadata)) {
        
        fitted_model <- model_entry$fitted_model
        spline_meta <- metadata$spline_config$spline_metadata[[spline_var]]
        
        # Get variable range from data
        model_structure <- metadata$model_structure
        model_data <- patient_data %>% filter(model == model_structure)
        var_values <- model_data[[spline_var]]
        var_values <- var_values[!is.na(var_values)]
        
        if (length(var_values) == 0) next
        
        var_range <- range(var_values, na.rm = TRUE)
        
        # Set reference value
        ref_value <- if (!is.null(reference_value)) {
          reference_value
        } else {
          median(var_values, na.rm = TRUE)
        }
        
        # Create evaluation grid
        eval_grid <- seq(var_range[1], var_range[2], length.out = n_points)
        
        # Calculate HRs at each point
        hr_data <- calculate_spline_hrs_at_points(
          fitted_model = fitted_model,
          spline_var = spline_var,
          eval_points = eval_grid,
          reference_value = ref_value,
          spline_meta = spline_meta
        )
        
        if (!is.null(hr_data)) {
          hr_data$model_name <- model_name
          hr_data$formula <- formula_name
          spline_results[[length(spline_results) + 1]] <- hr_data
        }
      }
    }
  }
  
  if (length(spline_results) == 0) {
    stop("No spline models found for variable: ", spline_var)
  }
  
  # Combine results
  plot_data <- do.call(rbind, spline_results)
  
  # Filter transitions if specified
  if (!is.null(transitions)) {
    plot_data <- plot_data %>% filter(transition %in% transitions)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = !!sym(spline_var), y = hazard_ratio)) +
    geom_line(aes(color = transition), size = 1) +
    geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper, fill = transition), 
                alpha = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_y_log10() +
    labs(
      title = paste("Hazard Ratio vs", tools::toTitleCase(spline_var)),
      x = tools::toTitleCase(spline_var),
      y = "Hazard Ratio (log scale)",
      caption = paste("Reference:", spline_var, "=", round(reference_value %||% 
                                                             median(patient_data[[spline_var]], na.rm = TRUE), 1))
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # Facet by model if multiple models
  if (length(unique(plot_data$model_name)) > 1) {
    p <- p + facet_wrap(~model_name)
  }
  
  return(p)
}

#' Calculate hazard ratios at specific points for spline variables
#' 
#' Internal function to calculate HRs relative to reference value
#' using the fitted spline coefficients.
#'
#' @param fitted_model MSM fitted model object
#' @param spline_var Name of spline variable
#' @param eval_points Vector of values to evaluate
#' @param reference_value Reference value (HR = 1)
#' @param spline_meta Spline metadata from model
#' @return Data frame with hazard ratios and CIs
calculate_spline_hrs_at_points <- function(fitted_model,
                                           spline_var,
                                           eval_points,
                                           reference_value,
                                           spline_meta) {
  
  tryCatch({
    # Get spline basis functions
    spline_type <- spline_meta$type
    spline_df <- spline_meta$df
    knots <- spline_meta$knots
    boundary_knots <- spline_meta$boundary_knots
    
    # Create spline basis for evaluation points
    if (spline_type == "ns") {
      eval_basis <- splines::ns(eval_points, df = spline_df, 
                                knots = knots, Boundary.knots = boundary_knots)
      ref_basis <- splines::ns(reference_value, df = spline_df,
                               knots = knots, Boundary.knots = boundary_knots)
    } else if (spline_type == "bs") {
      eval_basis <- splines::bs(eval_points, df = spline_df,
                                knots = knots, Boundary.knots = boundary_knots)
      ref_basis <- splines::bs(reference_value, df = spline_df,
                               knots = knots, Boundary.knots = boundary_knots)
    } else {
      warning("Unsupported spline type: ", spline_type)
      return(NULL)
    }
    
    # Get coefficient estimates and variance-covariance matrix
    coeffs <- fitted_model$estimates
    vcov_matrix <- fitted_model$covmat
    
    # Find spline coefficient indices
    spline_terms <- spline_meta$spline_terms
    spline_indices <- which(names(coeffs) %in% spline_terms)
    
    if (length(spline_indices) == 0) {
      warning("Spline coefficients not found")
      return(NULL)
    }
    
    # Calculate log hazard ratios
    results <- list()
    
    # Get transition names from qmatrix
    qmat <- fitted_model$qmodel$qmatrix
    transitions <- expand.grid(
      from = rownames(qmat),
      to = colnames(qmat)
    ) %>%
      filter(qmat[cbind(from, to)] != 0) %>%
      mutate(transition = paste(from, "->", to)) %>%
      pull(transition)
    
    # For each transition
    for (trans_idx in seq_along(transitions)) {
      transition <- transitions[trans_idx]
      
      # Get coefficients for this transition
      # MSM stores transition-specific coefficients
      trans_coeff_indices <- spline_indices + (trans_idx - 1) * length(coeffs) / length(transitions)
      
      # Calculate linear predictor difference (eval - ref)
      linear_pred_eval <- eval_basis %*% coeffs[trans_coeff_indices]
      linear_pred_ref <- ref_basis %*% coeffs[trans_coeff_indices]
      
      log_hr <- as.vector(linear_pred_eval - linear_pred_ref)
      hr <- exp(log_hr)
      
      # Calculate confidence intervals using delta method
      basis_diff <- eval_basis - ref_basis[rep(1, nrow(eval_basis)), ]
      se_log_hr <- sqrt(diag(basis_diff %*% vcov_matrix[trans_coeff_indices, trans_coeff_indices] %*% t(basis_diff)))
      
      hr_lower <- exp(log_hr - 1.96 * se_log_hr)
      hr_upper <- exp(log_hr + 1.96 * se_log_hr)
      
      # Store results
      result_df <- data.frame(
        age = eval_points,  # Use 'age' as column name for plotting
        hazard_ratio = hr,
        hr_lower = hr_lower,
        hr_upper = hr_upper,
        log_hr = log_hr,
        transition = transition,
        stringsAsFactors = FALSE
      )
      
      # Rename the variable column to match spline_var
      names(result_df)[1] <- spline_var
      
      results[[trans_idx]] <- result_df
    }
    
    return(do.call(rbind, results))
    
  }, error = function(e) {
    warning("Error calculating spline HRs: ", e$message)
    return(NULL)
  })
}

#' Simple wrapper for age vs hazard ratio plots
#' 
#' Convenience function specifically for age plots
#'
#' @param fitted_msm_models Your fitted models
#' @param patient_data Your patient data
#' @param reference_age Reference age (default: 65)
#' @param transitions Specific transitions to plot (optional)
#' @return ggplot object
plot_age_hazard_ratios <- function(fitted_msm_models,
                                   patient_data,
                                   reference_age = 65,
                                   transitions = NULL) {
  
  plot_spline_hazard_ratios(
    fitted_msm_models = fitted_msm_models,
    patient_data = patient_data,
    spline_var = "age",
    reference_value = reference_age,
    transitions = transitions,
    n_points = 50
  )
}

# Example usage:
#
# # Plot age vs HR for all transitions
# plot_age_hazard_ratios(fitted_univar_models, pt_stg)
#
# # Plot age vs HR for specific transitions
# plot_age_hazard_ratios(fitted_univar_models, pt_stg, 
#                       transitions = c("M -> D", "M1 -> D"))
#
# # Plot BMI vs HR
# plot_spline_hazard_ratios(fitted_univar_models, pt_stg, 
#                          spline_var = "BMI", reference_value = 25)