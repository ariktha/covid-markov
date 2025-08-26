
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

# Model Fit Functions -------------------------------------------------------

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

# Model Tidying and Extraction Functions --------------------------------------

#' Generate comprehensive model assessment report
#' @param fitted_models List of fitted models
#' @param patient_data Patient data
#' @param crude_rates List of crude rates for cross-validation
#' @param model_names Vector of model names to include in report
#' @param print_summary Logical, whether to print summary table
#' @param create_plots Logical, whether to generate plots
#' @param plot_patchwork Logical, whether to combine plots using patchwork
#' @return List containing all assessment results
generate_model_assessment_report <- function(fitted_models, patient_data, crude_rates,
                                             model_names = names(fitted_models),
                                             print_summary = TRUE,
                                             create_plots = TRUE, 
                                             plot_patchwork = FALSE) {
  
  cat("Generating comprehensive model assessment report...\n")
  
  report <- list()
  
  # Model fit statistics
  cat("1. Calculating model fit statistics...\n")
  report$model_fit <- create_model_summary_table(fitted_models[model_names])
  
  # Extract model formulas
  model_formulas <- map_dfr(model_names, function(name) {
    model <- fitted_models[[name]]
    if (is.null(model)) {
      return(data.frame(
        model = name,
        formula = "Failed to fit",
        covariates = "None"
      ))
    }
    
    # Extract main formula
    main_formula <- deparse(model$call$formula)
    
    # Extract covariates formula if present
    cov_formula <- if (!is.null(model$call$covariates)) {
      deparse(model$call$covariates)
    } else "None"
    
    data.frame(
      model = name,
      formula = main_formula,
      covariates = cov_formula
    )
  })
  
  # Predictive performance  
  cat("2. Calculating predictive performance...\n")
  report$predictive_performance <- calculate_predictive_performance(
    fitted_models[model_names], patient_data, crude_rates
  )
  
  # Combine model fit and predictive performance for summary
  if (!is.null(report$model_fit) && !is.null(report$predictive_performance)) {
    summary_table <- report$model_fit %>%
      left_join(model_formulas, by = "model") %>%
      left_join(
        report$predictive_performance %>%
          select(model, total_los_mae_cv, days_severe_mae_cv, death_auc_cv, severe_auc_cv),
        by = "model"
      ) %>%
      arrange(AIC) %>%
      select(model, formula, covariates, everything())
    
    report$summary_table <- summary_table
  }
  
  # State prevalence residuals
  cat("3. Calculating state prevalence residuals...\n")
  report$state_residuals <- map(model_names, function(name) {
    if (!is.null(fitted_models[[name]])) {
      model_data <- patient_data %>% filter(model == name)
      calculate_state_prevalence_residuals(fitted_models[[name]], model_data)
    } else NULL
  })
  names(report$state_residuals) <- model_names
  
  # Transition residuals
  cat("4. Calculating transition residuals...\n") 
  report$transition_residuals <- map(model_names, function(name) {
    if (!is.null(fitted_models[[name]])) {
      model_data <- patient_data %>% filter(model == name)
      calculate_transition_residuals(fitted_models[[name]], model_data)
    } else NULL
  })
  names(report$transition_residuals) <- model_names
  
  # Patient outliers
  cat("5. Identifying patient outliers...\n")
  report$patient_outliers <- map(model_names, function(name) {
    if (!is.null(fitted_models[[name]])) {
      model_data <- patient_data %>% filter(model == name)
      identify_outlier_patients(fitted_models[[name]], model_data)
    } else NULL
  })
  names(report$patient_outliers) <- model_names
  
  # Calibration
  cat("6. Calculating calibration...\n")
  report$calibration <- calculate_calibration_plots(
    fitted_models[model_names], patient_data
  )
  
  # Cumulative outcomes
  cat("7. Calculating cumulative outcomes...\n")
  report$cumulative_outcomes <- calculate_cumulative_outcomes(
    fitted_models[model_names], patient_data
  )
  
  # Sojourn times
  cat("8. Extracting sojourn times...\n")
  report$sojourn_times <- tidy_msm_sojourns(fitted_models[model_names])
  
  # Transition probabilities
  cat("9. Extracting transition probabilities...\n")
  report$transition_probs <- tidy_msm_pmats(fitted_models[model_names], 
                                            t_values = c(1, 5, 14, 30))
  
  # Transition intensities
  cat("10. Calculating transition intensities...\n")
  report$transition_intensities <- map(model_names, function(name) {
    if (!is.null(fitted_models[[name]])) {
      calculate_transition_intensities(fitted_models[[name]])
    } else NULL
  })
  names(report$transition_intensities) <- model_names
  
  # Create plots if requested
  if (create_plots) {
    cat("11. Generating plots...\n")
    report$plots <- list()
    
    # Calibration plots for each model
    calibration_plots <- map(model_names, function(name) {
      if (!is.null(fitted_models[[name]])) {
        tryCatch({
          plot_calibration(report$calibration, model_name = name)
        }, error = function(e) {
          warning(paste("Could not create calibration plot for", name, ":", e$message))
          NULL
        })
      } else NULL
    })
    names(calibration_plots) <- paste0(model_names, "_calibration")
    report$plots$calibration <- calibration_plots
    
    # State prevalence residual plots
    prevalence_plots <- map(model_names, function(name) {
      if (!is.null(report$state_residuals[[name]])) {
        tryCatch({
          plot_state_prevalence_residuals(report$state_residuals[[name]])
        }, error = function(e) {
          warning(paste("Could not create prevalence plot for", name, ":", e$message))
          NULL
        })
      } else NULL
    })
    names(prevalence_plots) <- paste0(model_names, "_prevalence")
    report$plots$state_prevalence <- prevalence_plots
    
    # Transition residual plots (time)
    transition_time_plots <- map(model_names, function(name) {
      if (!is.null(report$transition_residuals[[name]])) {
        tryCatch({
          plot_transition_residuals(report$transition_residuals[[name]], plot_type = "time")
        }, error = function(e) {
          warning(paste("Could not create transition time plot for", name, ":", e$message))
          NULL
        })
      } else NULL
    })
    names(transition_time_plots) <- paste0(model_names, "_transition_time")
    report$plots$transition_time <- transition_time_plots
    
    # Patient outlier plots
    outlier_plots <- map(model_names, function(name) {
      if (!is.null(report$patient_outliers[[name]])) {
        tryCatch({
          plot_patient_residuals(report$patient_outliers[[name]], plot_type = "outlier_score")
        }, error = function(e) {
          warning(paste("Could not create outlier plot for", name, ":", e$message))
          NULL
        })
      } else NULL
    })
    names(outlier_plots) <- paste0(model_names, "_outliers")
    report$plots$patient_outliers <- outlier_plots
    
    # Cumulative outcome plots
    if (!is.null(report$cumulative_outcomes)) {
      cumulative_plot <- tryCatch({
        library(ggplot2)
        
        report$cumulative_outcomes %>%
          select(model, time, prob_death, prob_recovery, observed_death_rate, observed_recovery_rate) %>%
          pivot_longer(cols = c(prob_death, prob_recovery, observed_death_rate, observed_recovery_rate),
                       names_to = "outcome_type", values_to = "probability") %>%
          mutate(
            outcome = case_when(
              grepl("death", outcome_type) ~ "Death",
              grepl("recovery", outcome_type) ~ "Recovery",
              TRUE ~ outcome_type
            ),
            data_type = ifelse(grepl("observed", outcome_type), "Observed", "Predicted")
          ) %>%
          filter(!is.na(probability)) %>%
          ggplot(aes(x = time, y = probability, color = model, linetype = data_type)) +
          geom_line(size = 1) +
          facet_wrap(~outcome) +
          scale_color_viridis_d(name = "Model") +
          scale_linetype_manual(values = c("Observed" = "dashed", "Predicted" = "solid"),
                                name = "Data Type") +
          labs(
            title = "Cumulative Outcome Probabilities",
            subtitle = "Predicted vs Observed Over Time",
            x = "Days Since Entry",
            y = "Probability"
          ) +
          theme_minimal()
      }, error = function(e) {
        warning("Could not create cumulative outcomes plot:", e$message)
        NULL
      })
      
      report$plots$cumulative_outcomes <- cumulative_plot
    }
    
    # Create patchwork if requested
    if (plot_patchwork) {
      cat("12. Creating patchwork plots...\n")
      
      # Check if patchwork is available
      if (!requireNamespace("patchwork", quietly = TRUE)) {
        warning("patchwork package not available. Install with: install.packages('patchwork')")
        report$patchwork <- NULL
      } else {
        library(patchwork)
        
        # Create patchwork for each model
        model_patchworks <- map(model_names, function(name) {
          plots_for_model <- list()
          
          # Collect plots for this model
          if (!is.null(report$plots$calibration[[paste0(name, "_calibration")]])) {
            plots_for_model$calibration <- report$plots$calibration[[paste0(name, "_calibration")]]
          }
          
          if (!is.null(report$plots$state_prevalence[[paste0(name, "_prevalence")]])) {
            plots_for_model$prevalence <- report$plots$state_prevalence[[paste0(name, "_prevalence")]]
          }
          
          # Handle transition time plots (might be a list with hospital_time and calendar_time)
          trans_time_plot <- report$plots$transition_time[[paste0(name, "_transition_time")]]
          if (is.list(trans_time_plot) && "hospital_time" %in% names(trans_time_plot)) {
            plots_for_model$transition_time <- trans_time_plot$hospital_time
          } else if (!is.null(trans_time_plot)) {
            plots_for_model$transition_time <- trans_time_plot
          }
          
          if (!is.null(report$plots$patient_outliers[[paste0(name, "_outliers")]])) {
            plots_for_model$outliers <- report$plots$patient_outliers[[paste0(name, "_outliers")]]
          }
          
          # Create patchwork if we have plots
          if (length(plots_for_model) >= 2) {
            tryCatch({
              # Arrange plots in 2x2 grid
              if (length(plots_for_model) == 4) {
                patchwork_plot <- (plots_for_model[[1]] + plots_for_model[[2]]) / 
                  (plots_for_model[[3]] + plots_for_model[[4]])
              } else if (length(plots_for_model) == 3) {
                patchwork_plot <- plots_for_model[[1]] / 
                  (plots_for_model[[2]] + plots_for_model[[3]])
              } else {
                patchwork_plot <- plots_for_model[[1]] + plots_for_model[[2]]
              }
              
              # Add title
              patchwork_plot <- patchwork_plot + 
                plot_annotation(title = paste("Model Assessment:", name),
                                theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
              
              return(patchwork_plot)
            }, error = function(e) {
              warning(paste("Could not create patchwork for", name, ":", e$message))
              return(NULL)
            })
          } else {
            return(NULL)
          }
        })
        names(model_patchworks) <- model_names
        
        # Overall summary patchwork
        summary_plots <- list()
        if (!is.null(report$plots$cumulative_outcomes)) {
          summary_plots$cumulative <- report$plots$cumulative_outcomes
        }
        
        # Create overall comparison patchwork if multiple models
        if (length(model_names) > 1 && length(summary_plots) > 0) {
          overall_patchwork <- tryCatch({
            summary_plots$cumulative + 
              plot_annotation(title = "Model Comparison: Cumulative Outcomes",
                              theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
          }, error = function(e) {
            warning("Could not create overall patchwork:", e$message)
            NULL
          })
          
          report$patchwork <- list(
            individual_models = model_patchworks,
            overall_comparison = overall_patchwork
          )
        } else {
          report$patchwork <- list(individual_models = model_patchworks)
        }
      }
    }
  }
  
  cat("Model assessment report completed!\n")
  
  # Add summary statistics
  report$summary <- list(
    n_models_assessed = length(model_names),
    n_patients = length(unique(patient_data$deid_enc_id)),
    n_observations = nrow(patient_data),
    assessment_date = Sys.Date()
  )
  
  # Print summary table if requested
  if (print_summary && !is.null(report$summary_table)) {
    cat("\n")
    cat("="*80, "\n")
    cat("MODEL ASSESSMENT SUMMARY\n")
    cat("="*80, "\n")
    cat("Assessment Date:", as.character(report$summary$assessment_date), "\n")
    cat("Number of Models:", report$summary$n_models_assessed, "\n")
    cat("Number of Patients:", report$summary$n_patients, "\n")
    cat("Number of Observations:", report$summary$n_observations, "\n\n")
    
    # Format summary table for printing
    print_table <- report$summary_table %>%
      mutate(
        # Format numeric columns
        AIC = round(AIC, 1),
        BIC = round(BIC, 1), 
        logLik = round(logLik, 1),
        total_los_mae_cv = round(total_los_mae_cv, 2),
        days_severe_mae_cv = round(days_severe_mae_cv, 2),
        death_auc_cv = round(death_auc_cv, 3),
        severe_auc_cv = round(severe_auc_cv, 3),
        
        # Shorten formulas for readability
        formula = ifelse(nchar(formula) > 50, 
                         paste0(substr(formula, 1, 47), "..."), 
                         formula),
        covariates = ifelse(nchar(covariates) > 30,
                            paste0(substr(covariates, 1, 27), "..."),
                            covariates)
      ) %>%
      rename(
        Model = model,
        Formula = formula,
        Covariates = covariates,
        Converged = converged,
        `N Params` = n_params,
        `N Subjects` = n_subjects,
        `LOS MAE` = total_los_mae_cv,
        `Severe MAE` = days_severe_mae_cv,
        `Death AUC` = death_auc_cv,
        `Severe AUC` = severe_auc_cv
      )
    
    cat("Model Performance Summary:\n")
    cat("-"*80, "\n")
    print(print_table)
    cat("\n")
    
    # Print model ranking
    cat("Model Rankings:\n")
    cat("-"*40, "\n")
    cat("Best AIC:", print_table$Model[which.min(print_table$AIC)], 
        "(AIC =", min(print_table$AIC, na.rm = TRUE), ")\n")
    cat("Best BIC:", print_table$Model[which.min(print_table$BIC)], 
        "(BIC =", min(print_table$BIC, na.rm = TRUE), ")\n")
    cat("Best Death AUC:", print_table$Model[which.max(print_table$`Death AUC`)], 
        "(AUC =", max(print_table$`Death AUC`, na.rm = TRUE), ")\n")
    cat("Best Severe AUC:", print_table$Model[which.max(print_table$`Severe AUC`)], 
        "(AUC =", max(print_table$`Severe AUC`, na.rm = TRUE), ")\n")
    cat("\n")
  }
  
  # Display patchwork plots if created
  if (create_plots && plot_patchwork && !is.null(report$patchwork)) {
    cat("Displaying patchwork plots...\n")
    
    # Print individual model plots
    for (name in model_names) {
      if (!is.null(report$patchwork$individual_models[[name]])) {
        print(report$patchwork$individual_models[[name]])
        cat("\nPress [Enter] to continue to next model plot...")
        readline()
      }
    }
    
    # Print overall comparison plot
    if (!is.null(report$patchwork$overall_comparison)) {
      print(report$patchwork$overall_comparison)
    }
  }
  
  return(report)
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

## State-level functions -------------------------------------------------------

### Sojourn times -------------------------------------------------------

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

### State prevalence residuals ------------------------------------------------

#' Calculate state prevalence residuals using prevalence.msm
#' @param fitted_model A fitted msm model
#' @param patient_data Patient data used to fit the model
#' @param time_points Vector of time points to evaluate
#' @return Data frame with prevalence residuals
calculate_state_prevalence_residuals <- function(fitted_model, patient_data, 
                                                 time_points = seq(1, 30, by = 1)) {
  
  if (is.null(fitted_model)) {
    warning("Model is NULL")
    return(NULL)
  }
  
  tryCatch({
    # Use prevalence.msm to get observed and expected prevalences
    prevalence_result <- prevalence.msm(fitted_model, times = time_points, ci = "normal")
    
    # Extract observed prevalences
    observed_prev <- as.data.frame(prevalence_result$"Observed percentages") %>%
      rownames_to_column("time_point") %>%
      mutate(time_point = as.numeric(time_point)) %>%
      pivot_longer(cols = -time_point, names_to = "state", values_to = "observed_prev") %>%
      mutate(observed_prev = observed_prev / 100)  # Convert from percentage
    
    # Extract expected prevalences
    expected_prev <- as.data.frame(prevalence_result$"Expected percentages") %>%
      rownames_to_column("time_point") %>%
      mutate(time_point = as.numeric(time_point)) %>%
      pivot_longer(cols = -time_point, names_to = "state", values_to = "expected_prev") %>%
      mutate(expected_prev = expected_prev / 100)  # Convert from percentage
    
    # Get sample sizes for each time point
    observed_n <- as.data.frame(prevalence_result$"Observed numbers") %>%
      rownames_to_column("time_point") %>%
      mutate(time_point = as.numeric(time_point)) %>%
      pivot_longer(cols = -time_point, names_to = "state", values_to = "observed_n")
    
    # Calculate total at risk for each time point
    total_at_risk <- observed_n %>%
      group_by(time_point) %>%
      summarise(total_at_risk = sum(observed_n, na.rm = TRUE), .groups = "drop")
    
    # Combine and calculate residuals
    prevalence_residuals <- observed_prev %>%
      left_join(expected_prev, by = c("time_point", "state")) %>%
      left_join(observed_n, by = c("time_point", "state")) %>%
      left_join(total_at_risk, by = "time_point") %>%
      mutate(
        # Raw residual
        residual = observed_prev - expected_prev,
        
        # Standardized residual (binomial standard error)
        standardized_residual = residual / sqrt(expected_prev * (1 - expected_prev) / total_at_risk),
        
        # Deviance residual for binomial data
        deviance_residual = ifelse(
          observed_prev > 0,
          sign(residual) * sqrt(2 * observed_n * (
            observed_prev * log(observed_prev / expected_prev) +
              (1 - observed_prev) * log((1 - observed_prev) / (1 - expected_prev))
          )),
          -sqrt(2 * observed_n * log(1 / (1 - expected_prev)))
        ),
        
        # Residual categorization
        residual_type = case_when(
          abs(standardized_residual) > 2.5 ~ "Very Large",
          abs(standardized_residual) > 2 ~ "Large",
          abs(standardized_residual) > 1.5 ~ "Moderate", 
          TRUE ~ "Small"
        ),
        
        # Chi-square contribution
        chi_sq_contrib = ifelse(expected_prev > 0 & expected_prev < 1,
                                (observed_n - total_at_risk * expected_prev)^2 / 
                                  (total_at_risk * expected_prev * (1 - expected_prev)),
                                NA_real_)
      ) %>%
      filter(!is.na(observed_prev) & !is.na(expected_prev))
    
    return(prevalence_residuals)
    
  }, error = function(e) {
    warning(paste("Error calculating prevalence residuals:", e$message))
    return(NULL)
  })
}

#' Plot state prevalence residuals
#' @param prevalence_residuals Output from calculate_state_prevalence_residuals
#' @param facet_by_state Logical, whether to facet by state
#' @return ggplot object
plot_state_prevalence_residuals <- function(prevalence_residuals, facet_by_state = TRUE) {
  
  if (is.null(prevalence_residuals) || nrow(prevalence_residuals) == 0) {
    warning("No prevalence residuals to plot")
    return(NULL)
  }
  
  library(ggplot2)
  
  p <- prevalence_residuals %>%
    filter(!is.na(residual)) %>%
    ggplot(aes(x = time_point, y = residual, color = state)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    scale_color_viridis_d() +
    labs(
      title = "State Prevalence Residuals",
      subtitle = "Observed - Expected Prevalence by Time",
      x = "Days Since Entry",
      y = "Prevalence Residual",
      color = "State"
    ) +
    theme_minimal()
  
  if (facet_by_state) {
    p <- p + facet_wrap(~state, scales = "free_y")
  }
  
  return(p)
}

## Transition-level functions -------------------------------------------------------

### Transition intensities -------------------------------------------------------

#' Calculate transition intensities with confidence intervals
#' @param fitted_model A fitted msm model
#' @param covariate_values Optional list of covariate values
#' @return Data frame with transition intensities
calculate_transition_intensities <- function(fitted_model, covariate_values = NULL) {
  
  if (is.null(fitted_model)) {
    warning("Model is NULL")
    return(NULL)
  }
  
  tryCatch({
    # Get Q-matrix (transition intensity matrix)
    if (!is.null(covariate_values)) {
      qmat <- qmatrix.msm(fitted_model, covariates = covariate_values, ci = "normal")
    } else {
      qmat <- qmatrix.msm(fitted_model, ci = "normal")
    }
    
    # Convert to long format
    intensities_df <- as.data.frame(qmat$estimates) %>%
      rownames_to_column("from_state") %>%
      pivot_longer(cols = -from_state, names_to = "to_state", values_to = "intensity") %>%
      left_join(
        as.data.frame(qmat$L) %>%
          rownames_to_column("from_state") %>%
          pivot_longer(cols = -from_state, names_to = "to_state", values_to = "intensity_L"),
        by = c("from_state", "to_state")
      ) %>%
      left_join(
        as.data.frame(qmat$U) %>%
          rownames_to_column("from_state") %>%
          pivot_longer(cols = -from_state, names_to = "to_state", values_to = "intensity_U"),
        by = c("from_state", "to_state")
      ) %>%
      mutate(
        is_diagonal = from_state == to_state,
        is_allowed_transition = intensity != 0 | intensity_L != 0 | intensity_U != 0,
        covariate_label = if(!is.null(covariate_values)) {
          paste(names(covariate_values), covariate_values, sep = "=", collapse = ", ")
        } else "Baseline"
      ) %>%
      filter(is_allowed_transition | !is_diagonal) %>%  # Remove impossible transitions but keep diagonal
      arrange(from_state, to_state)
    
    return(intensities_df)
    
  }, error = function(e) {
    warning(paste("Error calculating transition intensities:", e$message))
    return(NULL)
  })
}


### Tidy transition probability matrices -----------------------------------

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

#' Extract comprehensive hazard ratios from MSM models
#' @param fitted_model A fitted msm model
#' @param covariate_names Vector of covariate names in the model
#' @param hazard_scale Scale for hazard ratios (default: 1 unit change)
#' @return Data frame with hazard ratios for each transition and covariate
extract_hazard_ratios <- function(fitted_model, covariate_names = NULL, hazard_scale = 1) {
  
  if (is.null(fitted_model)) {
    warning("Model is NULL")
    return(NULL)
  }
  
  if (is.null(fitted_model$covariates)) {
    warning("Model has no covariates")
    return(NULL)
  }
  
  # Get covariate names from model if not provided
  if (is.null(covariate_names)) {
    covariate_names <- names(fitted_model$covariates)
  }
  
  hr_results <- list()
  
  for (cov_name in covariate_names) {
    tryCatch({
      # Calculate hazard ratios
      hr_result <- hazard.msm(fitted_model, hazard.scale = hazard_scale)
      
      # Extract results
      if (is.list(hr_result) && !is.null(hr_result[[cov_name]])) {
        hr_data <- hr_result[[cov_name]]
        
        # Convert to data frame
        hr_df <- as.data.frame(hr_data) %>%
          rownames_to_column("transition") %>%
          mutate(
            covariate = cov_name,
            hazard_scale = hazard_scale
          ) %>%
          rename(
            HR = 2,  # Usually the estimate column
            L = 3,   # Lower CI
            U = 4    # Upper CI
          )
        
        hr_results[[cov_name]] <- hr_df
        
      } else {
        # Alternative approach using coefficient estimates
        coef_data <- fitted_model$estimates
        se_data <- sqrt(diag(fitted_model$covmat))
        
        # Find coefficients for this covariate
        cov_indices <- grep(cov_name, names(coef_data))
        
        if (length(cov_indices) > 0) {
          hr_df <- data.frame(
            transition = names(coef_data)[cov_indices],
            covariate = cov_name,
            log_HR = coef_data[cov_indices] * hazard_scale,
            log_HR_se = se_data[cov_indices] * hazard_scale,
            hazard_scale = hazard_scale
          ) %>%
            mutate(
              HR = exp(log_HR),
              L = exp(log_HR - 1.96 * log_HR_se),
              U = exp(log_HR + 1.96 * log_HR_se),
              z_stat = log_HR / log_HR_se,
              p_value = 2 * (1 - pnorm(abs(z_stat)))
            )
          
          hr_results[[cov_name]] <- hr_df
        }
      }
    }, error = function(e) {
      warning(paste("Error calculating hazard ratios for", cov_name, ":", e$message))
      hr_results[[cov_name]] <- NULL
    })
  }
  
  # Combine results
  if (length(hr_results) > 0) {
    combined_results <- bind_rows(hr_results)
    return(combined_results)
  } else {
    return(NULL)
  }
}

# Predictive performance functions ----------------------------------------

#' Calculate comprehensive predictive outcomes using cross-validation
#' @param fitted_models List of fitted models
#' @param patient_data Patient data
#' @param crude_rates List of crude rates for model initialization
#' @param k_folds Number of CV folds
#' @param balanced_folds Logical, whether to stratify folds by outcome
#' @param prediction_times Vector of prediction time points
#' @param n_sim Number of simulations per patient
#' @return Data frame with detailed predictive metrics
calculate_predictive_performance <- function(fitted_models, patient_data, crude_rates,
                                             k_folds = 5, balanced_folds = TRUE,
                                             prediction_times = c(14, 30), n_sim = 100) {
  
  # Create fold assignments
  unique_patients <- patient_data %>%
    group_by(deid_enc_id) %>%
    summarise(
      final_state = last(state),
      max_day = max(DaysSinceEntry),
      model = first(model),
      .groups = "drop"
    )
  
  set.seed(123)
  
  if (balanced_folds) {
    # Stratify by final outcome for balanced folds
    fold_assignments <- unique_patients %>%
      group_by(final_state) %>%
      mutate(fold = sample(rep(1:k_folds, length.out = n()))) %>%
      ungroup() %>%
      select(deid_enc_id, fold)
  } else {
    # Simple random assignment
    fold_assignments <- unique_patients %>%
      mutate(fold = sample(rep(1:k_folds, length.out = n()))) %>%
      select(deid_enc_id, fold)
  }
  
  patient_data <- patient_data %>%
    left_join(fold_assignments, by = "deid_enc_id")
  
  predictive_results <- list()
  
  for (model_name in names(fitted_models)) {
    cat("Calculating predictive performance for:", model_name, "\n")
    
    model_data <- patient_data %>% filter(model == model_name)
    if (nrow(model_data) == 0) next
    
    fold_results <- map_dfr(1:k_folds, function(fold) {
      cat("  Fold", fold, "of", k_folds, "\n")
      
      # Split data
      train_data <- model_data %>% filter(fold != !!fold)
      test_data <- model_data %>% filter(fold == !!fold)
      
      if (nrow(train_data) == 0 || nrow(test_data) == 0) {
        return(data.frame(
          fold = fold,
          total_los_mae = NA,
          days_severe_mae = NA,
          death_auc = NA,
          severe_auc = NA,
          calibration_slope_death = NA,
          calibration_intercept_death = NA
        ))
      }
      
      # Use full model for prediction (in practice, would refit on training data)
      mod <- fitted_models[[model_name]]
      if (is.null(mod)) return(NULL)
      
      # Calculate observed outcomes for test patients
      test_outcomes <- test_data %>%
        group_by(deid_enc_id) %>%
        summarise(
          total_los = max(DaysSinceEntry, na.rm = TRUE),
          days_severe = sum(state %in% c("S", "S1", "S2"), na.rm = TRUE),
          died = any(state == "D"),
          ever_severe = any(state %in% c("S", "S1", "S2")),
          recovered = any(state == "R"),
          .groups = "drop"
        )
      
      # Simulate outcomes for test patients
      test_patients <- unique(test_data$deid_enc_id)
      simulated_outcomes <- map_dfr(test_patients, function(patient_id) {
        
        patient_covs <- test_data %>%
          filter(deid_enc_id == patient_id) %>%
          slice(1) %>%
          select(-deid_enc_id, -state, -state_num, -DaysSinceEntry, -date, 
                 -model, -fold) %>%
          as.list()
        
        # Simple simulation approach - get probability matrix
        tryCatch({
          # Simulate for longest prediction time
          max_time <- max(prediction_times)
          pmat <- pmatrix.msm(mod, t = max_time, 
                              covariates = if(length(patient_covs) > 0) patient_covs else NULL)
          
          # Starting from moderate state (state 1 typically)
          start_probs <- pmat[1, ]
          
          data.frame(
            deid_enc_id = patient_id,
            pred_prob_death = start_probs["D"] %||% 0,
            pred_prob_recovery = start_probs["R"] %||% 0,
            pred_prob_severe = sum(start_probs[grepl("S", names(start_probs))]) %||% 0,
            # Simple LOS prediction based on sojourn times
            pred_total_los = sum(sojourn.msm(mod, 
                                             covariates = if(length(patient_covs) > 0) patient_covs else NULL)$estimates),
            pred_days_severe = (sum(start_probs[grepl("S", names(start_probs))]) %||% 0) * 
              (sojourn.msm(mod, covariates = if(length(patient_covs) > 0) patient_covs else NULL)$estimates[grepl("S", rownames(sojourn.msm(mod)$estimates))] %||% c(0))[1]
          )
        }, error = function(e) {
          data.frame(
            deid_enc_id = patient_id,
            pred_prob_death = NA,
            pred_prob_recovery = NA,
            pred_prob_severe = NA,
            pred_total_los = NA,
            pred_days_severe = NA
          )
        })
      })
      
      # Combine observed and predicted
      fold_data <- test_outcomes %>%
        left_join(simulated_outcomes, by = "deid_enc_id")
      
      # Calculate performance metrics
      tryCatch({
        # MAE for continuous outcomes
        total_los_mae <- mean(abs(fold_data$total_los - fold_data$pred_total_los), na.rm = TRUE)
        days_severe_mae <- mean(abs(fold_data$days_severe - fold_data$pred_days_severe), na.rm = TRUE)
        
        # AUC for binary outcomes
        death_auc <- if(var(fold_data$died, na.rm = TRUE) > 0 && !all(is.na(fold_data$pred_prob_death))) {
          tryCatch({
            pROC::auc(fold_data$died, fold_data$pred_prob_death, quiet = TRUE)
          }, error = function(e) NA)
        } else NA
        
        severe_auc <- if(var(fold_data$ever_severe, na.rm = TRUE) > 0 && !all(is.na(fold_data$pred_prob_severe))) {
          tryCatch({
            pROC::auc(fold_data$ever_severe, fold_data$pred_prob_severe, quiet = TRUE)
          }, error = function(e) NA)
        } else NA
        
        # Calibration metrics (slope and intercept from logistic regression)
        calibration_death <- if(!all(is.na(fold_data$pred_prob_death)) && var(fold_data$died, na.rm = TRUE) > 0) {
          tryCatch({
            cal_model <- glm(died ~ I(qlogis(pmax(pmin(pred_prob_death, 0.999), 0.001))), 
                             data = fold_data, family = binomial())
            list(slope = coef(cal_model)[2], intercept = coef(cal_model)[1])
          }, error = function(e) list(slope = NA, intercept = NA))
        } else list(slope = NA, intercept = NA)
        
        data.frame(
          fold = fold,
          total_los_mae = total_los_mae,
          days_severe_mae = days_severe_mae,
          death_auc = as.numeric(death_auc),
          severe_auc = as.numeric(severe_auc),
          calibration_slope_death = calibration_death$slope,
          calibration_intercept_death = calibration_death$intercept
        )
      }, error = function(e) {
        data.frame(
          fold = fold,
          total_los_mae = NA,
          days_severe_mae = NA,
          death_auc = NA,
          severe_auc = NA,
          calibration_slope_death = NA,
          calibration_intercept_death = NA
        )
      })
    })
    
    # Aggregate across folds
    predictive_results[[model_name]] <- fold_results %>%
      summarise(
        model = model_name,
        total_los_mae_cv = mean(total_los_mae, na.rm = TRUE),
        total_los_mae_se = sd(total_los_mae, na.rm = TRUE) / sqrt(n()),
        days_severe_mae_cv = mean(days_severe_mae, na.rm = TRUE),
        days_severe_mae_se = sd(days_severe_mae, na.rm = TRUE) / sqrt(n()),
        death_auc_cv = mean(death_auc, na.rm = TRUE),
        death_auc_se = sd(death_auc, na.rm = TRUE) / sqrt(n()),
        severe_auc_cv = mean(severe_auc, na.rm = TRUE),
        severe_auc_se = sd(severe_auc, na.rm = TRUE) / sqrt(n()),
        calibration_slope_mean = mean(calibration_slope_death, na.rm = TRUE),
        calibration_slope_se = sd(calibration_slope_death, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
  }
  
  bind_rows(predictive_results)
}

# Model diagnostic functions --------------------------------------------

#' Calculate comprehensive transition residuals
#' @param fitted_model A fitted msm model
#' @param patient_data Patient data used to fit the model
#' @param residual_type Type of residuals ("deviance", "pearson", "raw")
#' @return Data frame with transition residuals
calculate_transition_residuals <- function(fitted_model, patient_data, 
                                           residual_type = "deviance") {
  
  if (is.null(fitted_model)) {
    warning("Model is NULL")
    return(NULL)
  }
  
  # Get observed transitions
  observed_transitions <- patient_data %>%
    arrange(deid_enc_id, DaysSinceEntry) %>%
    group_by(deid_enc_id) %>%
    mutate(
      from_state = state,
      to_state = lead(state),
      time_diff = lead(DaysSinceEntry) - DaysSinceEntry,
      transition_time = DaysSinceEntry,
      calendar_time = CalendarTime
    ) %>%
    filter(!is.na(to_state), time_diff > 0) %>%
    ungroup()
  
  # Calculate expected transitions using model
  transition_residuals <- observed_transitions %>%
    rowwise() %>%
    mutate(
      # Get patient covariates if model has covariates
      expected_prob = {
        if (!is.null(fitted_model$covariates)) {
          # Extract covariates for this patient/time
          patient_covs <- cur_data() %>%
            select(any_of(names(fitted_model$covariates))) %>%
            as.list()
          
          tryCatch({
            pmat <- pmatrix.msm(fitted_model, t = time_diff, covariates = patient_covs)
            if (from_state %in% rownames(pmat) && to_state %in% colnames(pmat)) {
              pmat[from_state, to_state]
            } else {
              NA_real_
            }
          }, error = function(e) NA_real_)
        } else {
          tryCatch({
            pmat <- pmatrix.msm(fitted_model, t = time_diff)
            if (from_state %in% rownames(pmat) && to_state %in% colnames(pmat)) {
              pmat[from_state, to_state]
            } else {
              NA_real_
            }
          }, error = function(e) NA_real_)
        }
      }
    ) %>%
    ungroup() %>%
    mutate(
      # Calculate residuals
      raw_residual = 1 - expected_prob,  # Observed is always 1 for actual transitions
      
      # Pearson residual
      pearson_residual = raw_residual / sqrt(expected_prob * (1 - expected_prob)),
      
      # Deviance residual  
      deviance_residual = sign(raw_residual) * sqrt(-2 * log(expected_prob)),
      
      # Transition identifier
      transition = paste(from_state, "->", to_state)
    ) %>%
    filter(!is.na(expected_prob))
  
  # Calculate standardized versions and categorization based on selected residual type
  if (residual_type == "deviance") {
    transition_residuals <- transition_residuals %>%
      mutate(
        residual = deviance_residual,
        residual_std = deviance_residual / sd(deviance_residual, na.rm = TRUE),
        residual_magnitude = abs(residual_std),
        residual_category = case_when(
          residual_magnitude > 2.5 ~ "Very Large",
          residual_magnitude > 2 ~ "Large", 
          residual_magnitude > 1.5 ~ "Moderate",
          TRUE ~ "Small"
        )
      )
  } else if (residual_type == "pearson") {
    transition_residuals <- transition_residuals %>%
      mutate(
        residual = pearson_residual,
        residual_std = pearson_residual / sd(pearson_residual, na.rm = TRUE),
        residual_magnitude = abs(residual_std),
        residual_category = case_when(
          residual_magnitude > 2.5 ~ "Very Large",
          residual_magnitude > 2 ~ "Large", 
          residual_magnitude > 1.5 ~ "Moderate",
          TRUE ~ "Small"
        )
      )
  } else {  # raw residuals
    transition_residuals <- transition_residuals %>%
      mutate(
        residual = raw_residual,
        residual_std = raw_residual / sd(raw_residual, na.rm = TRUE),
        residual_magnitude = abs(residual_std),
        residual_category = case_when(
          residual_magnitude > 2.5 ~ "Very Large",
          residual_magnitude > 2 ~ "Large", 
          residual_magnitude > 1.5 ~ "Moderate",
          TRUE ~ "Small"
        )
      )
  }
  
  return(transition_residuals)
}

#' Plot transition residuals over time and covariates
#' @param transition_residuals Output from calculate_transition_residuals
#' @param plot_type Type of plot ("time", "covariate", "qq")
#' @param covariate_name Name of covariate to plot against (for plot_type = "covariate")
#' @return ggplot object
plot_transition_residuals <- function(transition_residuals, plot_type = "time", 
                                      covariate_name = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  if (is.null(transition_residuals) || nrow(transition_residuals) == 0) {
    warning("No transition residuals to plot")
    return(NULL)
  }
  
  if (plot_type == "time") {
    # Residuals vs time
    p <- transition_residuals %>%
      ggplot(aes(x = transition_time, y = residual, color = transition)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
      scale_color_viridis_d(name = "Transition") +
      labs(
        title = "Transition Residuals vs Hospital Time",
        x = "Days Since Entry", 
        y = "Deviance Residual"
      ) +
      theme_minimal()
    
    # Add calendar time if available
    if ("calendar_time" %in% names(transition_residuals)) {
      p2 <- transition_residuals %>%
        ggplot(aes(x = calendar_time, y = residual, color = transition)) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
        scale_color_viridis_d(name = "Transition") +
        labs(
          title = "Transition Residuals vs Calendar Time",
          x = "Calendar Time", 
          y = "Deviance Residual"
        ) +
        theme_minimal()
      
      return(list(hospital_time = p, calendar_time = p2))
    }
    
  } else if (plot_type == "covariate" && !is.null(covariate_name)) {
    
    if (!covariate_name %in% names(transition_residuals)) {
      warning(paste("Covariate", covariate_name, "not found in residuals data"))
      return(NULL)
    }
    
    # Residuals vs covariate
    p <- transition_residuals %>%
      ggplot(aes_string(x = covariate_name, y = "residual", color = "transition")) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
      scale_color_viridis_d(name = "Transition") +
      labs(
        title = paste("Transition Residuals vs", covariate_name),
        x = covariate_name,
        y = "Deviance Residual"
      ) +
      theme_minimal()
    
  } else if (plot_type == "qq") {
    # Q-Q plot for residual normality
    p <- transition_residuals %>%
      ggplot(aes(sample = residual)) +
      geom_qq() +
      geom_qq_line() +
      facet_wrap(~transition, scales = "free") +
      labs(
        title = "Q-Q Plot of Transition Residuals",
        subtitle = "Checking normality assumption"
      ) +
      theme_minimal()
    
  } else {
    warning("Invalid plot_type. Choose 'time', 'covariate', or 'qq'")
    return(NULL)
  }
  
  return(p)
}


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

## Patient-Level Outlier Detection Functions -------------------------------

#' Identify outlier patients with large residuals
#' @param fitted_model A fitted msm model
#' @param patient_data Patient data used to fit the model  
#' @param threshold Threshold for defining outliers (standard deviations)
#' @return Data frame with patient outlier metrics
identify_outlier_patients <- function(fitted_model, patient_data, threshold = 2.5) {
  
  if (is.null(fitted_model)) {
    warning("Model is NULL")
    return(NULL)
  }
  
  # Calculate patient-level residuals
  patient_residuals <- calculate_patient_residuals(fitted_model, patient_data)
  
  if (is.null(patient_residuals)) {
    warning("Could not calculate patient residuals")
    return(NULL)
  }
  
  # Identify outliers based on various metrics
  outliers <- patient_residuals %>%
    mutate(
      # Overall outlier flags
      is_los_outlier = abs(los_residual_std) > threshold,
      is_transition_outlier = abs(mean_transition_residual) > threshold,
      is_prevalence_outlier = abs(mean_prevalence_residual) > threshold,
      
      # Combined outlier score
      outlier_score = sqrt(los_residual_std^2 + mean_transition_residual^2 + mean_prevalence_residual^2),
      is_outlier = outlier_score > threshold,
      
      # Outlier type
      outlier_type = case_when(
        is_los_outlier & is_transition_outlier & is_prevalence_outlier ~ "Multiple",
        is_los_outlier ~ "Length of Stay",
        is_transition_outlier ~ "Transition Pattern",
        is_prevalence_outlier ~ "State Occupancy",
        TRUE ~ "None"
      ),
      
      # Severity
      outlier_severity = case_when(
        outlier_score > threshold * 1.5 ~ "Severe",
        outlier_score > threshold ~ "Moderate", 
        outlier_score > threshold * 0.75 ~ "Mild",
        TRUE ~ "None"
      )
    ) %>%
    arrange(desc(outlier_score))
  
  return(outliers)
}

#' Calculate patient-level residuals and fit metrics
#' @param fitted_model A fitted msm model
#' @param patient_data Patient data used to fit the model
#' @return Data frame with patient-level residuals
calculate_patient_residuals <- function(fitted_model, patient_data) {
  
  if (is.null(fitted_model)) {
    warning("Model is NULL") 
    return(NULL)
  }
  
  # Get transition residuals
  transition_residuals <- calculate_transition_residuals(fitted_model, patient_data, "deviance")
  
  if (is.null(transition_residuals)) {
    warning("Could not calculate transition residuals")
    return(NULL)
  }
  
  # Calculate patient-level summaries
  patient_summaries <- patient_data %>%
    group_by(deid_enc_id) %>%
    summarise(
      observed_los = max(DaysSinceEntry, na.rm = TRUE),
      n_observations = n(),
      n_states = n_distinct(state),
      final_state = last(state),
      ever_severe = any(state %in% c("S", "S1", "S2")),
      days_severe = sum(state %in% c("S", "S1", "S2")),
      n_transitions = n() - 1,
      .groups = "drop"
    )
  
  # Calculate expected outcomes for each patient
  expected_outcomes <- map_dfr(unique(patient_data$deid_enc_id), function(pid) {
    
    patient_info <- patient_data %>%
      filter(deid_enc_id == pid) %>%
      slice(1)
    
    # Get covariates if model has them
    if (!is.null(fitted_model$covariates)) {
      patient_covs <- patient_info %>%
        select(any_of(names(fitted_model$covariates))) %>%
        as.list()
    } else {
      patient_covs <- NULL
    }
    
    tryCatch({
      # Predict LOS using sojourn times
      sojourn_times <- sojourn.msm(fitted_model, covariates = patient_covs)
      expected_los <- sum(sojourn_times$estimates)
      
      # Predict final state probabilities
      max_time <- patient_summaries$observed_los[patient_summaries$deid_enc_id == pid]
      pmat <- pmatrix.msm(fitted_model, t = max_time, covariates = patient_covs)
      
      # Assuming starting from first transient state
      start_state <- rownames(pmat)[1]
      final_probs <- pmat[start_state, ]
      
      data.frame(
        deid_enc_id = pid,
        expected_los = expected_los,
        expected_prob_death = final_probs["D"] %||% 0,
        expected_prob_recovery = final_probs["R"] %||% 0,
        expected_prob_severe = sum(final_probs[grepl("S", names(final_probs))]) %||% 0
      )
      
    }, error = function(e) {
      data.frame(
        deid_enc_id = pid,
        expected_los = NA,
        expected_prob_death = NA,
        expected_prob_recovery = NA,
        expected_prob_severe = NA
      )
    })
  })
  
  # Aggregate transition residuals by patient
  patient_transition_residuals <- transition_residuals %>%
    group_by(deid_enc_id) %>%
    summarise(
      mean_transition_residual = mean(abs(residual), na.rm = TRUE),
      max_transition_residual = max(abs(residual), na.rm = TRUE),
      n_large_residuals = sum(abs(residual) > 2, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate prevalence residuals by patient (simplified)
  patient_prevalence_residuals <- patient_data %>%
    group_by(deid_enc_id) %>%
    summarise(
      # This is a simplified approach - ideally would compare patient trajectory to expected
      mean_prevalence_residual = 0,  # Placeholder - would need more complex calculation
      .groups = "drop"
    )
  
  # Combine all patient-level metrics
  patient_residuals <- patient_summaries %>%
    left_join(expected_outcomes, by = "deid_enc_id") %>%
    left_join(patient_transition_residuals, by = "deid_enc_id") %>%
    left_join(patient_prevalence_residuals, by = "deid_enc_id") %>%
    mutate(
      # LOS residuals
      los_residual = observed_los - expected_los,
      los_residual_std = scale(los_residual)[, 1],
      
      # Death prediction residuals
      death_residual = as.numeric(final_state == "D") - expected_prob_death,
      
      # Recovery prediction residuals  
      recovery_residual = as.numeric(final_state == "R") - expected_prob_recovery,
      
      # Severe state prediction residuals
      severe_residual = as.numeric(ever_severe) - expected_prob_severe,
      
      # Replace NAs with 0 for residuals
      mean_transition_residual = replace_na(mean_transition_residual, 0),
      mean_prevalence_residual = replace_na(mean_prevalence_residual, 0)
    )
  
  return(patient_residuals)
}

#' Plot patient-level residuals and identify outliers
#' @param patient_residuals Output from calculate_patient_residuals or identify_outlier_patients
#' @param plot_type Type of plot ("outlier_score", "los_residual", "by_characteristic")
#' @param characteristic Name of patient characteristic to plot by
#' @return ggplot object
plot_patient_residuals <- function(patient_residuals, plot_type = "outlier_score", 
                                   characteristic = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  if (is.null(patient_residuals) || nrow(patient_residuals) == 0) {
    warning("No patient residuals to plot")
    return(NULL)
  }
  
  if (plot_type == "outlier_score") {
    # Plot overall outlier scores
    p <- patient_residuals %>%
      ggplot(aes(x = outlier_score)) +
      geom_histogram(bins = 50, alpha = 0.7, fill = "skyblue", color = "black") +
      geom_vline(xintercept = 2.5, linetype = "dashed", color = "red", 
                 alpha = 0.7, size = 1) +
      annotate("text", x = 2.5, y = Inf, label = "Outlier Threshold", 
               hjust = -0.1, vjust = 1.5, color = "red") +
      labs(
        title = "Distribution of Patient Outlier Scores",
        x = "Outlier Score",
        y = "Number of Patients"
      ) +
      theme_minimal()
    
  } else if (plot_type == "los_residual") {
    # Length of stay residuals
    p <- patient_residuals %>%
      ggplot(aes(x = observed_los, y = los_residual)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_point(aes(color = final_state), alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
      scale_color_viridis_d(name = "Final State") +
      labs(
        title = "Length of Stay Residuals",
        x = "Observed LOS (days)",
        y = "LOS Residual (Observed - Expected)"
      ) +
      theme_minimal()
    
  } else if (plot_type == "by_characteristic" && !is.null(characteristic)) {
    
    if (!characteristic %in% names(patient_residuals)) {
      warning(paste("Characteristic", characteristic, "not found in patient data"))
      return(NULL)
    }
    
    # Outlier scores by patient characteristic
    p <- patient_residuals %>%
      ggplot(aes_string(x = characteristic, y = "outlier_score")) +
      geom_boxplot(alpha = 0.7) +
      geom_hline(yintercept = 2.5, linetype = "dashed", color = "red", alpha = 0.7) +
      labs(
        title = paste("Outlier Scores by", characteristic),
        x = characteristic,
        y = "Outlier Score"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else {
    warning("Invalid plot_type or missing characteristic")
    return(NULL)
  }
  
  return(p)
}


## Calibration functions -----------------------------------------------------

#' Calculate calibration plots for specific outcomes and time points
#' @param fitted_models List of fitted models
#' @param patient_data Patient data
#' @param time_points Vector of time points (days) for calibration
#' @param outcomes Vector of outcomes to calibrate ("death", "recovery")
#' @param n_bins Number of risk bins for calibration
#' @return Data frame with calibration metrics
calculate_calibration_plots <- function(fitted_models, patient_data, 
                                        time_points = c(5, 10, 15, 20, 30),
                                        outcomes = c("death", "recovery"),
                                        n_bins = 10) {
  
  calibration_results <- list()
  
  for (model_name in names(fitted_models)) {
    cat("Calculating calibration for model:", model_name, "\n")
    
    model <- fitted_models[[model_name]]
    if (is.null(model)) next
    
    model_data <- patient_data %>% filter(model == model_name)
    if (nrow(model_data) == 0) next
    
    for (time_point in time_points) {
      for (outcome in outcomes) {
        
        # Get patients who have follow-up to this time point
        eligible_patients <- model_data %>%
          group_by(deid_enc_id) %>%
          filter(max(DaysSinceEntry) >= time_point) %>%
          slice(1) %>%
          ungroup()
        
        if (nrow(eligible_patients) == 0) next
        
        # Calculate observed outcomes at time point
        observed_outcomes <- map_dfr(unique(eligible_patients$deid_enc_id), function(pid) {
          patient_trajectory <- model_data %>%
            filter(deid_enc_id == pid, DaysSinceEntry <= time_point) %>%
            arrange(DaysSinceEntry)
          
          if (nrow(patient_trajectory) == 0) return(NULL)
          
          # Determine outcome by time_point
          outcome_occurred <- switch(outcome,
                                     "death" = any(patient_trajectory$state == "D"),
                                     "recovery" = any(patient_trajectory$state == "R"),
                                     FALSE
          )
          
          data.frame(
            deid_enc_id = pid,
            observed_outcome = as.numeric(outcome_occurred)
          )
        })
        
        # Calculate predicted probabilities
        predicted_probs <- map_dfr(unique(eligible_patients$deid_enc_id), function(pid) {
          patient_info <- eligible_patients %>% filter(deid_enc_id == pid)
          
          # Get covariates if model has them
          if (!is.null(model$covariates)) {
            patient_covs <- patient_info %>%
              select(any_of(names(model$covariates))) %>%
              as.list()
          } else {
            patient_covs <- NULL
          }
          
          tryCatch({
            pmat <- pmatrix.msm(model, t = time_point, covariates = patient_covs)
            
            # Starting from first transient state
            start_state <- rownames(pmat)[1]
            
            predicted_prob <- switch(outcome,
                                     "death" = pmat[start_state, "D"] %||% 0,
                                     "recovery" = pmat[start_state, "R"] %||% 0,
                                     0
            )
            
            data.frame(
              deid_enc_id = pid,
              predicted_prob = predicted_prob
            )
            
          }, error = function(e) {
            data.frame(
              deid_enc_id = pid,
              predicted_prob = NA
            )
          })
        })
        
        # Combine observed and predicted
        calibration_data <- observed_outcomes %>%
          inner_join(predicted_probs, by = "deid_enc_id") %>%
          filter(!is.na(predicted_prob))
        
        if (nrow(calibration_data) == 0) next
        
        # Create risk bins
        calibration_data <- calibration_data %>%
          mutate(
            risk_bin = cut(predicted_prob, 
                           breaks = quantile(predicted_prob, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = FALSE)
          )
        
        # Calculate calibration metrics by bin
        bin_metrics <- calibration_data %>%
          group_by(risk_bin) %>%
          summarise(
            n_patients = n(),
            observed_rate = mean(observed_outcome, na.rm = TRUE),
            predicted_rate = mean(predicted_prob, na.rm = TRUE),
            min_pred_prob = min(predicted_prob, na.rm = TRUE),
            max_pred_prob = max(predicted_prob, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(
            calibration_error = observed_rate - predicted_rate,
            model = model_name,
            time_point = time_point,
            outcome = outcome
          )
        
        # Overall calibration metrics
        overall_metrics <- calibration_data %>%
          summarise(
            model = model_name,
            time_point = time_point,
            outcome = outcome,
            n_patients = n(),
            
            # Calibration slope and intercept
            calibration_slope = {
              if (var(predicted_prob) > 0) {
                tryCatch({
                  cal_model <- glm(observed_outcome ~ I(qlogis(pmax(pmin(predicted_prob, 0.999), 0.001))), 
                                   family = binomial())
                  coef(cal_model)[2]
                }, error = function(e) NA)
              } else NA
            },
            
            calibration_intercept = {
              if (var(predicted_prob) > 0) {
                tryCatch({
                  cal_model <- glm(observed_outcome ~ I(qlogis(pmax(pmin(predicted_prob, 0.999), 0.001))), 
                                   family = binomial())
                  coef(cal_model)[1]
                }, error = function(e) NA)
              } else NA
            },
            
            # Hosmer-Lemeshow type metric
            hl_statistic = {
              if (length(unique(risk_bin)) > 1) {
                sum((bin_metrics$observed_rate - bin_metrics$predicted_rate)^2 * bin_metrics$n_patients, na.rm = TRUE)
              } else NA
            },
            
            # Mean absolute calibration error
            mace = mean(abs(bin_metrics$calibration_error), na.rm = TRUE),
            
            # Brier score
            brier_score = mean((observed_outcome - predicted_prob)^2, na.rm = TRUE),
            
            .groups = "drop"
          )
        
        calibration_results[[paste(model_name, time_point, outcome, sep = "_")]] <- list(
          overall = overall_metrics,
          by_bin = bin_metrics,
          raw_data = calibration_data
        )
      }
    }
  }
  
  return(calibration_results)
}

#' Calculate cumulative outcome probabilities
#' @param fitted_models List of fitted models
#' @param patient_data Patient data  
#' @param max_time Maximum follow-up time
#' @param time_step Time step for calculations
#' @return Data frame with cumulative probabilities over time
calculate_cumulative_outcomes <- function(fitted_models, patient_data, 
                                          max_time = 30, time_step = 1) {
  
  time_points <- seq(time_step, max_time, by = time_step)
  cumulative_results <- list()
  
  for (model_name in names(fitted_models)) {
    cat("Calculating cumulative outcomes for model:", model_name, "\n")
    
    model <- fitted_models[[model_name]]
    if (is.null(model)) next
    
    model_data <- patient_data %>% filter(model == model_name)
    if (nrow(model_data) == 0) next
    
    # Calculate for average patient (could extend to covariate-specific)
    cumulative_probs <- map_dfr(time_points, function(t) {
      
      tryCatch({
        pmat <- pmatrix.msm(model, t = t)
        
        # Starting from first transient state
        start_state <- rownames(pmat)[1]
        probs <- pmat[start_state, ]
        
        data.frame(
          model = model_name,
          time = t,
          prob_death = probs["D"] %||% 0,
          prob_recovery = probs["R"] %||% 0,
          prob_severe = sum(probs[grepl("S", names(probs))]) %||% 0,
          prob_moderate = sum(probs[grepl("M", names(probs))]) %||% 0
        )
        
      }, error = function(e) {
        data.frame(
          model = model_name,
          time = t,
          prob_death = NA,
          prob_recovery = NA,
          prob_severe = NA,
          prob_moderate = NA
        )
      })
    })
    
    # Calculate observed cumulative frequencies for comparison
    observed_cumulative <- map_dfr(time_points, function(t) {
      
      outcomes_by_time <- model_data %>%
        group_by(deid_enc_id) %>%
        filter(DaysSinceEntry <= t) %>%
        summarise(
          died = any(state == "D"),
          recovered = any(state == "R"),
          ever_severe = any(state %in% c("S", "S1", "S2")),
          .groups = "drop"
        )
      
      data.frame(
        model = model_name,
        time = t,
        observed_death_rate = mean(outcomes_by_time$died, na.rm = TRUE),
        observed_recovery_rate = mean(outcomes_by_time$recovered, na.rm = TRUE),
        observed_severe_rate = mean(outcomes_by_time$ever_severe, na.rm = TRUE),
        n_patients = nrow(outcomes_by_time)
      )
    })
    
    # Combine predicted and observed
    combined_cumulative <- cumulative_probs %>%
      left_join(observed_cumulative, by = c("model", "time"))
    
    cumulative_results[[model_name]] <- combined_cumulative
  }
  
  return(bind_rows(cumulative_results))
}

#' Plot calibration results
#' @param calibration_results Output from calculate_calibration_plots
#' @param model_name Specific model to plot (optional)
#' @param outcome Specific outcome to plot (optional)
#' @return ggplot object
plot_calibration <- function(calibration_results, model_name = NULL, outcome = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  if (length(calibration_results) == 0) {
    warning("No calibration results to plot")
    return(NULL)
  }
  
  # Extract bin-level data
  plot_data <- map_dfr(names(calibration_results), function(result_name) {
    result <- calibration_results[[result_name]]
    if (!is.null(result$by_bin)) {
      result$by_bin
    } else {
      NULL
    }
  })
  
  if (nrow(plot_data) == 0) {
    warning("No bin-level data found for plotting")
    return(NULL)
  }
  
  # Filter if specific model/outcome requested
  if (!is.null(model_name)) {
    plot_data <- plot_data %>% filter(model == model_name)
  }
  if (!is.null(outcome)) {
    plot_data <- plot_data %>% filter(outcome == outcome)
  }
  
  # Create calibration plot
  p <- plot_data %>%
    ggplot(aes(x = predicted_rate, y = observed_rate)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5, color = "red") +
    geom_point(aes(size = n_patients, color = factor(time_point)), alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.3, color = "blue") +
    scale_size_continuous(name = "N Patients", range = c(1, 4)) +
    scale_color_viridis_d(name = "Time Point") +
    labs(
      title = "Calibration Plot",
      subtitle = "Observed vs Predicted Event Rates",
      x = "Predicted Probability",
      y = "Observed Rate",
      caption = "Red dashed line = perfect calibration"
    ) +
    theme_minimal() +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1)
  
  # Add faceting if multiple models/outcomes
  if (length(unique(plot_data$model)) > 1 && length(unique(plot_data$outcome)) > 1) {
    p <- p + facet_grid(outcome ~ model)
  } else if (length(unique(plot_data$model)) > 1) {
    p <- p + facet_wrap(~model)
  } else if (length(unique(plot_data$outcome)) > 1) {
    p <- p + facet_wrap(~outcome)
  }
  
  return(p)
}



# Covariate functions -----------------------------------------------------


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

