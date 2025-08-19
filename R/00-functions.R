
date_fn <- function(x) as.Date(as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S"))

`%nin%` <- Negate(`%in%`)

tbl_fn <- function(x){
  tbl <- table(x)
  res <- cbind(tbl, round(prop.table(tbl)*100, 2))
  colnames(res) <- c('Count', 'Percentage')
  res
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

fit_msm_models <- function(patient_data, crude_rates, covariates = NULL) {
  fitted_models <- list()
  
  # Create formula string for covariates if provided
  covariate_formula <- if (!is.null(covariates)) {
    paste("~", paste(covariates, collapse = " + "))
  } else {
    NULL
  }
  
  for (modelname in unique(patient_data$model)) {
    model_data <- patient_data[which(patient_data$model == modelname), ]
    crude_result <- crude_rates[[modelname]]
    
    if (is.null(crude_result)) {
      warning(paste("Crude rates missing for", modelname, "- skipping model fitting"))
      next
    }
    
    fitted_model <- tryCatch({
      if (!is.null(covariate_formula)) {
        msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, 
            qmatrix = crude_result$qmat, covariates = as.formula(covariate_formula),
            control = list(fnscale = 10000, maxit = 1000)
        )
      } else {
        msm(state_num ~ DaysSinceEntry, subject = deid_enc_id, data = model_data, 
            qmatrix = crude_result$qmat, control = list(fnscale = 10000, maxit = 1000)
        )
      }
    }, error = function(e) {
      warning(paste("Error fitting msm for", modelname, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(fitted_model)) {
      fitted_models[[modelname]] <- fitted_model
    }
  }
  
  return(fitted_models)
}

tidy_msm_models <- function(fitted_msm_models) {
  tidied_models <- data.frame()
  for (modelname in names(fitted_msm_models)) {
    fitted_model <- fitted_msm_models[[modelname]]
    
    if (is.null(fitted_model)) {
      warning(paste("Skipping tidying for", modelname, "- model is NULL"))
      next
    }
    
    model_tidy <- tryCatch({
      tidy(fitted_model) %>%
        mutate(model = modelname)  # Add the model name as a column
    }, error = function(e) {
      warning(paste("Error tidying model for", modelname, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(model_tidy)) {
      tidied_models <- bind_rows(tidied_models, model_tidy)
    }
  }
  
  return(tidied_models)
}

tidy_msm_pmats <- function(fitted_msm_models) {
  tidied_pmats <- data.frame()
  for (modelname in names(fitted_msm_models)) {
    fitted_model <- fitted_msm_models[[modelname]]
    
    if (is.null(fitted_model)) {
      warning(paste("Skipping tidying for", modelname, "- model is NULL"))
      next
    }
    
    model_tidy <- tryCatch({
      pmat_to_tib(pmatrix.msm(fitted_model, ci = "normal"), model = modelname, cov_name = NA, cov_value = NA)
    }, error = function(e) {
      warning(paste("Error tidying model for", modelname, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(model_tidy)) {
      tidied_pmats <- bind_rows(tidied_pmats, model_tidy)
    }
  }
  
  return(tidied_pmats)
}

# Create function to convert P matrix to tibble
pmat_to_tib <- function(pmatrix, model = NA, cov_name = NA, cov_value = NA) {
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
      cov_value = cov_value
    )
  
  return(tibble_data)
}

# Parallelize this outer function
fit_univariate_models <- function(patient_data, crude_rates, covariates) {
  # Parallel processing of covariates
  univariate_models <- mclapply(covariates, function(cov) {
    # Fit model with single covariate
    models <- fit_msm_models(patient_data, crude_rates, covariates = c(cov))
    return(models)
  }, mc.cores = n.cores)
  
  # Name the list elements
  names(univariate_models) <- covariates
  
  return(univariate_models)
}

tidy_msm_models_univ <- function(fitted_msm_models) {
  tidied_models <- data.frame()
  
  # Handle nested structure from fit_univariate_models
  for (covariate in names(fitted_msm_models)) {
    for (modelname in names(fitted_msm_models[[covariate]])) {
      fitted_model <- fitted_msm_models[[covariate]][[modelname]]
      if (is.null(fitted_model)) {
        warning(paste("Skipping tidying for", covariate, modelname, "- model is NULL"))
        next
      }
      
      model_tidy <- tryCatch({
        tidy(fitted_model) %>%
          mutate(
            model = modelname,
            covariate = covariate
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
  return(tidied_models)
}

tidy_msm_pmats_univ <- function(fitted_msm_models, cov_value) {
  tidied_pmats <- data.frame()
  for (covariate in names(fitted_msm_models)) {
    for (modelname in names(fitted_msm_models[[covariate]])) {
      fitted_model <- fitted_msm_models[[covariate]][[modelname]]
      if (is.null(fitted_model)) {
        warning(paste("Skipping tidying for", covariate, modelname, "- model is NULL"))
        next
      }
      
      model_tidy <- tryCatch({
        pmat_to_tib(pmatrix.msm(fitted_model, ci = "normal"), model = modelname, cov_name = covariate, cov_value = cov_value)
      }, error = function(e) {
        warning(paste("Error tidying model for", modelname, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(model_tidy)) {
        tidied_pmats <- bind_rows(tidied_pmats, model_tidy)
      }
    }
  }
  return(tidied_pmats)
}

tidy_msm_hr_univ <- function(fitted_msm_models, hr_scale = 1) {
  tidied_pmats <- data.frame()
  for (covariate in names(fitted_msm_models)) {
    for (modelname in names(fitted_msm_models[[covariate]])) {
      fitted_model <- fitted_msm_models[[covariate]][[modelname]]
      if (is.null(fitted_model)) {
        warning(paste("Skipping tidying for", covariate, modelname, "- model is NULL"))
        next
      }
      
      model_tidy <- tryCatch({
        model_tidy <- as.data.frame(hazard.msm(fitted_model, hazard.scale = 10, cl = 0.95)) %>%
          rename_with(~ sub("^[^.]*\\.(.*)$", "\\1", .x), .cols = everything()) %>% 
          rownames_to_column(var = "transition") %>%
          mutate(
            model = modelname,
            covariate = covariate,
            hr_scale = hr_scale
          )
      }, error = function(e) {
        warning(paste("Error tidying model for", modelname, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(model_tidy)) {
        tidied_pmats <- bind_rows(tidied_pmats, model_tidy)
      }
    }
  }
  return(tidied_pmats)
}
