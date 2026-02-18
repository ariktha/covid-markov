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

