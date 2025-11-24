date_fn <- function(x) as.Date(as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S"))

`%nin%` <- Negate(`%in%`)

tbl_fn <- function(x){
  tbl <- table(x)
  res <- cbind(tbl, round(prop.table(tbl)*100, 2))
  colnames(res) <- c('Count', 'Percentage')
  res
}

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
